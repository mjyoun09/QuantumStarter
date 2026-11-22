"""Real closed-shell singlet TDA using native Coulomb, XC and hybrid exchange response."""
import csv
import numpy as np
from scipy.constants import c, h, electron_volt
from scipy.linalg import eigh
from scipy.sparse.linalg import LinearOperator, eigsh

from .native_response import FockResponse
from .pyscf_engine import HARTREE_TO_EV


def singlet_spectrum(mf, requested, directory):
    occupied = np.flatnonzero(mf.mo_occ == 2)
    virtual = np.flatnonzero(mf.mo_occ == 0)
    co, cv = mf.mo_coeff[:, occupied], mf.mo_coeff[:, virtual]
    gaps = mf.mo_energy[virtual] - mf.mo_energy[occupied, None]
    dimension = gaps.size
    if not dimension:
        raise ValueError('TDA에는 점유 궤도와 가상 궤도가 각각 하나 이상 필요합니다.')
    count = min(requested, dimension)
    response = FockResponse(mf, mf.make_rdm1())

    def multiply(vector):
        amplitudes = vector.reshape(gaps.shape)
        transition = 2 * cv @ amplitudes.T @ co.T
        induced = response(transition, transition=True)
        return (gaps * amplitudes + (cv.T @ induced @ co).T).ravel()

    if dimension <= max(32, count+1):
        matrix = np.column_stack([multiply(v) for v in np.eye(dimension)])
        values, vectors = eigh((matrix + matrix.T)/2, subset_by_index=(0, count-1))
    else:
        operator = LinearOperator((dimension, dimension), matvec=multiply, dtype=float)
        values, vectors = eigsh(operator, k=count, which='SA', tol=1e-10,
                                v0=np.random.default_rng(2026).normal(size=dimension))
    if values[0] <= 0:
        raise ValueError('비양의 TDA 들뜸 에너지가 있습니다. 바닥상태 SCF의 안정성을 확인하세요.')
    integrals = mf.mol.intor_symmetric('int1e_r', comp=3)
    dipoles = np.einsum('pi,xpq,qa->xia', co, integrals, cv, optimize=True)
    peaks = []
    for state, (omega, vector) in enumerate(zip(values, vectors.T), 1):
        x = vector.reshape(gaps.shape)
        mu = np.sqrt(2.) * np.einsum('ia,kia->k', x, dipoles)
        energy_ev = float(omega * HARTREE_TO_EV)
        leading = []
        for index in np.argsort(vector**2)[::-1][:3]:
            i, a = np.unravel_index(index, x.shape)
            leading.append(dict(occupied=int(occupied[i])+1, virtual=int(virtual[a])+1, weight=float(vector[index]**2)))
        peaks.append(dict(state=state, energy_eV=energy_ev, wavelength_nm=h*c*1e9/(energy_ev*electron_volt),
                          oscillator_strength=float((2/3)*omega*(mu@mu)), transition_dipole_au=mu.tolist(),
                          leading_transitions=leading))
    with (directory/'uvvis.csv').open('w', newline='', encoding='utf-8') as stream:
        writer = csv.writer(stream, lineterminator='\n')
        writer.writerow(['state','energy_eV','wavelength_nm','oscillator_strength','leading_transitions (MO:MO:weight)'])
        for row in peaks:
            labels = ' '.join(f"{p['occupied']}:{p['virtual']}:{p['weight']:.6f}" for p in row['leading_transitions'])
            writer.writerow([row[k] for k in ('state','energy_eV','wavelength_nm','oscillator_strength')] + [labels])
    return dict(method='native closed-shell singlet Tamm-Dancoff approximation', gauge='length',
                requested_states=requested, excitation_dimension=dimension, peaks=peaks, table_file='uvvis.csv')
