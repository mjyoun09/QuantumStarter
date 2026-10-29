"""DFT-D4 dispersion from the reference dftd4 library (Grimme group).

D4: E. Caldeweyher et al., J. Chem. Phys. 150, 154122 (2019). Charge-dependent
C6 from EEQ partial charges, rational (Becke-Johnson) damping and the
Axilrod-Teller-Muto three-body term (s9). The correction depends only on
atomic numbers, positions and the total charge, never on the SCF density,
so it is added to the converged SCF energy as a separate component.
"""
import numpy as np


def d4_energy(atomic_numbers, coords_bohr, charge, xc):
    try:
        import dftd4
        from dftd4.interface import DampingParam, DispersionModel
        from dftd4.parameters import get_damping_param
    except ImportError as exc:
        raise RuntimeError('D4 분산 보정에는 dftd4 패키지가 필요합니다. '
                           'wsl --exec .venv-dft/bin/python -m pip install -r requirements-dispersion.txt') from exc
    # Use the published parameters explicitly so the recorded values are the ones applied.
    parameters = {key: float(value) for key, value in get_damping_param(xc).items()}
    coords = np.asarray(coords_bohr, dtype=float)
    model = DispersionModel(np.asarray(atomic_numbers, dtype=int), coords, charge=float(charge))
    # dftd4 truncates pairs beyond 60 bohr (two-body) / 40 bohr (three-body) by default, a
    # speed approximation meant for large/periodic systems. An isolated molecule has no
    # images, so include every pair and triple; the coordination-number cutoff stays default.
    span = float(np.max(np.linalg.norm(coords[:, None] - coords[None], axis=-1)))
    cutoff = max(60.0, span + 1.0)
    model.set_realspace_cutoff(disp2=cutoff, disp3=cutoff, cn=30.0)
    energy = float(model.get_dispersion(DampingParam(**parameters), grad=False)['energy'])
    properties = model.get_properties()
    return {'model': 'D4', 'energy_hartree': energy, 'functional_parameters_for': xc,
            'damping': 'rational (Becke-Johnson) two-body + Axilrod-Teller-Muto three-body',
            'parameters': parameters, 'charge_model': 'EEQ', 'total_charge': int(charge),
            'eeq_partial_charges': [float(q) for q in properties['partial charges']],
            'coordination_numbers': [float(c) for c in properties['coordination numbers']],
            'realspace_cutoff_bohr': {'two_body': cutoff, 'three_body': cutoff, 'coordination_number': 30.0},
            'pair_truncation': 'none: cutoff exceeds the largest interatomic distance',
            'library': 'dftd4', 'library_version': dftd4.__version__,
            'reference': 'E. Caldeweyher et al., J. Chem. Phys. 150, 154122 (2019)'}
