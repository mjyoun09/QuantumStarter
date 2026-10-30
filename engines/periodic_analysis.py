"""Post-SCF analysis of periodic Kohn-Sham results: density of states and band structure.

DOS: Gaussian-broadened k-mesh eigenvalues, weighted by the k-point weights, in
states/eV per input cell and spin channel (a restricted orbital counts twice).
Bands: non-self-consistent eigenvalues on a high-symmetry path from ASE's
bandpath (2D path for slabs), evaluated with the converged density (PySCF
get_bands). Each analysis records its own consistency check.
"""
import numpy as np


def _channels(values, unrestricted):
    array = np.asarray(values, dtype=float)
    return array if unrestricted else array[None]


def frontier(mo_energy, mo_occ, unrestricted):
    """Highest occupied / lowest unoccupied eigenvalue over all k (Hartree)."""
    energies, occupations = _channels(mo_energy, unrestricted), _channels(mo_occ, unrestricted)
    half = 0.5 if unrestricted else 1.0
    occupied, empty = energies[occupations > half], energies[occupations <= half]
    return (float(occupied.max()) if occupied.size else None, float(empty.min()) if empty.size else None)


def density_of_states(mo_energy, mo_occ, weights, unrestricted, sigma_ev, hartree_to_ev):
    energies = _channels(mo_energy, unrestricted) * hartree_to_ev  # (spin, k, band)
    weights = np.asarray(weights, dtype=float)
    degeneracy = 1.0 if unrestricted else 2.0
    grid = np.arange(energies.min() - 5 * sigma_ev, energies.max() + 5 * sigma_ev, sigma_ev / 5)
    dos = []
    for channel in energies:
        gauss = np.exp(-0.5 * ((grid[:, None, None] - channel[None]) / sigma_ev)**2) / (sigma_ev * np.sqrt(2 * np.pi))
        dos.append(degeneracy * np.einsum('ekn,k->e', gauss, weights))
    dos = np.asarray(dos)
    step = grid[1] - grid[0]
    homo, lumo = frontier(mo_energy, mo_occ, unrestricted)
    occupation_sum = float(np.einsum('skn,k->', _channels(mo_occ, unrestricted), weights))
    return {
        'method': 'Gaussian-broadened k-mesh eigenvalues', 'sigma_ev': sigma_ev, 'units': 'states/eV per input cell',
        'energy_grid_ev': grid.tolist(), 'dos_by_spin': dos.tolist(),
        'homo_ev': None if homo is None else homo * hartree_to_ev,
        'lumo_ev': None if lumo is None else lumo * hartree_to_ev,
        'mesh_gap_ev': None if homo is None or lumo is None else max(lumo - homo, 0.0) * hartree_to_ev,
        # Checks: the integral equals the number of states; weighted occupations equal the electron count.
        'integrated_states': float(dos.sum() * step),
        'expected_states': float(degeneracy * energies.shape[0] * energies.shape[2]),
        'occupied_electrons_from_weights': occupation_sum}


def band_structure(mf, cell, config, unrestricted, hartree_to_ev):
    from ase.cell import Cell
    pbc = (True, True, False) if config.boundary == 'slab' else (True, True, True)
    path = None if config.band_path == 'auto' else config.band_path
    bandpath = Cell(np.asarray(config.lattice)).bandpath(path, npoints=config.band_points, pbc=pbc)
    absolute = cell.get_abs_kpts(np.asarray(bandpath.kpts))
    energies = _channels(mf.get_bands(absolute)[0], unrestricted)  # (spin, k, band), Hartree
    x, special_x, labels = bandpath.get_linear_kpoint_axis()
    # Consistency: bands evaluated at the SCF mesh must reproduce the SCF eigenvalues.
    mesh = _channels(mf.get_bands(np.asarray(mf.kpts))[0], unrestricted)
    consistency = float(np.abs(mesh - _channels(mf.mo_energy, unrestricted)).max()) * hartree_to_ev
    result = {'method': 'non-self-consistent eigenvalues on an ASE high-symmetry path (PySCF get_bands)',
              'path': bandpath.path, 'pbc': list(pbc), 'kpoints_scaled': np.asarray(bandpath.kpts).tolist(),
              'axis_inverse_angstrom': np.asarray(x).tolist(), 'special_axis': np.asarray(special_x).tolist(),
              'special_labels': list(labels), 'bands_ev_by_spin': (energies * hartree_to_ev).tolist(),
              'mesh_consistency_max_ev': consistency}
    if not config.smearing_ev:  # Integer occupations: band index of the highest occupied state is known.
        counts = cell.nelec if unrestricted else (cell.nelectron // 2,)
        gaps = []
        for channel, nocc in zip(energies, counts):
            if 0 < nocc < channel.shape[1]:
                vbm, cbm = channel[:, nocc - 1].max(), channel[:, nocc].min()
                gaps.append({'vbm_ev': vbm * hartree_to_ev, 'cbm_ev': cbm * hartree_to_ev,
                             'gap_ev': (cbm - vbm) * hartree_to_ev,
                             'direct_gap_ev': float((channel[:, nocc] - channel[:, nocc - 1]).min()) * hartree_to_ev,
                             'vbm_k_index': int(channel[:, nocc - 1].argmax()),
                             'cbm_k_index': int(channel[:, nocc].argmin())})
        result['gaps_by_spin'] = gaps
    return result
