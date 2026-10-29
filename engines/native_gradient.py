"""Analytic nuclear gradient of the native Kohn-Sham energy; no PySCF gradient driver.

dE/dR_A = sum D dh/dR_A + dE_J/dR_A + dE_K/dR_A + dE_xc/dR_A - sum W dS/dR_A + dV_nn/dR_A

External primitives: PySCF/libcint derivative integrals (int1e_ipovlp, int1e_ipkin,
int1e_ipnuc, int1e_iprinv, int2e_ip1), AO values/derivatives on the grid, Libxc
functional derivatives. PySCF 'ip' integrals differentiate the bra function with
respect to the electron coordinate; a function centred on atom A obeys
d(phi)/dR_A = -grad_r(phi). Grid points and weights stay fixed (no grid-response
term), which is also PySCF's default; the result is then the exact derivative of
the energy evaluated on a frozen grid.
"""
import numpy as np
from scf_control import check_cancelled

# PySCF AO second-derivative component for (d/dx_i)(d/dx_j): xx xy xz yy yz zz -> 4..9
SECOND = ((4, 5, 6), (5, 7, 8), (6, 8, 9))


def native_gradient(mf):
    """Gradient components (natm, 3) in Hartree/bohr for a converged NativeKS object."""
    mol, n = mf.mol, mf.mol.nao_nr()
    slices = mol.aoslice_by_atom()
    ao_atom = np.empty(n, dtype=int)
    for atom, (_, _, p0, p1) in enumerate(slices):
        ao_atom[p0:p1] = atom
    shell_sizes = np.diff(mol.ao_loc_nr())
    if 3 * 3 * int(shell_sizes.max()) * n**3 * 8 + 64 * 1024**2 > mf.config.memory_mb * 1024**2:
        raise ValueError('자체 기울기의 2전자 미분 블록이 메모리 목표를 넘습니다. memory_mb를 늘리거나 작은 기저를 쓰세요.')
    dm = mf._dm
    spin_dms = dm if mf.unrestricted else dm[None]
    total = dm.sum(axis=0) if mf.unrestricted else dm
    _, fock, _ = mf.evaluate(dm)
    # Energy-weighted density for idempotent converged densities:
    # restricted D = 2 C C^T -> W = 2 C eps C^T = D F D / 2; unrestricted W = sum_s D_s F_s D_s.
    weighted = (sum(d @ f @ d for d, f in zip(dm, fock)) if mf.unrestricted else total @ fock @ total / 2)
    parts = {key: np.zeros((mol.natm, 3)) for key in
             ('one_electron', 'coulomb', 'exact_exchange', 'xc', 'overlap_pulay', 'nuclear_repulsion')}

    # One-electron: basis-function motion (rows of A) plus motion of nucleus A in -Z_A/|r-R_A|.
    ip_overlap = mol.intor('int1e_ipovlp', comp=3)
    basis_motion = -(mol.intor('int1e_ipkin', comp=3) + mol.intor('int1e_ipnuc', comp=3))
    for atom, (_, _, p0, p1) in enumerate(slices):
        with mol.with_rinv_at_nucleus(atom):
            half = -mol.atom_charge(atom) * mol.intor('int1e_iprinv', comp=3)
        half[:, p0:p1] += basis_motion[:, p0:p1]
        parts['one_electron'][atom] = np.einsum('xij,ij->x', half + half.transpose(0, 2, 1), total)
        # -sum W dS/dR_A with dS/dR_A = -(ip rows of A) - transpose.
        parts['overlap_pulay'][atom] = 2 * np.einsum('xij,ij->x', ip_overlap[:, p0:p1], weighted[p0:p1])

    # Two-electron: ((d p) q | r s) one bra shell at a time keeps memory at 3 n_shell n^3.
    # E_J = 1/2 D D (pq|rs); E_K = -c hyb sum_s D_s (pr|qs) D_s with c = 1/4 (RKS) or 1/2 (UKS).
    # The four index positions contribute equally, so dE/dR_A = -4 x (bra-derivative term).
    exchange_scale = .5 if mf.unrestricted else .25
    ao_loc = mol.ao_loc_nr()
    for shell in range(mol.nbas):
        check_cancelled(mf.directory)
        p0, p1 = ao_loc[shell], ao_loc[shell + 1]
        atom = mol.bas_atom(shell)
        block = mol.intor('int2e_ip1', comp=3, aosym='s1',
                          shls_slice=(shell, shell + 1, 0, mol.nbas, 0, mol.nbas, 0, mol.nbas))
        coulomb = np.tensordot(block, total, axes=([3, 4], [0, 1]))  # x, p, q
        parts['coulomb'][atom] -= 2 * np.einsum('xpq,pq->x', coulomb, total[p0:p1])
        if mf.hybrid:
            for matrix in spin_dms:
                exchange = np.tensordot(block, matrix, axes=([2, 4], [0, 1]))  # x, p, q
                parts['exact_exchange'][atom] += (4 * exchange_scale * mf.hybrid
                                                  * np.einsum('xpq,pq->x', exchange, matrix[p0:p1]))

    # XC on the fixed grid: rho_s = sum D_s phi phi, so d rho_s/dR_A = -2 sum_{mu in A} D_s (d phi_mu) phi.
    from pyscf.dft import numint, libxc
    gga = mf.xctype == 'GGA'
    per_ao = np.zeros((3, n))
    for start in range(0, len(mf.grids.weights), 1024):
        check_cancelled(mf.directory)
        coords = mf.grids.coords[start:start + 1024]
        weights = mf.grids.weights[start:start + 1024]
        ao = numint.eval_ao(mol, coords, deriv=2 if gga else 1)
        values = ao[0]
        products = [values @ matrix for matrix in spin_dms]  # (g, mu): sum_nu D_mu,nu phi_nu
        rho = []
        for matrix, product in zip(spin_dms, products):
            density = np.einsum('gi,gi->g', product, values)
            rho.append(np.vstack((density, 2 * np.einsum('gi,xgi->xg', product, ao[1:4]))) if gga else density)
        rho = np.asarray(rho)
        _, vxc, _, _ = libxc.eval_xc(mf.xc, rho if mf.unrestricted else rho[0], spin=int(mf.unrestricted), deriv=1)
        vrho = vxc[0].T if mf.unrestricted else vxc[0][None]
        for spin, (matrix, product) in enumerate(zip(spin_dms, products)):
            potential = weights * vrho[spin]
            for x in range(3):
                per_ao[x] += np.einsum('gm,gm->m', ao[1 + x] * potential[:, None], product)
            if gga:
                # w_s = dE/d(grad rho_s): 2 v_ss grad rho_s + v_ab grad rho_other (RKS: 2 v grad rho).
                if mf.unrestricted:
                    vsigma = vxc[1].T
                    w = 2 * vsigma[0 if spin == 0 else 2] * rho[spin, 1:4] + vsigma[1] * rho[1 - spin, 1:4]
                else:
                    w = 2 * vxc[1] * rho[0, 1:4]
                w = w * weights
                gradient_products = [ao[1 + y] @ matrix for y in range(3)]
                weighted_gradient = sum(w[y][:, None] * gradient_products[y] for y in range(3))
                for x in range(3):
                    second = sum(w[y][:, None] * ao[SECOND[x][y]] for y in range(3))
                    per_ao[x] += np.einsum('gm,gm->m', second, product)
                    per_ao[x] += np.einsum('gm,gm->m', ao[1 + x], weighted_gradient)
    for atom in range(mol.natm):
        parts['xc'][atom] = -2 * per_ao[:, ao_atom == atom].sum(axis=1)

    charges, coords = mol.atom_charges(), mol.atom_coords()
    for a in range(mol.natm):
        for b in range(mol.natm):
            if a != b:
                r = coords[a] - coords[b]
                parts['nuclear_repulsion'][a] -= charges[a] * charges[b] * r / np.linalg.norm(r)**3
    return parts
