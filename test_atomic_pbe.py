"""PBE pointwise oracle, variational derivative, SCF, and CLI regression checks."""
import json
from decimal import Decimal, localcontext
from pathlib import Path
import tempfile
import unittest

import numpy as np

from atomic_pbe import pbe_pointwise, pbe_radial, radial_gradient, radial_gradient_transpose
from radial_atom import AtomConfig, log_grid, run_atom, load_atom_config, save_atom_result
from test_radial_atom import cli

try:
    from pyscf.dft import libxc
except ImportError:
    libxc = None


def sample(spin, seed=917, extreme=False):
    rng = np.random.default_rng(seed)
    n = 10.0 ** rng.uniform(-6, 3, (1 + spin, 180))
    if spin and not extreme:
        total = n[0].copy()
        fraction = rng.uniform(0.01, 0.99, total.shape)
        n = np.array([total * fraction, total * (1 - fraction)])
    g = rng.normal(size=(1 + spin, 3, 180)) * n[:, None, :] ** (4 / 3)
    g[:, :, :10] = 0
    sigma = np.sum(g * g, axis=1) if not spin else np.array(
        [np.sum(g[0] ** 2, axis=0), np.sum(g[0] * g[1], axis=0), np.sum(g[1] ** 2, axis=0)])
    return n, g, sigma


class PbeLocalTests(unittest.TestCase):
    def test_extreme_spin_partials_against_60_digit_energy_differences(self):
        # Libxc 7 uses 1 +/- zeta in these regions: at n_up/n_down~1e-9
        # subtraction loses ~7 digits in the small spin fraction. An independent
        # high-precision evaluation of the ORIGINAL eps+H formula decides which
        # implementation is accurate without relaxing the ordinary oracle test.
        n, _, sigma = sample(1, extreme=True)
        terms = pbe_pointwise(n, sigma)
        ratio = np.min(n, axis=0) / np.max(n, axis=0)
        indices = [i for i in np.argsort(ratio) if i >= 10][:6]
        with localcontext() as context:
            context.prec = 60
            D = Decimal
            pi = D('3.14159265358979323846264338327950288419716939937510582097494459')
            one, two, three = D(1), D(2), D(3)
            beta = D('0.06672455060314922')
            gamma = (one - two.ln()) / pi ** 2
            mu, kappa = beta * pi ** 2 / 3, D('.804')

            def energy(v, term):
                up, down, suu, sud, sdd = v
                def exchange(density, s):
                    t = s / (4 * (3 * pi ** 2) ** (two / three) * density ** (D(8) / 3))
                    return -D('.75') * (3 / pi) ** (one / three) * density ** (D(4) / 3) * (1 + kappa - kappa / (1 + mu / kappa * t))
                if term == 'x':
                    return (exchange(2 * up, 4 * suu) + exchange(2 * down, 4 * sdd)) / 2
                total, z = up + down, (up - down) / (up + down)
                rs = (3 / (4 * pi * total)) ** (one / three)
                def g(a, alpha, b1, b2, b3, b4):
                    a, alpha, b1, b2, b3, b4 = map(D, (a, alpha, b1, b2, b3, b4))
                    q = 2 * a * (b1 * rs.sqrt() + b2 * rs + b3 * rs ** D('1.5') + b4 * rs ** 2)
                    return -2 * a * (1 + alpha * rs) * (1 + 1 / q).ln()
                g0 = g('.0310907', '.21370', '7.5957', '3.5876', '1.6382', '.49294')
                g1 = g('.01554535', '.20548', '14.1189', '6.1977', '3.3662', '.62517')
                g2 = g('.0168869', '.11125', '10.357', '3.6231', '.88026', '.49671')
                f = ((1 + z) ** (D(4) / 3) + (1 - z) ** (D(4) / 3) - 2) / (two ** (D(4) / 3) - 2)
                f20 = D('1.709920934161365617563962776245')
                ec = g0 + (g1 - g0) * f * z ** 4 - g2 * f * (1 - z ** 4) / f20
                phi = ((1 + z) ** (two / three) + (1 - z) ** (two / three)) / 2
                t2 = (suu + 2 * sud + sdd) / (16 / pi * (3 * pi ** 2) ** (one / three) * phi ** 2 * total ** (D(7) / 3))
                a = beta / gamma / ((-ec / (gamma * phi ** 3)).exp() - 1)
                h = gamma * phi ** 3 * (1 + beta / gamma * t2 * (1 + a * t2) / (1 + a * t2 + a ** 2 * t2 ** 2)).ln()
                return total * (ec + h)

            for index in indices:
                variables = list(map(lambda x: D(str(x)), np.concatenate((n[:, index], sigma[:, index]))))
                for term in ('x', 'c'):
                    self.assertAlmostEqual(float(energy(variables, term) / sum(variables[:2])),
                                           terms['eps_' + term][index], delta=1e-12)
                    derivatives = np.concatenate((terms['vrho_' + term][:, index], terms['vsigma_' + term][:, index]))
                    for axis, expected in enumerate(derivatives):
                        step = abs(variables[axis]) * D('1e-18')
                        high, low = variables.copy(), variables.copy()
                        high[axis] += step
                        low[axis] -= step
                        actual = float((energy(high, term) - energy(low, term)) / (2 * step))
                        np.testing.assert_allclose(expected, actual, rtol=3e-12, atol=1e-14)

    @unittest.skipIf(libxc is None, 'Libxc oracle requires WSL .venv-dft')
    def test_libxc_energy_and_all_partials(self):
        for spin in (0, 1):
            n, gradient, sigma = sample(spin)
            terms = pbe_pointwise(n, sigma)
            rho = np.concatenate((n[:, None, :], gradient), axis=1)
            for name, code in (('x', 'GGA_X_PBE,'), ('c', ',GGA_C_PBE'), ('xc', 'PBE')):
                exc, vxc = libxc.eval_xc(code, rho if spin else rho[0], spin=spin, deriv=1)[:2]
                with self.subTest(spin=spin, term=name):
                    np.testing.assert_allclose(terms['eps_' + name], exc, rtol=3e-12, atol=2e-15)
                    np.testing.assert_allclose(terms['vrho_' + name], vxc[0].T if spin else vxc[0][None],
                                               rtol=3e-12, atol=2e-14)
                    # PBE's uniform unpolarized x and c gradient coefficients
                    # cancel. Scale error by the terms, not their near-zero sum.
                    ref = vxc[1].T if spin else vxc[1][None]
                    scale = np.abs(terms['vsigma_x']) + np.abs(terms['vsigma_c'])
                    self.assertTrue(np.all(np.abs(terms['vsigma_' + name] - ref) <= 3e-12 * scale + 1e-15))

    @unittest.skipIf(libxc is None, 'Libxc oracle requires WSL .venv-dft')
    def test_zero_gradient_correlation_is_modified_pw92(self):
        for spin in (0, 1):
            n, _, sigma = sample(spin)
            terms = pbe_pointwise(n, sigma * 0)
            exc, vxc = libxc.eval_xc(',LDA_C_PW_MOD', n if spin else n[0], spin=spin, deriv=1)[:2]
            np.testing.assert_allclose(terms['eps_c'], exc, rtol=3e-13)
            np.testing.assert_allclose(terms['vrho_c'], vxc[0].T if spin else vxc[0][None], rtol=3e-12)

    def test_zero_density_polarized_limits_and_spin_symmetry(self):
        n = np.array([[0, 1e-35, 1e-29, 1e-8, 1, 1000.]])
        s = (2 * n) ** 2
        plain, split = pbe_pointwise(n, s), pbe_pointwise(np.repeat(n / 2, 2, axis=0), np.repeat(s / 4, 3, axis=0))
        np.testing.assert_allclose(plain['eps_xc'], split['eps_xc'], rtol=1e-13, atol=1e-25)
        np.testing.assert_allclose(plain['vrho_xc'][0], split['vrho_xc'][0], rtol=1e-12, atol=1e-25)
        polarized = pbe_pointwise(np.concatenate((n, n * 0)), np.concatenate((s, s * 0, s * 0)))
        swapped = pbe_pointwise(np.concatenate((n * 0, n)), np.concatenate((s * 0, s * 0, s)))
        for key, value in polarized.items():
            self.assertTrue(np.all(np.isfinite(value)), key)
        np.testing.assert_array_equal(polarized['eps_xc'], swapped['eps_xc'])
        np.testing.assert_array_equal(polarized['vrho_xc'], swapped['vrho_xc'][::-1])
        # A squared minority gradient can underflow while the cross product is
        # still representable (Cr's far tail). The physical Gram matrix is valid.
        tiny = pbe_pointwise([[1e-10], [1e-200]], [[1e-20], [1e-210], [0]])
        self.assertTrue(np.all(np.isfinite(tiny['vrho_xc'])))
        with self.assertRaises(ValueError):
            pbe_pointwise([[1.0], [1.0]], [[1.0], [2.0], [1.0]])

    def test_pointwise_directional_derivative(self):
        for spin in (0, 1):
            n, gradient, sigma = sample(spin)
            rng = np.random.default_rng(981)
            dn = n * rng.normal(size=n.shape)
            dg = gradient * rng.normal(size=gradient.shape)
            terms = pbe_pointwise(n, sigma)
            ds = 2 * np.sum(gradient * dg, axis=1) if not spin else np.array(
                [2 * np.sum(gradient[0] * dg[0], axis=0),
                 np.sum(gradient[0] * dg[1] + gradient[1] * dg[0], axis=0),
                 2 * np.sum(gradient[1] * dg[1], axis=0)])
            analytic = np.sum(terms['vrho_xc'] * dn, axis=0) + np.sum(terms['vsigma_xc'] * ds, axis=0)
            energies = []
            for sign in (1, -1):
                h = sign * 1e-5
                nh, gh = n + h * dn, gradient + h * dg
                sh = np.sum(gh * gh, axis=1) if not spin else np.array(
                    [np.sum(gh[0] ** 2, axis=0), np.sum(gh[0] * gh[1], axis=0), np.sum(gh[1] ** 2, axis=0)])
                energies.append(nh.sum(axis=0) * pbe_pointwise(nh, sh)['eps_xc'])
            np.testing.assert_allclose((energies[0] - energies[1]) / 2e-5, analytic, rtol=2e-7, atol=1e-10)


class PbeGridTests(unittest.TestCase):
    def test_gradient_transpose_including_boundaries(self):
        grid = log_grid(AtomConfig('He', grid_step=0.01))
        rng = np.random.default_rng(23)
        a, b = rng.normal(size=(2, 2, grid.r.size))
        np.testing.assert_allclose(np.sum(a * radial_gradient(b, grid)),
                                   np.sum(radial_gradient_transpose(a, grid) * b), rtol=1e-14)

    def test_discrete_potential_is_energy_derivative(self):
        grid = log_grid(AtomConfig('He', grid_step=0.007))
        r = grid.r
        for polarized in (False, True):
            n = np.array([np.exp(-2 * r)])
            if polarized:
                n = np.concatenate((n, 0.7 * np.exp(-1.5 * r)[None]))
            # Exercise both smooth perturbations and arbitrary grid perturbations;
            # include endpoints so that incorrect boundary transposes fail.
            rng = np.random.default_rng(111)
            for direction in (n * np.cos(r)[None], n * rng.normal(size=n.shape)):
                weights = 4 * np.pi * r ** 2 * grid.weights
                terms = pbe_radial(n, grid)
                analytic = np.sum(weights * terms['v_xc'] * direction)
                h = 1e-6
                def energy(d):
                    return np.sum(weights * d.sum(axis=0) * pbe_radial(d, grid)['eps_xc'])
                finite = (energy(n + h * direction) - energy(n - h * direction)) / (2 * h)
                self.assertAlmostEqual(analytic, finite, delta=2e-7 * max(abs(analytic), 1))


class PbeScfTests(unittest.TestCase):
    def test_he_second_order_and_polarized_equivalence(self):
        runs = [run_atom(AtomConfig('He', functional='pbe', grid_step=h)) for h in (0.004, 0.002, 0.001)]
        for result in runs:
            self.assertTrue(result.converged)
            self.assertAlmostEqual(result.electron_count, 2, delta=1e-10)
            self.assertAlmostEqual(result.components['total_direct'], result.total_energy, delta=1e-11)
        ratio = (runs[0].total_energy - runs[1].total_energy) / (runs[1].total_energy - runs[2].total_energy)
        self.assertAlmostEqual(ratio, 4, delta=0.15)
        extrapolated = (4 * runs[-1].total_energy - runs[-2].total_energy) / 3
        # External PySCF RKS (basis ratio 1.3), docs/VALIDATION_KO.md 0.9.7.
        self.assertAlmostEqual(extrapolated, -2.892934865842369, delta=1e-6)
        polarized = run_atom(AtomConfig('He', functional='pbe', spin_polarized=True))
        self.assertTrue(polarized.converged)
        self.assertAlmostEqual(polarized.total_energy, runs[-1].total_energy, delta=1e-11)

    def test_open_p_and_filled_d_shell_reference_energies(self):
        # B: fractional PySCF UKS. Kr: independent Gaussian Galerkin RKS.
        # The references use different orbitals, quadratures and XC code.
        for symbol, polarized, reference in (('B', True, -24.60538130018798),
                                             ('Kr', False, -2753.416108819394)):
            with self.subTest(element=symbol):
                coarse, fine = [run_atom(AtomConfig(symbol, functional='pbe', spin_polarized=polarized,
                                                    grid_step=h)) for h in (.002, .001)]
                self.assertTrue(coarse.converged and fine.converged)
                self.assertAlmostEqual((4 * fine.total_energy - coarse.total_energy) / 3,
                                       reference, delta=1e-6)

    def test_nonconvergence_and_configuration(self):
        with self.assertRaises(ValueError):
            AtomConfig(functional='PBE0')
        result = run_atom(AtomConfig('N', functional='pbe', spin_polarized=True, max_iterations=2))
        self.assertFalse(result.converged)
        with tempfile.TemporaryDirectory() as temp:
            saved = save_atom_result(result, temp)
            self.assertEqual(load_atom_config(Path(temp) / 'settings.json'), result.config)
            legacy = Path(temp) / 'legacy.json'
            legacy.write_text('{"model": "atom", "element": "He"}', encoding='utf-8')
            self.assertEqual(load_atom_config(legacy), AtomConfig('He'))
        self.assertEqual(saved['xc']['correlation']['libxc_id'], 130)
        self.assertEqual(saved['status'], 'not_converged')

    def test_cli_selection_roundtrip_and_conflicts(self):
        with tempfile.TemporaryDirectory() as temp:
            done = cli('--model', 'atom', '--element', 'He', '--functional', 'pbe', '--no-show', '--output', temp)
            self.assertEqual(done.returncode, 0, done.stderr)
            self.assertIn('PBE', done.stdout)
            (run,) = Path(temp).iterdir()
            summary = json.loads((run / 'results.json').read_text(encoding='utf-8'))
            self.assertEqual(summary['config']['functional'], 'pbe')
            self.assertEqual(summary['xc']['exchange']['libxc_id'], 101)
            self.assertGreater((run / 'atom.png').stat().st_size, 1000)
            again = cli('--model', 'atom', '--config', str(run / 'settings.json'), '--no-show', '--output', temp)
            self.assertEqual(again.returncode, 0, again.stderr)
            repeated = json.loads((sorted(Path(temp).iterdir())[-1] / 'results.json').read_text(encoding='utf-8'))
            self.assertEqual(repeated['total_energy_hartree'], summary['total_energy_hartree'])
            stopped = cli('--model', 'atom', '--element', 'He', '--functional', 'pbe', '--max-iterations', '1',
                          '--no-show', '--output', temp)
            self.assertEqual(stopped.returncode, 2, stopped.stderr)
        for args in (('--model', 'he-lda'), ('--hydrogen-check',), ()):
            self.assertNotEqual(cli(*args, '--functional', 'pbe', '--no-show').returncode, 0)


if __name__ == '__main__':
    unittest.main()
