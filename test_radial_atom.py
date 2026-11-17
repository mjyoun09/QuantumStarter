"""로그 방사 격자 원자 LDA/LSDA(radial_atom.py) 검증. 실행: python -m unittest test_radial_atom -v

Windows .venv와 WSL .venv-dft에서 모두 실행한다. Libxc 대조는 PySCF 번들 Libxc가 있을 때만(WSL).
NIST 기준값은 data/nist_atomic_reference.json(원 페이지에서 파싱, scripts/build_nist_atomic_reference.py).
"""
import csv
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

import numpy as np

from atomic_dft import lda_xc
from radial_atom import (GROUND_CONFIGURATIONS, SYMBOLS, AtomConfig, hartree_potential, load_atom_config, log_grid,
                         lsda_xc, parse_configuration, run_atom, save_atom_config, save_atom_result, solve_channel,
                         spin_occupations)

ROOT = Path(__file__).resolve().parent
NIST = json.loads((ROOT / 'data' / 'nist_atomic_reference.json').read_text(encoding='utf-8'))
# 0.9.5 독립 기준(scripts/validate_he_lda.py, PySCF RKS 짝수 템퍼드 s 39개).
PYSCF_HE_LDA = -2.834835624033

try:
    from pyscf.dft import libxc
except ImportError:
    libxc = None


def cli(*arguments):
    return subprocess.run([sys.executable, str(ROOT / 'main.py'), *arguments], cwd=ROOT, capture_output=True,
                          encoding='utf-8', errors='replace', env=dict(os.environ, PYTHONIOENCODING='utf-8'),
                          timeout=600)


def nist_value(symbol, column, orbital):
    entry = NIST['elements'][symbol][column]['eigenvalues'][orbital['label']]
    return entry if column == 'LDA' else entry[orbital['spin']]


class ConfigurationTests(unittest.TestCase):
    def test_parse_cores_and_notation(self):
        argon = parse_configuration('[Ne] 3s2 3p6')
        self.assertEqual([s.label for s in argon], ['1s', '2s', '2p', '3s', '3p'])
        self.assertEqual(sum(s.occupation for s in argon), 18)
        boron = parse_configuration('1s^2 2s^2 2p^1')
        self.assertEqual([(s.n, s.l, s.occupation) for s in boron], [(1, 0, 2.0), (2, 0, 2.0), (2, 1, 1.0)])
        self.assertEqual(parse_configuration('2p0.5 1s2')[0].label, '1s')
        for bad in ('2d1', '1s3', '1s2 1s1', '[Xe] 6s1', 'abc', ''):
            with self.subTest(bad=bad), self.assertRaises(ValueError):
                parse_configuration(bad)

    def test_ground_configurations_match_nist_headers(self):
        self.assertEqual(len(NIST['elements']), 36)
        for z, symbol in enumerate(SYMBOLS, start=1):
            ours = [(s.n, s.l, s.occupation) for s in parse_configuration(GROUND_CONFIGURATIONS[z])]
            theirs = [(s.n, s.l, s.occupation) for s in parse_configuration(NIST['elements'][symbol]['configuration'])]
            self.assertEqual(ours, theirs, symbol)
            self.assertEqual(sum(o for _, _, o in ours), z)
            self.assertEqual(len(NIST['elements'][symbol]['page_sha256']), 64)

    def test_hund_spin_occupations(self):
        cases = {'2p3': (3, 0), '2p4': (3, 1), '3d6': (5, 1), '4s1': (1, 0), '3d10': (5, 5)}
        for text, expected in cases.items():
            (subshell,) = parse_configuration(text)
            self.assertEqual(tuple(q for _, q in spin_occupations(subshell, True)), expected)
            self.assertEqual(spin_occupations(subshell, False), [('unpolarized', subshell.occupation)])

    def test_invalid_config(self):
        for changes in ({'element': 'Xx'}, {'element': 'H', 'configuration': '1s2'}, {'grid_step': 0},
                        {'mixing': 1.5}, {'spin_polarized': 'yes'}, {'max_iterations': 0}, {'r_max_bohr': 0.5}):
            with self.subTest(**changes), self.assertRaises(ValueError):
                AtomConfig(**changes)
        self.assertEqual(AtomConfig('ar').element, 'Ar')
        with tempfile.TemporaryDirectory() as temp:
            path = Path(temp) / 'settings.json'
            config = AtomConfig('O', spin_polarized=True, mixing=0.4)
            save_atom_config(config, path)
            self.assertEqual(load_atom_config(path), config)
            path.write_text(json.dumps({'model': 'he-lda'}), encoding='utf-8')
            with self.assertRaises(ValueError):
                load_atom_config(path)


class SpinPolarizedXcTests(unittest.TestCase):
    rng = np.random.default_rng(20261117)
    totals = 10.0 ** rng.uniform(-6, 2, 200)
    zetas = np.concatenate([rng.uniform(-1, 1, 196), [-1.0, 1.0, 0.0, 0.5]])
    up, down = totals * (1 + zetas) / 2, totals * (1 - zetas) / 2

    def test_unpolarized_limit_equals_lda(self):
        n = self.totals
        spin, plain = lsda_xc(n / 2, n / 2), lda_xc(n)
        np.testing.assert_allclose(spin['eps_xc'], plain['eps_xc'], rtol=1e-13)
        np.testing.assert_allclose(spin['v_xc'][0], plain['v_xc'], rtol=1e-13)
        np.testing.assert_allclose(spin['v_xc'][1], plain['v_xc'], rtol=1e-13)

    def test_potentials_are_spin_density_derivatives(self):
        energy = lambda a, b: (a + b) * lsda_xc(a, b)['eps_xc']  # noqa: E731
        terms = lsda_xc(self.up, self.down)
        for index, (base, other) in enumerate(((self.up, self.down), (self.down, self.up))):
            step = 1e-6 * (base + other)
            plus = (base + step, other) if index == 0 else (other, base + step)
            minus = (np.maximum(base - step, 0.0), other) if index == 0 else (other, np.maximum(base - step, 0.0))
            width = base + step - np.maximum(base - step, 0.0)
            numeric = (energy(*plus) - energy(*minus)) / width
            interior = base > 2 * step     # 한쪽 스핀 밀도가 0인 경계는 한쪽 차분이라 제외
            np.testing.assert_allclose(terms['v_xc'][index][interior], numeric[interior], rtol=1e-6)

    @unittest.skipIf(libxc is None, 'Libxc 대조에는 PySCF(번들 Libxc)가 필요합니다(WSL .venv-dft).')
    def test_matches_libxc_spin_polarized(self):
        # Libxc 7은 1e-15 bohr^-3보다 작은 스핀 밀도를 1e-15로 바꿔 계산한다(2026-11-17 실측: n↑ = 1e-6에서
        # Libxc의 n↓ = 1e-16 값이 자체 식의 n↓ = 1e-15 값과 10자리 일치). 그런 점(여기서는 ζ = ±1의 0 밀도)은
        # 느슨한 기준으로 따로 비교한다. 자체 식은 그 점에서도 해석적이며 한쪽 차분과 일치한다.
        ours = lsda_xc(self.up, self.down)
        partial = np.minimum(self.up, self.down) >= 1e-15
        for name, key in (('LDA_X,', 'x'), (',LDA_C_VWN', 'c'), ('LDA_X,LDA_C_VWN', 'xc')):
            with self.subTest(xc=name):
                exc, vxc = libxc.eval_xc(name, (self.up, self.down), spin=1, deriv=1)[:2]
                for mask, rtol in ((partial, 1e-11), (~partial, 1e-9)):
                    np.testing.assert_allclose(ours['eps_' + key][mask], exc[mask], rtol=rtol, atol=1e-16)
                    for spin, density in enumerate((self.up, self.down)):
                        keep = mask & (density > 1e-14)
                        np.testing.assert_allclose(ours['v_' + key][spin][keep], vxc[0][keep, spin], rtol=rtol)


class LogGridTests(unittest.TestCase):
    def test_hydrogenic_eigenvalues_converge_second_order(self):
        z = 18
        errors = []
        for h in (0.004, 0.002, 0.001):
            grid = log_grid(AtomConfig('Ar', grid_step=h))
            for ell in (0, 1, 2):
                eigenvalues, raw, kinetic, phi, nodes = solve_channel(grid, -z / grid.r, ell, 3 - (ell == 2))
                exact = -z ** 2 / (2.0 * np.arange(ell + 1, ell + 1 + eigenvalues.size) ** 2)
                if ell == 0:
                    errors.append(abs(eigenvalues[0] - exact[0]))
                np.testing.assert_allclose(eigenvalues, exact, rtol=5e-5 * (h / 0.001) ** 2)
                self.assertEqual(nodes, list(range(eigenvalues.size)))
                # 비리얼: 순수 쿨롱 고유상태는 T = -E.
                np.testing.assert_allclose(kinetic, -eigenvalues, rtol=1e-6 * (h / 0.001) ** 2 * 50)
                u = np.sqrt(grid.jac)[:, None] * phi
                np.testing.assert_allclose(np.sum(u ** 2 * grid.weights[:, None], axis=0), 1.0, rtol=1e-12)
        self.assertAlmostEqual(errors[0] / errors[1], 4.0, delta=0.1)
        self.assertAlmostEqual(errors[1] / errors[2], 4.0, delta=0.1)

    def test_hartree_matches_hydrogenic_analytic_potential(self):
        zeta = 2.0
        grid = log_grid(AtomConfig('He'))
        u2 = 4 * zeta ** 3 * grid.r ** 2 * np.exp(-2 * zeta * grid.r)
        v_h = hartree_potential(u2, grid)
        exact = 1.0 / grid.r - (zeta + 1.0 / grid.r) * np.exp(-2 * zeta * grid.r)
        self.assertLess(np.max(np.abs(v_h - exact)), 2e-6)
        self.assertAlmostEqual(0.5 * float(np.sum(v_h * u2 * grid.weights)), 5 * zeta / 16, delta=1e-6)


class AtomScfTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.helium = run_atom(AtomConfig('He'))
        cls.nitrogen = run_atom(AtomConfig('N', spin_polarized=True))

    def test_helium_matches_pyscf_and_uniform_grid_engine(self):
        result = self.helium
        self.assertTrue(result.converged)
        self.assertLess(abs(result.total_energy - PYSCF_HE_LDA), 1e-7)
        polarized = run_atom(AtomConfig('He', spin_polarized=True))
        self.assertAlmostEqual(polarized.total_energy, result.total_energy, delta=1e-11)

    def test_nitrogen_lsd_matches_nist(self):
        result, nist = self.nitrogen, NIST['elements']['N']['LSD']
        self.assertTrue(result.converged)
        self.assertLess(abs(result.total_energy - nist['Etot']), 1e-5)
        for orbital in result.orbitals:
            self.assertLess(abs(orbital['eigenvalue_hartree'] - nist_value('N', 'LSD', orbital)), 2e-6, orbital)
        components = result.components
        for term, key in (('Ekin', 'kinetic'), ('Ecoul', 'hartree'), ('Eenuc', 'electron_nuclear'),
                          ('Exc', 'exchange_correlation')):
            self.assertLess(abs(components[key] - nist[term]), 1e-5, term)

    def test_richardson_reaches_nist_rounding(self):
        coarse = run_atom(AtomConfig('N', spin_polarized=True, grid_step=0.002)).total_energy
        extrapolated = self.nitrogen.total_energy + (self.nitrogen.total_energy - coarse) / 3.0
        self.assertLess(abs(extrapolated - NIST['elements']['N']['LSD']['Etot']), 2e-6)

    def test_heavier_closed_shell_lda_matches_nist(self):
        for symbol, tolerance in (('Ne', 1e-5), ('Ar', 1e-4)):
            with self.subTest(element=symbol):
                result = run_atom(AtomConfig(symbol))
                self.assertTrue(result.converged)
                self.assertLess(abs(result.total_energy - NIST['elements'][symbol]['LDA']['Etot']), tolerance)
                for orbital in result.orbitals:
                    self.assertLess(abs(orbital['eigenvalue_hartree'] - nist_value(symbol, 'LDA', orbital)), 1e-5)
                    self.assertEqual(orbital['nodes'], orbital['expected_nodes'])

    def test_electron_count_and_energy_expression(self):
        for result, electrons in ((self.helium, 2), (self.nitrogen, 7)):
            self.assertAlmostEqual(result.electron_count, electrons, delta=1e-10)
            c = result.components
            self.assertAlmostEqual(c['total_direct'], result.total_energy, delta=1e-10)
            self.assertAlmostEqual(c['kinetic'] + c['electron_nuclear'] + c['hartree'] + c['exchange_correlation'],
                                   result.total_energy, delta=1e-10)
        majority = sum(o['occupation'] for o in self.nitrogen.orbitals if o['spin'] == 'majority')
        self.assertEqual(majority, 5)
        self.assertAlmostEqual(float(np.sum(self.nitrogen.spin_radial[0] * self.nitrogen.grid.weights)), 5, delta=1e-10)

    def test_small_iteration_limit_is_not_converged(self):
        result = run_atom(AtomConfig('Ne', max_iterations=3))
        self.assertFalse(result.converged)
        self.assertEqual(len(result.history), 3)
        with tempfile.TemporaryDirectory() as temp:
            summary = save_atom_result(result, temp)
            saved = json.loads((Path(temp) / 'results.json').read_text(encoding='utf-8'))
        self.assertEqual(saved, summary)
        self.assertEqual(saved['status'], 'not_converged')
        self.assertGreater(saved['convergence']['final_density_residual_electrons'], 1e-3)
        self.assertTrue(any('미수렴' in warning for warning in saved['warnings']))

    def test_saved_files(self):
        with tempfile.TemporaryDirectory() as temp:
            summary = save_atom_result(self.nitrogen, temp, elapsed_seconds=1.0)
            directory = Path(temp)
            self.assertEqual(load_atom_config(directory / 'settings.json'), self.nitrogen.config)
            with (directory / 'orbitals.csv').open(encoding='utf-8-sig') as stream:
                orbitals = list(csv.DictReader(stream))
            with (directory / 'radial.csv').open(encoding='utf-8') as stream:
                header = stream.readline().strip().split(',')
        self.assertEqual(len(orbitals), 6)
        self.assertIn('u_2p_majority', header)
        self.assertEqual(summary['system']['spin_moment_electrons'], 3.0)
        self.assertEqual(summary['xc']['correlation']['libxc_id'], 7)


class AtomCommandLineTests(unittest.TestCase):
    def test_atom_cli_converged_and_not_converged(self):
        with tempfile.TemporaryDirectory() as temp:
            done = cli('--model', 'atom', '--element', 'Li', '--spin-polarized', '--no-show', '--output', temp)
            self.assertEqual(done.returncode, 0, done.stderr)
            self.assertIn('총에너지', done.stdout)
            (run,) = Path(temp).iterdir()
            summary = json.loads((run / 'results.json').read_text(encoding='utf-8'))
            self.assertEqual(summary['status'], 'converged')
            self.assertLess(abs(summary['total_energy_hartree'] - NIST['elements']['Li']['LSD']['Etot']), 1e-5)
            self.assertTrue((run / 'atom.png').stat().st_size > 1000)
            stopped = cli('--model', 'atom', '--element', 'Li', '--max-iterations', '2', '--no-show', '--output', temp)
            self.assertEqual(stopped.returncode, 2, stopped.stderr)
            again = cli('--model', 'atom', '--config', str(run / 'settings.json'), '--no-show', '--output', temp)
            self.assertEqual(again.returncode, 0, again.stderr)
            repeated = json.loads((sorted(Path(temp).iterdir())[-1] / 'results.json').read_text(encoding='utf-8'))
            self.assertEqual(repeated['total_energy_hartree'], summary['total_energy_hartree'])

    def test_option_conflicts(self):
        for arguments in (('--model', 'atom'), ('--model', 'atom', '--element', 'He', '--grid-points', '999'),
                          ('--model', 'he-lda', '--element', 'He'), ('--element', 'He'),
                          ('--model', 'atom', '--element', 'He', '--xc', 'pbe')):
            with self.subTest(arguments=arguments):
                done = cli(*arguments, '--no-show')
                self.assertNotEqual(done.returncode, 0)
                self.assertTrue(done.stderr.strip())


if __name__ == '__main__':
    unittest.main()
