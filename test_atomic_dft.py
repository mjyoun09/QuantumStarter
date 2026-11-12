"""중성 He 방사 격자 Kohn–Sham LDA-SCF 검증. 실행: python -m unittest test_atomic_dft -v

Windows .venv와 WSL .venv-dft에서 모두 실행한다. Libxc 대조는 PySCF에 들어 있는 Libxc가 있을 때만
실행한다(WSL). PySCF는 여기서도 비교 검증에만 쓰고 계산 엔진에는 쓰지 않는다.
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

from atomic_dft import (ELECTRONS, XC_FUNCTIONALS, AtomicDFTConfig, density_from_orbital, electron_count,
                        hartree_potential, initial_orbital, lda_xc, load_atomic_config, run_scf,
                        save_atomic_config, save_atomic_result, slater_exchange, vwn5_correlation)
from solver import radial_grid

ROOT = Path(__file__).resolve().parent
# 독립 기준값: scripts/validate_he_lda.py가 PySCF 2.14 RKS(LDA_X,LDA_C_VWN, 짝수 템퍼드 s 기저
# 39개 0.02~1e8 bohr^-2 β=1.8, 방사 300점, small_rho_cutoff=0)로 구한 총에너지. 기저·격자 변화 ≤3e-9 Ha.
# 회귀 검사용 고정값이며, 재현: wsl --exec .venv-dft/bin/python scripts/validate_he_lda.py
PYSCF_RKS_TOTAL_ENERGY = -2.834835624033
TARGET_ERROR_HARTREE = 1e-4

try:
    from pyscf.dft import libxc
except ImportError:
    libxc = None


def cli(*arguments):
    return subprocess.run([sys.executable, str(ROOT / 'main.py'), *arguments], cwd=ROOT, capture_output=True,
                          encoding='utf-8', errors='replace', env=dict(os.environ, PYTHONIOENCODING='utf-8'),
                          timeout=300)


class LdaFormulaValidation(unittest.TestCase):
    densities = np.logspace(-10, 4, 141)

    def test_slater_exchange_closed_form(self):
        eps, v = slater_exchange(self.densities)
        np.testing.assert_allclose(eps, -0.75 * (3.0 * self.densities / np.pi) ** (1.0 / 3.0), rtol=1e-14)
        np.testing.assert_allclose(v, -(3.0 * self.densities / np.pi) ** (1.0 / 3.0), rtol=1e-14)

    def test_potentials_are_density_derivatives(self):
        # v = d(n eps)/dn, 상대 간격 1e-5 중앙차분.
        for function in (slater_exchange, vwn5_correlation):
            with self.subTest(function=function.__name__):
                n = self.densities
                h = 1e-5 * n
                numeric = ((n + h) * function(n + h)[0] - (n - h) * function(n - h)[0]) / (2 * h)
                np.testing.assert_allclose(function(n)[1], numeric, rtol=1e-8)

    def test_zero_density_is_finite(self):
        with np.errstate(all='raise'):
            terms = lda_xc(np.array([0.0, 1e-40, 1.0]))
        self.assertEqual(float(terms['eps_xc'][0]), 0.0)
        self.assertEqual(float(terms['v_xc'][1]), 0.0)
        self.assertLess(float(terms['eps_xc'][2]), 0.0)

    @unittest.skipIf(libxc is None, 'Libxc 대조에는 PySCF(번들 Libxc)가 필요합니다(WSL .venv-dft).')
    def test_matches_libxc_lda_x_and_lda_c_vwn(self):
        _, ((x_id, _), (c_id, _)) = libxc.parse_xc(XC_FUNCTIONALS['pyscf_xc_string'])
        self.assertEqual((int(x_id), int(c_id)), (XC_FUNCTIONALS['exchange']['libxc_id'],
                                                  XC_FUNCTIONALS['correlation']['libxc_id']))
        ours = lda_xc(self.densities)
        for name, key in (('LDA_X,', 'x'), (',LDA_C_VWN', 'c'), ('LDA_X,LDA_C_VWN', 'xc')):
            with self.subTest(xc=name):
                exc, vxc = libxc.eval_xc(name, self.densities, spin=0, deriv=1)[:2]
                np.testing.assert_allclose(ours['eps_' + key], exc, rtol=1e-12, atol=0)
                np.testing.assert_allclose(ours['v_' + key], vxc[0], rtol=1e-12, atol=0)


class HartreeValidation(unittest.TestCase):
    zeta = 27.0 / 16.0

    def hydrogenic(self, grid_points, r_max=20.0):
        dr, r = radial_grid(r_max, grid_points)
        u, _ = initial_orbital(r, dr, self.zeta)
        radial, _ = density_from_orbital(u, r, occupation=1.0)
        return dr, r, radial, hartree_potential(radial, r, dr)

    def test_hydrogenic_density_matches_analytic_potential(self):
        # 한 전자 1s 밀도: v_H = 1/r - (zeta + 1/r) exp(-2 zeta r), E_H = 5 zeta / 16.
        errors = []
        for points in (4999, 9999, 19999):
            dr, r, radial, v_h = self.hydrogenic(points)
            exact = 1.0 / r - (self.zeta + 1.0 / r) * np.exp(-2.0 * self.zeta * r)
            errors.append(np.max(np.abs(v_h - exact)))
            self.assertAlmostEqual(0.5 * dr * float(np.dot(v_h, radial)), 5.0 * self.zeta / 16.0, delta=2e-6)
        self.assertLess(errors[-1], 2e-6)
        self.assertGreater(errors[0] / errors[1], 3.8)
        self.assertGreater(errors[1] / errors[2], 3.8)

    def test_inner_and_outer_charge_contributions(self):
        dr, r, radial, v_h = self.hydrogenic(19999)
        # 모든 전하가 안쪽인 바깥 영역은 Q/r, 원점 근처는 바깥 구껍질 ∫ rho/r dr = zeta가 지배한다.
        self.assertAlmostEqual(float(v_h[-1] * r[-1]), 1.0, places=12)
        self.assertAlmostEqual(float(v_h[0]), self.zeta, delta=1e-5)


class HeliumScfValidation(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.result = run_scf(AtomicDFTConfig())

    def test_default_run_converges_with_both_criteria(self):
        result, last = self.result, self.result.history[-1]
        self.assertTrue(result.converged)
        self.assertEqual(result.status, 'converged')
        self.assertLess(len(result.history), result.config.max_iterations)
        self.assertLess(abs(last['delta_energy_hartree']), result.config.energy_tol_hartree)
        self.assertLess(last['density_residual_electrons'], result.config.density_tol_electrons)
        self.assertIsNone(result.history[0]['delta_energy_hartree'])
        # 출력 오비탈의 Kohn–Sham 에너지는 변분적이라 수렴 과정에서 위에서부터 내려온다.
        energies = [h['total_energy_hartree'] for h in result.history]
        self.assertTrue(all(b <= a + 1e-12 for a, b in zip(energies, energies[1:])))

    def test_electron_count_is_two(self):
        result = self.result
        self.assertAlmostEqual(result.electron_count, ELECTRONS, delta=1e-10)
        self.assertAlmostEqual(result.initial_guess['electron_count'], ELECTRONS, delta=1e-10)
        for record in result.history:
            self.assertAlmostEqual(record['electrons_in'], ELECTRONS, delta=1e-10)
        # 전자 수 적분이 4 pi r^2 dr 부피 요소를 쓰는지: 밀도 n에서 직접 다시 적분한다.
        self.assertAlmostEqual(electron_count(result.density, result.r_bohr, result.dr), 2.0, delta=1e-10)
        self.assertAlmostEqual(float(result.dr * np.dot(result.u, result.u)), 1.0, delta=1e-12)

    def test_double_counting_total_equals_direct_functional(self):
        c, total = self.result.components, self.result.total_energy
        direct = c['kinetic'] + c['electron_nuclear'] + c['hartree'] + c['exchange_correlation']
        self.assertAlmostEqual(direct, total, delta=1e-10)
        self.assertAlmostEqual(c['eigenvalue_sum'] + c['double_counting_hartree'] + c['double_counting_xc'],
                               total, delta=1e-12)
        self.assertAlmostEqual(c['eigenvalue_sum'], 2.0 * self.result.eigenvalue_1s, delta=1e-14)
        # 수렴 밀도에서 하트리 이중 계산 보정은 -E_H이다.
        self.assertAlmostEqual(c['double_counting_hartree'], -c['hartree'], delta=1e-6)

    def test_total_energy_matches_independent_pyscf_reference(self):
        self.assertLess(abs(self.result.total_energy - PYSCF_RKS_TOTAL_ENERGY), TARGET_ERROR_HARTREE)
        self.assertTrue(-0.58 < self.result.eigenvalue_1s < -0.56)

    def test_grid_refinement_is_second_order_and_extrapolates_to_reference(self):
        energies = [run_scf(AtomicDFTConfig(r_max_bohr=20.0, grid_points=n)).total_energy
                    for n in (4999, 9999, 19999)]
        first, second = energies[1] - energies[0], energies[2] - energies[1]
        self.assertLess(second, 0)
        self.assertAlmostEqual(first / second, 4.0, delta=0.05)
        richardson = energies[2] + second / 3.0
        self.assertLess(abs(richardson - PYSCF_RKS_TOTAL_ENERGY), 1e-8)

    def test_domain_enlargement_converges(self):
        energies = {r_max: run_scf(AtomicDFTConfig(r_max_bohr=r_max, grid_points=int(r_max / 0.002) - 1)).total_energy
                    for r_max in (5.0, 10.0, 20.0)}
        self.assertGreater(energies[5.0] - energies[20.0], 1e-4)
        self.assertLess(abs(energies[10.0] - energies[20.0]), 1e-7)

    def test_small_iteration_limit_is_reported_as_not_converged(self):
        result = run_scf(AtomicDFTConfig(max_iterations=3))
        self.assertFalse(result.converged)
        self.assertEqual(result.status, 'not_converged')
        self.assertEqual(len(result.history), 3)
        with tempfile.TemporaryDirectory() as temp:
            summary = save_atomic_result(result, temp)
            saved = json.loads((Path(temp) / 'results.json').read_text(encoding='utf-8'))
        self.assertEqual(saved, summary)
        self.assertEqual(saved['status'], 'not_converged')
        self.assertFalse(saved['converged'])
        self.assertGreater(saved['convergence']['final_density_residual_electrons'], 1e-3)
        self.assertEqual(saved['convergence']['final_density_residual_electrons'],
                         result.history[-1]['density_residual_electrons'])
        self.assertFalse(saved['convergence']['density_criterion_met'])
        self.assertTrue(any('미수렴' in warning for warning in saved['warnings']))

    def test_saved_files(self):
        with tempfile.TemporaryDirectory() as temp:
            summary = save_atomic_result(self.result, temp, elapsed_seconds=1.0)
            directory = Path(temp)
            self.assertEqual(load_atomic_config(directory / 'settings.json'), self.result.config)
            with (directory / 'scf_history.csv').open(encoding='utf-8-sig') as stream:
                rows = list(csv.DictReader(stream))
            radial = np.loadtxt(directory / 'radial.csv', delimiter=',', skiprows=1)
        self.assertEqual(len(rows), summary['iterations'])
        self.assertEqual(float(rows[-1]['total_energy_hartree']), self.result.total_energy)
        self.assertEqual(radial.shape, (self.result.config.grid_points, 8))
        np.testing.assert_array_equal(radial[:, 1], self.result.u)
        for key in ('total_energy_hartree', 'total_energy_ev', 'eigenvalue_1s_hartree', 'electron_count',
                    'scf_history', 'xc', 'grid', 'status'):
            self.assertIn(key, summary)
        self.assertEqual(summary['grid']['dr_bohr'], self.result.dr)
        self.assertEqual(summary['xc']['correlation']['libxc_name'], 'LDA_C_VWN')

    def test_invalid_config(self):
        for changes in ({'mixing': 0}, {'mixing': 1.5}, {'mixing': float('nan')}, {'grid_points': 1.5},
                        {'grid_points': True}, {'grid_points': 50}, {'max_iterations': 0},
                        {'energy_tol_hartree': -1e-9}, {'r_max_bohr': 0}):
            with self.subTest(**changes), self.assertRaises(ValueError):
                AtomicDFTConfig(**changes)
        with tempfile.TemporaryDirectory() as temp:
            path = Path(temp) / 'settings.json'
            save_atomic_config(AtomicDFTConfig(mixing=0.7), path)
            self.assertEqual(load_atomic_config(path).mixing, 0.7)
            path.write_text(json.dumps({'model': 'hydrogen'}), encoding='utf-8')
            with self.assertRaises(ValueError):
                load_atomic_config(path)
            path.write_text(json.dumps({'unknown': 1}), encoding='utf-8')
            with self.assertRaises(ValueError):
                load_atomic_config(path)


class CommandLineValidation(unittest.TestCase):
    def test_he_lda_cli_converged_and_not_converged(self):
        with tempfile.TemporaryDirectory() as temp:
            done = cli('--model', 'he-lda', '--no-show', '--output', temp, '--r-max', '20', '--grid-points', '9999')
            self.assertEqual(done.returncode, 0, done.stderr)
            self.assertIn('총에너지', done.stdout)
            self.assertIn('Ha', done.stdout)
            (run,) = Path(temp).iterdir()
            summary = json.loads((run / 'results.json').read_text(encoding='utf-8'))
            self.assertEqual(summary['status'], 'converged')
            self.assertEqual(summary['config']['grid_points'], 9999)
            self.assertTrue((run / 'he_lda.png').stat().st_size > 1000)
            stopped = cli('--model', 'he-lda', '--no-show', '--output', temp, '--max-iterations', '2')
            self.assertEqual(stopped.returncode, 2, stopped.stderr)
            self.assertIn('미수렴', stopped.stdout)
            again = cli('--model', 'he-lda', '--no-show', '--output', temp, '--config', str(run / 'settings.json'))
            self.assertEqual(again.returncode, 0, again.stderr)
            runs = sorted(Path(temp).iterdir())
            self.assertEqual(len(runs), 3)
            repeated = json.loads((runs[-1] / 'results.json').read_text(encoding='utf-8'))
            self.assertEqual(repeated['total_energy_hartree'], summary['total_energy_hartree'])

    def test_he_lda_rejects_molecular_options(self):
        for arguments in (('--xc', 'pbe'), ('--hydrogen-check',), ('--charge', '1')):
            with self.subTest(arguments=arguments):
                done = cli('--model', 'he-lda', '--no-show', *arguments)
                self.assertNotEqual(done.returncode, 0)
                self.assertIn('he-lda', done.stderr)
        done = cli('--mixing', '0.3', '--no-show')
        self.assertNotEqual(done.returncode, 0)
        self.assertIn('he-lda', done.stderr)

    def test_hydrogen_check_cli_is_unchanged(self):
        with tempfile.TemporaryDirectory() as temp:
            done = cli('--hydrogen-check', '--no-show', '--output', temp)
            self.assertEqual(done.returncode, 0, done.stderr)
            summary = json.loads((Path(temp) / 'results.json').read_text(encoding='utf-8'))
        self.assertEqual([s['label'] for s in summary['states']], ['1s', '2s', '3s'])
        np.testing.assert_allclose([s['energy_hartree'] for s in summary['states']],
                                   [-0.5, -0.125, -1.0 / 18.0], atol=2e-5, rtol=0)


if __name__ == '__main__':
    unittest.main()
