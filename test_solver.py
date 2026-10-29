"""물리적으로 알려진 해와 수치 수렴 검증. 실행: python -m unittest -v."""
import unittest

import numpy as np

from solver import Config, solve


class HydrogenValidation(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.result = solve(Config())

    def test_hydrogen_1s_2s_3s_energies(self):
        # Hartree 원자 단위, 무한 핵질량 근사에서 알려진 에너지.
        np.testing.assert_allclose(
            self.result.energies_hartree, [-0.5, -0.125, -1.0 / 18.0],
            atol=2e-5, rtol=0,
        )

    def test_normalization_and_orthogonality(self):
        # 같은 l에서 서로 다른 상태는 직교하고 각 상태의 총 확률은 1.
        overlap = self.result.dr * (self.result.u.T @ self.result.u)
        np.testing.assert_allclose(overlap, np.eye(3), atol=1e-11, rtol=0)

    def test_1s_wavefunction_and_mean_radius(self):
        result = self.result
        overlap = result.dr * np.dot(result.u[:, 0], result.u_exact[:, 0])
        self.assertGreater(overlap, 0.999999)
        # 수소 1s의 <r> = 3*a0/2.
        self.assertAlmostEqual(result.mean_r_bohr[0], 1.5, delta=1e-4)

    def test_grid_refinement_reduces_ground_state_error(self):
        # 계산영역을 고정하고 간격을 반으로 줄이면 2차 차분 오차가 약 1/4.
        errors = []
        for count in (399, 799, 1599):
            result = solve(Config(num_states=1, r_max_bohr=40, grid_points=count))
            errors.append(abs(float(result.energies_hartree[0]) + 0.5))
        self.assertGreater(errors[0] / errors[1], 3.8)
        self.assertGreater(errors[1] / errors[2], 3.8)
        self.assertLess(errors[0] / errors[1], 4.2)

    def test_l_and_nuclear_charge(self):
        # He+ 한 전자의 2p. 전자 두 개인 중성 He 계산이 아니다.
        result = solve(Config(nuclear_charge=2, angular_momentum_l=1,
                              num_states=1, r_max_bohr=30, grid_points=6000))
        self.assertEqual(int(result.principal_n[0]), 2)
        self.assertAlmostEqual(float(result.energies_hartree[0]), -0.5, delta=1e-5)

    def test_box_size_error(self):
        # 격자 간격이 비슷해도 영역이 너무 작으면 들뜬 상태가 틀어진다.
        small = solve(Config(num_states=2, r_max_bohr=8, grid_points=799))
        large = solve(Config(num_states=2, r_max_bohr=40, grid_points=3999))
        small_error = abs(float(small.energies_hartree[1]) + 0.125)
        large_error = abs(float(large.energies_hartree[1]) + 0.125)
        self.assertGreater(small_error, 100 * large_error)

    def test_invalid_config(self):
        for value in (True, 1.5, 0, -1, "1"):
            with self.subTest(nuclear_charge=value):
                with self.assertRaises(ValueError):
                    Config(nuclear_charge=value)
        with self.assertRaises(ValueError):
            Config(r_max_bohr=float("nan"))


if __name__ == "__main__":
    unittest.main()
