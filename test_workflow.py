"""사용 화면과 CLI가 공유하는 재현성·저장·단위 변환 검증."""
from dataclasses import asdict
import json
from pathlib import Path
import tempfile
import unittest

import numpy as np

from plotting import create_figure
from solver import BOHR_TO_ANGSTROM, Config, solve
from workflow import (config_from_fields, list_runs, load_config, new_run_dir,
                      recommended_config, save_config, save_result)


class WorkflowValidation(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.result = solve(Config(num_states=2))

    def test_settings_roundtrip(self):
        with tempfile.TemporaryDirectory() as temp:
            path = Path(temp) / '설정.json'
            cfg = Config(2, 1, 2, 50, 8000)
            save_config(cfg, path)
            self.assertEqual(load_config(path), cfg)

    def test_presets_follow_atomic_size(self):
        hydrogen = recommended_config(1, 0, 3)
        helium_ion = recommended_config(2, 0, 3)
        self.assertEqual(hydrogen, Config())
        self.assertEqual(hydrogen.r_max_bohr, 2 * helium_ion.r_max_bohr)
        self.assertGreater(recommended_config(1, 1, 3).r_max_bohr, hydrogen.r_max_bohr)
        self.assertGreater(recommended_config(precision='정밀').grid_points, hydrogen.grid_points)

    def test_invalid_fields_are_explained_in_korean(self):
        fields = {k: str(v) for k, v in asdict(Config()).items()}
        fields['grid_points'] = 'abc'
        with self.assertRaisesRegex(ValueError, '격자점 수'):
            config_from_fields(fields)
        fields['grid_points'] = '6000'
        fields['r_max_bohr'] = 'nan'
        with self.assertRaisesRegex(ValueError, '계산영역 반지름'):
            config_from_fields(fields)

    def test_unique_runs_and_raw_data_reproducibility(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            first, second = new_run_dir(root), new_run_dir(root)
            self.assertNotEqual(first, second)
            summary = save_result(self.result, first)
            self.assertEqual(load_config(first/'settings.json'), self.result.config)
            loaded = json.loads((first/'results.json').read_text(encoding='utf-8'))
            self.assertEqual(summary, loaded)
            raw = np.loadtxt(first/'radial.csv', skiprows=1, delimiter=',')
            np.testing.assert_allclose(raw[:, 1], self.result.u[:, 0], rtol=0, atol=1e-15)
            self.assertEqual(len(list_runs(root)), 1)
            (second/'results.json').write_text('{broken', encoding='utf-8')
            self.assertEqual(len(list_runs(root)), 1)

    def test_angstrom_plot_preserves_probability(self):
        figure = create_figure(self.result, unit='angstrom')
        probability_line = figure.axes[1].lines[0]
        distance, density = probability_line.get_data()
        self.assertAlmostEqual(float(distance[1]), self.result.dr * BOHR_TO_ANGSTROM)
        norm = np.sum(density) * (distance[1]-distance[0])
        self.assertAlmostEqual(float(norm), 1.0, places=10)


if __name__ == '__main__':
    unittest.main()
