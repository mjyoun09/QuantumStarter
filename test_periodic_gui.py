"""The desktop uses the real WSL bridge for a small PBC calculation."""
import os
from pathlib import Path
import tempfile
import time
import tkinter as tk
import unittest
from unittest.mock import patch
from dataclasses import replace
from ase import Atoms
from periodic import PeriodicConfig, load_calculation_config
from periodic_app import PeriodicApp

HIDDEN = os.environ.get('QUANTUMSTARTER_TK_HIDDEN') == '1'


class PeriodicGuiTests(unittest.TestCase):
    def setUp(self):
        if os.name == 'nt':
            import ctypes
            ctypes.windll.shcore.SetProcessDpiAwareness(1)
        try:
            self.root = tk.Tk()
        except tk.TclError:
            self.skipTest('Tk display unavailable')
        if HIDDEN:
            self.root.withdraw()
            create_toplevel = tk.Toplevel.__init__
            def hidden_toplevel(window, *args, **kwargs):
                create_toplevel(window, *args, **kwargs)
                window.withdraw()
            patcher = patch.object(tk.Toplevel, '__init__', hidden_toplevel)
            patcher.start()
            self.addCleanup(patcher.stop)
        self.errors = []
        self.root.report_callback_exception = lambda *args: self.errors.append(args)
        self.app = PeriodicApp(self.root)
        self.root.update()

    def tearDown(self):
        self.app.close()

    def test_settings_roundtrip_and_validation(self):
        original = self.app.read_config()
        self.app.apply_config(original)
        self.assertEqual(self.app.read_config(), original)
        self.app.vars['kmesh'].set('1 0 1')
        with patch('periodic_app.messagebox.showerror') as error:
            self.app.start_calculation()
            error.assert_called_once()
        self.assertFalse(self.app.busy)
        self.assertFalse(self.errors)
        if not HIDDEN:
            self.assertLess(self.app.run_button.winfo_rooty() + self.app.run_button.winfo_height(),
                            self.root.winfo_rooty() + self.root.winfo_height())

    def test_cell_settings_and_unexposed_fields_roundtrip(self):
        config = replace(self.app.read_config(), integral_scheme='fftdf', fft_mesh=(21, 23, 25),
                         calculation='cell_optimize', cell_relaxation='isotropic', cell_max_steps=13,
                         cell_force_tol=.015, cell_stress_tol_gpa=.05, cell_energy_tol_ev=1e-6,
                         band_path='GX', band_points=17, dos_sigma_ev=.23)
        self.app.apply_config(config)
        self.assertEqual(self.app.read_config(), config)
        with patch('periodic_app.messagebox.showerror') as error:
            self.app.mesh_sequence.set('21 21 21;19 19 19;25 25 25')
            self.app.start_mesh_check()
            error.assert_called_once()
        self.assertFalse(self.app.busy)

    def wait_for_job(self, seconds=180):
        deadline = time.monotonic() + seconds
        while self.app.busy and time.monotonic() < deadline:
            self.root.update()
            time.sleep(.03)
        self.assertFalse(self.app.busy)
        self.assertFalse(self.errors, self.errors)

    @unittest.skipUnless(os.environ.get('QS_PERIODIC_GUI_DFT') == '1', 'Enable real WSL GUI test explicitly')
    def test_real_cell_step_limit_result(self):
        cfg = replace(load_calculation_config('examples/silicon_cell.json'), kmesh=(1, 1, 1),
                      fft_mesh=(17, 17, 17), cell_max_steps=1)
        self.app.apply_config(cfg)
        with patch('periodic_app.messagebox.showerror') as error:
            self.app.start_calculation()
            self.wait_for_job()
            error.assert_not_called()
        self.assertEqual(self.app.summary['status'], 'not_converged')
        self.assertIn('최대 원자 힘', self.app.result.get('1.0', 'end'))
        self.assertEqual(self.app.summary['cell_optimization']['steps'], 1)
        self.assertTrue((self.app.current_dir / 'final_settings.json').is_file())
        self.assertTrue((self.app.current_dir / 'periodic.png').is_file())

    @unittest.skipUnless(os.environ.get('QS_PERIODIC_GUI_DFT') == '1', 'Enable real WSL GUI test explicitly')
    def test_real_periodic_cancel_before_worker(self):
        self.app.apply_config(load_calculation_config('examples/silicon_cell.json'))
        with patch('periodic_app.messagebox.showerror') as error:
            self.app.start_calculation()
            self.app.cancel()
            self.wait_for_job()
            error.assert_not_called()
        self.assertEqual(self.app.summary['status'], 'cancelled')
        self.assertFalse(self.app.summary['converged'])
        self.assertTrue(self.app.cancel_button.instate(['disabled']))

    def test_mesh_report_display(self):
        self.app.apply_config(load_calculation_config('examples/silicon_cell.json'))
        report = {'status': 'not_converged', 'objective': 'free_energy_ev', 'selected_mesh': None,
                  'rows': [{'mesh': [21, 21, 21], 'passed': False}]}
        with tempfile.TemporaryDirectory() as temp, patch('periodic_convergence.run_mesh_convergence',
                                                          return_value=(report, Path(temp))):
            self.app.start_mesh_check()
            self.wait_for_job()
        self.assertIn('FFT mesh 수렴 확인: not_converged', self.app.result.get('1.0', 'end'))
        self.assertIsNone(self.app.control)

    @unittest.skipUnless(os.environ.get('QS_PERIODIC_GUI_DFT') == '1', 'Enable real WSL GUI test explicitly')
    def test_real_periodic_worker_result_and_plot(self):
        config = PeriodicConfig.from_ase(Atoms('H', positions=[[1, 1, 1]], cell=[4, 4, 4], pbc=True),
                     spin=1, kmesh=(2, 1, 1), basis='gth-dzvp', grid_level=0, precision=1e-5, conv_tol=1e-6, threads=2)
        self.app.apply_config(config)
        self.app.start_calculation()
        deadline = time.monotonic() + 120
        while self.app.busy and time.monotonic() < deadline:
            self.root.update()
            time.sleep(.03)
        self.assertFalse(self.app.busy)
        self.assertFalse(self.errors, self.errors)
        self.assertIn('converged', self.app.result.get('1.0', 'end'))
        self.assertTrue(self.app.canvas.figure.axes)


if __name__ == '__main__':
    unittest.main()
