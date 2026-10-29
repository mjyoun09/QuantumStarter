"""The desktop uses the real WSL bridge for a small PBC calculation."""
import os
from pathlib import Path
import tempfile
import time
import tkinter as tk
import unittest
from unittest.mock import patch
from ase import Atoms
from periodic import PeriodicConfig
from periodic_app import PeriodicApp


class PeriodicGuiTests(unittest.TestCase):
    def setUp(self):
        if os.name == 'nt':
            import ctypes
            ctypes.windll.shcore.SetProcessDpiAwareness(1)
        self.root = tk.Tk()
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
        self.assertLess(self.app.run_button.winfo_rooty() + self.app.run_button.winfo_height(),
                        self.root.winfo_rooty() + self.root.winfo_height())

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
