"""Actual Tk workflow, including Windows -> WSL computation when enabled."""
from dataclasses import replace
import os
from pathlib import Path
import tempfile
import time
import tkinter as tk
import unittest
from unittest.mock import patch
import zipfile

from molecular import Atom, DFTConfig, load_dft_config
from molecular_app import MolecularApp


class MolecularGuiTests(unittest.TestCase):
    def setUp(self):
        try:
            self.root = tk.Tk()
        except tk.TclError:
            self.skipTest('Tk display unavailable')
        self.temp = tempfile.TemporaryDirectory()
        self.base = Path(self.temp.name)
        self.errors = []
        self.root.report_callback_exception = lambda *args: self.errors.append(args)
        self.app = MolecularApp(self.root, self.base / 'runs')
        self.root.update()

    def tearDown(self):
        if hasattr(self, 'app'):
            if not self.app.closed:
                self.app.close()
            self.temp.cleanup()

    def wait_finished(self):
        deadline = time.monotonic() + 120
        while self.app.busy and time.monotonic() < deadline:
            self.root.update()
            time.sleep(.02)
        self.assertFalse(self.app.busy)
        self.assertFalse(self.errors, self.errors)
        self.root.update()

    def test_input_validation_xyz_and_settings(self):
        app = self.app
        # Essential structure controls must be visible without scrolling at the default size.
        for key in ('charge', 'spin'):
            control = next(w for w, _ in app.controls if isinstance(w, tk.ttk.Entry)
                           and str(w.cget('textvariable')) == str(app.vars[key]))
            self.assertLess(control.winfo_rooty() + control.winfo_height(), app.run_button.winfo_rooty())
        app.vars['spin'].set('1')
        with patch('molecular_app.messagebox.showerror') as error:
            app.start_calculation()
            error.assert_called_once()
        self.assertFalse(app.busy)
        app.apply_config(DFTConfig())
        xyz = self.base / 'test.xyz'
        xyz.write_text('2\n\nH 0 0 0\nH 0 0 .74\n', encoding='utf-8')
        with patch('molecular_app.filedialog.askopenfilename', return_value=str(xyz)):
            app.load_xyz()
        self.assertEqual(app.read_config().formula, 'H2')
        settings = self.base / 'settings.json'
        with patch('molecular_app.filedialog.asksaveasfilename', return_value=str(settings)):
            app.save_settings()
        self.assertEqual(load_dft_config(settings), app.read_config())

    def test_per_element_basis_entry(self):
        app = self.app
        app.apply_config(DFTConfig(basis='sto-3g', basis_overrides={'O': 'def2-svp'}))
        self.assertEqual(app.vars['basis_overrides'].get(), 'O=def2-svp')
        app.vars['basis_overrides'].set('o=def2-tzvp; H=sto-3g')
        cfg = app.read_config()
        self.assertEqual(cfg.basis_overrides, {'O': 'def2-tzvp'})
        self.assertEqual(cfg.basis_by_element, {'H': 'sto-3g', 'O': 'def2-tzvp'})
        app.vars['basis_overrides'].set('')
        self.assertEqual(app.read_config().basis_overrides, {})
        app.vars['basis_overrides'].set('O=not-a-basis')
        with patch('molecular_app.messagebox.showerror') as error:
            app.start_calculation()
            error.assert_called_once()
        self.assertFalse(app.busy)
        self.assertFalse(self.errors, self.errors)

    def test_period4_examples_and_incompatible_basis_blocked(self):
        app = self.app
        app.example.set('코발트 이온 · Co2+ (2S=3 예제)')
        app.use_example()
        self.root.update()
        cfg = app.read_config()
        self.assertEqual((cfg.formula, cfg.electrons, cfg.spin), ('Co', 25, 3))
        self.assertEqual(cfg.scf_solver, 'newton')
        self.assertNotIn('6-31+g*', app.basis_box['values'])
        app.vars['basis'].set('6-31+g*')
        with patch('molecular_app.messagebox.showerror') as error:
            app.start_calculation()
            error.assert_called_once()
        self.assertFalse(app.busy)
        self.assertEqual(app.vars['basis'].get(), '6-31+g*')
        app.example.set('염화칼륨 · KCl (예제)')
        app.use_example()
        self.assertNotIn('cc-pvdz', app.basis_box['values'])
        self.assertEqual(app.read_config().electrons, 36)
        app.example.set('물 · H2O (예제)')
        app.use_example()
        self.assertIn('6-31+g*', app.basis_box['values'])

    def test_ecp_preview_and_pairing_validation(self):
        app = self.app
        app.example.set('지르코늄 이온 · Zr4+ (ECP 예제)')
        app.use_example()
        cfg = app.read_config()
        self.assertEqual((cfg.electrons, cfg.explicit_electrons, cfg.ecp_core_electrons), (36, 8, 28))
        self.assertIn('전체 전자 36개', app.preview.get())
        self.assertIn('계산 전자 8개', app.preview.get())
        self.assertIn('def2-ECP: Zr', app.ecp_hint.get())
        self.assertEqual(tuple(app.basis_box['values']), ('sto-3g', 'def2-svp', 'def2-tzvp', 'dyall-v2z', 'dyall-v3z'))
        app.vars['ecp'].set('none')
        with patch('molecular_app.messagebox.showerror') as error:
            app.start_calculation()
            error.assert_called_once()
        self.assertFalse(app.busy)
        self.assertEqual(app.vars['ecp'].get(), 'none')
        app.vars['basis'].set('sto-3g')
        self.assertEqual(app.read_config().explicit_electrons, 36)
        self.assertIn('전체전자', app.ecp_hint.get())
        app.example.set('물 · H2O (예제)')
        app.use_example()
        self.assertNotIn('코어 28', app.preview.get())

    def test_all_element_example_support_table_and_save(self):
        app = self.app
        app.example.set('오가네손 · Og (118번 · X2C)')
        app.use_example()
        cfg = app.read_config()
        self.assertEqual(cfg.electrons, 118)
        self.assertEqual(cfg.resolved_relativity, 'x2c')
        self.assertIn('x2c', app.ecp_hint.get())
        self.assertEqual(tuple(app.basis_box['values']), ('dyall-v2z', 'dyall-v3z'))
        window, tree = app.show_elements()
        self.root.update()
        self.assertEqual(len(tree.get_children()), 118)
        self.assertEqual(tree.item('Og', 'values')[-1], '미지원')
        self.assertEqual(tree.item('Au', 'values')[3], '60')
        window.destroy()
        settings = self.base / 'og.json'
        with patch('molecular_app.filedialog.asksaveasfilename', return_value=str(settings)):
            app.save_settings()
        self.assertEqual(load_dft_config(settings), cfg)
        self.assertFalse(self.errors, self.errors)

    def test_native_selection_preset_and_geometry(self):
        app = self.app
        app.example.set('자체 DFT · 물 (LDA 검증)')
        app.use_example()
        self.assertEqual(app.read_config().engine, 'native')
        self.assertEqual(app.read_config().xc, 'lda,vwn')
        app.use_preset('빠른 확인')
        self.assertEqual(app.read_config().grid_level, 1)
        self.assertIn('conv_tol=1e-06', app.status.get())
        window = app.preview_geometry()
        self.root.update()
        self.assertTrue(window.winfo_exists())
        window.destroy()
        analysis_window = app.analyze_adsorption()
        self.root.update()
        button = next(w for w in analysis_window.winfo_children() if isinstance(w, tk.ttk.Button))
        self.assertTrue(button.winfo_ismapped())
        self.assertLessEqual(button.winfo_y() + button.winfo_height(), analysis_window.winfo_height())
        analysis_window.destroy()
        self.assertFalse(self.errors, self.errors)

    @unittest.skipUnless(os.environ.get('QUANTUMSTARTER_GUI_DFT_TEST') == '1', 'Enable real GUI DFT integration explicitly')
    def test_real_native_cancel_and_restart(self):
        app = self.app
        cfg = DFTConfig(engine='native', basis='def2-tzvp', xc='lda,vwn', cube_points=12, grid_level=3, threads=1,
                        native_damping=.6)
        app.apply_config(cfg)
        app.start_calculation()
        def cancel_after_checkpoint():
            folder = app.run_control.directory
            if folder is not None and (folder / 'checkpoint.npz').exists():
                app.cancel_button.invoke()
            elif app.busy:
                self.root.after(30, cancel_after_checkpoint)
        self.root.after(30, cancel_after_checkpoint)
        self.wait_finished()
        self.assertEqual(app.history[0]['summary']['status'], 'cancelled')
        self.assertTrue((app.history[0]['directory'] / 'checkpoint.npz').exists())
        app.history_tree.selection_set('0')
        app.restart_history()
        self.wait_finished()
        self.assertTrue(app.summary['converged'])
        self.assertEqual(app.summary['engine'], 'quantumstarter-native')
        self.assertTrue(app.summary['provenance']['restarted'])
        self.assertIn('자체', app.metrics[2].get())
        self.assertTrue(app.scf_rows)
        self.assertTrue((app.current_dir / 'scf_convergence.png').exists())

    @unittest.skipUnless(os.environ.get('QUANTUMSTARTER_GUI_DFT_TEST') == '1', 'Enable real GUI DFT integration explicitly')
    def test_real_x2c_run_and_metadata(self):
        app = self.app
        cfg = DFTConfig(atoms=(Atom('H', (0, 0, 0)),), spin=1, basis='dyall-v2z',
                        cube_points=16, grid_level=2, threads=1, scf_solver='newton')
        app.apply_config(cfg)
        with patch('molecular_app.messagebox.showerror') as error:
            app.start_calculation()
            self.wait_finished()
            error.assert_not_called()
        self.assertTrue(app.summary['converged'])
        self.assertEqual(app.summary['relativity'], 'x2c')
        self.assertEqual(app.summary['nuclear_model'], 'gaussian')
        self.assertEqual(app.summary['integration_grid']['scheme'], 'universal')
        self.assertIn('X2C', app.figure.axes[1].get_title())
        self.assertAlmostEqual(app.summary['electron_count_from_density'], 1, places=8)
        self.assertEqual(load_dft_config(app.current_dir / 'settings.json'), cfg)

    @unittest.skipUnless(os.environ.get('QUANTUMSTARTER_GUI_DFT_TEST') == '1', 'Enable real GUI DFT integration explicitly')
    def test_real_native_gradient_run(self):
        app = self.app
        cfg = DFTConfig(basis='sto-3g', xc='pbe', engine='native', calculation='gradient',
                        cube_points=12, grid_level=2, threads=1)
        app.apply_config(cfg)
        self.assertEqual(app.vars['calculation'].get(), 'gradient')
        with patch('molecular_app.messagebox.showerror') as error:
            app.start_calculation()
            self.wait_finished()
            error.assert_not_called()
        gradient = app.summary['gradient']
        self.assertEqual(gradient['driver'], 'QuantumStarter')
        self.assertEqual(len(gradient['forces_ev_per_angstrom']), 3)
        self.assertIn('최대 힘', app.metrics[1].get())
        self.assertEqual(load_dft_config(app.current_dir / 'settings.json'), cfg)

    @unittest.skipUnless(os.environ.get('QUANTUMSTARTER_GUI_DFT_TEST') == '1', 'Enable real GUI DFT integration explicitly')
    def test_real_native_optimization_run(self):
        app = self.app
        cfg = DFTConfig(atoms=(Atom('H', (0, 0, 0)), Atom('H', (0, 0, .9))), basis='sto-3g', xc='pbe',
                        engine='native', calculation='optimize', cube_points=12, grid_level=3, threads=1)
        app.apply_config(cfg)
        self.assertEqual(app.read_config(), cfg)
        with patch('molecular_app.messagebox.showerror') as error:
            app.start_calculation()
            self.wait_finished()
            error.assert_not_called()
        self.assertTrue(app.summary['optimization']['converged'])
        self.assertIn('최적화 수렴', app.metrics[1].get())
        self.assertNotIn('OPT_JSON', '\n'.join(app.log_lines))
        self.assertTrue((app.current_dir / 'optimized.xyz').is_file())

    @unittest.skipUnless(os.environ.get('QUANTUMSTARTER_GUI_DFT_TEST') == '1', 'Enable real GUI DFT integration explicitly')
    def test_real_ecp_run_and_serialized_metadata(self):
        app = self.app
        app.example.set('지르코늄 이온 · Zr4+ (ECP 예제)')
        app.use_example()
        cfg = replace(app.read_config(), cube_points=16, grid_level=2, threads=1)
        app.apply_config(cfg)
        with patch('molecular_app.messagebox.showerror') as error:
            app.start_calculation()
            self.wait_finished()
            error.assert_not_called()
        self.assertTrue(app.summary['converged'])
        self.assertEqual(app.summary['electrons'], 36)
        self.assertEqual(app.summary['explicit_electrons'], 8)
        self.assertEqual(app.summary['atomic_numbers'], [40])
        self.assertAlmostEqual(app.summary['electron_count_from_density'], 8, places=8)
        self.assertIn('ECP', app.figure.axes[1].get_title())
        self.assertEqual(load_dft_config(app.current_dir / 'settings.json'), cfg)
        app.apply_config(DFTConfig())
        app.history_tree.selection_set('0')
        app.load_history_settings()
        self.assertEqual(app.read_config(), cfg)

    @unittest.skipUnless(os.environ.get('QUANTUMSTARTER_GUI_DFT_TEST') == '1', 'Enable real GUI DFT integration explicitly')
    def test_real_cobalt_newton_run(self):
        app = self.app
        app.example.set('코발트 이온 · Co2+ (2S=3 예제)')
        app.use_example()
        cfg = replace(app.read_config(), cube_points=16, grid_level=2, threads=1)
        app.apply_config(cfg)
        with patch('molecular_app.messagebox.showerror') as error:
            app.start_calculation()
            self.wait_finished()
            error.assert_not_called()
        self.assertTrue(app.summary['converged'])
        self.assertEqual(app.summary['scf_solver'], 'newton')
        self.assertIn('외부 반복', app.metrics[1].get())
        self.assertAlmostEqual(app.summary['spin_from_density'], 3, places=8)
        self.assertEqual(load_dft_config(app.current_dir / 'settings.json'), cfg)
        self.assertIn('spin_density.cube', app.field_box['values'])
        self.assertEqual(len(app.table.get_children()), len(app.summary['orbitals']))

    @unittest.skipUnless(os.environ.get('QUANTUMSTARTER_GUI_DFT_TEST') == '1', 'Enable real GUI DFT integration explicitly')
    def test_real_run_nonconvergence_export_history_and_error_recovery(self):
        app = self.app
        config = DFTConfig(basis='sto-3g', cube_points=16, grid_level=1, threads=1)
        app.apply_config(config)
        with patch('molecular_app.messagebox.showerror') as error:
            app.start_calculation()
            self.assertTrue(app.busy)
            self.wait_finished()
            error.assert_not_called()
        self.assertTrue(app.summary['converged'])
        self.assertEqual(len(app.table.get_children()), len(app.summary['orbitals']))
        self.assertAlmostEqual(app.figure.bbox.width, app.canvas._tkphoto.width(), delta=1)
        self.assertTrue((app.current_dir / 'molecular.png').exists())
        app.field.set('restricted_homo.cube')
        app.plane.set('yz')
        app.redraw()
        image = self.base / 'plot.png'
        with patch('molecular_app.filedialog.asksaveasfilename', return_value=str(image)):
            app.export_plot()
        self.assertGreater(image.stat().st_size, 1000)
        archive = self.base / 'result.zip'
        with patch('app.filedialog.asksaveasfilename', return_value=str(archive)):
            app.export_results()
        with zipfile.ZipFile(archive) as z:
            self.assertIsNone(z.testzip())
            self.assertIn('density.cube', z.namelist())
        first_dir = app.current_dir
        app.apply_config(replace(config, max_cycle=1, conv_tol=1e-12))
        self.assertIn('설정 변경됨', app.changed.get())
        app.start_calculation()
        self.wait_finished()
        self.assertFalse(app.summary['converged'])
        self.assertIn('미수렴', app.metrics[1].get())
        self.assertNotEqual(first_dir, app.current_dir)
        app.history_tree.selection_set('1')
        app.load_history_settings()
        self.assertEqual(app.read_config(), config)
        with patch('molecular_app.run_calculation', side_effect=RuntimeError('엔진 오류 검증')), \
                patch('molecular_app.messagebox.showerror') as error:
            app.start_calculation()
            self.wait_finished()
            error.assert_called_once()
        self.assertIn('이전 실행', app.changed.get())
        self.assertEqual(str(app.run_button['state']), 'normal')


if __name__ == '__main__':
    unittest.main()
