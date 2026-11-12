"""실제 Tk 창에서 입력→계산→저장→복원 검증. 화면이 없는 환경에서는 건너뛴다."""
from dataclasses import asdict
import hashlib
import json
import os
from pathlib import Path
import tempfile
import time
import tkinter as tk
import unittest
from unittest.mock import patch
import zipfile

from app import QuantumApp
from solver import Config

# QUANTUMSTARTER_TK_HIDDEN=1: 창을 매핑 전에 숨겨 사용자 화면에 아무것도 띄우지 않는다
# (사용자 요청 2026-10-29, test_molecular_gui와 같은 방식). 화면 픽셀이 필요한 검사만 생략한다.
HIDDEN = os.environ.get('QUANTUMSTARTER_TK_HIDDEN') == '1'


class GuiValidation(unittest.TestCase):
    def setUp(self):
        try:
            self.root = tk.Tk()
        except tk.TclError:
            self.skipTest('GUI 검증에는 그래픽 화면이 필요합니다.')
        if HIDDEN:
            self.root.withdraw()  # 첫 이벤트 처리 전에 숨겨 한 번도 매핑되지 않게 한다.
            create_toplevel = tk.Toplevel.__init__

            def hidden_toplevel(window, *args, **kwargs):
                create_toplevel(window, *args, **kwargs)
                window.withdraw()
            patcher = patch.object(tk.Toplevel, '__init__', hidden_toplevel)
            patcher.start()
            self.addCleanup(patcher.stop)
        self.temp = tempfile.TemporaryDirectory()
        self.base = Path(self.temp.name)
        self.callback_errors = []
        self.root.report_callback_exception = lambda *args: self.callback_errors.append(args)
        self.app = QuantumApp(self.root, output_root=self.base/'runs')
        self.root.update()

    def tearDown(self):
        if hasattr(self, 'app'):
            if not self.app.closed:
                self.app.close()
            self.temp.cleanup()

    def wait_finished(self):
        deadline = time.monotonic() + 30
        while self.app.busy and not self.app.closed and time.monotonic() < deadline:
            self.root.update()
            time.sleep(0.02)
        self.assertFalse(self.app.busy, 'GUI 계산이 완료되지 않았습니다.')
        self.assertFalse(self.callback_errors, self.callback_errors)
        if not self.app.closed:
            self.root.update()

    def test_calculate_export_restore_and_recover(self):
        app = self.app
        app.vars['grid_points'].set('abc')
        with patch('app.messagebox.showerror') as error:
            app.start_calculation()
            error.assert_called_once()
        self.assertFalse(app.busy)
        app.apply_config(Config(), initial=True)
        app.run_button.invoke()
        self.assertTrue(app.busy)
        self.assertIn('disabled', str(app.run_button['state']))
        self.wait_finished()
        self.assertEqual(app.result.config, Config())
        self.assertEqual(len(app.table.get_children()), 3)
        # 화면 픽셀보다 도표가 크게 그려져 축/범례가 잘리는 DPI 회귀 방지(화면 배치가 필요해 숨김 모드에서는 생략).
        if not HIDDEN:
            self.assertAlmostEqual(app.figure.bbox.width, app.canvas._tkphoto.width(), delta=1)
            self.assertAlmostEqual(app.figure.bbox.height, app.canvas._tkphoto.height(), delta=1)
        first = app.current_dir
        self.assertTrue((first/'wavefunctions.png').is_file())
        self.assertTrue((first/'execution.log').is_file())
        original_hash = hashlib.sha256((first/'results.json').read_bytes()).hexdigest()

        # 선택한 상태와 단위를 적용한 내장 그래프를 실제로 저장한다.
        app.state_filter.set('2s')
        app.unit.set('Å')
        app.compare.set(False)
        app.redraw()
        self.assertEqual(len(app.figure.axes[1].lines), 1)
        plot = self.base/'선택_그래프.png'
        with patch('app.filedialog.asksaveasfilename', return_value=str(plot)):
            app.export_plot()
        self.assertGreater(plot.stat().st_size, 1000)
        settings = self.base/'계산_설정.json'
        with patch('app.filedialog.asksaveasfilename', return_value=str(settings)):
            app.save_settings()
        self.assertEqual(json.loads(settings.read_text()), asdict(Config()))

        # 두 번째 실행은 새 폴더에 저장하며 첫 결과를 바꾸지 않는다.
        app.target.set('헬륨 이온 He+')
        app.precision.set('권장')
        app._target_changed()
        self.assertIn('설정 변경됨', app.changed.get())
        app.run_button.invoke()
        self.wait_finished()
        self.assertEqual(app.result.config.nuclear_charge, 2)
        self.assertNotEqual(app.current_dir, first)
        self.assertEqual(hashlib.sha256((first/'results.json').read_bytes()).hexdigest(), original_hash)
        self.assertEqual(len(app.history), 2)
        archive = self.base/'결과.zip'
        with patch('app.filedialog.asksaveasfilename', return_value=str(archive)):
            app.export_results()
        with zipfile.ZipFile(archive) as bundle:
            self.assertIsNone(bundle.testzip())
            self.assertTrue({'results.json', 'radial.csv', 'wavefunctions.png', 'settings.json', 'execution.log'} <= set(bundle.namelist()))

        # 저장한 설정을 다시 읽어도 이전 결과와 섞이지 않는다.
        with patch('app.filedialog.askopenfilename', return_value=str(settings)):
            app.load_settings()
        self.assertEqual(app.read_config(), Config())
        self.assertEqual(app.result.config.nuclear_charge, 2)
        self.assertIn('설정 변경됨', app.changed.get())
        app.history_tree.selection_set('1')
        app.load_history_settings()
        self.assertEqual(app.read_config(), Config())

        # 계산 실패 후에도 다시 실행할 수 있고 오류 원인을 확인할 수 있다.
        with patch('app.solve', side_effect=ValueError('의도한 오류 검증')), patch('app.messagebox.showerror') as error:
            app.start_calculation()
            self.wait_finished()
            error.assert_called_once()
        self.assertIn('normal', str(app.run_button['state']))
        self.assertTrue(any('의도한 오류' in line for line in app.log_lines))

    def test_close_waits_for_saving(self):
        self.app.start_calculation()
        self.app.close()
        self.assertTrue(self.app.closing)
        self.wait_finished()
        self.assertTrue(self.app.closed)
        self.assertTrue((self.app.current_dir/'wavefunctions.png').is_file())


if __name__ == '__main__':
    unittest.main()
