"""Desktop periodic calculations and a GPU-screening/DFT-verification workflow."""
from dataclasses import replace
from datetime import datetime
import json
import os
from pathlib import Path
import queue
import subprocess
import sys
import threading
import tkinter as tk
from tkinter import filedialog, messagebox, scrolledtext, ttk

from app import AppCanvas, VERSION
from dft_workflow import ROOT, RunControl, _wsl_path, process_options, run_calculation
from molecular import save_dft_config
from periodic import PERIODIC_BASES, PeriodicAtom, PeriodicConfig, load_calculation_config, read_periodic_structure
from periodic_plotting import create_periodic_figure


class PeriodicApp:
    def __init__(self, root):
        self.root, self.busy, self.closing = root, False, False
        self.messages = queue.Queue()
        self.canvas = None
        self.screen_path = None
        self.control = None
        self.summary, self.current_dir = None, None
        root.title(f'QuantumStarter {VERSION} | 고체·표면 · DFT / 힘 / 셀 최적화')
        root.geometry('1320x900')
        root.minsize(1100, 760)
        self.status = tk.StringVar(value='주기계 단일점 PBE · 격자와 k점을 확인하세요.')
        ttk.Label(root, text='고체·표면 계산', font=('맑은 고딕', 20, 'bold'), padding=12).pack(anchor='w')
        ttk.Label(root, textvariable=self.status, padding=8).pack(side='bottom', fill='x')
        self.tabs = ttk.Notebook(root)
        self.tabs.pack(fill='both', expand=True, padx=12, pady=8)
        page, hybrid, logs = [ttk.Frame(self.tabs, padding=12) for _ in range(3)]
        for frame, title in ((page, '주기적 DFT'), (hybrid, '표면 H 탐색 → DFT 검증'), (logs, '실행 로그')):
            self.tabs.add(frame, text=title)
        self.log = scrolledtext.ScrolledText(logs, wrap='word')
        self.log.pack(fill='both', expand=True)
        sidebar = ttk.Frame(page, width=420)
        sidebar.pack(side='left', fill='y', padx=(0, 12))
        self.run_button = ttk.Button(sidebar, text='주기적 DFT 계산 시작', command=self.start_calculation)
        self.run_button.pack(side='bottom', fill='x', pady=10)
        self.cancel_button = ttk.Button(sidebar, text='계산 중단 요청', command=self.cancel, state='disabled')
        self.cancel_button.pack(side='bottom', fill='x')
        ttk.Label(sidebar, text='PBE / GTH-PBE · 중성 셀 · FFTDF bulk: 힘·응력·셀 완화\nslab: GDF 단일점 · SOC/분산 미포함', wraplength=395).pack(side='bottom', anchor='w')
        scroll_frame = ttk.Frame(sidebar)
        scroll_frame.pack(fill='both', expand=True)
        self.settings_canvas = tk.Canvas(scroll_frame, width=420, highlightthickness=0)
        scrollbar = ttk.Scrollbar(scroll_frame, orient='vertical', command=self.settings_canvas.yview)
        scrollbar.pack(side='right', fill='y')
        self.settings_canvas.configure(yscrollcommand=scrollbar.set)
        self.settings_canvas.pack(side='left', fill='both', expand=True)
        left = ttk.Frame(self.settings_canvas)
        content_id = self.settings_canvas.create_window((0, 0), window=left, anchor='nw')
        left.bind('<Configure>', lambda _: self.settings_canvas.configure(scrollregion=self.settings_canvas.bbox('all')))
        self.settings_canvas.bind('<Configure>', lambda e: self.settings_canvas.itemconfigure(content_id, width=e.width))
        right = ttk.Frame(page)
        right.pack(side='left', fill='both', expand=True)
        actions = ttk.Frame(left)
        actions.pack(fill='x')
        for text, command in [('구조 불러오기', self.import_structure), ('설정 열기', self.load_settings), ('설정 저장', self.save_settings)]:
            ttk.Button(actions, text=text, command=command).pack(side='left', padx=2)
        examples = ttk.Frame(left)
        examples.pack(fill='x', pady=7)
        for name, label in (('silicon_bulk', 'Si 벌크'), ('pt111_small', 'Pt(111)'), ('silicon_cell', 'Si 셀 완화')):
            ttk.Button(examples, text=label,
                       command=lambda n=name: self.apply_config(load_calculation_config(ROOT / 'examples' / (n + '.json')))).pack(side='left', padx=2)
        ttk.Label(left, text='원소 x y z · Cartesian Å').pack(anchor='w')
        self.coordinates = scrolledtext.ScrolledText(left, height=8, width=44, wrap='none')
        self.coordinates.pack(fill='x', pady=4)
        ttk.Label(left, text='격자 벡터 · 각 행 a₁, a₂, a₃ · Å').pack(anchor='w')
        self.lattice = tk.Text(left, height=3, width=44)
        self.lattice.pack(fill='x', pady=4)
        fields = ttk.Frame(left)
        fields.pack(fill='x')
        self.vars = {}
        definitions = [('boundary', '경계', 'bulk', ['bulk', 'slab']),
                       ('calculation', '계산 종류', 'single_point', ['single_point', 'gradient', 'cell_optimize']),
                       ('integral_scheme', '적분 방식', 'gdf', ['gdf', 'fftdf']),
                       ('fft_mesh', 'FFT mesh (빈칸=자동)', '', None),
                       ('basis', '기저', PERIODIC_BASES[0], PERIODIC_BASES),
                       ('kmesh', 'k점 x y z', '2 2 2', None), ('spin', '셀당 스핀 2S', '0', None),
                       ('smearing_ev', 'Smearing [eV]', '.1', None), ('grid_level', '적분 격자', '3', None),
                       ('conv_tol', 'SCF 기준 [Ha]', '1e-8', None), ('precision', '적분 정밀도', '1e-7', None),
                       ('max_cycle', '최대 반복', '100', None), ('threads', 'CPU 스레드', '4', None),
                       ('memory_mb', '메모리 목표 [MB]', '4000', None),
                       ('cell_relaxation', '셀 완화 범위', 'full', ['full', 'isotropic']),
                       ('cell_max_steps', '셀 최대 단계', '50', None),
                       ('cell_force_tol', '셀 원자 힘 [eV/Å]', '.02', None),
                       ('cell_stress_tol_gpa', '셀 응력 잔차 [GPa]', '.1', None),
                       ('cell_energy_tol_ev', '셀 ΔE/ΔF [eV]', '1e-5', None),
                       ('band_path', '밴드 경로 (빈칸=생략)', '', None),
                       ('band_points', '밴드 점 수', '60', None)]
        for i, (key, label, default, choices) in enumerate(definitions):
            self.vars[key] = tk.StringVar(value=default)
            ttk.Label(fields, text=label).grid(row=i, column=0, sticky='w', pady=2)
            widget = (ttk.Combobox(fields, textvariable=self.vars[key], values=choices, state='readonly', width=24)
                      if choices else ttk.Entry(fields, textvariable=self.vars[key], width=26))
            widget.grid(row=i, column=1, sticky='ew', pady=2)
        ttk.Label(left, text='FFT 수렴 검사 mesh (세미콜론으로 구분)').pack(anchor='w', pady=(10, 0))
        self.mesh_sequence = tk.StringVar(value='21 21 21;25 25 25;29 29 29')
        ttk.Entry(left, textvariable=self.mesh_sequence).pack(fill='x')
        ttk.Label(left, text='고정 구조 · 마지막 두 구간 ΔE/F ≤ 1 meV/atom,\nΔ힘 ≤ 0.01 eV/Å, Δ응력 ≤ 0.1 GPa', wraplength=395).pack(anchor='w')
        self.mesh_button = ttk.Button(left, text='FFT mesh 수렴 확인', command=self.start_mesh_check)
        self.mesh_button.pack(fill='x', pady=5)
        self.result = scrolledtext.ScrolledText(right, height=10, wrap='word')
        self.result.pack(fill='x')
        self.result.insert('end', '구조와 격자·k점을 확인한 뒤 계산을 시작하세요.\nSi / Pt 소형 예제는 작동 확인용으로 수치 설정을 낮춘 예시입니다.\n정량 비교에는 기저·k점·진공·슬래브 두께 수렴이 필요합니다.')
        self.chart = ttk.Frame(right)
        self.chart.pack(fill='both', expand=True)
        self._build_hybrid(hybrid)
        self.apply_config(load_calculation_config(ROOT / 'examples' / 'silicon_bulk.json'))
        root.protocol('WM_DELETE_WINDOW', self.close)
        self.after_id = root.after(100, self._poll)

    def _build_hybrid(self, frame):
        ttk.Label(frame, text='후보 생성 → GPU에서 H만 완화 → 선택 후보를 PBE로 재계산', font=('맑은 고딕', 14, 'bold')).pack(anchor='w', pady=8)
        ttk.Label(frame, text='표면은 고정합니다. MLIP 상위 후보 3개와 무작위 탐색 후보 1개를 검증 대상으로 선택합니다.\n'
                  'MACE-MP medium / float64 사용. H 위치의 중복을 제거한 뒤 후보를 고릅니다. DFT는 위 탭의 수치 조건을 사용합니다.\n'
                  'DFT 결과는 선택된 구조의 상대 사이트 에너지입니다. 절대 흡착에너지·전체 사이트의 최적해를 보증하지 않습니다.', wraplength=1100).pack(anchor='w', pady=8)
        row = ttk.Frame(frame)
        row.pack(fill='x', pady=12)
        self.material = tk.StringVar(value='HEA')
        self.size = tk.StringVar(value='3 3 4')
        self.seed = tk.StringVar(value='0')
        self.acceleration = tk.BooleanVar(value=True)
        for label, variable, choices in [('표면', self.material, ['HEA', 'Pt', 'Ag', 'Au', 'Cu', 'Pd']),
                                          ('크기 x y 층', self.size, None), ('무작위 seed', self.seed, None)]:
            ttk.Label(row, text=label).pack(side='left', padx=6)
            (ttk.Combobox(row, textvariable=variable, values=choices, width=9, state='readonly') if choices else
             ttk.Entry(row, textvariable=variable, width=10)).pack(side='left')
        ttk.Checkbutton(frame, text='cuEquivariance GPU 가속 (동일 MACE 모델·float64)', variable=self.acceleration).pack(anchor='w', pady=6)
        self.screen_button = ttk.Button(frame, text='1. 후보 생성 + GPU 탐색', command=self.start_screen)
        self.screen_button.pack(anchor='w', pady=10)
        ttk.Button(frame, text='기존 screen.json 불러오기', command=self.load_screen).pack(anchor='w', pady=4)
        self.verify_button = ttk.Button(frame, text='2. 선택 후보 DFT 검증', command=self.start_verify)
        self.verify_button.pack(anchor='w', pady=10)
        self.hybrid_text = scrolledtext.ScrolledText(frame, height=22, wrap='word')
        self.hybrid_text.pack(fill='both', expand=True, pady=8)

    def apply_config(self, config):
        if self.busy:
            return
        if not isinstance(config, PeriodicConfig):
            raise ValueError('주기계 JSON 설정을 선택하세요.')
        self.base_config = config
        for widget, text in [(self.coordinates, '\n'.join(a.element + ' ' + ' '.join(f'{x:.12g}' for x in a.position) for a in config.atoms)),
                             (self.lattice, '\n'.join(' '.join(f'{x:.12g}' for x in row) for row in config.lattice))]:
            widget.delete('1.0', 'end')
            widget.insert('end', text)
        for key, variable in self.vars.items():
            value = getattr(config, key)
            variable.set(' '.join(map(str, value)) if key in ('kmesh', 'fft_mesh') else str(value))
        self._plot(dict(config=config.to_dict(), formula=config.formula, status='Input structure'))

    def numeric_settings(self):
        values = {key: variable.get().strip() for key, variable in self.vars.items()}
        values['kmesh'] = tuple(int(x) for x in values['kmesh'].split())
        values['fft_mesh'] = tuple(int(x) for x in values['fft_mesh'].split())
        for key in ('spin', 'grid_level', 'max_cycle', 'threads', 'memory_mb', 'cell_max_steps', 'band_points'):
            values[key] = int(values[key])
        for key in ('smearing_ev', 'conv_tol', 'precision', 'cell_force_tol', 'cell_stress_tol_gpa', 'cell_energy_tol_ev'):
            values[key] = float(values[key])
        return values

    def read_config(self):
        atoms = []
        for line in self.coordinates.get('1.0', 'end').splitlines():
            if line.strip():
                element, *coords = line.split()
                atoms.append(PeriodicAtom(element, tuple(float(x) for x in coords)))
        lattice = tuple(tuple(float(x) for x in row.split()) for row in self.lattice.get('1.0', 'end').splitlines() if row.strip())
        return replace(self.base_config, atoms=tuple(atoms), lattice=lattice, **self.numeric_settings())

    def import_structure(self):
        if self.busy: return
        path = filedialog.askopenfilename(parent=self.root, title='POSCAR / CIF / extxyz 선택')
        if path:
            try: self.apply_config(read_periodic_structure(path, boundary=self.vars['boundary'].get()))
            except Exception as exc: messagebox.showerror('구조 오류', str(exc), parent=self.root)

    def load_settings(self):
        if self.busy: return
        path = filedialog.askopenfilename(parent=self.root, filetypes=[('JSON', '*.json')])
        if path:
            try: self.apply_config(load_calculation_config(path))
            except Exception as exc: messagebox.showerror('설정 오류', str(exc), parent=self.root)

    def save_settings(self):
        try:
            config = self.read_config()
            path = filedialog.asksaveasfilename(parent=self.root, defaultextension='.json')
            if path: save_dft_config(config, path)
        except Exception as exc: messagebox.showerror('설정 오류', str(exc), parent=self.root)

    def _launch(self, operation, control=None):
        if self.busy: return
        self.busy = True
        self.control = control
        self.cancel_button.state(['!disabled'] if control is not None else ['disabled'])
        self.status.set('계산 중 · 실행 로그에서 진행 상황을 확인하세요.')
        for button in (self.run_button, self.screen_button, self.verify_button, self.mesh_button): button.state(['disabled'])
        def worker():
            try: self.messages.put(('done', operation()))
            except Exception as exc: self.messages.put(('error', str(exc)))
        threading.Thread(target=worker, daemon=True).start()

    def start_calculation(self):
        if self.busy: return
        try: config = self.read_config()
        except Exception as exc:
            messagebox.showerror('입력 오류', str(exc), parent=self.root)
            return
        control = RunControl()
        self._launch(lambda: ('dft', run_calculation(config, progress=lambda msg: self.messages.put(('log', msg)),
                                                  control=control)), control)

    def start_mesh_check(self):
        if self.busy: return
        from periodic_convergence import parse_meshes, run_mesh_convergence, validate_mesh_request
        try:
            config, meshes = self.read_config(), parse_meshes(self.mesh_sequence.get())
            validate_mesh_request(config, meshes, (.001, .01, .1))
        except Exception as exc:
            messagebox.showerror('입력 오류', str(exc), parent=self.root)
            return
        control = RunControl()
        self._launch(lambda: ('mesh', run_mesh_convergence(config, meshes,
                     progress=lambda msg: self.messages.put(('log', msg)), control=control)), control)

    def cancel(self):
        if self.control is not None:
            self.control.cancel()
            self.cancel_button.state(['disabled'])
            self.status.set('중단 요청 · 현재 SCF/힘 계산의 안전한 경계에서 종료합니다.')

    def start_screen(self):
        from surface_workflow import generate_candidates
        try:
            size, seed, material = tuple(map(int, self.size.get().split())), int(self.seed.get()), self.material.get()
            acceleration = 'cueq' if self.acceleration.get() else 'none'
            directory = ROOT / 'outputs' / 'surfaces' / datetime.now().strftime('%Y%m%d_%H%M%S_%f')
            manifest = generate_candidates(directory, material, size, seed)
        except Exception as exc:
            messagebox.showerror('입력 오류', str(exc), parent=self.root)
            return
        def operation():
            if sys.platform == 'win32':
                root = _wsl_path(ROOT)
                command = ['wsl', '--exec', root + '/.venv-mlip/bin/python', '-u', root + '/surface_workflow.py', 'screen', _wsl_path(manifest)]
            else:
                command = [str(ROOT / '.venv-mlip/bin/python'), '-u', str(ROOT / 'surface_workflow.py'), 'screen', str(manifest)]
            command += ['--acceleration', acceleration, '--seed', str(seed)]
            with (directory / 'screen_execution.log').open('w', encoding='utf-8') as log:
                with subprocess.Popen(command, cwd=ROOT, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                                      encoding='utf-8', errors='replace', **process_options()) as process:
                    for line in process.stdout:
                        log.write(line); log.flush()
                        self.messages.put(('log', line.rstrip()))
                    code = process.wait()
            outputs = sorted(directory.glob('screen_*/screen.json'))
            if code not in (0, 2) or not outputs:
                raise RuntimeError(f'GPU 탐색 실패. {directory / "screen_execution.log"} 확인')
            return 'screen', outputs[-1]
        self._launch(operation)

    def load_screen(self):
        if self.busy: return
        path = filedialog.askopenfilename(parent=self.root, title='screen.json 선택', filetypes=[('JSON', '*.json')])
        if path:
            try: self._show_screen(Path(path))
            except (ValueError, KeyError, OSError) as exc:
                messagebox.showerror('탐색 결과 오류', str(exc), parent=self.root)

    def _show_screen(self, path):
        data = json.loads(Path(path).read_text(encoding='utf-8'))
        if data.get('status') not in ('complete', 'failed') or 'selected_ids' not in data:
            raise ValueError('완료된 screen.json 파일이 필요합니다.')
        self.screen_path = path
        self.hybrid_text.delete('1.0', 'end')
        lines = [f"모델: {data['model'].get('model')} · {data['model'].get('hardware', '')}",
                 f"탐색 시간: {data['elapsed_seconds']:.1f}초 (모델 로딩 별도)",
                 f"구별된 H 위치: {data.get('distinct_minima', '미분류')} · 가속: {data['model'].get('acceleration', 'none')}",
                 'DFT 선택 후보: ' + ', '.join(data['selected_ids']), '',
                 '후보                     상태             상대에너지[eV]   자유 원자 최대힘[eV/Å]']
        for row in data['rows']:
            energy = f"{row['relative_energy_ev']:.5f}" if row.get('relative_energy_ev') is not None else '—'
            force = f"{row['max_free_force_ev_a']:.5f}" if row.get('max_free_force_ev_a') is not None else '—'
            lines.append(f"{row['id']:<24} {row['status']:<16} {energy:>12} {force:>18}")
            if row.get('error'): lines.append('  오류: ' + row['error'])
        lines += ['', '시작 사이트 이름은 완화 후의 사이트 종류를 보증하지 않습니다.', f'원시 결과: {path}']
        self.hybrid_text.insert('end', '\n'.join(lines))

    def start_verify(self):
        from surface_workflow import verify_selected
        if not self.screen_path:
            messagebox.showerror('탐색 결과 필요', '먼저 GPU 탐색 또는 screen.json 불러오기를 진행하세요.', parent=self.root)
            return
        try:
            settings = self.numeric_settings()
            if settings['integral_scheme'] != 'gdf' or settings['calculation'] != 'single_point':
                raise ValueError('표면 검증에는 주기계 설정의 gdf / single_point를 선택하세요.')
            # A surface job always uses 2D boundaries; preserve in-plane k points.
            settings.update(boundary='slab', kmesh=(*settings['kmesh'][:2], 1))
            settings.pop('spin')  # template parity minimum; recorded per candidate, no ground-state claim
            source = self.screen_path
        except Exception as exc:
            messagebox.showerror('설정 오류', str(exc), parent=self.root)
            return
        self._launch(lambda: ('comparison', verify_selected(source, settings=settings,
                         progress=lambda msg: self.messages.put(('log', msg)))))

    def _plot(self, summary):
        if self.canvas:
            self.canvas.get_tk_widget().destroy()
            self.canvas.figure.clear()
        figure = create_periodic_figure(summary)
        figure.set_dpi(float(self.root.winfo_fpixels('1i')))
        self.canvas = AppCanvas(figure, master=self.chart)
        self.canvas.draw()
        self.canvas.get_tk_widget().pack(fill='both', expand=True)

    def _poll(self):
        while not self.messages.empty():
            kind, value = self.messages.get()
            if kind == 'log':
                self.log.insert('end', value + '\n'); self.log.see('end')
                self.status.set(value)
                continue
            self.busy = False
            self.control = None
            self.cancel_button.state(['disabled'])
            for button in (self.run_button, self.screen_button, self.verify_button, self.mesh_button): button.state(['!disabled'])
            if kind == 'error':
                self.status.set('실행 실패 · 로그 확인')
                self.log.insert('end', value + '\n')
                messagebox.showerror('계산 실패', value, parent=self.root)
            else:
                task, data = value
                if task == 'dft':
                    summary, directory = data
                    self.summary, self.current_dir = summary, directory
                    self.result.delete('1.0', 'end')
                    lines = [f"{summary['formula']} · {summary['status']}"]
                    if 'total_energy_ev' in summary:
                        lines += [f"E = {summary['total_energy_ev']:.9f} eV / cell",
                                  f"F = {summary['free_energy_ev']:.9f} eV / cell",
                                  f"(E+F)/2 = {summary['extrapolated_energy_ev']:.9f} eV / cell"]
                    if summary.get('gradient'):
                        gradient = summary['gradient']
                        lines += [f"최대 원자 힘 = {gradient['max_atomic_force_ev_per_angstrom']:.6g} eV/Å",
                                  f"압력 = {gradient['pressure_gpa']:.6g} GPa", '응력 [GPa]:']
                        lines += ['  '.join(f'{x: .5g}' for x in row) for row in gradient['stress_gpa']]
                    if summary.get('cell_optimization'):
                        block = summary['cell_optimization']
                        lines += [f"셀 완화: {block['mode']} · {block['steps']}단계 · {block['objective']}",
                                  f"수렴 기준 충족: {block['final_criteria_met']}"]
                    lines += summary.get('warnings', []) + [f'저장: {directory}']
                    self.result.insert('end', '\n'.join(lines))
                    self._plot(summary)
                    self.canvas.figure.savefig(directory / 'periodic.png', dpi=140)
                    self.status.set(f"{summary['status']} · {directory}")
                elif task == 'mesh':
                    report, directory = data
                    self.summary, self.current_dir = report, directory
                    lines = [f"FFT mesh 수렴 확인: {report['status']} · 비교량 {report['objective']}",
                             'mesh | ΔE/F [eV/atom] | Δ힘 [eV/Å] | Δ응력 [GPa] | 통과']
                    for row in report['rows']:
                        changes = ' | '.join(f'{row[key]:.4g}' if key in row else '—' for key in
                                            ('delta_energy_ev_per_atom', 'max_force_difference_ev_per_angstrom',
                                             'max_stress_difference_gpa'))
                        lines.append(f"{row['mesh']} | {changes} | {row['passed']}")
                    lines += [f"선택 mesh: {report['selected_mesh']}",
                              '입력 구조에서의 검사입니다. 셀 완화 후 최종 구조와 k점 수렴을 별도로 확인하세요.', f'저장: {directory}']
                    self.result.delete('1.0', 'end')
                    self.result.insert('end', '\n'.join(lines))
                    self.status.set(f"FFT 검사 {report['status']} · {directory}")
                elif task == 'screen':
                    self._show_screen(data)
                    self.status.set(f'탐색 종료 · 수렴/미수렴 상태와 선택 후보를 확인하세요. {data}')
                else:
                    report = json.loads(Path(data).read_text(encoding='utf-8'))
                    comparison = report['comparison']
                    lines = [f"선택 후보 검증: {comparison['status']}",
                             f"DFT 최저 후보: {comparison.get('dft_best_id', '비교 불가')}",
                             f"MLIP 최저 후보: {comparison.get('ml_best_id', '비교 불가')}",
                             f"상대에너지 MAE: {comparison.get('relative_energy_mae_ev', '—')} eV", '',
                             '후보                     DFT 상태          DFT 외삽에너지[eV/cell]']
                    for row in report['rows']:
                        energy = f"{row['dft_energy_ev']:.8f}" if 'dft_energy_ev' in row else '—'
                        lines.append(f"{row['id']:<24} {row['status']:<18} {energy}")
                        if row.get('error'): lines.append('  오류: ' + row['error'])
                    lines += ['', report['note'], '기저·k점·진공·슬래브 두께 수렴 검증은 별도로 필요합니다.', f'원시 결과: {data}']
                    self.hybrid_text.delete('1.0', 'end')
                    self.hybrid_text.insert('end', '\n'.join(lines))
                    self.status.set(f'DFT 검증 종료 · 개별 수렴 상태 확인: {data}')
        if self.closing and not self.busy:
            self._destroy()
        else:
            self.after_id = self.root.after(100, self._poll)

    def _destroy(self):
        if self.after_id:
            self.root.after_cancel(self.after_id)
        if self.canvas:
            pending = getattr(self.canvas, '_idle_draw_id', None)
            if pending: self.canvas.get_tk_widget().after_cancel(pending)
        self.root.destroy()

    def close(self):
        if self.busy:
            self.closing = True
            self.status.set('실행 결과를 저장한 뒤 창을 닫습니다.')
        else:
            self._destroy()
