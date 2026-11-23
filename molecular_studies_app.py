"""Desktop views of gas-mixture properties, basis/grid studies and isotope shifts."""
from pathlib import Path
from queue import Empty, Queue
from threading import Thread
import tkinter as tk
from tkinter import filedialog, messagebox, ttk

from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg

from app import open_directory
from dft_workflow import RunControl
from gas_mixture import DEFAULT_DATA, DEFAULT_MIX, mixture_figure, run_mixture
from engines.infrared import DEFAULT_FWHM_CM1
from molecular_convergence import DEFAULT_BASES, DEFAULT_GRIDS, OPERATIONS, run_comparison
from molecular_isotopes import WATER_IR_EXAMPLE, isotope_spectrum, run_isotopes
from molecular_plotting import infrared_figure


class MolecularStudies:
    def __init__(self, root, molecular_app):
        self.root, self.parent = root, molecular_app
        self.busy = self.closing = False
        self.files = (DEFAULT_DATA,)
        self.events = Queue()
        self.canvas = None
        self.root.title('혼합기체와 분자 계산 비교')
        self.root.geometry('960x640')
        book = ttk.Notebook(root)
        book.pack(fill='both', expand=True)
        mix = ttk.Frame(book, padding=14)
        conv = ttk.Frame(book, padding=14)
        book.add(mix, text='혼합기체 비열·음속')
        book.add(conv, text='기저·격자 비교')
        self.build_isotope_tab(book)

        bar = ttk.Frame(mix)
        bar.pack(fill='x')
        ttk.Label(bar, text='몰분율 (합 1)').pack(side='left')
        self.mix = tk.StringVar(value=DEFAULT_MIX)
        ttk.Entry(bar, textvariable=self.mix, width=33).pack(side='left', padx=6)
        ttk.Button(bar, text='계수 파일 선택', command=self.choose_files).pack(side='left')
        ttk.Button(bar, text='계산', command=self.calculate_mix).pack(side='left', padx=6)
        self.coeff_label = tk.StringVar(value='기본 N₂·O₂·CO₂ NASA7 · 1 bar')
        ttk.Label(mix, textvariable=self.coeff_label).pack(anchor='w', pady=6)
        mix_views = ttk.Notebook(mix)
        mix_views.pack(fill='both', expand=True)
        self.chart = ttk.Frame(mix_views)
        mix_views.add(self.chart, text='온도별 그래프')
        self.properties = self.table(mix_views, '온도별 수치', ('T (K)', 'Cp (J/kg/K)', 'Cv (J/kg/K)', 'γ',
                                  '음속 (m/s)', 'γ=1.4 (m/s)', '편차 (%)'))
        self.mix_info = tk.StringVar(value='298.15–2,000 K · 고정 몰분율 이상기체')
        ttk.Label(mix, textvariable=self.mix_info).pack(anchor='w')
        self.mix_folder = ttk.Button(mix, text='결과 폴더 열기', state='disabled',
                                     command=lambda: open_directory(self.mix_dir))
        self.mix_folder.pack(anchor='e')

        ttk.Label(conv, text='분자 화면의 현재 구조와 계산 종류를 사용합니다. UV–Vis, 진동수, 스핀 해리, 기체 열물성을 선택할 수 있습니다.').pack(anchor='w')
        bar = ttk.Frame(conv)
        bar.pack(fill='x', pady=8)
        self.bases = tk.StringVar(value=' '.join(DEFAULT_BASES))
        self.grids = tk.StringVar(value=' '.join(map(str, DEFAULT_GRIDS)))
        for label, var, width in [('기저', self.bases, 26), ('격자 수준', self.grids, 8)]:
            ttk.Label(bar, text=label).pack(side='left')
            ttk.Entry(bar, textvariable=var, width=width).pack(side='left', padx=6)
        self.start = ttk.Button(bar, text='비교 시작', command=self.start_comparison)
        self.start.pack(side='left')
        self.stop = ttk.Button(bar, text='중단', state='disabled', command=lambda: self.control.cancel())
        self.stop.pack(side='left', padx=6)
        self.conv_folder = ttk.Button(bar, text='결과 폴더 열기', state='disabled',
                                      command=lambda: open_directory(self.conv_dir))
        self.conv_folder.pack(side='right')
        self.status = tk.StringVar(value='기저와 격자 값은 공백으로 구분하세요.')
        ttk.Label(conv, textvariable=self.status, wraplength=880).pack(anchor='w', pady=4)
        tables = ttk.Notebook(conv)
        tables.pack(fill='both', expand=True)
        self.runs = self.table(tables, '실행 시간', ('기저', '격자', '상태', '경과 시간 (s)'))
        self.diffs = self.table(tables, '조건별 변화량', ('변경 축', '이전 조건', '다음 조건', '물성', '순번 / K', '단위', '차이'))
        self.root.protocol('WM_DELETE_WINDOW', self.close)

    def build_isotope_tab(self, book):
        page = ttk.Frame(book, padding=14)
        book.add(page, text='동위원소 진동수 · IR')
        first = ttk.Frame(page)
        first.pack(fill='x')
        ttk.Label(first, text='기준 진동수 파일').pack(side='left')
        source = WATER_IR_EXAMPLE
        if self.parent.current_dir and (self.parent.summary or {}).get('calculation') == 'frequencies':
            source = Path(self.parent.current_dir) / 'results.json'
        self.isotope_source = tk.StringVar(value=str(source))
        ttk.Entry(first, textvariable=self.isotope_source).pack(side='left', fill='x', expand=True, padx=6)
        ttk.Button(first, text='파일 선택', command=self.choose_isotope_source).pack(side='left')
        ttk.Button(first, text='물 예제', command=lambda: self.isotope_source.set(str(WATER_IR_EXAMPLE))).pack(side='left', padx=4)
        second = ttk.Frame(page)
        second.pack(fill='x', pady=8)
        ttk.Label(second, text='원자 순번=질량(u)').pack(side='left')
        self.isotope_masses = tk.StringVar()
        ttk.Entry(second, textvariable=self.isotope_masses, width=20).pack(side='left', padx=6)
        ttk.Label(second, text='빈칸은 H 전체 → D   예시  2=2.0141').pack(side='left')
        ttk.Label(second, text='IR 선폭 FWHM (cm⁻¹)').pack(side='left', padx=(14, 0))
        self.ir_width = tk.StringVar(value=f'{DEFAULT_FWHM_CM1:g}')
        ttk.Entry(second, textvariable=self.ir_width, width=6).pack(side='left', padx=5)
        ttk.Button(second, text='비교', command=self.calculate_isotopes).pack(side='right')
        self.isotope_info = tk.StringVar(value='기준과 치환 구조의 진동수 차이와 영점에너지(ZPE)를 계산합니다.')
        ttk.Label(page, textvariable=self.isotope_info, wraplength=850).pack(anchor='w', pady=6)
        self.isotope_folder = ttk.Button(page, text='결과 폴더 열기', state='disabled',
                                         command=lambda: open_directory(self.isotope_dir))
        self.isotope_folder.pack(side='bottom', anchor='e')
        self.isotope_views = ttk.Notebook(page)
        self.isotope_views.pack(fill='both', expand=True)
        self.isotope_modes = self.table(self.isotope_views, '모드별 수치', ('기준 모드', '치환 모드', '기준 (cm⁻¹)',
            '치환 (cm⁻¹)', '차이 (cm⁻¹)', '모드 겹침', '기준 IR (km/mol)', '치환 IR (km/mol)'))
        self.ir_chart = ttk.Frame(self.isotope_views)
        self.isotope_views.add(self.ir_chart, text='IR 스펙트럼', state='disabled')
        self.ir_canvas = None

    def choose_isotope_source(self):
        filename = filedialog.askopenfilename(parent=self.root, title='분자 진동수 결과 선택',
                                              filetypes=[('진동수 결과', '*.json')])
        if filename:
            self.isotope_source.set(filename)

    def calculate_isotopes(self):
        assignments = self.isotope_masses.get().split()
        try:
            report, directory = run_isotopes(self.isotope_source.get(), self.parent.output_root / 'isotopes',
                                             deuterate=not assignments, assignments=assignments, fwhm=self.ir_width.get())
        except (ValueError, OSError) as exc:
            messagebox.showerror('동위원소 비교', str(exc), parent=self.root)
            return
        self.isotope_modes.delete(*self.isotope_modes.get_children())
        for row in report['modes']:
            strengths = ['—' if row.get(key) is None else f'{row[key]:.3f}'
                         for key in ('reference_ir_km_mol', 'substituted_ir_km_mol')]
            self.isotope_modes.insert('', 'end', values=(row['reference_mode'], row['substituted_mode'],
                f"{row['reference_cm1']:.2f}", f"{row['substituted_cm1']:.2f}", f"{row['shift_cm1']:+.2f}",
                f"{row['absolute_overlap']:.3f}", *strengths))
        if self.ir_canvas:
            self.ir_canvas.get_tk_widget().destroy()
            self.ir_canvas = None
        if report['ir_labels']:
            grid, curves = isotope_spectrum(report)
            figure = infrared_figure(grid, curves, report['ir_labels'], report['ir_fwhm_cm1'])
            self.ir_canvas = FigureCanvasTkAgg(figure, master=self.ir_chart)
            self.isotope_views.tab(self.ir_chart, state='normal')
            self.isotope_views.select(self.ir_chart)
            self.ir_canvas.get_tk_widget().pack(fill='both', expand=True)
            self.ir_canvas.draw()
        else:
            self.isotope_views.select(0)
            self.isotope_views.tab(self.ir_chart, state='disabled')
        change = 'H 전체 → D' if not assignments else ', '.join(assignments) + ' u'
        energy = (f"{report['formula'] or '기준 구조'} · {change} · ZPE {report['reference_zpe_kj_mol']:.3f} → "
                  f"{report['substituted_zpe_kj_mol']:.3f} kJ/mol (Δ {report['delta_zpe_kj_mol']:+.3f} kJ/mol)")
        self.isotope_info.set('\n'.join([energy, *report['warnings']]))
        self.isotope_dir = directory
        self.isotope_folder.configure(state='normal')

    @staticmethod
    def table(book, label, headings):
        frame = ttk.Frame(book)
        if label:
            book.add(frame, text=label)
        else:
            frame.pack(fill='both', expand=True)
        keys = tuple(str(i) for i in range(len(headings)))
        tree = ttk.Treeview(frame, columns=keys, show='headings')
        for key, heading in zip(keys, headings):
            tree.heading(key, text=heading)
            tree.column(key, width=108, anchor='center')
        scroll = ttk.Scrollbar(frame, command=tree.yview)
        tree.configure(yscrollcommand=scroll.set)
        scroll.pack(side='right', fill='y')
        tree.pack(fill='both', expand=True)
        return tree

    def choose_files(self):
        files = filedialog.askopenfilenames(parent=self.root, title='NASA7 계수 선택',
                                            filetypes=[('NASA7', '*.yaml'), ('JSON', '*.json')])
        if files:
            self.files = tuple(map(Path, files))
            self.coeff_label.set(' · '.join(str(path) for path in self.files))

    def calculate_mix(self):
        try:
            report, directory = run_mixture(self.mix.get(), self.parent.output_root/'gas_mixtures', self.files)
        except (ValueError, KeyError, TypeError, OSError) as exc:
            messagebox.showerror('혼합기체 계산', str(exc), parent=self.root)
            return
        if self.canvas:
            self.canvas.get_tk_widget().destroy()
        self.canvas = FigureCanvasTkAgg(mixture_figure(report), master=self.chart)
        self.canvas.draw()
        self.canvas.get_tk_widget().pack(fill='both', expand=True)
        self.properties.delete(*self.properties.get_children())
        for r in report['rows']:
            self.properties.insert('', 'end', values=(f"{r['temperature_K']:g}", f"{r['cp_J_kg_K']:.3f}",
                f"{r['cv_J_kg_K']:.3f}", f"{r['gamma']:.5f}", f"{r['sound_speed_m_s']:.3f}",
                f"{r['sound_speed_gamma14_m_s']:.3f}", f"{r['gamma14_error_percent']:+.3f}"))
        self.mix_dir = directory
        self.mix_folder.configure(state='normal')
        row = report['rows'][-1]
        self.mix_info.set(f"2,000 K에서 γ = {row['gamma']:.5f}, 음속 {row['sound_speed_m_s']:.2f} m/s, "
                          f"γ=1.4 가정의 상대 편차 {row['gamma14_error_percent']:+.2f}%")

    def start_comparison(self):
        if self.busy:
            return
        try:
            cfg = self.parent.read_config()
            if cfg.calculation not in OPERATIONS:
                raise ValueError('분자 화면에서 UV–Vis, 진동수, 스핀 해리 또는 기체 열물성 계산을 먼저 선택하세요.')
            grids = tuple(map(int, self.grids.get().split()))
            bases = self.bases.get().lower().split()
            backend = self.parent.backend.get()
        except ValueError as exc:
            messagebox.showerror('비교 설정', str(exc), parent=self.root)
            return
        self.control = RunControl()
        self.busy = True
        self.start.configure(state='disabled')
        self.stop.configure(state='normal')
        self.conv_folder.configure(state='disabled')
        self.runs.delete(*self.runs.get_children())
        self.diffs.delete(*self.diffs.get_children())

        def work():
            try:
                report, directory = run_comparison(cfg, bases, grids, self.parent.output_root/'basis_grid',
                    backend, lambda text: self.events.put(('progress', text)), self.control)
                self.events.put(('done', report, directory))
            except Exception as exc:
                self.events.put(('error', str(exc)))

        Thread(target=work, daemon=True).start()
        self.root.after(75, self.poll)

    def poll(self):
        try:
            while True:
                event, *data = self.events.get_nowait()
                if event == 'progress':
                    self.status.set(data[0])
                    continue
                self.busy = False
                self.start.configure(state='normal')
                self.stop.configure(state='disabled')
                if event == 'error':
                    self.status.set(data[0])
                else:
                    report, self.conv_dir = data
                    self.conv_folder.configure(state='normal')
                    self.status.set(f"{report['status']} · {self.conv_dir}")
                    for r in report['runs']:
                        self.runs.insert('', 'end', values=(r['basis'], r['grid'], r['status'], f"{r['wall_seconds']:.2f}"))
                    for r in report['differences']:
                        value = f"{r['delta']:+.6g}" if r['delta'] is not None else r['status']
                        self.diffs.insert('', 'end', values=(r['axis'], f"{r['from_basis']} / {r['from_grid']}",
                            f"{r['to_basis']} / {r['to_grid']}", r['observable'], r['rank'], r['unit'], value))
                if self.closing:
                    self.root.destroy()
                return
        except Empty:
            self.root.after(75, self.poll)

    def close(self):
        if self.busy:
            self.closing = True
            self.control.cancel()
            self.status.set('현재 계산을 중단하고 기록을 저장하는 중입니다.')
        else:
            self.root.destroy()
