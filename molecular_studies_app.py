"""Desktop views of gas-mixture properties and molecular basis/grid comparisons."""
from pathlib import Path
from queue import Empty, Queue
from threading import Thread
import tkinter as tk
from tkinter import filedialog, messagebox, ttk

from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg

from app import open_directory
from dft_workflow import RunControl
from gas_mixture import DEFAULT_DATA, DEFAULT_MIX, mixture_figure, run_mixture
from molecular_convergence import DEFAULT_BASES, DEFAULT_GRIDS, OPERATIONS, run_comparison


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

    @staticmethod
    def table(book, label, headings):
        frame = ttk.Frame(book)
        book.add(frame, text=label)
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
