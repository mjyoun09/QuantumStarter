"""분자 DFT 기본 진입점과 기존 수소 검증 화면. python app.py [--hydrogen-check]."""
from dataclasses import asdict
from datetime import datetime
import json
import os
from pathlib import Path
import queue
import subprocess
import sys
import threading
import time
import traceback
import zipfile
import tkinter as tk
from tkinter import filedialog, messagebox, scrolledtext, ttk

from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk

from plotting import create_figure
from solver import BOHR_TO_ANGSTROM, HARTREE_TO_EV, Config, solve
from workflow import (ORBITALS, PRESETS, VERSION, config_from_fields, list_runs,
                      load_config, new_run_dir, recommended_config, save_config, save_result)

ROOT = Path(__file__).resolve().parent
BG, WHITE, INK, MUTED, BLUE = '#f2f5fb', '#ffffff', '#17253d', '#607089', '#2563eb'
TARGETS = {'수소 H': 1, '헬륨 이온 He+': 2, '리튬 이온 Li2+': 3}

HELP = """처음 실행하는 순서

1. 왼쪽에서 수소 H를 선택합니다.
2. 정밀도는 권장으로 두고 계산 시작을 누릅니다.
3. 결과 탭에서 에너지, 오차, 그래프를 확인합니다.
4. 결과 폴더 열기를 누르면 그래프·CSV·설정 파일을 볼 수 있습니다.

설정을 바꾸고 다시 계산하기

방위 양자수 l을 s에서 p로 바꾸면 2p부터 계산합니다. 상태 수는 해당 l에서 낮은 순서대로 구할 상태 수입니다.
정밀도를 정밀로 바꾸면 격자점 수가 늘어납니다. 고급 탭에서 반지름과 격자점을 직접 수정하면 직접 설정으로 표시됩니다.
설정을 바꿔도 이전 결과를 새 조건의 결과처럼 표시하지 않습니다. 계산 시작을 다시 누르면 반영됩니다.

그래프 읽기

왼쪽은 u(r)=rR(r), 오른쪽은 P(r)=|u(r)|²입니다. P(r)dr은 구껍질에서 전자를 발견할 확률입니다.
3차원 부피당 확률밀도 |psi|²와 다릅니다. 실선은 수치해, 점선은 동일한 모델의 해석해입니다.
상태 선택으로 한 상태만 볼 수 있습니다. 거리 단위를 Å로 바꾸면 함수의 세로축 단위도 함께 변환됩니다.
영문 그래프 도구막대에서 확대·이동·초기화가 가능합니다. 현재 그래프 저장 버튼으로 화면에 선택한 그래프를 저장합니다.

결과 보관과 재현

매번 outputs/runs 안에 새 폴더를 만듭니다. 기존 실행을 덮어쓰지 않습니다.
각 폴더의 settings.json에는 재실행에 필요한 모든 계산 설정이 있습니다.
실행 기록 탭에서 과거 조건 불러오기를 누르고 다시 계산할 수 있습니다.
결과 ZIP 내보내기는 마지막 실행의 원시 데이터와 기본 그래프 전체를 묶습니다.
그래프 필터 변경은 원시 데이터나 저장된 기본 그래프를 바꾸지 않습니다.

현재 계산 범위

수소와 한 전자 이온만 계산합니다. He+는 전자 한 개이며 중성 헬륨이 아닙니다.
이 버전은 DFT 이전의 수치 양자역학 모델입니다. 촉매, 분자, 전자 간 상호작용, SCF는 아직 포함하지 않습니다.
사용법은 간단하게 만들되 물리적 가정과 계산 정보는 설정과 원리 탭에 그대로 제공합니다.
"""

THEORY = """현재 모델과 계산 원리

한 전자 원자의 방사형 슈뢰딩거 방정식을 풉니다. 아직 DFT/SCF가 아닙니다.
핵은 고정된 무한 질량 점핵이며, 전자는 비상대론적으로 취급합니다.
스핀, 미세구조, 유한 핵질량 효과, 분자 및 주기적 표면은 포함하지 않습니다.

방정식, 원자 단위계

[-1/2 d²/dr² + l(l+1)/(2r²) - Z/r] u(r) = E u(r)
u(r) = r R(r),   u(0) = u(r_max) = 0

내부 격자 간격 dr = r_max/(N+1)입니다. r=0은 경계이므로 핵 퍼텐셜의 0 나눗셈을 피합니다.
2차 중앙차분으로 삼중대각 해밀토니안을 구성합니다.
대각 성분 = 1/dr² + l(l+1)/(2r²) - Z/r
인접 비대각 성분 = -1/(2dr²)
SciPy 고유값 계산기로 낮은 상태를 구합니다. 에너지 해석식으로 수치 결과를 생성하지 않습니다.

검증 기준

E_n = -Z²/(2n²) hartree이며, 고정된 l에서 첫 주양자수는 n=l+1입니다.
해석적 파동함수는 그래프 비교와 부호 정렬에만 사용합니다.
정규화는 sum(|u|²)dr=1입니다. 서로 다른 n의 방사형 상태는 직교합니다.
1 hartree는 약 27.2114 eV, 1 bohr는 약 0.529177 Å입니다.

정확도 해석

표의 오차는 동일한 무한 핵질량 모델의 해석해와의 절대 에너지 차이입니다.
실험값에 대한 오차, 촉매 성능의 정확도, 이후 구현할 DFT의 정확도가 아닙니다.
권장·정밀은 격자 선택의 출발점이며 정확도를 보증하는 등급이 아닙니다.
격자 수렴: 반지름을 고정하고 격자점을 두 배로 늘려 비교합니다.
경계 수렴: 반지름과 격자점을 함께 두 배로 늘려 간격을 유지하며 비교합니다.
들뜬 상태일수록 큰 영역이 필요합니다. 양의 에너지를 원자의 결합상태로 해석하지 마세요.

원시 데이터

results.json: 모든 조건, 상수, 모델 가정, 버전, 에너지, 정규화, 평균 반지름.
energies.csv: 상태별 에너지, 오차, 정규화, 평균 반지름.
radial.csv: r, 수치 u(r), P(r), 해석 u(r). 기본 단위는 bohr입니다.
settings.json: 재실행 설정. execution.log: 실행 단계와 완료 상태.
wavefunctions.png: 모든 상태와 해석해를 표시한 기본 그래프.
"""


class AppCanvas(FigureCanvasTkAgg):
    """Tk가 제공하는 화면 픽셀과 도표 픽셀을 일치시켜 축/범례 잘림을 방지한다."""

    def _update_device_pixel_ratio(self, event=None):
        # 글꼴 크기는 Figure DPI로 반영한다. Tk 위젯에 배율을 이중 적용하지 않는다.
        return


def open_directory(path):
    path = Path(path).resolve()
    if not path.is_dir():
        raise OSError('결과 폴더를 찾을 수 없습니다.')
    if sys.platform == 'win32':
        os.startfile(str(path))
    elif sys.platform == 'darwin':
        subprocess.Popen(['open', str(path)])
    else:
        subprocess.Popen(['xdg-open', str(path)])


class QuantumApp:
    def __init__(self, root, output_root=None):
        self.root = root
        self.output_root = Path(output_root) if output_root else ROOT / 'outputs' / 'runs'
        self.busy = False
        self.closing = False
        self.closed = False
        self.applying = True
        self.messages = queue.Queue()
        self.result = self.summary = self.current_dir = None
        self.figure = self.canvas = self.toolbar = None
        self.controls = []
        self.history = []
        self.log_lines = []
        self.vars = {k: tk.StringVar(value=str(v)) for k, v in asdict(Config()).items()}
        self.target = tk.StringVar(value='수소 H')
        self.orbital = tk.StringVar(value='s · l=0')
        self.precision = tk.StringVar(value='권장')
        self.status = tk.StringVar(value='준비됨 · 수소 H / 권장으로 첫 계산을 시작하세요.')
        self.preview = tk.StringVar()
        self.changed = tk.StringVar(value='계산 결과가 아직 없습니다.')
        self.state_filter = tk.StringVar(value='전체 상태')
        self.unit = tk.StringVar(value='bohr')
        self.compare = tk.BooleanVar(value=True)
        self._setup_styles()
        self._build()
        for variable in self.vars.values():
            variable.trace_add('write', self._settings_changed)
        self.applying = False
        try:
            if (ROOT / 'config.json').exists():
                cfg = load_config(ROOT / 'config.json')
                self.apply_config(cfg, initial=True)
        except (ValueError, TypeError, OSError) as exc:
            self.status.set('기본 설정으로 시작했습니다. config.json을 확인하세요.')
            self._log(f'설정 파일 읽기 실패, 기본값 적용: {exc}')
        self._settings_changed()
        self.refresh_history()
        root.protocol('WM_DELETE_WINDOW', self.close)
        root.after(80, self._poll)

    def _setup_styles(self):
        root = self.root
        root.title(f'QuantumStarter {VERSION} | 양자 계산 워크벤치')
        root.geometry('1360x920')
        root.minsize(1140, 780)
        root.configure(bg=BG)
        import tkinter.font as tkfont
        available = set(tkfont.families(root))
        family = next((f for f in ('맑은 고딕', 'Malgun Gothic', 'Noto Sans CJK KR',
                                   'Noto Sans CJK JP', 'Apple SD Gothic Neo') if f in available),
                      tkfont.nametofont('TkDefaultFont').actual('family'))
        self.font = family
        for name in ('TkDefaultFont', 'TkTextFont', 'TkMenuFont'):
            tkfont.nametofont(name).configure(family=family, size=10)
        style = ttk.Style(root)
        style.theme_use('clam')
        style.configure('.', font=(family, 10), foreground=INK)
        style.configure('TFrame', background=BG)
        style.configure('Card.TFrame', background=WHITE)
        style.configure('TLabel', background=BG, foreground=INK)
        style.configure('Card.TLabel', background=WHITE)
        style.configure('Muted.TLabel', background=WHITE, foreground=MUTED, font=(family, 9))
        style.configure('Title.TLabel', background=WHITE, foreground=INK, font=(family, 20, 'bold'))
        style.configure('Section.TLabel', background=WHITE, font=(family, 12, 'bold'))
        style.configure('Big.TLabel', background=WHITE, font=(family, 16, 'bold'), foreground=INK)
        style.configure('TButton', padding=(10, 7), background='#eaf0fa', borderwidth=0)
        style.map('TButton', background=[('active', '#dae6f8'), ('disabled', '#edf1f7')])
        style.configure('Primary.TButton', background=BLUE, foreground=WHITE, font=(family, 12, 'bold'), padding=12)
        style.map('Primary.TButton', background=[('disabled', '#94a3b8'), ('active', '#1d4ed8')],
                  foreground=[('disabled', WHITE)])
        style.configure('TNotebook', background=BG, borderwidth=0)
        style.configure('TNotebook.Tab', padding=(14, 10), background='#e8eef8')
        style.map('TNotebook.Tab', background=[('selected', WHITE)], foreground=[('selected', BLUE)])
        style.configure('TEntry', fieldbackground=WHITE, padding=6)
        style.configure('TCombobox', fieldbackground=WHITE, padding=6)
        style.configure('Treeview', background=WHITE, fieldbackground=WHITE, rowheight=27, borderwidth=0)
        style.configure('Treeview.Heading', background='#edf2fa', font=(family, 9, 'bold'), padding=7)
        style.configure('TCheckbutton', background=WHITE)

    def _button(self, parent, text, command, tracked=False, **kwargs):
        widget = ttk.Button(parent, text=text, command=command, **kwargs)
        if tracked:
            self.controls.append((widget, 'normal'))
        return widget

    def _build(self):
        header = ttk.Frame(self.root, style='Card.TFrame', padding=(24, 14))
        header.pack(fill='x')
        ttk.Label(header, text='QuantumStarter', style='Title.TLabel').pack(side='left')
        subtitle = ttk.Frame(header, style='Card.TFrame', padding=(20, 0))
        subtitle.pack(side='left')
        ttk.Label(subtitle, text='양자 계산 워크벤치', style='Section.TLabel').pack(anchor='w')
        ttk.Label(subtitle, text='설정부터 결과 해석까지, 한 화면에서', style='Muted.TLabel').pack(anchor='w', pady=(3, 0))
        ttk.Label(header, text='현재 모델  ·  수소 / 한 전자 이온', style='Muted.TLabel').pack(side='right')
        body = ttk.Frame(self.root, padding=(18, 16))
        body.pack(fill='both', expand=True)
        sidebar = ttk.Frame(body, style='Card.TFrame', padding=18, width=300)
        sidebar.pack(side='left', fill='y', padx=(0, 16))
        sidebar.pack_propagate(False)
        self._build_settings(sidebar)
        self.tabs = ttk.Notebook(body)
        self.tabs.pack(side='left', fill='both', expand=True)
        self.result_tab = ttk.Frame(self.tabs, style='Card.TFrame', padding=16)
        self.detail_tab = ttk.Frame(self.tabs, style='Card.TFrame', padding=16)
        self.history_tab = ttk.Frame(self.tabs, style='Card.TFrame', padding=16)
        self.log_tab = ttk.Frame(self.tabs, style='Card.TFrame', padding=16)
        self.help_tab = ttk.Frame(self.tabs, style='Card.TFrame', padding=16)
        for frame, label in ((self.result_tab, '계산 결과'), (self.detail_tab, '설정과 원리'),
                             (self.history_tab, '실행 기록'), (self.log_tab, '실행 로그'), (self.help_tab, '사용 안내')):
            self.tabs.add(frame, text=label)
        self._build_results()
        self.detail_text = self._text(self.detail_tab, THEORY)
        self.log_text = self._text(self.log_tab, '실행 단계와 오류 원인이 여기에 표시됩니다.\n')
        self._text(self.help_tab, HELP)
        self._build_history()
        footer = ttk.Frame(self.root, padding=(22, 7))
        footer.pack(side='bottom', fill='x', before=body)
        ttk.Label(footer, textvariable=self.status, foreground=MUTED).pack(side='left')
        ttk.Label(footer, text='에너지 eV · 거리 bohr / Å', foreground=MUTED).pack(side='right')

    def _settings_page(self, title):
        outer = ttk.Frame(self.settings_tabs, style='Card.TFrame')
        self.settings_tabs.add(outer, text=title)
        canvas = tk.Canvas(outer, bg=WHITE, borderwidth=0, highlightthickness=0, height=1)
        scroll = ttk.Scrollbar(outer, orient='vertical', command=canvas.yview)
        scroll.pack(side='right', fill='y')
        canvas.pack(side='left', fill='both', expand=True)
        canvas.configure(yscrollcommand=scroll.set)
        inner = ttk.Frame(canvas, style='Card.TFrame', padding=(8, 12))
        window = canvas.create_window((0, 0), window=inner, anchor='nw')
        inner.bind('<Configure>', lambda event: canvas.configure(scrollregion=canvas.bbox('all')))
        canvas.bind('<Configure>', lambda event: canvas.itemconfigure(window, width=event.width))
        return inner

    def _build_settings(self, parent):
        ttk.Label(parent, text='01  계산 설정', style='Section.TLabel').pack(anchor='w')
        ttk.Label(parent, text='처음에는 권장 설정으로 시작하세요.', style='Muted.TLabel').pack(anchor='w', pady=(5, 12))
        # 실행 버튼을 먼저 아래쪽에 고정한다. 작은 창에서도 버튼이 밀려나지 않는다.
        bottom = ttk.Frame(parent, style='Card.TFrame')
        bottom.pack(side='bottom', fill='x')
        ttk.Separator(bottom).pack(fill='x', pady=(10, 8))
        ttk.Label(bottom, textvariable=self.preview, style='Muted.TLabel', wraplength=252,
                  justify='left').pack(anchor='w', pady=(0, 8))
        row = ttk.Frame(bottom, style='Card.TFrame')
        row.pack(fill='x', pady=(0, 8))
        self._button(row, '설정 불러오기', self.load_settings, True).pack(side='left', expand=True, fill='x', padx=(0, 5))
        self._button(row, '설정 저장', self.save_settings, True).pack(side='left', expand=True, fill='x')
        self.run_button = self._button(bottom, '계산 시작', self.start_calculation, style='Primary.TButton')
        self.run_button.pack(fill='x')
        self.progress = ttk.Progressbar(bottom, mode='indeterminate')
        self.progress.pack(fill='x', pady=(8, 0))
        self.settings_tabs = ttk.Notebook(parent)
        self.settings_tabs.pack(fill='both', expand=True)
        basic = self._settings_page('기본')
        advanced = self._settings_page('고급')
        self._combo(basic, '계산 대상', self.target, [*TARGETS, '직접 지정'],
                    '모든 대상의 전자 수는 한 개입니다.', self._target_changed)
        self._combo(basic, '오비탈 종류', self.orbital,
                    [f'{name} · l={i}' for i, name in enumerate(ORBITALS)],
                    's는 1s부터, p는 2p부터 계산합니다.', self._orbital_changed)
        ttk.Label(basic, text='계산할 상태 수', style='Card.TLabel').pack(anchor='w', pady=(10, 5))
        count = ttk.Spinbox(basic, from_=1, to=6, textvariable=self.vars['num_states'], command=self._suggest)
        count.pack(fill='x')
        count.bind('<FocusOut>', self._suggest)
        count.bind('<Return>', self._suggest)
        self.controls.append((count, 'normal'))
        ttk.Label(basic, text='같은 종류에서 낮은 순서대로 1~6개.', style='Muted.TLabel').pack(anchor='w', pady=(5, 0))
        self._combo(basic, '정밀도', self.precision, [*PRESETS, '직접 설정'],
                    '빠른 확인 2천 · 권장 6천 · 정밀 1.2만 점', self._suggest)
        ttk.Label(basic, text='현재는 한 전자 모델입니다.\n분자·촉매 DFT는 다음 개발 단계예요.',
                  style='Muted.TLabel', wraplength=220, justify='left').pack(anchor='w', pady=(18, 0))
        self._entry(advanced, '핵전하 Z', 'nuclear_charge', '직접 지정: Z=1~20, 전자는 항상 한 개.', self._z_changed)
        self._entry(advanced, '계산영역 반지름 [bohr]', 'r_max_bohr', '바깥 경계입니다. 들뜬 상태는 큰 영역이 필요해요.')
        self._entry(advanced, '내부 격자점 수', 'grid_points', '100~100,000. 늘릴수록 간격이 작아집니다.')
        ttk.Label(advanced, text='직접 바꾼 설정도 모두 결과와 함께 저장됩니다.\n권장값으로 돌아가려면 기본 탭에서 정밀도를 다시 선택하세요.',
                  style='Muted.TLabel', wraplength=220, justify='left').pack(anchor='w', pady=(15, 0))

    def _combo(self, parent, label, variable, values, hint, callback):
        ttk.Label(parent, text=label, style='Card.TLabel').pack(anchor='w', pady=(10, 5))
        control = ttk.Combobox(parent, textvariable=variable, values=values, state='readonly')
        control.pack(fill='x')
        control.bind('<<ComboboxSelected>>', callback)
        self.controls.append((control, 'readonly'))
        ttk.Label(parent, text=hint, style='Muted.TLabel', wraplength=235).pack(anchor='w', pady=(5, 0))

    def _entry(self, parent, label, key, hint, callback=None):
        ttk.Label(parent, text=label, style='Card.TLabel').pack(anchor='w', pady=(10, 5))
        entry = ttk.Entry(parent, textvariable=self.vars[key])
        entry.pack(fill='x')
        entry.bind('<FocusOut>', callback or self._manual)
        entry.bind('<Return>', callback or self._manual)
        self.controls.append((entry, 'normal'))
        ttk.Label(parent, text=hint, style='Muted.TLabel', wraplength=235).pack(anchor='w', pady=(5, 4))

    def _build_results(self):
        parent = self.result_tab
        parent.columnconfigure(0, weight=1)
        parent.rowconfigure(3, weight=1)
        ttk.Label(parent, textvariable=self.changed, style='Muted.TLabel', wraplength=920).grid(row=0, column=0, sticky='w', pady=(0, 12))
        cards = ttk.Frame(parent, style='Card.TFrame')
        cards.grid(row=1, column=0, sticky='ew', pady=(0, 12))
        self.metrics = []
        for i, label in enumerate(('가장 낮은 상태의 에너지', '최대 절대오차 · 해석해 기준', '계산 + 데이터 저장 시간')):
            cards.columnconfigure(i, weight=1)
            card = tk.Frame(cards, bg='#f2f6fc', padx=14, pady=10)
            card.grid(row=0, column=i, sticky='ew', padx=(0, 8 if i < 2 else 0))
            tk.Label(card, text=label, bg='#f2f6fc', fg=MUTED, font=(self.font, 9)).pack(anchor='w')
            var = tk.StringVar(value='—')
            tk.Label(card, textvariable=var, bg='#f2f6fc', fg=INK, font=(self.font, 16, 'bold')).pack(anchor='w', pady=(5, 0))
            self.metrics.append(var)
        filters = ttk.Frame(parent, style='Card.TFrame')
        filters.grid(row=2, column=0, sticky='ew', pady=(0, 5))
        ttk.Label(filters, text='그래프', style='Card.TLabel').pack(side='left', padx=(0, 10))
        self.filter_box = ttk.Combobox(filters, textvariable=self.state_filter, values=['전체 상태'], state='readonly', width=12)
        self.filter_box.pack(side='left', padx=(0, 12))
        self.filter_box.bind('<<ComboboxSelected>>', self.redraw)
        ttk.Label(filters, text='거리 단위', style='Muted.TLabel').pack(side='left', padx=(0, 6))
        units = ttk.Combobox(filters, textvariable=self.unit, values=['bohr', 'Å'], state='readonly', width=7)
        units.pack(side='left', padx=(0, 12))
        units.bind('<<ComboboxSelected>>', self.redraw)
        ttk.Checkbutton(filters, text='해석해와 비교', variable=self.compare, command=self.redraw).pack(side='left')
        self._button(filters, '현재 그래프 저장', self.export_plot).pack(side='right')
        self.chart_frame = ttk.Frame(parent, style='Card.TFrame')
        self.chart_frame.grid(row=3, column=0, sticky='nsew')
        self.empty = ttk.Frame(self.chart_frame, style='Card.TFrame', padding=30)
        self.empty.pack(expand=True)
        ttk.Label(self.empty, text='첫 계산을 시작해 보세요', style='Title.TLabel').pack(pady=(25, 12))
        ttk.Label(self.empty, text='왼쪽에서 대상과 정밀도를 고른 뒤 계산 시작을 누르면\n에너지, 파동함수와 오차를 한 번에 확인할 수 있어요.',
                  style='Muted.TLabel', justify='center').pack()
        ttk.Label(parent, text='P(r)는 구껍질 확률을 나타내는 방사확률밀도입니다. 원시 데이터는 bohr 단위로 저장됩니다.',
                  style='Muted.TLabel', wraplength=920).grid(row=4, column=0, sticky='w', pady=(7, 8))
        table_frame = ttk.Frame(parent, style='Card.TFrame')
        table_frame.grid(row=5, column=0, sticky='ew')
        columns = ('state', 'energy', 'exact', 'error', 'norm', 'radius')
        self.table = ttk.Treeview(table_frame, columns=columns, show='headings', height=3)
        names = ('상태', '에너지 [eV]', '해석해 [eV]', '절대오차 [eV]', '정규화', '평균 반지름 [bohr]')
        widths = (45, 125, 125, 120, 90, 140)
        for col, name, width in zip(columns, names, widths):
            self.table.heading(col, text=name)
            self.table.column(col, width=width, minwidth=50, anchor='center')
        self.table.pack(side='left', fill='x', expand=True)
        table_scroll = ttk.Scrollbar(table_frame, orient='vertical', command=self.table.yview)
        table_scroll.pack(side='right', fill='y')
        self.table.configure(yscrollcommand=table_scroll.set)
        self.quality = ttk.Label(parent, text='계산이 완료되면 정확도 확인 사항이 표시됩니다.', style='Muted.TLabel', wraplength=780)
        self.quality.grid(row=6, column=0, sticky='w', pady=(10, 8))
        actions = ttk.Frame(parent, style='Card.TFrame')
        actions.grid(row=7, column=0, sticky='ew')
        self._button(actions, '결과 폴더 열기', self.open_current).pack(side='left', padx=(0, 8))
        self._button(actions, '결과 ZIP 내보내기', self.export_results).pack(side='left')
        self._button(actions, '처음 사용하는 방법', lambda: self.tabs.select(self.help_tab)).pack(side='right')
        parent.bind('<Configure>', self._resize_labels)

    def _resize_labels(self, event):
        width = max(300, event.width - 32)
        for child in self.result_tab.winfo_children():
            if isinstance(child, ttk.Label):
                child.configure(wraplength=width)

    def _text(self, parent, content):
        widget = scrolledtext.ScrolledText(parent, wrap='word', relief='flat', bg=WHITE, fg=INK,
                                         font=(self.font, 10), padx=12, pady=10, spacing1=3, spacing3=5)
        widget.pack(fill='both', expand=True)
        self._set_text(widget, content)
        return widget

    @staticmethod
    def _set_text(widget, content):
        widget.configure(state='normal')
        widget.delete('1.0', 'end')
        widget.insert('1.0', content)
        widget.configure(state='disabled')

    def _build_history(self):
        ttk.Label(self.history_tab, text='실행할 때마다 조건과 결과를 따로 보관합니다.', style='Section.TLabel').pack(anchor='w', pady=(0, 12))
        self.history_tree = ttk.Treeview(self.history_tab, columns=('time', 'system', 'grid', 'energy'), show='headings', height=8)
        for key, title, width in (('time', '실행 시각', 190), ('system', '조건', 170), ('grid', '격자점', 100), ('energy', '첫 상태 에너지 [eV]', 180)):
            self.history_tree.heading(key, text=title)
            self.history_tree.column(key, width=width, anchor='center')
        self.history_tree.pack(fill='x')
        self.history_tree.bind('<<TreeviewSelect>>', self._history_selected)
        row = ttk.Frame(self.history_tab, style='Card.TFrame')
        row.pack(fill='x', pady=12)
        self._button(row, '선택한 조건 불러오기', self.load_history_settings, True).pack(side='left', padx=(0, 8))
        self._button(row, '선택한 결과 폴더 열기', self.open_history).pack(side='left', padx=(0, 8))
        self._button(row, '새로고침', self.refresh_history).pack(side='left')
        self.history_text = self._text(self.history_tab, '실행 기록을 선택하면 전체 설정과 수치 결과를 볼 수 있습니다.')

    def _settings_changed(self, *_):
        if self.applying:
            return
        try:
            cfg = self.read_config()
            labels = ', '.join(f'{n}{ORBITALS[cfg.angular_momentum_l]}' for n in range(cfg.angular_momentum_l + 1, cfg.angular_momentum_l + cfg.num_states + 1))
            self.preview.set(f'전자 1개 · Z={cfg.nuclear_charge}\n{labels}\n격자 간격 {cfg.r_max_bohr/(cfg.grid_points+1):.6g} bohr')
            dirty = self.result is not None and cfg != self.result.config
        except ValueError:
            self.preview.set('설정값을 입력한 뒤 계산 시작을 눌러 주세요.')
            dirty = self.result is not None
        if self.result is not None:
            cfg = self.result.config
            prefix = '설정 변경됨 · 다시 계산해야 반영됩니다.  ' if dirty else '마지막 실행 결과  ·  '
            self.changed.set(prefix + f'아래 결과: Z={cfg.nuclear_charge}, l={cfg.angular_momentum_l}, '
                             f'N={cfg.grid_points}, 반지름={cfg.r_max_bohr:g} bohr')

    def read_config(self):
        return config_from_fields({key: var.get() for key, var in self.vars.items()})

    def _target_changed(self, *_):
        if self.target.get() in TARGETS:
            self.vars['nuclear_charge'].set(str(TARGETS[self.target.get()]))
            self._suggest()
        else:
            self.settings_tabs.select(1)

    def _orbital_changed(self, *_):
        self.vars['angular_momentum_l'].set(str(ORBITALS.index(self.orbital.get().split()[0])))
        self._suggest()

    def _z_changed(self, *_):
        z = self.vars['nuclear_charge'].get().strip()
        self.target.set(next((label for label, charge in TARGETS.items() if str(charge) == z), '직접 지정'))
        self._suggest()

    def _manual(self, *_):
        self.precision.set('직접 설정')
        self._settings_changed()

    def _suggest(self, *_):
        if self.precision.get() == '직접 설정':
            self.settings_tabs.select(1)
            self._settings_changed()
            return
        try:
            cfg = recommended_config(int(self.vars['nuclear_charge'].get()),
                                     int(self.vars['angular_momentum_l'].get()),
                                     int(self.vars['num_states'].get()), self.precision.get())
        except (ValueError, KeyError, ZeroDivisionError, OverflowError):
            return
        self.vars['r_max_bohr'].set(f'{cfg.r_max_bohr:.10g}')
        self.vars['grid_points'].set(str(cfg.grid_points))
        self._settings_changed()

    def apply_config(self, cfg, initial=False):
        self.applying = True
        for key, value in asdict(cfg).items():
            self.vars[key].set(str(value))
        self.target.set(next((label for label, z in TARGETS.items() if z == cfg.nuclear_charge), '직접 지정'))
        self.orbital.set(f'{ORBITALS[cfg.angular_momentum_l]} · l={cfg.angular_momentum_l}')
        if initial and cfg == recommended_config(cfg.nuclear_charge, cfg.angular_momentum_l, cfg.num_states):
            self.precision.set('권장')
        else:
            self.precision.set('직접 설정')
            self.settings_tabs.select(1)
        self.applying = False
        self._settings_changed()

    def load_settings(self):
        name = filedialog.askopenfilename(title='계산 설정 불러오기', filetypes=[('JSON 설정', '*.json')])
        if name:
            try:
                self.apply_config(load_config(Path(name)))
                self.status.set('설정을 불러왔습니다. 계산 시작을 눌러 실행하세요.')
            except (ValueError, TypeError, OSError) as exc:
                messagebox.showerror('설정을 읽을 수 없습니다', str(exc), parent=self.root)

    def save_settings(self):
        try:
            cfg = self.read_config()
        except ValueError as exc:
            messagebox.showerror('설정 확인', str(exc), parent=self.root)
            return
        name = filedialog.asksaveasfilename(title='계산 설정 저장', initialfile='settings.json', defaultextension='.json', filetypes=[('JSON 설정', '*.json')])
        if name:
            try:
                save_config(cfg, Path(name))
                self.status.set('설정 파일을 저장했습니다.')
            except OSError as exc:
                messagebox.showerror('저장 실패', str(exc), parent=self.root)

    def _set_busy(self, busy):
        self.busy = busy
        self.run_button.configure(state='disabled' if busy else 'normal', text='계산 중…' if busy else '계산 시작')
        for widget, normal in self.controls:
            widget.configure(state='disabled' if busy else normal)
        if busy:
            self.progress.start(12)
        else:
            self.progress.stop()

    def _log(self, message):
        self.log_lines.append(f'[{datetime.now().strftime("%H:%M:%S")}] {message}')
        if hasattr(self, 'log_text'):
            self._set_text(self.log_text, '\n'.join(self.log_lines))
            self.log_text.see('end')

    def start_calculation(self):
        if self.busy:
            return
        try:
            config = self.read_config()
        except ValueError as exc:
            self.status.set(str(exc))
            messagebox.showerror('계산 설정을 확인해 주세요', str(exc), parent=self.root)
            return
        self._set_busy(True)
        self.tabs.select(self.result_tab)
        self.status.set('계산 중 · 해밀토니안을 구성하고 에너지 상태를 구합니다.')
        self.log_lines = []
        self._log('계산 시작 · ' + json.dumps(asdict(config), ensure_ascii=False))
        # 워커는 Tk나 Matplotlib에 접근하지 않는다. 결과만 큐로 주 스레드에 전달한다.
        threading.Thread(target=self._worker, args=(config,), daemon=True).start()

    def _worker(self, config):
        directory = None
        try:
            started = time.perf_counter()
            result = solve(config)
            self.messages.put(('stage', '고유상태 계산 완료 · 원시 데이터와 설정을 저장합니다.'))
            directory = new_run_dir(self.output_root)
            summary = save_result(result, directory, time.perf_counter() - started)
            elapsed = time.perf_counter() - started
            # 표기 시간에는 수치 계산과 원시 데이터 저장을 포함한다.
            summary['elapsed_seconds'] = elapsed
            (directory / 'results.json').write_text(json.dumps(summary, ensure_ascii=False, indent=2), encoding='utf-8')
            self.messages.put(('done', result, summary, directory))
        except Exception as exc:
            detail = traceback.format_exc()
            if directory is not None:
                try:
                    (directory / 'error.log').write_text(detail, encoding='utf-8')
                except OSError:
                    pass
            self.messages.put(('error', str(exc), detail))

    def _poll(self):
        try:
            while True:
                event = self.messages.get_nowait()
                if event[0] == 'stage':
                    self.status.set(event[1])
                    self._log(event[1])
                elif event[0] == 'done':
                    self._complete(*event[1:])
                elif event[0] == 'error':
                    self._set_busy(False)
                    self.status.set('계산을 완료하지 못했습니다. 설정 또는 실행 로그를 확인하세요.')
                    self._log(event[2])
                    if self.closing:
                        self._destroy()
                        return
                    messagebox.showerror('계산을 완료하지 못했습니다',
                                         event[1] + '\n\n권장 설정으로 다시 시도하거나 실행 로그에서 원인을 확인해 주세요.', parent=self.root)
                if self.closed:
                    return
        except queue.Empty:
            pass
        if not self.closed:
            self.root.after(80, self._poll)

    def _complete(self, result, summary, directory):
        self.result, self.summary, self.current_dir = result, summary, directory
        self.filter_box.configure(values=['전체 상태'] + [state['label'] for state in summary['states']])
        self.state_filter.set('전체 상태')
        self.unit.set('bohr')
        self.compare.set(True)
        self.metrics[0].set(f"{summary['states'][0]['energy_ev']:.6f} eV")
        self.metrics[1].set(f"{max(state['abs_error_ev'] for state in summary['states']):.3g} eV")
        self.metrics[2].set(f"{summary['elapsed_seconds']:.2f} 초")
        self.table.delete(*self.table.get_children())
        for state in summary['states']:
            self.table.insert('', 'end', values=(state['label'], f"{state['energy_ev']:.8f}",
                f"{state['exact_energy_hartree']*HARTREE_TO_EV:.8f}", f"{state['abs_error_ev']:.6g}",
                f"{state['radial_normalization']:.8f}", f"{state['mean_r_bohr']:.6f}"))
        if summary['warnings']:
            self.quality.configure(text='확인 필요 · ' + ' '.join(summary['warnings']), foreground='#b45309')
        else:
            self.quality.configure(text='기본 격자·영역 점검 통과 · 정밀도를 바꿔 수렴을 확인하면 더 신뢰할 수 있어요.', foreground='#0f766e')
        self._set_text(self.detail_text, THEORY + '\n\n마지막 실행의 전체 설정과 결과\n\n' + json.dumps(summary, ensure_ascii=False, indent=2))
        try:
            self.redraw()
            # 고해상도 파일 저장이 내장 캔버스의 DPI/레이아웃을 바꾸지 않도록 분리한다.
            create_figure(result).savefig(directory / 'wavefunctions.png', dpi=170)
            self._log('계산 완료 · 그래프, 원시 데이터, 설정 저장 완료')
            self.status.set(f"완료 · {len(summary['states'])}개 상태 계산 / 실행 기록에 저장했습니다.")
        except Exception:
            self._log(traceback.format_exc())
            self.status.set('수치 데이터는 저장했지만 그래프를 표시하지 못했습니다. 실행 로그를 확인하세요.')
        try:
            (directory / 'execution.log').write_text('\n'.join(self.log_lines), encoding='utf-8')
        except OSError as exc:
            self._log('실행 로그 파일 저장 실패: ' + str(exc))
        self._set_busy(False)
        self._settings_changed()
        self.refresh_history()
        if self.closing:
            self._destroy()

    def redraw(self, *_):
        if self.result is None:
            return
        labels = [state['label'] for state in self.summary['states']]
        selected = None if self.state_filter.get() == '전체 상태' else labels.index(self.state_filter.get())
        figure = create_figure(self.result, selected, 'angstrom' if self.unit.get() == 'Å' else 'bohr', self.compare.get())
        if self.toolbar is not None:
            self.toolbar.destroy()
        if self.canvas is not None:
            pending = getattr(self.canvas, '_idle_draw_id', None)
            if pending is not None:
                self.canvas.get_tk_widget().after_cancel(pending)
            self.canvas.get_tk_widget().destroy()
        self.empty.pack_forget()
        self.figure = figure
        figure.set_dpi(float(self.root.winfo_fpixels('1i')))
        self.canvas = AppCanvas(figure, master=self.chart_frame)
        self.toolbar = NavigationToolbar2Tk(self.canvas, self.chart_frame, pack_toolbar=False)
        self.toolbar.update()
        self.toolbar.pack(side='bottom', fill='x')
        self.canvas.get_tk_widget().pack(side='top', fill='both', expand=True)
        self.canvas.draw_idle()

    def export_plot(self):
        if self.figure is None:
            messagebox.showinfo('계산 후 사용할 수 있어요', '먼저 계산 시작을 눌러 주세요.', parent=self.root)
            return
        name = filedialog.asksaveasfilename(title='현재 선택한 그래프 저장', initialfile='wavefunctions_view.png', defaultextension='.png', filetypes=[('PNG 그림', '*.png')])
        if name:
            try:
                labels = [state['label'] for state in self.summary['states']]
                selected = None if self.state_filter.get() == '전체 상태' else labels.index(self.state_filter.get())
                exported = create_figure(self.result, selected,
                                         'angstrom' if self.unit.get() == 'Å' else 'bohr', self.compare.get())
                for source, target in zip(self.figure.axes, exported.axes):
                    target.set_xlim(source.get_xlim())
                    target.set_ylim(source.get_ylim())
                exported.savefig(name, dpi=200)
                self.status.set('현재 그래프를 저장했습니다.')
            except OSError as exc:
                messagebox.showerror('저장 실패', str(exc), parent=self.root)

    def export_results(self):
        if self.current_dir is None:
            messagebox.showinfo('계산 후 사용할 수 있어요', '먼저 계산 시작을 눌러 주세요.', parent=self.root)
            return
        name = filedialog.asksaveasfilename(title='마지막 실행 결과 ZIP 내보내기', initialfile=self.current_dir.name + '.zip', defaultextension='.zip', filetypes=[('ZIP 압축', '*.zip')])
        if name:
            try:
                target = Path(name).resolve()
                files = [p for p in self.current_dir.iterdir() if p.is_file() and p.resolve() != target]
                with zipfile.ZipFile(target, 'w', zipfile.ZIP_DEFLATED) as archive:
                    for path in files:
                        archive.write(path, path.name)
                self.status.set('설정, 데이터, 기본 그래프를 ZIP으로 내보냈습니다.')
            except OSError as exc:
                messagebox.showerror('내보내기 실패', str(exc), parent=self.root)

    def open_current(self):
        if self.current_dir is None:
            messagebox.showinfo('계산 후 사용할 수 있어요', '먼저 계산 시작을 눌러 주세요.', parent=self.root)
            return
        self._open(self.current_dir)

    def _open(self, path):
        try:
            open_directory(path)
        except OSError as exc:
            messagebox.showerror('폴더를 열 수 없습니다', str(exc), parent=self.root)

    def refresh_history(self):
        self.history = list_runs(self.output_root)
        self.history_tree.delete(*self.history_tree.get_children())
        for i, entry in enumerate(self.history):
            summary = entry['summary']
            cfg = summary['config']
            stamp = summary.get('created_at_utc', '')
            try:
                stamp = datetime.fromisoformat(stamp).astimezone().strftime('%Y-%m-%d %H:%M:%S')
            except ValueError:
                pass
            self.history_tree.insert('', 'end', iid=str(i), values=(stamp, f"Z={cfg['nuclear_charge']} / l={cfg['angular_momentum_l']}",
                cfg['grid_points'], f"{summary['states'][0]['energy_ev']:.8f}"))

    def _selected_history(self):
        selected = self.history_tree.selection()
        return self.history[int(selected[0])] if selected else None

    def _history_selected(self, *_):
        entry = self._selected_history()
        if entry:
            self._set_text(self.history_text, json.dumps(entry['summary'], ensure_ascii=False, indent=2))

    def load_history_settings(self):
        entry = self._selected_history()
        if entry:
            self.apply_config(Config(**entry['summary']['config']))
            self.status.set('과거 조건을 불러왔습니다. 계산 시작을 눌러 재실행하세요.')
            self.tabs.select(self.result_tab)
        else:
            messagebox.showinfo('실행 기록 선택', '불러올 실행 기록을 먼저 선택하세요.', parent=self.root)

    def open_history(self):
        entry = self._selected_history()
        if entry:
            self._open(entry['directory'])
        else:
            messagebox.showinfo('실행 기록 선택', '열 실행 기록을 먼저 선택하세요.', parent=self.root)

    def _destroy(self):
        self.closed = True
        # Tk/Matplotlib의 예약된 그리기와 폴링을 해제한 뒤 창을 닫는다.
        for job in self.root.tk.call('after', 'info'):
            try:
                # Tcl 명령의 소유 위젯은 destroy가 정리하도록 예약만 취소한다.
                self.root.tk.call('after', 'cancel', job)
            except tk.TclError:
                pass
        self.root.destroy()

    def close(self):
        if self.busy:
            self.closing = True
            self.status.set('실행 중인 계산과 저장을 마무리한 뒤 창을 닫습니다.')
        else:
            self._destroy()


def main():
    try:
        if sys.platform == 'win32':
            try:
                import ctypes
                ctypes.windll.shcore.SetProcessDpiAwareness(1)
            except (AttributeError, OSError):
                pass
        root = tk.Tk()
    except tk.TclError as exc:
        print('사용 화면을 열 수 없습니다. 그래픽 환경과 Python의 Tkinter 설치를 확인하세요.\n'
              '창 없이 계산하려면 python main.py --no-show를 실행하세요.\n' + str(exc), file=sys.stderr)
        return 1
    if '--hydrogen-check' in sys.argv:
        QuantumApp(root)
    elif '--periodic' in sys.argv:
        from periodic_app import PeriodicApp
        PeriodicApp(root)
    else:
        from molecular_app import MolecularApp
        MolecularApp(root)
    root.mainloop()
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
