"""General molecular DFT desktop; legacy QuantumApp supplies shared UI utilities."""
from dataclasses import fields
import json
from pathlib import Path
import queue
import threading
import time
import traceback
import tkinter as tk
from tkinter import filedialog, messagebox, scrolledtext, ttk

from matplotlib.backends.backend_tkagg import NavigationToolbar2Tk
from matplotlib.figure import Figure

from app import AppCanvas, QuantumApp, ROOT, WHITE, INK, MUTED
from dft_workflow import CalculationError, RunControl, list_dft_runs, run_calculation
from element_catalog import PERIODIC_CAPABILITIES
from molecular import (BASES, CALCULATIONS, DISPERSIONS, ELEMENTS, OPT_CONVERGENCE, FUNCTIONALS, DEF2_CORE_ELECTRONS, Atom, DFTConfig, available_bases, coordinates_text, load_dft_config,
                       describe_bases, format_basis_overrides, parse_basis_overrides, parse_coordinates, read_xyz,
                       save_dft_config)
from molecular_plotting import create_molecular_figure
from geometry_tools import PRESETS, apply_preset, measure, BOHR_ANGSTROM

HELP = '''분자 DFT 사용 안내

엔진 선택: native는 QuantumStarter 자체 밀도/Fock/DIIS/SCF/에너지 계산,
pyscf는 외부 PySCF SCF입니다. native도 기저·적분·격자에는 PySCF, XC 함수에는 Libxc를 사용합니다.
자체 엔진은 H~Kr, 비상대론 전체전자 점핵, RKS/UKS, 지정 LDA/GGA/global hybrid, 최대160 AO입니다.
에너지 성분과 출처를 결과에 저장합니다. 기존118원소/ECP/X2C는 PySCF 엔진에서 제공합니다.
중단은 안전한 SCF/적분 단계에서 반영되며 즉시 프로세스를 강제 종료하지 않습니다.
체크포인트는 밀도만 저장하므로 재시작 시 DIIS 이력은 새로 구성합니다.

1. 원소와 x y z 좌표를 입력하거나 XYZ 파일을 불러옵니다.
   XYZ는 한 구조만 읽으며 좌표 단위는 Å입니다. 파일의 주석은 실행하지 않습니다.
2. 전체 전하와 스핀 2S = Nα − Nβ를 지정합니다.
   전하 0은 중성, +1은 전자 한 개 제거, −1은 전자 한 개 추가입니다.
   2S=0은 singlet, 1은 doublet, 2는 triplet입니다. 다중도는 2S+1입니다.
   전자 수와 2S의 홀짝이 같아야 합니다. 열린껍질 singlet은 UKS를 직접 선택할 수 있으나,
   깨진 대칭 초기 추정과 상태 탐색은 현재 자동 제공하지 않습니다.
3. 계산 방법 탭에서 기저함수와 교환·상관 근사를 선택합니다.
   ECP auto는 무거운 원소의 def2 기저에 대응하는 def2-ECP를 적용합니다.
   코어 전자 수는 원소에 따라 28/46/60개이며, 전체 전자와 명시적 계산 전자를 구분합니다.
   Dyall 기저는 H~Og 전체전자 기저이며 auto에서 스칼라 X2C와 Gaussian 유한 핵을 사용합니다.
   스핀-궤도 결합(SOC)은 포함하지 않습니다. ECP와 X2C를 동시에 적용할 수 없습니다.
   auto는 2S=0에서 RKS, 그 외에는 UKS를 사용합니다.
   SCF 반복 방식은 기본 diis 또는 2차 newton을 선택할 수 있습니다.
   diis가 수렴하지 않을 때 newton으로 다시 시도할 수 있으며 자동으로 바꾸지는 않습니다.
   newton의 반복 수는 바깥 반복(macroiteration) 수입니다. 내부 반복까지 합친 횟수는 아닙니다.
4. 계산 시작을 누르면 SCF 반복 상태가 표시됩니다.
5. 총에너지, 수렴 상태, 오비탈 점유/에너지와 밀도 단면을 확인합니다.

현재 지원 범위

H~Og (원자번호 1~118), 고립 원자·분자·이온, 고정 핵에서 단일점 KS-DFT 에너지.
Dyall-v2z/v3z로 모든 원소를 입력할 수 있으며, 원소별 지원표에서 다른 기저도 확인할 수 있습니다.
기저함수 800개 이하 / 원자 100개 이하. 크기가 커지면 시간과 메모리가 급격히 증가합니다.
이 분자 화면은 임의 ECP, 스핀-궤도 결합, 구조 최적화, 들뜬 상태를 지원하지 않습니다.
고체·촉매 표면은 상단 '고체·표면' 버튼의 별도 주기계 엔진을 사용하세요. Ag·Au·Pd·Pt도 입력할 수 있습니다.
6-31G는 H~Zn, 6-31G*는 H~Kr, 6-31+G*는 H~Ca, cc-pVDZ는 K를 제외한 H~Kr를 지원합니다.
STO-3G는 H~I, def2-SVP/TZVP는 Ce~Lu를 제외한 H~Rn의 72개 원소를 지원합니다.
Dyall의 auto 격자는 Gauss-Chebyshev 방사 격자와 각도 축소 없는 Becke 분할을 사용합니다.
Z>36은 해당 수준 기본 방사 점 수의 3배를 쓰며 실제 점 수를 결과 JSON에 기록합니다.
격자 수준을 높여 에너지를 비교하세요. 악티늄족·초중원소는 특히 상대론·전자상태 검증이 필요합니다.
전이금속은 여러 스핀 상태·전자배치가 가능하므로 지정한 2S의 수렴이 바닥상태의 보증은 아닙니다.
범용 입력은 모든 물질에 대한 정확도 보증을 뜻하지 않습니다.

결과 해석

총에너지는 전자 에너지와 핵 간 반발을 포함합니다. 오비탈 에너지의 합과 다릅니다.
오비탈 에너지를 실험적 여기 에너지나 밴드갭과 동일시하지 마세요.
수렴은 해당 SCF 조건을 만족한다는 뜻이며 기저/격자 수렴이나 올바른 전자상태를 보증하지 않습니다.
미수렴 결과는 마지막 반복의 진단 데이터이며 주황색으로 표시합니다.
밀도 그림은 핵 좌표 중심에 가까운 격자면의 단면입니다. 원자 전체가 그 평면 위에 있지는 않습니다.
그림은 보간한 시각화이며 SCF 적분 격자와 다릅니다. 전자 수는 Tr(DS)로 확인합니다.
밀도 단위는 electron/bohr³, 오비탈 진폭은 bohr^(-3/2), 그림의 좌표는 Å입니다.
ECP 계산의 그림·CUBE에는 코어 전자밀도가 없으며 명시적 전자의 의사밀도를 표시합니다.
X2C 계산의 그림·CUBE는 X2C 표현이며 물리적 사성분 밀도의 picture-change 보정은 없습니다.
ECP와 전체전자 모델의 절대 총에너지는 직접 비교하지 마세요.

저장과 재현

실행마다 outputs/molecular 아래 새 폴더를 만듭니다.
settings.json: 원소·좌표·전하·스핀·계산 방법과 수치 설정.
results.json: 수렴 여부, 총에너지, 오비탈, 경고, 버전과 단위.
density.cube / spin_density.cube: 총 전자밀도 / UKS 스핀밀도.
*_homo.cube / *_lumo.cube: 존재하는 각 스핀 채널의 HOMO/LUMO.
electronic_structure.npz: 모든 MO 계수·에너지·점유수, AO 라벨, 밀도/겹침 행렬.
density_grid.npz: 밀도 그림 원시 격자. orbitals.csv / scf.csv: 오비탈/SCF 표.
execution.log: 엔진 로그. 실패한 실행에는 error.json과 설정을 남깁니다.
실행 기록에서 과거 설정을 불러와 재계산할 수 있습니다.

Windows 실행

화면은 Windows Python, 계산은 프로젝트 .venv-dft 안의 WSL PySCF를 사용합니다.
setup_dft_windows.bat로 엔진을 설치하고 start_windows.bat로 화면을 여세요.
Linux/macOS에서는 requirements-dft.txt 설치 후 같은 화면/명령행을 실행합니다.
수소 한 전자 검증 화면: python app.py --hydrogen-check
'''

EXAMPLES = {
    '자체 DFT · 물 (LDA 검증)': ('O 0 0 0\nH 0 -0.757 0.587\nH 0 0.757 0.587', 0, 0,
                                {'engine': 'native', 'basis': 'sto-3g', 'xc': 'lda,vwn'}),
    '물 · H2O (예제)': ('O 0 0 0\nH 0 -0.757 0.587\nH 0 0.757 0.587', 0, 0),
    '수소 분자 · H2 (검증)': ('H 0 0 0\nH 0 0 0.74', 0, 0),
    '헬륨 · He (검증)': ('He 0 0 0', 0, 0),
    '수소 원자 · H (검증)': ('H 0 0 0', 0, 1),
    'OH 라디칼 (예제)': ('O 0 0 0\nH 0 0 0.97', 0, 1),
    '암모늄 · NH4+ (예제)': ('N 0 0 0\nH 0.59 0.59 0.59\nH -0.59 -0.59 0.59\nH -0.59 0.59 -0.59\nH 0.59 -0.59 -0.59', 1, 0),
    '염화칼륨 · KCl (예제)': ('K 0 0 0\nCl 0 0 2.67', 0, 0),
    '철 이온 · Fe3+ (2S=5 예제)': ('Fe 0 0 0', 3, 5),
    '코발트 이온 · Co2+ (2S=3 예제)': ('Co 0 0 0', 2, 3, {'scf_solver': 'newton'}),
    '구리 이온 · Cu+ (예제)': ('Cu 0 0 0', 1, 0),
    '크립톤 · Kr (예제)': ('Kr 0 0 0', 0, 0),
    '염화루비듐 · RbCl (ECP 예제)': ('Rb 0 0 0\nCl 0 0 2.8', 0, 0),
    '지르코늄 이온 · Zr4+ (ECP 예제)': ('Zr 0 0 0', 4, 0),
    '금 이온 · Au+ (ECP 예제)': ('Au 0 0 0', 1, 0, {'scf_solver': 'newton'}),
    '오가네손 · Og (118번 · X2C)': ('Og 0 0 0', 0, 0, {'basis': 'dyall-v2z', 'scf_solver': 'newton'}),
}


class MolecularApp(QuantumApp):
    def __init__(self, root, output_root=None):
        self.root = root
        self.output_root = Path(output_root) if output_root else ROOT / 'outputs' / 'molecular'
        self.busy = self.closing = self.closed = False
        self.applying = True
        self.messages = queue.Queue()
        self.result = self.summary = self.current_dir = None
        self.figure = self.canvas = self.toolbar = None
        self.controls, self.history, self.log_lines = [], [], []
        self.vars = {f.name: tk.StringVar(value=str(getattr(DFTConfig(), f.name)))
                     for f in fields(DFTConfig) if f.name != 'atoms'}
        self.backend = tk.StringVar(value='auto')
        self.example = tk.StringVar(value='물 · H2O (예제)')
        self.status = tk.StringVar(value='준비됨 · 원소와 좌표를 입력하고 전하·스핀을 확인하세요.')
        self.preview = tk.StringVar()
        self.ecp_hint = tk.StringVar()
        self.changed = tk.StringVar(value='아직 계산한 결과가 없습니다.')
        self.field = tk.StringVar(value='density.cube')
        self.plane = tk.StringVar(value='yz')
        self._setup_styles()
        root.title('QuantumStarter 0.8 | 자체 / 외부 분자 DFT')
        root.geometry('1440x960')
        root.minsize(1220, 820)
        self._build()
        for variable in self.vars.values():
            variable.trace_add('write', self._settings_changed)
        self.apply_config(DFTConfig())
        if (ROOT / 'molecule.json').exists():
            try:
                self.apply_config(load_dft_config(ROOT / 'molecule.json'))
            except (OSError, ValueError) as exc:
                self._log(f'기본 분자 설정 오류: {exc}')
        self.refresh_history()
        root.protocol('WM_DELETE_WINDOW', self.close)
        root.bind('<F5>', lambda event: self.start_calculation())
        root.after(80, self._poll)

    def _build(self):
        header = ttk.Frame(self.root, style='Card.TFrame', padding=(24, 16))
        header.pack(fill='x')
        ttk.Label(header, text='QuantumStarter', style='Title.TLabel').pack(side='left')
        ttk.Label(header, text='  /  범용 분자 DFT', style='Section.TLabel').pack(side='left')
        ttk.Label(header, text='H~Og · 118개 원소', style='Muted.TLabel').pack(side='right')
        ttk.Button(header, text='원소별 지원표', command=self.show_elements).pack(side='right', padx=6)
        ttk.Button(header, text='고체·표면 / GPU 탐색', command=self.open_periodic).pack(side='right', padx=14)
        footer = ttk.Frame(self.root, padding=(20, 8))
        footer.pack(side='bottom', fill='x')
        ttk.Label(footer, textvariable=self.status, wraplength=1300).pack(anchor='w')
        body = ttk.Frame(self.root, padding=16)
        body.pack(fill='both', expand=True)
        sidebar = ttk.Frame(body, style='Card.TFrame', padding=14, width=370)
        sidebar.pack(side='left', fill='y', padx=(0, 14))
        sidebar.pack_propagate(False)
        ttk.Label(sidebar, text='01  구조와 계산 조건', style='Section.TLabel').pack(anchor='w', pady=(0, 12))
        bottom = ttk.Frame(sidebar, style='Card.TFrame')
        bottom.pack(side='bottom', fill='x')
        state_row = ttk.Frame(bottom, style='Card.TFrame')
        state_row.pack(fill='x', pady=(10, 0))
        for i, (key, label) in enumerate((('unit', '좌표 단위'), ('charge', '전체 전하'), ('spin', '스핀 2S'))):
            state_row.columnconfigure(i, weight=1)
            ttk.Label(state_row, text=label, style='Card.TLabel').grid(row=0, column=i, sticky='w', pady=(0, 4))
            if key == 'unit':
                widget = ttk.Combobox(state_row, textvariable=self.vars[key], values=['angstrom', 'bohr'], state='readonly', width=9)
                normal = 'readonly'
            else:
                widget = ttk.Entry(state_row, textvariable=self.vars[key], width=6)
                normal = 'normal'
            widget.grid(row=1, column=i, sticky='ew', padx=(0, 5))
            self.controls.append((widget, normal))
        ttk.Label(bottom, textvariable=self.preview, style='Muted.TLabel', wraplength=330).pack(fill='x', pady=10)
        row = ttk.Frame(bottom, style='Card.TFrame')
        row.pack(fill='x', pady=(0, 8))
        for label, command in (('설정 불러오기', self.load_settings), ('설정 저장', self.save_settings)):
            self._button(row, label, command, True).pack(side='left', expand=True, fill='x', padx=2)
        self.run_button = self._button(bottom, '계산 시작', self.start_calculation, style='Primary.TButton')
        self.run_button.pack(fill='x')
        self.cancel_button = ttk.Button(bottom, text='계산 중단 · 체크포인트 보존', command=self.cancel_calculation, state='disabled')
        self.cancel_button.pack(fill='x', pady=(5, 0))
        self.progress = ttk.Progressbar(bottom, mode='indeterminate')
        self.progress.pack(fill='x', pady=(8, 0))
        self.settings_tabs = ttk.Notebook(sidebar)
        self.settings_tabs.pack(fill='both', expand=True)
        structure = self._settings_page('구조')
        method = self._settings_page('계산 방법')
        advanced = self._settings_page('수치 설정')
        self._combo(structure, '시작 예제 (좌표는 자유롭게 수정)', self.example, list(EXAMPLES),
                    '예제 이름은 계산 대상의 제한이 아닙니다.', self.use_example)
        self._button(structure, 'XYZ 파일 불러오기', self.load_xyz, True).pack(fill='x', pady=10)
        self._button(structure, '구조 미리보기 · 거리/각도', self.preview_geometry, True).pack(fill='x', pady=(0, 7))
        ttk.Label(structure, text='원소   x   y   z  ·  한 줄에 원자 하나', style='Card.TLabel', wraplength=280).pack(anchor='w')
        self.coordinates = scrolledtext.ScrolledText(structure, height=6, width=26, wrap='none',
                                                    font=('Consolas', 11), relief='solid', borderwidth=1)
        self.coordinates.pack(fill='x', pady=6)
        self.controls.append((self.coordinates, 'normal'))
        self.coordinates.bind('<<Modified>>', self._coordinates_changed)
        ttk.Label(structure, text='전하 0: 중성 · 스핀 0/1/2: singlet/doublet/triplet\nXYZ 단위는 Å. 전하와 스핀은 직접 확인하세요.',
                  style='Muted.TLabel', wraplength=280).pack(anchor='w', pady=5)
        self._combo(method, '계산 엔진', self.vars['engine'], ['pyscf', 'native'],
                    'pyscf: 외부 SCF / 118원소. native: 자체 SCF / H~Kr·전체전자·점핵.', self._settings_changed)
        self._combo(method, '기저함수', self.vars['basis'], BASES, '크기가 다른 기저로 수렴을 확인하세요.', self._settings_changed)
        self.basis_box = self.controls[-1][0]
        self._entry(method, '원소별 기저 (선택)', 'basis_overrides',
                    '예: O=def2-tzvp, H=sto-3g. 비우면 모든 원소에 위 기저. '
                    '반응·흡착에너지의 모든 항에서 원소별 기저를 같게 유지하세요.', self._settings_changed)
        self._combo(method, '유효 코어 퍼텐셜 (ECP)', self.vars['ecp'], ['auto', 'none', 'def2'],
                    'auto: 무거운 원소의 def2 기저에 대응 ECP 적용.', self._settings_changed)
        ttk.Label(method, textvariable=self.ecp_hint, style='Muted.TLabel', wraplength=280).pack(anchor='w', pady=5)
        self._combo(method, '상대론 보정', self.vars['relativity'], ['auto', 'none', 'x2c'],
                    'auto: Dyall 기저에 스칼라 X2C 적용. SOC 제외.', self._settings_changed)
        self._combo(method, '핵 모형', self.vars['nuclear_model'], ['auto', 'point', 'gaussian'],
                    'auto: Dyall은 Gaussian 유한 핵, 그 외 점핵.', self._settings_changed)
        self._combo(method, '교환·상관 함수', self.vars['xc'], FUNCTIONALS, 'LDA / PBE / B3LYP / PBE0', self._settings_changed)
        self._combo(method, '계산 종류', self.vars['calculation'], list(CALCULATIONS),
                    'gradient: 에너지와 원자 힘. optimize: geomeTRIC 구조 최적화(격자 수준 3 이상 권장).',
                    self._settings_changed)
        self._combo(method, '분산 보정', self.vars['dispersion'], list(DISPERSIONS),
                    'd4: Grimme DFT-D4 (dftd4 라이브러리). PBE/PBE0/B3LYP만, LDA 불가. '
                    'SCF 후 총에너지에 더하는 기하 보정입니다.', self._settings_changed)
        self._combo(method, 'Kohn–Sham 방법', self.vars['method'], ['auto', 'rks', 'uks'],
                    'auto: 2S=0이면 RKS, 그 외 UKS.', self._settings_changed)
        self._combo(method, 'SCF 반복 방식', self.vars['scf_solver'], ['diis', 'newton'],
                    '기본 diis · 수렴이 어려우면 2차 newton으로 재시도.', self._settings_changed)
        self._combo(method, '실행 환경', self.backend, ['auto', 'local', 'wsl'],
                    'Windows는 WSL 엔진을 자동 사용합니다.', lambda *_: None)
        ttk.Label(method, text='현재 지원: H~Og의 고립 분자·원자·이온·착물.\n\n'
                  '무거운 원소는 기저·격자·스핀 상태를 비교해 확인하세요.\n'
                  '고체·표면은 상단 버튼을 사용하세요. 주기계와 ML 모델의 원소 지원 범위는 별도입니다.',
                  style='Muted.TLabel', wraplength=280).pack(anchor='w', pady=20)
        self._combo(advanced, 'DFT 적분 격자', self.vars['grid_scheme'], ['auto', 'treutler', 'universal'],
                    'auto: Dyall은 universal. 원소별 반경 표에 의존하지 않는 격자.', self._settings_changed)
        for name, values in PRESETS.items():
            self._button(advanced, name + ' 적용', lambda key=name: self.use_preset(key), True).pack(fill='x', pady=3)
            ttk.Label(advanced, text=', '.join(f'{key}={value}' for key, value in values.items()),
                      wraplength=280, style='Muted.TLabel').pack(anchor='w')
        self._combo(advanced, '구조 최적화 수렴 기준', self.vars['opt_convergence'], list(OPT_CONVERGENCE),
                    'default: geomeTRIC 기본(최대 기울기 4.5e-4 Ha/bohr). tight: 10배 엄격.', self._settings_changed)
        self._entry(advanced, '구조 최적화 최대 단계', 'opt_max_steps', '1~1000. 기본 100.', self._settings_changed)
        self._combo(advanced, '초기 전자 추정', self.vars['initial_guess'], ['auto', 'core', 'perturbed'],
                    'auto: 자체 core / PySCF minao. perturbed는 자체 엔진의 재현 가능한 궤도 혼합.', self._settings_changed)
        for key, label, hint in (
            ('native_density_tol', '자체 SCF 밀도 수렴', 'ΔD Frobenius norm. 기본 1e-7.'),
            ('native_residual_tol', '자체 SCF 잔차 수렴', '직교 기저의 FDS−SDF norm. 기본 1e-6.'),
            ('native_damping', '자체 밀도 혼합 비율', '이전 밀도 비율 0~0.8. 기본 0.15.'),
            ('native_diis_space', '자체 DIIS 이력 수', '2~16. 기본 8. PySCF에는 적용하지 않습니다.'),
            ('max_cycle', 'SCF 최대 반복', '1~1000. 초과하면 미수렴으로 저장합니다.'),
            ('conv_tol', 'SCF 에너지 수렴 기준 [Ha]', '1e-12~1e-3. 기본 1e-9.'),
            ('grid_level', 'DFT 적분 격자 수준', '0~7. 기본 3. 에너지 계산에 사용합니다.'),
            ('cube_points', '밀도·오비탈 격자 (축당)', '12~100. 시각화용, SCF 격자와 별개입니다.'),
            ('cube_margin_bohr', '시각화 상자 여백 [bohr]', '1~20. 기본 4.'),
            ('threads', 'CPU 스레드', '1~64. 기본 2.'),
            ('memory_mb', '엔진 메모리 목표 [MB]', '256~64000. PySCF의 메모리 목표값입니다.')):
            self._entry(advanced, label, key, hint, self._settings_changed)
        self.tabs = ttk.Notebook(body)
        self.tabs.pack(fill='both', expand=True)
        self.result_tab = ttk.Frame(self.tabs, style='Card.TFrame', padding=16)
        self.orbital_tab = ttk.Frame(self.tabs, style='Card.TFrame', padding=16)
        self.detail_tab = ttk.Frame(self.tabs, style='Card.TFrame', padding=12)
        self.history_tab = ttk.Frame(self.tabs, style='Card.TFrame', padding=16)
        self.log_tab = ttk.Frame(self.tabs, style='Card.TFrame', padding=12)
        self.help_tab = ttk.Frame(self.tabs, style='Card.TFrame', padding=12)
        for frame, title in ((self.result_tab, '계산 결과'), (self.orbital_tab, '오비탈'), (self.detail_tab, '설정과 원시 결과'),
                             (self.history_tab, '실행 기록'), (self.log_tab, '실행 로그'), (self.help_tab, '사용 안내')):
            self.tabs.add(frame, text=title)
        self._build_results()
        self._build_scf_view()
        self.detail_text = self._text(self.detail_tab, '계산하면 재현 가능한 설정과 결과 JSON을 표시합니다.')
        self.log_text = self._text(self.log_tab, 'SCF 반복과 실행 오류를 표시합니다.')
        self._text(self.help_tab, HELP)
        self._build_history()

    def open_periodic(self):
        from periodic_app import PeriodicApp
        if getattr(self, 'periodic_window', None) is not None and self.periodic_window.winfo_exists():
            self.periodic_window.lift()
            return
        self.periodic_window = tk.Toplevel(self.root)
        self.periodic_app = PeriodicApp(self.periodic_window)

    def close(self):
        child = getattr(self, 'periodic_app', None)
        if child is not None and child.busy:
            self.status.set('고체·표면 계산 결과를 저장한 뒤 창을 닫습니다.')
            self.root.after(500, self.close)
            return
        super().close()

    def _build_results(self):
        parent = self.result_tab
        parent.columnconfigure(0, weight=1)
        parent.rowconfigure(3, weight=1)
        ttk.Label(parent, textvariable=self.changed, style='Muted.TLabel', wraplength=850).grid(row=0, column=0, sticky='w', pady=(0, 10))
        cards = ttk.Frame(parent, style='Card.TFrame')
        cards.grid(row=1, column=0, sticky='ew')
        self.metrics = []
        for i, label in enumerate(('총에너지 · 선택한 계산 모델', 'SCF 수렴 상태', '계산 방법')):
            cards.columnconfigure(i, weight=1)
            card = tk.Frame(cards, bg='#f2f6fc', padx=12, pady=12)
            card.grid(row=0, column=i, sticky='ew', padx=3)
            tk.Label(card, text=label, bg='#f2f6fc', fg=MUTED, font=(self.font, 9)).pack(anchor='w')
            value = tk.StringVar(value='—')
            tk.Label(card, textvariable=value, bg='#f2f6fc', fg=INK, font=(self.font, 13, 'bold')).pack(anchor='w', pady=(6, 0))
            self.metrics.append(value)
        filters = ttk.Frame(parent, style='Card.TFrame')
        filters.grid(row=2, column=0, sticky='ew', pady=10)
        self.field_box = ttk.Combobox(filters, textvariable=self.field, values=['density.cube'], state='readonly', width=23)
        self.field_box.pack(side='left', padx=(0, 8))
        self.field_box.bind('<<ComboboxSelected>>', self.redraw)
        plane = ttk.Combobox(filters, textvariable=self.plane, values=['xy', 'xz', 'yz'], state='readonly', width=5)
        plane.pack(side='left')
        plane.bind('<<ComboboxSelected>>', self.redraw)
        self._button(filters, '그래프 저장', self.export_plot).pack(side='right')
        self.chart_frame = ttk.Frame(parent, style='Card.TFrame')
        self.chart_frame.grid(row=3, column=0, sticky='nsew')
        self.empty = ttk.Label(self.chart_frame, text='원자 좌표로 시작하는 분자 DFT\n\n'
                               '왼쪽에서 구조·전하·스핀을 입력한 뒤 계산 시작을 누르세요.',
                               style='Section.TLabel', justify='center')
        self.empty.pack(expand=True)
        ttk.Label(parent, text='원자 중심에 가까운 격자면의 단면 · 전자밀도 e/bohr³ · 좌표 Å',
                  style='Muted.TLabel').grid(row=4, column=0, sticky='w', pady=5)
        self.orbital_context = ttk.Label(self.orbital_tab, text='계산하면 모든 오비탈의 점유수와 에너지를 표시합니다.',
                                         style='Muted.TLabel', wraplength=820)
        self.orbital_context.pack(anchor='w', pady=(0, 12))
        table_frame = ttk.Frame(self.orbital_tab, style='Card.TFrame')
        table_frame.pack(fill='both', expand=True)
        self.table = ttk.Treeview(table_frame, columns=('channel', 'index', 'label', 'occ', 'ha', 'ev'), show='headings', height=14)
        for key, title, width in (('channel', '채널', 85), ('index', 'MO', 45), ('label', '경계 오비탈', 95),
                                  ('occ', '점유수', 70), ('ha', '에너지 [Ha]', 140), ('ev', '에너지 [eV]', 140)):
            self.table.heading(key, text=title)
            self.table.column(key, width=width, anchor='center')
        self.table.pack(side='left', fill='both', expand=True)
        scroll = ttk.Scrollbar(table_frame, command=self.table.yview)
        scroll.pack(side='right', fill='y')
        self.table.configure(yscrollcommand=scroll.set)
        ttk.Label(self.orbital_tab, text='오비탈 에너지의 합은 총에너지가 아닙니다. '
                  'Kohn–Sham 에너지 차이는 실험적 여기 에너지와 다릅니다.', style='Muted.TLabel',
                  wraplength=820).pack(anchor='w', pady=12)
        self.quality = ttk.Label(parent, text='SCF 수렴 여부와 계산의 한계를 함께 표시합니다.', style='Muted.TLabel', wraplength=850)
        self.quality.grid(row=6, column=0, sticky='w', pady=10)
        actions = ttk.Frame(parent, style='Card.TFrame')
        actions.grid(row=7, column=0, sticky='ew')
        self._button(actions, '결과 폴더 열기', self.open_current).pack(side='left', padx=(0, 8))
        self._button(actions, '결과 ZIP 내보내기', self.export_results).pack(side='left')
        parent.bind('<Configure>', self._resize_labels)

    def _coordinates_changed(self, *_):
        if self.coordinates.edit_modified():
            self.coordinates.edit_modified(False)
            self._settings_changed()

    def read_config(self):
        data = {key: var.get().strip() for key, var in self.vars.items()}
        for key in ('charge', 'spin', 'max_cycle', 'grid_level', 'cube_points', 'threads', 'memory_mb', 'native_diis_space',
                    'opt_max_steps'):
            try:
                data[key] = int(data[key])
            except ValueError:
                raise ValueError(f'{key}: 정수를 입력하세요.') from None
        for key in ('conv_tol', 'cube_margin_bohr', 'native_density_tol', 'native_residual_tol', 'native_damping'):
            try:
                data[key] = float(data[key])
            except ValueError:
                raise ValueError(f'{key}: 숫자를 입력하세요.') from None
        data['basis_overrides'] = parse_basis_overrides(data['basis_overrides'])
        return DFTConfig(atoms=parse_coordinates(self.coordinates.get('1.0', 'end')), **data)

    def apply_config(self, cfg, initial=False):
        self.applying = True
        for key, var in self.vars.items():
            value = getattr(cfg, key)
            var.set(format_basis_overrides(value) if key == 'basis_overrides' else str(value))
        self.coordinates.delete('1.0', 'end')
        self.coordinates.insert('1.0', coordinates_text(cfg.atoms))
        self.coordinates.edit_modified(False)
        self.applying = False
        self._settings_changed()

    def _settings_changed(self, *_):
        if self.applying:
            return
        # Filter using only the geometry, so an invalid charge/spin does not hide
        # the compatible basis choices. Keep the user's selection until changed.
        try:
            atoms = parse_coordinates(self.coordinates.get('1.0', 'end'))
            self.basis_box.configure(values=available_bases(atoms))
        except ValueError:
            self.basis_box.configure(values=BASES)
        try:
            cfg = self.read_config()
            ecp_info = (f'ECP 코어 {cfg.ecp_core_electrons}개 제외 · 계산 전자 {cfg.explicit_electrons}개\n'
                        if cfg.ecp_core_electrons else '')
            self.preview.set(f'{cfg.formula} · 원자 {len(cfg.atoms)}개 · 전체 전자 {cfg.electrons}개\n' + ecp_info +
                             f'전하 {cfg.charge:+d} · 2S={cfg.spin} · {cfg.resolved_method.upper()} · {cfg.engine}')
            self.ecp_hint.set('def2-ECP: ' + ', '.join(cfg.resolved_ecp) + '\n저장 밀도에는 코어 전자가 없습니다.'
                              if cfg.resolved_ecp else f'전체전자 · 상대론 {cfg.resolved_relativity} · 핵 {cfg.resolved_nuclear_model}')
            matches = self.summary and cfg == DFTConfig.from_dict(self.summary['config'])
        except ValueError as exc:
            self.preview.set(str(exc))
            self.ecp_hint.set('기저·ECP·전하·스핀 설정을 확인하세요.')
            matches = False
        if self.summary:
            self.changed.set('표시된 결과는 현재 입력 조건의 계산입니다.' if matches else
                             '설정 변경됨 · 아래 결과는 이전 조건입니다. 계산 시작을 눌러 갱신하세요.')
            status = '수렴' if self.summary['converged'] else '미수렴 · 진단값'
            self.orbital_context.configure(text=f'{self.summary["formula"]} · {self.summary["method"]} / '
                f'{self.summary["xc"]} / '
                f'{describe_bases(self.summary["basis"], self.summary.get("basis_by_element", {}))} · {status}\n'
                + self.changed.get())

    def show_elements(self):
        window = tk.Toplevel(self.root)
        window.title('원소별 계산 지원 · H~Og')
        window.geometry('1100x720')
        ttk.Label(window, text='분자: Dyall-v2z/v3z로 118개 원소 지원 · 기저 데이터와 계산 정확도는 별개입니다.\n'
                  'def2 코어는 대체되는 전자 수입니다. 주기계 열은 GTH-PBE 기저 조합이며 ML 학습 범위와 다릅니다.',
                  padding=12, wraplength=1000).pack(fill='x')
        body = ttk.Frame(window, padding=10)
        body.pack(fill='both', expand=True)
        columns = ('z', 'symbol', 'molecular', 'core', 'periodic')
        tree = ttk.Treeview(body, columns=columns, show='headings')
        for key, label, width in zip(columns, ('번호', '원소', '분자 기저', 'def2 코어', '주기계 기저'), (55, 60, 460, 85, 330)):
            tree.heading(key, text=label)
            tree.column(key, width=width, minwidth=width, stretch=False)
        vertical = ttk.Scrollbar(body, orient='vertical', command=tree.yview)
        horizontal = ttk.Scrollbar(body, orient='horizontal', command=tree.xview)
        tree.configure(yscrollcommand=vertical.set, xscrollcommand=horizontal.set)
        tree.grid(row=0, column=0, sticky='nsew')
        vertical.grid(row=0, column=1, sticky='ns')
        horizontal.grid(row=1, column=0, sticky='ew')
        body.rowconfigure(0, weight=1)
        body.columnconfigure(0, weight=1)
        for z, symbol in enumerate(ELEMENTS, 1):
            periodic = ', '.join(b for b, info in PERIODIC_CAPABILITIES.items() if symbol in info['valence_electrons'])
            tree.insert('', 'end', iid=symbol, values=(z, symbol,
                        ', '.join(available_bases((Atom(symbol, (0, 0, 0)),))),
                        DEF2_CORE_ELECTRONS.get(symbol, '—'), periodic or '미지원'))
        return window, tree

    def use_example(self, *_):
        text, charge, spin, *options = EXAMPLES[self.example.get()]
        self.apply_config(DFTConfig(atoms=parse_coordinates(text), charge=charge, spin=spin,
                                   **(options[0] if options else {})))

    def use_preset(self, name):
        try:
            self.apply_config(apply_preset(self.read_config(), name))
            self.status.set(name + ': ' + ', '.join(f'{k}={v}' for k, v in PRESETS[name].items()))
        except ValueError as exc:
            messagebox.showerror('입력 확인', str(exc), parent=self.root)

    def preview_geometry(self):
        try:
            atoms = parse_coordinates(self.coordinates.get('1.0', 'end'))
            unit = self.vars['unit'].get()
            if unit not in ('angstrom', 'bohr'):
                raise ValueError('좌표 단위를 확인하세요.')
        except ValueError as exc:
            messagebox.showerror('구조 입력 확인', str(exc), parent=self.root)
            return
        window = tk.Toplevel(self.root)
        window.title('입력 구조 · 계산 전 미리보기')
        window.geometry('820x720')
        figure = Figure(figsize=(7, 5), layout='constrained')
        axis = figure.add_subplot(111, projection='3d')
        factor = BOHR_ANGSTROM if unit == 'bohr' else 1
        points = [tuple(x * factor for x in atom.position) for atom in atoms]
        axis.scatter(*zip(*points), s=100)
        for i, (atom, point) in enumerate(zip(atoms, points), 1):
            axis.text(*point, f' {i}:{atom.element}')
        axis.set(xlabel='x [Å]', ylabel='y [Å]', zlabel='z [Å]')
        canvas = AppCanvas(figure, master=window)
        canvas.get_tk_widget().pack(fill='both', expand=True)
        canvas.draw_idle()
        row = ttk.Frame(window, padding=10)
        row.pack(fill='x')
        ttk.Label(row, text='원자 번호 (거리: 1 2 / 각도: 1 2 3, 가운데가 꼭짓점)').pack(anchor='w')
        indices = tk.StringVar(value='1 2' if len(atoms) > 1 else '1')
        ttk.Entry(row, textvariable=indices).pack(fill='x')
        result = tk.StringVar(value='구조의 원자 번호를 입력하세요.')
        def calculate_measurement():
            try:
                value, label = measure(atoms, [int(i) for i in indices.get().split()], unit)
                result.set(f'{value:.8f} {label}')
            except ValueError as exc:
                result.set(str(exc))
        ttk.Button(row, text='거리/각도 확인', command=calculate_measurement).pack(anchor='w', pady=6)
        ttk.Label(row, textvariable=result).pack(anchor='w')
        # Keep references alive for Tk callbacks without modifying the result view.
        window.preview_canvas = canvas
        return window

    def load_xyz(self):
        name = filedialog.askopenfilename(title='단일 구조 XYZ 불러오기', filetypes=[('XYZ 좌표', '*.xyz')])
        if name:
            try:
                atoms = read_xyz(name)
                self.coordinates.delete('1.0', 'end')
                self.coordinates.insert('1.0', coordinates_text(atoms))
                self.vars['unit'].set('angstrom')
                self._settings_changed()
                self.status.set('XYZ 좌표를 불러왔습니다. 전체 전하와 스핀은 직접 확인하세요.')
            except (ValueError, OSError) as exc:
                messagebox.showerror('XYZ 입력 오류', str(exc), parent=self.root)

    def load_settings(self):
        name = filedialog.askopenfilename(title='분자 설정 불러오기', filetypes=[('JSON 설정', '*.json')])
        if name:
            try:
                self.apply_config(load_dft_config(name))
            except (ValueError, OSError) as exc:
                messagebox.showerror('설정 입력 오류', str(exc), parent=self.root)

    def save_settings(self):
        try:
            cfg = self.read_config()
            name = filedialog.asksaveasfilename(initialfile='molecule.json', defaultextension='.json', filetypes=[('JSON 설정', '*.json')])
            if name:
                save_dft_config(cfg, name)
        except (ValueError, OSError) as exc:
            messagebox.showerror('설정 저장 오류', str(exc), parent=self.root)

    def _build_scf_view(self):
        self.scf_tab = ttk.Frame(self.tabs, padding=12)
        self.tabs.add(self.scf_tab, text='SCF 실시간')
        ttk.Label(self.scf_tab, text='잔차: 자체 엔진은 직교 기저 commutator, PySCF는 orbital gradient. '
                  '서로 다른 잔차를 같은 기준으로 비교하지 마세요.', wraplength=850).pack(fill='x')
        self.scf_figure = Figure(figsize=(8, 5), layout='constrained')
        self.scf_axes = self.scf_figure.subplots(2, 1)
        self.scf_canvas = AppCanvas(self.scf_figure, master=self.scf_tab)
        self.scf_canvas.get_tk_widget().pack(fill='both', expand=True)
        self.scf_rows = []

    def _plot_scf(self):
        for axis in self.scf_axes:
            axis.clear()
        if self.scf_rows:
            x = [r['cycle'] for r in self.scf_rows]
            self.scf_axes[0].plot(x, [r['energy_hartree'] for r in self.scf_rows], marker='.')
            for key, label in (('gradient_norm', 'Residual'), ('density_change_norm', 'Density change')):
                rows = [r for r in self.scf_rows if r.get(key) is not None]
                if rows:
                    self.scf_axes[1].semilogy([r['cycle'] for r in rows], [max(r[key], 1e-16) for r in rows], label=label)
            if self.scf_axes[1].lines:
                self.scf_axes[1].legend()
        self.scf_axes[0].set_ylabel('Total energy [Ha]')
        self.scf_axes[1].set(xlabel='SCF iteration', ylabel='Norm')
        self.scf_canvas.draw_idle()

    def _set_busy(self, busy):
        super()._set_busy(busy)
        self.cancel_button.configure(state='normal' if busy else 'disabled')

    def cancel_calculation(self):
        if self.busy and hasattr(self, 'run_control'):
            self.run_control.cancel()
            self.status.set('중단 요청됨 · 현재 적분/SCF 단계가 끝나면 체크포인트를 보존합니다.')
            self.cancel_button.configure(state='disabled')

    def restart_history(self):
        entry = self._selected_history()
        if not entry:
            messagebox.showinfo('실행 기록 선택', '재시작할 실행을 먼저 선택하세요.', parent=self.root)
            return
        checkpoint = entry['directory'] / 'checkpoint.npz'
        if not checkpoint.exists():
            messagebox.showerror('체크포인트 없음', '이 기록에는 재시작할 밀도 체크포인트가 없습니다.', parent=self.root)
            return
        self.apply_config(DFTConfig.from_dict(entry['summary']['config']))
        self.start_calculation(restart_from=checkpoint)

    def start_calculation(self, restart_from=None):
        if self.busy:
            return
        try:
            config = self.read_config()
        except ValueError as exc:
            messagebox.showerror('계산 조건 확인', str(exc), parent=self.root)
            return
        self._set_busy(True)
        self.run_control = RunControl()
        self.scf_rows = []
        self._plot_scf()
        self.tabs.select(self.scf_tab)
        self.status.set('분자 DFT 엔진을 시작합니다.')
        self.log_lines = []
        self._log('DFT 계산 시작')
        threading.Thread(target=self._worker, args=(config, self.backend.get(), restart_from), daemon=True).start()

    def _worker(self, config, backend, restart_from=None):
        try:
            summary, directory = run_calculation(config, self.output_root, backend,
                lambda text: self.messages.put(('stage', text)), control=self.run_control, restart_from=restart_from)
            self.messages.put(('done', summary, directory))
        except Exception as exc:
            self.messages.put(('error', str(exc), traceback.format_exc(), getattr(exc, 'directory', None)))

    def _poll(self):
        deadline = time.monotonic() + .03
        redraw_scf = False
        try:
            while time.monotonic() < deadline:
                event = self.messages.get_nowait()
                if event[0] == 'stage':
                    if event[1].startswith('SCF_JSON '):
                        self.scf_rows.append(json.loads(event[1][9:]))
                        redraw_scf = True
                    elif event[1].startswith('OPT_JSON '):
                        pass  # Machine-readable twin of the optimization line already logged.
                    else:
                        self._log(event[1])
                        self.status.set(event[1][-180:])
                elif event[0] == 'done':
                    if redraw_scf:
                        self._plot_scf()
                        redraw_scf = False
                    self._complete(*event[1:])
                elif event[0] == 'error':
                    self._log(event[2])
                    self._set_busy(False)
                    self.refresh_history()
                    self.status.set('계산 실패 · 실행 로그와 기록에서 원인을 확인하세요.')
                    if self.summary:
                        self.changed.set('최근 계산은 실패했습니다. 아래에는 이전 실행 결과가 남아 있습니다.')
                    if self.closing:
                        self._destroy()
                        return
                    messagebox.showerror('계산 실패', event[1] + f'\n\n기록: {event[3]}', parent=self.root)
                if self.closed:
                    return
        except queue.Empty:
            pass
        if not self.closed:
            if redraw_scf:
                self._plot_scf()
            self.root.after(80, self._poll)

    def _complete(self, summary, directory):
        if summary['status'] == 'cancelled':
            self._set_busy(False)
            self.refresh_history()
            self.status.set(f'계산 중단됨 · 체크포인트: {summary.get("checkpoint_file") or "첫 저장 전 중단"} · {directory}')
            self._set_text(self.detail_text, json.dumps(summary, ensure_ascii=False, indent=2))
            if self.summary:
                self.changed.set('최근 실행은 중단되었습니다. 표시된 계산 결과는 이전 실행입니다.')
            if self.closing:
                self._destroy()
            return
        self.summary = self.result = summary
        self.current_dir = directory
        self.plane.set('yz')
        dispersion = summary.get('dispersion')
        self.metrics[0].set(f'{summary["total_energy_hartree"]:.9f} Ha'
                            + (f' (D4 {dispersion["energy_hartree"]:+.6f} 포함)' if dispersion else ''))
        iteration = '외부 반복' if summary.get('scf_solver') == 'newton' else '반복'
        gradient = summary.get('gradient')
        self.metrics[1].set(f'{"수렴" if summary["converged"] else "미수렴"} / {iteration} {summary["scf_cycles"]}회'
                            + (f' · 최대 힘 {gradient["max_atomic_force_ev_per_angstrom"]:.4f} eV/Å' if gradient else '')
                            + (f' · 최적화 {"수렴" if summary["optimization"]["converged"] else "미수렴"} '
                               f'{summary["optimization"]["steps"]}단계' if summary.get('optimization') else ''))
        engine_label = '자체' if summary['engine'] == 'quantumstarter-native' else 'PySCF'
        self.metrics[2].set(f'{engine_label} · {summary["method"]} / {summary["xc"].upper()}'
                            + ('-D4' if dispersion else ''))
        self.field_box.configure(values=summary['cube_files'])
        self.field.set('density.cube')
        self.table.delete(*self.table.get_children())
        for mo in summary['orbitals']:
            self.table.insert('', 'end', values=(mo['channel'], mo['index'], mo['label'],
                f'{mo["occupation"]:.1f}', f'{mo["energy_hartree"]:.7f}', f'{mo["energy_ev"]:.6f}'))
        self.quality.configure(text=' '.join(summary['warnings']), foreground='#b45309')
        self._set_text(self.detail_text, json.dumps(summary, ensure_ascii=False, indent=2))
        self.scf_figure.savefig(directory / 'scf_convergence.png', dpi=140)
        status = '수렴 완료' if summary['converged'] else 'SCF 미수렴 · 성공한 계산이 아닙니다'
        try:
            self.redraw()
            create_molecular_figure(summary, directory).savefig(directory / 'molecular.png', dpi=150)
        except Exception as exc:
            self._log(traceback.format_exc())
            status += f' · 원시 결과 저장됨 / 그래프 표시 실패: {exc}'
        self._set_busy(False)
        self._settings_changed()
        self.refresh_history()
        self.status.set(status)
        if self.closing:
            self._destroy()

    def redraw(self, *_):
        if self.summary is None:
            return
        figure = create_molecular_figure(self.summary, self.current_dir, self.field.get(), self.plane.get())
        if self.toolbar is not None:
            self.toolbar.destroy()
        if self.canvas is not None:
            pending = getattr(self.canvas, '_idle_draw_id', None)
            if pending:
                self.canvas.get_tk_widget().after_cancel(pending)
            self.canvas.get_tk_widget().destroy()
        self.empty.pack_forget()
        self.figure = figure
        figure.set_dpi(float(self.root.winfo_fpixels('1i')))
        self.canvas = AppCanvas(figure, master=self.chart_frame)
        self.toolbar = NavigationToolbar2Tk(self.canvas, self.chart_frame, pack_toolbar=False)
        self.toolbar.pack(side='bottom', fill='x')
        self.canvas.get_tk_widget().pack(fill='both', expand=True)
        self.canvas.draw_idle()

    def export_plot(self):
        if self.summary is None:
            return
        name = filedialog.asksaveasfilename(initialfile='molecular_view.png', defaultextension='.png', filetypes=[('PNG 그림', '*.png')])
        if name:
            try:
                create_molecular_figure(self.summary, self.current_dir, self.field.get(), self.plane.get()).savefig(name, dpi=180)
            except (ValueError, OSError) as exc:
                messagebox.showerror('그림 저장 실패', str(exc), parent=self.root)

    def _build_history(self):
        ttk.Label(self.history_tab, text='수렴·미수렴·중단·실패와 엔진 출처를 구분합니다.', style='Section.TLabel').pack(anchor='w', pady=(0, 12))
        self.history_tree = ttk.Treeview(self.history_tab, columns=('time', 'system', 'status', 'energy'), show='headings', height=8)
        for key, title, width in (('time', '실행 폴더', 230), ('system', '분자 / 방법', 210),
                                  ('status', '상태', 110), ('energy', '총에너지 [Ha]', 160)):
            self.history_tree.heading(key, text=title)
            self.history_tree.column(key, width=width)
        self.history_tree.pack(fill='x')
        self.history_tree.bind('<<TreeviewSelect>>', self._history_selected)
        row = ttk.Frame(self.history_tab, style='Card.TFrame')
        row.pack(fill='x', pady=10)
        self._button(row, '과거 조건 불러오기', self.load_history_settings, True).pack(side='left', padx=(0, 8))
        self._button(row, '체크포인트 재시작', self.restart_history, True).pack(side='left', padx=(0, 8))
        self._button(row, '결과 폴더 열기', self.open_history).pack(side='left')
        actions = ttk.Frame(self.history_tab)
        actions.pack(fill='x', pady=5)
        self._button(actions, '선택한 기록 비교', self.compare_history).pack(side='left', padx=(0, 8))
        self._button(actions, '흡착에너지 분석', self.analyze_adsorption).pack(side='left')
        self.history_text = self._text(self.history_tab, '실행을 선택하면 전체 기록을 표시합니다.')

    def compare_history(self):
        from energy_analysis import compare_runs, save_analysis
        try:
            paths = [self.history[int(i)]['directory'] for i in self.history_tree.selection()]
            result = compare_runs(paths)
            path = save_analysis(result, ROOT / 'outputs' / 'analysis')
            self._set_text(self.history_text, json.dumps(result, ensure_ascii=False, indent=2))
            self.status.set(f'비교 기록 저장: {path}')
        except (ValueError, OSError, KeyError) as exc:
            messagebox.showerror('기록 비교', str(exc), parent=self.root)

    def analyze_adsorption(self):
        from energy_analysis import adsorption_energy, save_analysis, ENERGY_FIELDS
        window = tk.Toplevel(self.root)
        window.title('흡착에너지 · 기준과 계산 조건 확인')
        window.geometry('800x580')
        ttk.Label(window, text='E_ads = E(촉매+흡착종) − E(촉매) − 계수 × E(흡착종 기준)\n'
                  '음수: 선택 기준에서 전자적으로 유리한 흡착. 활성·선택성의 증거는 아닙니다.\n'
                  '수렴·엔진·기저·XC·격자·원소/전하 수지를 검사합니다. 주기계는 셀/k점도 같아야 합니다.',
                  padding=12, wraplength=760).pack(fill='x')
        paths = []
        for title in ('촉매+흡착종', '촉매', '흡착종 기준'):
            row = ttk.Frame(window, padding=8)
            row.pack(fill='x')
            ttk.Label(row, text=title, width=16).pack(side='left')
            value = tk.StringVar()
            paths.append(value)
            ttk.Entry(row, textvariable=value).pack(side='left', fill='x', expand=True)
            def choose(var=value):
                path = filedialog.askopenfilename(title='results.json 선택', filetypes=[('계산 결과', '*.json')], parent=window)
                if path:
                    var.set(path)
            ttk.Button(row, text='결과 선택', command=choose).pack(side='right')
        factor = tk.StringVar(value='1')
        ttk.Label(window, text='흡착종 기준 계수 (예: H 기준으로 H₂를 사용할 때 0.5)').pack(anchor='w', padx=12)
        ttk.Entry(window, textvariable=factor).pack(fill='x', padx=12)
        energy = tk.StringVar(value=ENERGY_FIELDS[0])
        ttk.Combobox(window, textvariable=energy, values=ENERGY_FIELDS, state='readonly').pack(fill='x', padx=12, pady=8)
        button = ttk.Button(window, text='조건 검사 및 흡착에너지 계산', command=lambda: analyze())
        button.pack(side='bottom', pady=10)
        output = self._text(window, '결과 파일을 지정한 뒤 분석하세요. ZPE/열역학 보정은 자동 포함하지 않습니다.')
        def analyze():
            try:
                result = adsorption_energy(*(v.get() for v in paths), float(factor.get()), energy.get())
                path = save_analysis(result, ROOT / 'outputs' / 'analysis')
                self._set_text(output, f"{result['definition']}\n{result['energy_ev']:.9f} eV\n"
                               + result['interpretation'] + f'\n저장: {path}')
            except (ValueError, OSError, KeyError) as exc:
                self._set_text(output, '분석 중단: ' + str(exc))
        return window

    def refresh_history(self):
        self.history = list_dft_runs(self.output_root)
        self.history_tree.delete(*self.history_tree.get_children())
        for i, entry in enumerate(self.history):
            data = entry['summary']
            energy = data.get('total_energy_hartree')
            self.history_tree.insert('', 'end', iid=str(i), values=(entry['directory'].name,
                f'{data["formula"]} / {data["config"].get("engine", "pyscf")} / {data["config"]["xc"]}',
                {'converged': '수렴', 'not_converged': '미수렴', 'cancelled': '중단', 'error': '실패'}[data['status']],
                f'{energy:.9f}' if energy is not None else '—'))

    def load_history_settings(self):
        entry = self._selected_history()
        if entry:
            self.apply_config(DFTConfig.from_dict(entry['summary']['config']))
            self.status.set('과거 조건을 불러왔습니다. 계산 시작을 눌러 재실행하세요.')
