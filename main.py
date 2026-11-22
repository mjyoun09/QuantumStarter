"""명령행 실행 유지: python main.py --no-show. 사용 화면: python app.py.

중성 He 방사 격자 LDA-SCF: python main.py --model he-lda --no-show
원자 LDA/LSDA 및 PBE: python main.py --model atom --element Ar --functional pbe --no-show
"""
import argparse
from pathlib import Path
import sys
import time
import numpy as np
from plotting import save_plot
from solver import HARTREE_TO_EV, Config, solve
from workflow import VERSION, load_config, new_run_dir, save_result

ROOT = Path(__file__).resolve().parent
# CLI 옵션 이름 -> atomic_dft.AtomicDFTConfig 필드.
ATOMIC_OPTIONS = {'mixing': 'mixing', 'max_iterations': 'max_iterations', 'energy_tol': 'energy_tol_hartree',
                  'density_tol': 'density_tol_electrons', 'r_max': 'r_max_bohr', 'grid_points': 'grid_points'}
# CLI 옵션 이름 -> radial_atom.AtomConfig 필드 (--model atom).
ATOM_OPTIONS = {'element': 'element', 'configuration': 'configuration', 'spin_polarized': 'spin_polarized',
                'mixing': 'mixing', 'max_iterations': 'max_iterations', 'energy_tol': 'energy_tol_hartree',
                'density_tol': 'density_tol_electrons', 'r_max': 'r_max_bohr', 'grid_step': 'grid_step',
                'grid_scale': 'grid_scale', 'functional': 'functional', 'mixing_method': 'mixing_method'}
ATOM_ONLY_OPTIONS = ('element', 'configuration', 'spin_polarized', 'grid_step', 'grid_scale', 'functional', 'mixing_method')
# --model he-lda / atom과 함께 쓸 수 없는 분자·주기계 옵션.
DFT_ONLY_OPTIONS = ('charge', 'spin', 'basis', 'xc', 'dispersion', 'calculation', 'opt_convergence', 'opt_max_steps',
                    'freq_step_bohr', 'thermo_temperature_k', 'thermo_pressure_pa', 'symmetry_number',
                    'native_stability', 'opt_constraints', 'scan_spec', 'neb_images', 'scf_solver', 'ecp',
                    'relativity', 'nuclear_model', 'grid_scheme', 'engine', 'initial_guess', 'xyz', 'structure',
                    'restart', 'basis_element', 'neb_product')


def run_he_lda(config, output_root: Path, show: bool):
    """중성 He LDA-SCF 실행·저장·출력. 종료 코드: 수렴 0, 미수렴 2."""
    from atomic_dft import XC_FUNCTIONALS, run_scf, save_atomic_result
    from atomic_plotting import create_atomic_figure
    directory = new_run_dir(output_root)
    xc = XC_FUNCTIONALS
    print(f"QuantumStarter {VERSION} | 중성 He Kohn–Sham LDA-SCF (구대칭 방사 격자, 자체 구현)")
    print(f"XC: {xc['exchange']['libxc_name']} (Libxc ID {xc['exchange']['libxc_id']}) + "
          f"{xc['correlation']['libxc_name']} (Libxc ID {xc['correlation']['libxc_id']}, VWN5), 스핀 비분극 | "
          f"Z = 2, 전자 2개, 1s²")
    print(f"격자: dr = {config.dr_bohr:.6g} bohr, r_max = {config.r_max_bohr:g} bohr, 내부 점 {config.grid_points} | "
          f"혼합 {config.mixing:g}, 반복 상한 {config.max_iterations}, "
          f"|ΔE| < {config.energy_tol_hartree:.1e} Ha, 밀도 잔차 < {config.density_tol_electrons:.1e} e")
    print(" 반복    E_total [Ha]        ΔE [Ha]     밀도 잔차 [e]   ε_1s [Ha]")

    def report(record):
        delta = record['delta_energy_hartree']
        print(f"{record['iteration']:5d}  {record['total_energy_hartree']:17.12f}  "
              f"{'-' if delta is None else f'{delta:.3e}':>10s}  {record['density_residual_electrons']:12.3e}  "
              f"{record['eigenvalue_1s_hartree']:16.12f}")

    started = time.perf_counter()
    result = run_scf(config, progress=report)
    summary = save_atomic_result(result, directory, time.perf_counter() - started)
    conv = summary['convergence']
    delta = conv['final_delta_energy_hartree']
    detail = (f"|ΔE| = {'-' if delta is None else f'{abs(delta):.2e}'} Ha, "
              f"밀도 잔차 = {conv['final_density_residual_electrons']:.2e} e")
    if result.converged:
        print(f"SCF 수렴: {summary['iterations']}회 ({detail})")
        label = ''
    else:
        print(f"SCF 미수렴: 반복 상한 {config.max_iterations}회 도달 (마지막 {detail}; "
              f"기준 {config.energy_tol_hartree:.1e} Ha / {config.density_tol_electrons:.1e} e)")
        label = ' [미수렴 마지막 값]'
    print(f"총에너지   E_total = {summary['total_energy_hartree']:.12f} Ha ({summary['total_energy_ev']:.8f} eV){label}")
    print(f"1s 고윳값  ε_1s    = {summary['eigenvalue_1s_hartree']:.12f} Ha ({summary['eigenvalue_1s_ev']:.8f} eV){label}")
    print(f"적분 전자 수 N = {summary['electron_count']:.12f}")
    for warning in summary['warnings']:
        print('확인:', warning)
    figure = create_atomic_figure(result)
    plot_path = directory / 'he_lda.png'
    figure.savefig(plot_path, dpi=150)
    print(f"상태: {summary['status']} | 결과 폴더: {directory}")
    if show:
        import matplotlib.pyplot as plt
        import matplotlib.image as mpimg
        plt.figure(figsize=figure.get_size_inches())
        plt.imshow(mpimg.imread(plot_path))
        plt.axis('off')
        plt.show()
    return 0 if result.converged else 2


def run_atom_cli(config, output_root: Path, show: bool):
    """구면 평균 원자 SCF를 실행하고 저장한다. 수렴 시 0, 미수렴 시 2를 반환한다."""
    from radial_atom import functional_label, functional_metadata, log_grid, run_atom, save_atom_result
    from atomic_plotting import create_atom_figure
    directory = new_run_dir(output_root)
    grid = log_grid(config)
    xc = functional_metadata(config)
    kind = functional_label(config)
    print(f"QuantumStarter {VERSION} | 원자 Kohn–Sham {kind}, 구면 평균 밀도, 로그 방사 격자 (자체 구현)")
    print(f"원소 {config.element} (Z = {config.nuclear_charge}) | 배치 {config.configuration_text} | "
          f"XC: {xc['exchange']['libxc_name']} (Libxc ID {xc['exchange']['libxc_id']}) + "
          f"{xc['correlation']['libxc_name']} (Libxc ID {xc['correlation']['libxc_id']})")
    print(f"격자: r = a(e^x − 1), h = {grid.h:g}, a = {grid.a:.3g} bohr, 점 {grid.r.size}, r_end = {grid.r_end:.4g} bohr | "
          f"혼합 {config.mixing_method} {config.mixing:g}, 반복 상한 {config.max_iterations}, |ΔE| < {config.energy_tol_hartree:.1e} Ha, "
          f"밀도 잔차 < {config.density_tol_electrons:.1e} e")
    print(" 반복    E_total [Ha]            ΔE [Ha]     밀도 잔차 [e]")

    def report(record):
        delta = record['delta_energy_hartree']
        print(f"{record['iteration']:5d}  {record['total_energy_hartree']:20.12f}  "
              f"{'-' if delta is None else f'{delta:.3e}':>10s}  {record['density_residual_electrons']:12.3e}")

    started = time.perf_counter()
    result = run_atom(config, progress=report)
    summary = save_atom_result(result, directory, time.perf_counter() - started)
    conv = summary['convergence']
    delta = conv['final_delta_energy_hartree']
    detail = (f"|ΔE| = {'-' if delta is None else f'{abs(delta):.2e}'} Ha, "
              f"밀도 잔차 = {conv['final_density_residual_electrons']:.2e} e")
    label = '' if result.converged else ' [미수렴 마지막 값]'
    if result.converged:
        print(f"SCF 수렴: {summary['iterations']}회 ({detail})")
    else:
        print(f"SCF 미수렴: 반복 상한 {config.max_iterations}회 도달 (마지막 {detail}; "
              f"기준 {config.energy_tol_hartree:.1e} Ha / {config.density_tol_electrons:.1e} e)")
    print(f"{'오비탈':6s} {'스핀':12s} {'점유':>6s}  {'고윳값 [Ha]':>18s}  {'고윳값 [eV]':>14s}")
    for orbital in summary['orbitals']:
        print(f"{orbital['label']:6s} {orbital['spin']:12s} {orbital['occupation']:6.3g}  "
              f"{orbital['eigenvalue_hartree']:18.9f}  {orbital['eigenvalue_hartree'] * HARTREE_TO_EV:14.6f}")
    print(f"총에너지 E_total = {summary['total_energy_hartree']:.10f} Ha ({summary['total_energy_ev']:.6f} eV){label}")
    parts = summary['energy_components_hartree']
    print(f"  운동 {parts['kinetic']:.8f}, 전자-핵 {parts['electron_nuclear']:.8f}, 하트리 {parts['hartree']:.8f}, "
          f"교환·상관 {parts['exchange_correlation']:.8f} Ha")
    print(f"적분 전자 수 N = {summary['electron_count']:.12f}")
    for warning in summary['warnings']:
        print('확인:', warning)
    figure = create_atom_figure(result)
    plot_path = directory / 'atom.png'
    figure.savefig(plot_path, dpi=150)
    print(f"상태: {summary['status']} | 결과 폴더: {directory}")
    if show:
        import matplotlib.pyplot as plt
        import matplotlib.image as mpimg
        plt.figure(figsize=figure.get_size_inches())
        plt.imshow(mpimg.imread(plot_path))
        plt.axis('off')
        plt.show()
    return 0 if result.converged else 2


def run(config: Config, output_dir: Path, show: bool):
    started = time.perf_counter()
    result = solve(config)
    summary = save_result(result, output_dir, time.perf_counter() - started)
    print(f"QuantumStarter {VERSION} | 한 전자 원자의 방사형 슈뢰딩거 방정식")
    print("수소 검증 모드: 한 전자 방사방정식. 분자 DFT는 --hydrogen-check 없이 실행하세요.")
    print("state       numerical [eV]       exact [eV]       abs. error [eV]")
    for state in summary["states"]:
        print(f"{state['label']:5s}    {state['energy_ev']:16.8f}"
              f"    {state['exact_energy_hartree'] * HARTREE_TO_EV:14.8f}"
              f"    {state['abs_error_ev']:.8f}")
    for warning in summary["warnings"]:
        print("확인:", warning)
    print(f"결과 폴더: {Path(output_dir).resolve()}")
    save_plot(result, output_dir / "wavefunctions.png", show=show)
    return result, summary


def main():
    parser = argparse.ArgumentParser(description="원자 좌표·전하·스핀으로 계산하는 분자 DFT. 화면: app.py")
    parser.add_argument("--config", type=Path, help="분자 JSON 설정 (기본 molecule.json)")
    parser.add_argument("--xyz", type=Path, help="단일 구조 XYZ, 좌표 단위 Å")
    parser.add_argument('--structure', type=Path, help='주기 구조 POSCAR / CIF / extxyz')
    parser.add_argument('--boundary', choices=['bulk', 'slab'])
    parser.add_argument('--kmesh', nargs=3, type=int, metavar=('KX', 'KY', 'KZ'))
    parser.add_argument('--smearing-ev', type=float)
    parser.add_argument('--band-path', help="주기계 밴드 경로: auto 또는 GXWKGL 같은 고대칭점 (없으면 계산 안 함)")
    parser.add_argument('--band-points', type=int, help='밴드 경로 점 수 (기본 60)')
    parser.add_argument('--integral-scheme', choices=['gdf', 'fftdf'],
                        help='주기계 적분: gdf(기본, slab 가능) 또는 fftdf(bulk 전용, 힘·응력 가능)')
    parser.add_argument('--fft-mesh', nargs=3, type=int, metavar=('NX', 'NY', 'NZ'), help='FFTDF FFT 격자 (생략 시 자동)')
    parser.add_argument('--cell-relaxation', choices=['full', 'isotropic'], help='셀 모양+체적 또는 모양 고정 체적 완화')
    parser.add_argument('--cell-max-steps', type=int)
    parser.add_argument('--cell-force-tol', type=float, help='최대 원자 힘 [eV/Å]')
    parser.add_argument('--cell-stress-tol-gpa', type=float, help='셀 응력 잔차 [GPa]')
    parser.add_argument('--cell-energy-tol-ev', type=float, help='연속 단계 E/F 변화 [eV/cell]')
    parser.add_argument("--charge", type=int, help="전체 전하")
    parser.add_argument("--spin", type=int, help="2S = N_alpha - N_beta; 다중도 아님")
    parser.add_argument("--basis", help="기저함수 (원소별 지정이 없는 원소의 기본값)")
    parser.add_argument('--basis-element', action='append', metavar='EL=BASIS',
                        help='원소별 기저, 반복 가능: --basis-element O=def2-tzvp --basis-element H=sto-3g. '
                             '설정 파일의 원소별 지정을 대체합니다.')
    parser.add_argument("--xc", help="교환·상관 함수")
    parser.add_argument('--calculation', choices=['single_point', 'gradient', 'optimize', 'frequencies', 'scan',
                                                  'transition_state', 'neb', 'irc', 'cell_optimize', 'uvvis', 'gas_thermo', 'spin_scan'],
                        help='gradient: 에너지+원자 힘. optimize: geomeTRIC 구조 최적화(geometric 필요). '
                             'frequencies: 조화 진동수·ZPE')
    parser.add_argument('--freq-step', dest='freq_step_bohr', type=float, help='진동수 유한차분 변위 [bohr] (기본 0.005)')
    parser.add_argument('--states', dest='excited_states', type=int, help='uvvis 계산의 가장 낮은 단일항 들뜬 상태 수 (기본 6)')
    parser.add_argument('--temperature', dest='thermo_temperature_k', type=float, help='열역학 온도 [K]')
    parser.add_argument('--pressure', dest='thermo_pressure_pa', type=float, help='열역학 압력 [Pa]. 기본 101325, gas_thermo 표준 상태 100000.')
    parser.add_argument('--symmetry-number', dest='symmetry_number', type=int, help='회전 대칭수 (0=자동 판정)')
    parser.add_argument('--opt-convergence', choices=['default', 'tight'], help='구조 최적화 수렴 기준')
    parser.add_argument('--opt-max-steps', type=int, help='구조 최적화 최대 단계 (기본 100)')
    parser.add_argument('--constraints', dest='opt_constraints', help='geomeTRIC 제약. 예: "$freeze;xyz 1" 또는 "$set;angle 2 1 3 100"')
    parser.add_argument('--neb-product', type=Path, help='NEB 생성물 XYZ(반응물과 같은 원자 순서, Å)')
    parser.add_argument('--neb-images', dest='neb_images', type=int, help='NEB 이미지 수(양 끝 포함, 기본 9)')
    parser.add_argument('--scan', dest='scan_spec', help='scan 또는 spin_scan의 좌표. 예: "distance 1 2 0.7 1.2 6" (거리 Å, 원자 번호 1부터)')
    parser.add_argument('--dispersion', choices=['none', 'd4'],
                        help='분산 보정. d4: Grimme DFT-D4 (dftd4 필요, PBE/PBE0/B3LYP)')
    parser.add_argument("--scf-solver", choices=['diis', 'newton'], help="SCF 반복 방식 (기본 diis)")
    parser.add_argument("--ecp", choices=['auto', 'none', 'def2'], help="ECP 선택; auto는 무거운 원소의 def2 기저에 대응 ECP 적용")
    parser.add_argument('--relativity', choices=['auto', 'none', 'x2c'])
    parser.add_argument('--nuclear-model', choices=['auto', 'point', 'gaussian'])
    parser.add_argument('--grid-scheme', choices=['auto', 'treutler', 'universal'])
    parser.add_argument('--engine', choices=['pyscf', 'native'], help='SCF 구현: 외부 PySCF 또는 자체 엔진')
    parser.add_argument('--initial-guess', choices=['auto', 'core', 'minao', 'perturbed'])
    parser.add_argument('--stability', dest='native_stability', choices=['none', 'check', 'follow'],
                        help='자체 엔진 SCF 안정성 분석(실수 내부 회전). follow: 불안정하면 더 낮은 해 추적')
    parser.add_argument('--restart', type=Path, help='기존 실행 폴더 또는 checkpoint.npz')
    parser.add_argument("--backend", choices=['auto', 'local', 'wsl'], default='auto')
    parser.add_argument("--hydrogen-check", action='store_true', help="기존 수소형 한 전자 검증 계산")
    parser.add_argument('--model', choices=['he-lda', 'atom'],
                        help='he-lda: 중성 He LDA-SCF (균일 방사 격자, 0.9.5). atom: H–Kr 원자 LDA/LSDA '
                             '또는 PBE (로그 방사 격자, 구면 평균). 자체 XC 구현. --config는 그 모델의 settings.json')
    parser.add_argument('--element', help='[atom] 원소 기호 H–Kr')
    parser.add_argument('--functional', choices=['lda', 'pbe'], help='[atom] XC 함수 (기본 lda: Slater+VWN5)')
    parser.add_argument('--configuration', help='[atom] 전자 배치, 예: "[Ne] 3s2 3p6" (생략 시 NIST 바닥 배치)')
    parser.add_argument('--spin-polarized', action='store_const', const=True,
                        help='[atom] 스핀 분극(껍질별 다수 스핀 우선 점유). 생략 시 스핀 비분극')
    parser.add_argument('--grid-step', type=float, help='[atom] 로그 격자 x 간격 h (기본 0.001)')
    parser.add_argument('--grid-scale', type=float, help='[atom] 격자 척도 a·Z [bohr] (기본 1e-4)')
    parser.add_argument('--mixing', type=float, help='[he-lda/atom] 밀도 혼합 감쇠 계수 0<a≤1 (기본 0.5 / 0.3)')
    parser.add_argument('--mixing-method', choices=('linear', 'anderson'), help='[atom] 밀도 혼합 방법 (기본 linear)')
    parser.add_argument('--max-iterations', type=int, help='[he-lda/atom] SCF 반복 상한 (기본 100 / 300)')
    parser.add_argument('--energy-tol', type=float, help='[he-lda/atom] |ΔE| 수렴 기준 [Ha] (기본 1e-10)')
    parser.add_argument('--density-tol', type=float, help='[he-lda/atom] 혼합 전 밀도 잔차 기준 [전자] (기본 1e-7)')
    parser.add_argument('--r-max', type=float, help='[he-lda/atom] 계산영역 반지름 [bohr] (기본 30 / 50)')
    parser.add_argument('--grid-points', type=int, help='[he-lda] 경계 제외 방사 격자점 수 (기본 29999, dr=0.001)')
    parser.add_argument("--output", type=Path, help="실행 폴더들을 저장할 상위 경로")
    parser.add_argument("--no-show", action="store_true", help="그래프 창 없이 파일만 저장")
    args = parser.parse_args()
    if args.excited_states is not None and (args.model or args.hydrogen_check):
        parser.error('--states는 분자 uvvis 계산에서 사용합니다.')
    periodic_keys = ('boundary', 'kmesh', 'smearing_ev', 'band_path', 'band_points', 'integral_scheme', 'fft_mesh',
                     'cell_relaxation', 'cell_max_steps', 'cell_force_tol', 'cell_stress_tol_gpa', 'cell_energy_tol_ev')
    has_periodic_options = any(getattr(args, key) is not None for key in periodic_keys)
    try:
        if args.model in ('he-lda', 'atom'):
            from dataclasses import replace
            used = [key for key in DFT_ONLY_OPTIONS + periodic_keys if getattr(args, key) not in (None, False)]
            if args.hydrogen_check or used or args.backend != 'auto':
                parser.error(f'--model {args.model}에는 --hydrogen-check·분자·주기계·--backend 옵션을 함께 사용할 수 없습니다: '
                             + ', '.join(used or ['hydrogen_check' if args.hydrogen_check else 'backend']))
            if args.model == 'he-lda':
                from atomic_dft import AtomicDFTConfig, load_atomic_config
                wrong = [key for key in ATOM_ONLY_OPTIONS if getattr(args, key) is not None]
                if wrong:
                    parser.error('--model he-lda에서는 쓸 수 없는 원자 옵션: ' + ', '.join(wrong)
                                 + ' (임의 원자는 --model atom)')
                config = load_atomic_config(args.config) if args.config else AtomicDFTConfig()
                config = replace(config, **{field: getattr(args, key) for key, field in ATOMIC_OPTIONS.items()
                                            if getattr(args, key) is not None})
                return run_he_lda(config, args.output or ROOT / 'outputs' / 'he_lda', show=not args.no_show)
            from radial_atom import AtomConfig, load_atom_config
            if args.grid_points is not None:
                parser.error('--model atom은 --grid-points 대신 --grid-step(로그 격자 간격)을 씁니다.')
            if not args.config and args.element is None:
                parser.error('--model atom에는 --element 또는 --config가 필요합니다.')
            config = load_atom_config(args.config) if args.config else AtomConfig(element=args.element)
            config = replace(config, **{field: getattr(args, key) for key, field in ATOM_OPTIONS.items()
                                        if getattr(args, key) is not None})
            return run_atom_cli(config, args.output or ROOT / 'outputs' / 'atom', show=not args.no_show)
        if any(getattr(args, key) is not None for key in ATOMIC_OPTIONS) or \
                any(getattr(args, key) is not None for key in ATOM_ONLY_OPTIONS):
            parser.error('--mixing·--max-iterations·--energy-tol·--density-tol·--r-max·--grid-points와 '
                         '--element·--configuration·--spin-polarized·--grid-step·--grid-scale·--functional·--mixing-method는 '
                         '--model he-lda / atom 전용입니다.')
        if args.hydrogen_check:
            if has_periodic_options:
                parser.error('--hydrogen-check에는 주기계 옵션을 함께 사용할 수 없습니다.')
            if args.xyz or args.structure or args.boundary or args.kmesh or args.smearing_ev is not None or args.engine or args.initial_guess or args.restart or args.basis_element or any(getattr(args, key) is not None for key in ('charge', 'spin', 'basis', 'xc', 'dispersion', 'calculation', 'opt_convergence', 'opt_max_steps', 'freq_step_bohr', 'thermo_temperature_k', 'thermo_pressure_pa', 'symmetry_number', 'native_stability', 'opt_constraints', 'scan_spec', 'neb_images', 'scf_solver', 'ecp', 'relativity', 'nuclear_model', 'grid_scheme')):
                parser.error('--hydrogen-check에는 분자 DFT 옵션을 함께 사용할 수 없습니다.')
            run(load_config(args.config or ROOT / 'config.json'), args.output or ROOT / 'outputs' / 'hydrogen', show=not args.no_show)
        else:
            from dataclasses import replace
            from molecular import coordinates_text, load_dft_config, parse_basis_overrides, read_xyz
            from dft_workflow import run_calculation
            from molecular_plotting import create_molecular_figure
            from periodic import load_calculation_config, read_periodic_structure
            if args.structure and (args.config or args.xyz):
                parser.error('--structure는 --config / --xyz와 함께 사용할 수 없습니다.')
            config = (read_periodic_structure(args.structure, **({'boundary': args.boundary} if args.boundary else {}))
                      if args.structure else load_calculation_config(args.config or ROOT / 'molecule.json'))
            overrides = {key: getattr(args, key) for key in ('charge', 'spin', 'basis', 'xc', 'dispersion', 'calculation', 'opt_convergence', 'opt_max_steps', 'freq_step_bohr', 'thermo_temperature_k', 'thermo_pressure_pa', 'symmetry_number', 'opt_constraints', 'scan_spec', 'neb_images', 'scf_solver', 'ecp', 'relativity', 'nuclear_model', 'grid_scheme', 'engine', 'initial_guess', 'native_stability') if getattr(args, key) is not None}
            if config.system == 'periodic':
                if args.xyz or args.ecp or args.scf_solver or args.relativity or args.nuclear_model or args.grid_scheme or args.engine or args.initial_guess or args.restart or args.basis_element or args.dispersion or args.calculation not in (None, 'single_point', 'gradient', 'cell_optimize') or args.opt_convergence or args.opt_max_steps or args.freq_step_bohr is not None or args.thermo_temperature_k is not None or args.thermo_pressure_pa is not None or args.symmetry_number is not None or args.native_stability or args.opt_constraints or args.scan_spec or args.neb_product or args.neb_images is not None:
                    parser.error('주기계는 single_point·gradient·cell_optimize만 지원하며 분자 전용 옵션은 사용할 수 없습니다.')
                overrides.update({key: getattr(args, key) for key in periodic_keys if getattr(args, key) is not None})
            elif has_periodic_options or args.calculation == 'cell_optimize':
                parser.error('주기계 옵션에는 --structure 또는 주기계 --config가 필요합니다.')
            if args.xyz:
                overrides.update(atoms=read_xyz(args.xyz), unit='angstrom')
            if args.basis_element:
                overrides['basis_overrides'] = parse_basis_overrides(','.join(args.basis_element))
            if args.neb_product:
                overrides['neb_product'] = coordinates_text(read_xyz(args.neb_product))
            if args.excited_states is not None:
                if config.system != 'molecule' or (args.calculation or config.calculation) != 'uvvis':
                    parser.error('--states에는 uvvis 계산이 필요합니다.')
                overrides['excited_states'] = args.excited_states
            if args.calculation in ('uvvis', 'spin_scan', 'gas_thermo') and not args.engine and not args.config:
                overrides['engine'] = 'native'
            if args.calculation == 'gas_thermo' and args.thermo_pressure_pa is None and not args.config:
                overrides['thermo_pressure_pa'] = 100000
            config = replace(config, **overrides)
            summary, directory = run_calculation(config, args.output, args.backend, print, restart_from=args.restart)
            if summary['status'] == 'cancelled':
                print(f'계산 중단 · 결과 폴더: {directory}')
                return 3
            if config.system == 'periodic':
                from periodic_plotting import create_periodic_figure
                figure = create_periodic_figure(summary)
            else:
                figure = create_molecular_figure(summary, directory)
            plot_path = directory / ('periodic.png' if config.system == 'periodic' else 'molecular.png')
            figure.savefig(plot_path, dpi=150)
            print(f"총에너지: {summary['total_energy_hartree']:.12f} Ha")
            if 'uvvis' in summary:
                print('상태    에너지 (eV)    파장 (nm)    진동자 세기    주요 전이 (MO)')
                for p in summary['uvvis']['peaks']:
                    mo = p['leading_transitions'][0]
                    print(f"{p['state']:4d} {p['energy_eV']:13.5f} {p['wavelength_nm']:12.3f} "
                          f"{p['oscillator_strength']:12.5g}    {mo['occupied']} → {mo['virtual']}")
            if 'spin_scan' in summary:
                print('스핀 곡률 부호 전환 구간 (Å):', summary['spin_scan']['brackets_angstrom'])
            if 'gas_thermo' in summary:
                print('열물성 표·NASA7 계수:', directory / 'gas_thermo.csv', directory / 'nasa7.yaml')
            print(f"상태: {summary['status']} | 결과 폴더: {directory}")
            for warning in summary['warnings']:
                print('확인:', warning)
            if not args.no_show:
                # A pyplot manager is only needed for an interactive CLI window.
                import matplotlib.pyplot as plt
                import matplotlib.image as mpimg
                plt.figure(figsize=figure.get_size_inches())
                plt.imshow(mpimg.imread(plot_path))
                plt.axis('off')
                plt.show()
            return 0 if summary['converged'] else 2
    except (ValueError, OSError, RuntimeError, np.linalg.LinAlgError) as exc:
        print(f"실행 오류: {exc}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    for stream in (sys.stdout, sys.stderr):
        if hasattr(stream, 'reconfigure'):
            stream.reconfigure(encoding='utf-8', errors='replace')
    raise SystemExit(main())
