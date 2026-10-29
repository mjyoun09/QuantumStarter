"""명령행 실행 유지: python main.py --no-show. 사용 화면: python app.py."""
import argparse
from pathlib import Path
import sys
import time
import numpy as np
from plotting import save_plot
from solver import HARTREE_TO_EV, Config, solve
from workflow import VERSION, load_config, save_result

ROOT = Path(__file__).resolve().parent


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
    parser.add_argument("--charge", type=int, help="전체 전하")
    parser.add_argument("--spin", type=int, help="2S = N_alpha - N_beta; 다중도 아님")
    parser.add_argument("--basis", help="기저함수 (원소별 지정이 없는 원소의 기본값)")
    parser.add_argument('--basis-element', action='append', metavar='EL=BASIS',
                        help='원소별 기저, 반복 가능: --basis-element O=def2-tzvp --basis-element H=sto-3g. '
                             '설정 파일의 원소별 지정을 대체합니다.')
    parser.add_argument("--xc", help="교환·상관 함수")
    parser.add_argument('--calculation', choices=['single_point', 'gradient', 'optimize', 'frequencies'],
                        help='gradient: 에너지+원자 힘. optimize: geomeTRIC 구조 최적화(geometric 필요). '
                             'frequencies: 조화 진동수·ZPE')
    parser.add_argument('--freq-step', dest='freq_step_bohr', type=float, help='진동수 유한차분 변위 [bohr] (기본 0.005)')
    parser.add_argument('--temperature', dest='thermo_temperature_k', type=float, help='열역학 온도 [K]')
    parser.add_argument('--pressure', dest='thermo_pressure_pa', type=float, help='열역학 압력 [Pa] (기본 101325)')
    parser.add_argument('--symmetry-number', dest='symmetry_number', type=int, help='회전 대칭수 (0=자동 판정)')
    parser.add_argument('--opt-convergence', choices=['default', 'tight'], help='구조 최적화 수렴 기준')
    parser.add_argument('--opt-max-steps', type=int, help='구조 최적화 최대 단계 (기본 100)')
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
    parser.add_argument("--output", type=Path, help="실행 폴더들을 저장할 상위 경로")
    parser.add_argument("--no-show", action="store_true", help="그래프 창 없이 파일만 저장")
    args = parser.parse_args()
    try:
        if args.hydrogen_check:
            if args.xyz or args.structure or args.boundary or args.kmesh or args.smearing_ev is not None or args.engine or args.initial_guess or args.restart or args.basis_element or any(getattr(args, key) is not None for key in ('charge', 'spin', 'basis', 'xc', 'dispersion', 'calculation', 'opt_convergence', 'opt_max_steps', 'freq_step_bohr', 'thermo_temperature_k', 'thermo_pressure_pa', 'symmetry_number', 'native_stability', 'scf_solver', 'ecp', 'relativity', 'nuclear_model', 'grid_scheme')):
                parser.error('--hydrogen-check에는 분자 DFT 옵션을 함께 사용할 수 없습니다.')
            run(load_config(args.config or ROOT / 'config.json'), args.output or ROOT / 'outputs' / 'hydrogen', show=not args.no_show)
        else:
            from dataclasses import replace
            from molecular import load_dft_config, parse_basis_overrides, read_xyz
            from dft_workflow import run_calculation
            from molecular_plotting import create_molecular_figure
            from periodic import load_calculation_config, read_periodic_structure
            if args.structure and (args.config or args.xyz):
                parser.error('--structure는 --config / --xyz와 함께 사용할 수 없습니다.')
            config = (read_periodic_structure(args.structure, **({'boundary': args.boundary} if args.boundary else {}))
                      if args.structure else load_calculation_config(args.config or ROOT / 'molecule.json'))
            overrides = {key: getattr(args, key) for key in ('charge', 'spin', 'basis', 'xc', 'dispersion', 'calculation', 'opt_convergence', 'opt_max_steps', 'freq_step_bohr', 'thermo_temperature_k', 'thermo_pressure_pa', 'symmetry_number', 'scf_solver', 'ecp', 'relativity', 'nuclear_model', 'grid_scheme', 'engine', 'initial_guess', 'native_stability') if getattr(args, key) is not None}
            if config.system == 'periodic':
                if args.xyz or args.ecp or args.scf_solver or args.relativity or args.nuclear_model or args.grid_scheme or args.engine or args.initial_guess or args.restart or args.basis_element or args.dispersion or args.calculation or args.opt_convergence or args.opt_max_steps or args.freq_step_bohr is not None or args.thermo_temperature_k is not None or args.thermo_pressure_pa is not None or args.symmetry_number is not None or args.native_stability:
                    parser.error('주기계는 --xyz / --ecp / --scf-solver / --relativity / --nuclear-model / --grid-scheme / --engine / --initial-guess / --restart / --basis-element / --dispersion / --calculation을 지원하지 않습니다.')
                overrides.update({key: getattr(args, key) for key in ('boundary', 'kmesh', 'smearing_ev') if getattr(args, key) is not None})
            elif args.boundary or args.kmesh or args.smearing_ev is not None:
                parser.error('주기계 옵션에는 --structure 또는 주기계 --config가 필요합니다.')
            if args.xyz:
                overrides.update(atoms=read_xyz(args.xyz), unit='angstrom')
            if args.basis_element:
                overrides['basis_overrides'] = parse_basis_overrides(','.join(args.basis_element))
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
