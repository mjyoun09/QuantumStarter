"""Process boundary for native Linux/macOS and Windows -> WSL execution."""
import argparse
import sys

from engines import calculate
from periodic import load_calculation_config
from scf_control import CalculationCancelled


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--config', required=True)
    parser.add_argument('--output', required=True)
    args = parser.parse_args()
    try:
        summary = calculate(load_calculation_config(args.config), args.output,
                            lambda message: print(message, flush=True))
        return 0 if summary['converged'] else 2
    except CalculationCancelled as exc:
        import json
        from pathlib import Path
        config = load_calculation_config(args.config)
        directory = Path(args.output)
        result = {'status': 'cancelled', 'converged': False, 'config': config.to_dict(),
                  'formula': config.formula, 'engine': config.engine if config.system == 'molecule' else 'pyscf-pbc',
                  'message': str(exc),
                  'checkpoint_file': 'checkpoint.npz' if (directory / 'checkpoint.npz').exists() else None}
        (directory / 'results.json').write_text(json.dumps(result, ensure_ascii=False, indent=2), encoding='utf-8')
        print(str(exc), flush=True)
        return 3
    except Exception as exc:
        print(f'계산 오류: {exc}', file=sys.stderr, flush=True)
        import traceback
        traceback.print_exc()
        return 1


if __name__ == '__main__':
    raise SystemExit(main())
