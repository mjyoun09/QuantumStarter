"""Run explicitly selected unittest modules, preserving complete output and exit code."""
import argparse
import json
import os
from pathlib import Path
import subprocess
import sys
import time

ROOT = Path(__file__).resolve().parents[1]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--label', required=True)
    parser.add_argument('tests', nargs='+')
    args = parser.parse_args()
    directory = ROOT / 'outputs' / 'cell_validation' / args.label
    directory.mkdir(parents=True, exist_ok=False)
    command = [sys.executable, '-u', '-m', 'unittest', *args.tests, '-v']
    started = time.perf_counter()
    with (directory / 'tests.log').open('w', encoding='utf-8', newline='\n') as stream:
        result = subprocess.run(command, cwd=ROOT, stdout=stream, stderr=subprocess.STDOUT,
                                env=dict(os.environ, PYTHONIOENCODING='utf-8'))
    report = {'command': command, 'exit_code': result.returncode,
              'elapsed_seconds': time.perf_counter() - started,
              'environment': {key: os.environ.get(key) for key in ('QUANTUMSTARTER_TK_HIDDEN', 'QS_PERIODIC_GUI_DFT')}}
    (directory / 'tests.json').write_text(json.dumps(report, indent=2), encoding='utf-8')
    print(json.dumps(report), flush=True)
    print((directory / 'tests.log').read_text(encoding='utf-8')[-4000:], flush=True)
    return result.returncode


if __name__ == '__main__':
    raise SystemExit(main())
