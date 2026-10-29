"""Backend transport, immutable run directories and common DFT result storage."""
from datetime import datetime, timezone
import importlib.util
import json
import os
from pathlib import Path
import subprocess
import sys
import shutil
import threading
from tempfile import mkdtemp

from molecular import DFTConfig, save_dft_config
from periodic import PeriodicConfig, config_from_dict

ROOT = Path(__file__).resolve().parent


class RunControl:
    """Thread-safe, file-based cancellation also works across the WSL boundary."""
    def __init__(self):
        self.lock = threading.Lock()
        self.directory = None
        self.requested = False

    def bind(self, directory):
        with self.lock:
            self.directory = Path(directory)
            if self.requested:
                (self.directory / 'cancel.request').touch()

    def cancel(self):
        with self.lock:
            self.requested = True
            if self.directory is not None:
                (self.directory / 'cancel.request').touch()


class CalculationError(RuntimeError):
    def __init__(self, message, directory):
        super().__init__(message)
        self.directory = directory


def process_options():
    return {'creationflags': subprocess.CREATE_NO_WINDOW} if sys.platform == 'win32' else {}


def _wsl_path(path):
    result = subprocess.run(['wsl', '--exec', 'wslpath', '-a', str(Path(path).resolve())],
                            capture_output=True, encoding='utf-8', errors='replace',
                            timeout=30, **process_options())
    if result.returncode:
        raise RuntimeError('WSL에 접근할 수 없습니다. START_HERE_KO.md의 설치 안내를 확인하세요.\n'
                           + result.stderr.replace('\x00', '').strip())
    return result.stdout.strip()


def worker_command(config_path, directory, backend='auto'):
    if backend not in ('auto', 'local', 'wsl'):
        raise ValueError('실행 환경은 auto, local, wsl 중에서 선택하세요.')
    if backend == 'auto':
        backend = 'local' if importlib.util.find_spec('pyscf') else ('wsl' if sys.platform == 'win32' else 'local')
    if backend == 'local':
        if not importlib.util.find_spec('pyscf'):
            raise RuntimeError('현재 Python에 PySCF가 없습니다. pip install -r requirements-dft.txt를 실행하세요.')
        return [sys.executable, '-u', str(ROOT / 'dft_worker.py'), '--config', str(config_path),
                '--output', str(directory)]
    if sys.platform != 'win32':
        raise RuntimeError('WSL 실행은 Windows에서만 사용하세요. Linux/macOS에서는 local을 선택하세요.')
    root_linux = _wsl_path(ROOT)
    # Only invoke an executable and arguments, never interpolate input into a shell command.
    python = root_linux + '/.venv-dft/bin/python'
    if not (ROOT / '.venv-dft' / 'bin').is_dir():
        raise RuntimeError('DFT용 WSL 환경이 없습니다. setup_dft_windows.bat를 먼저 실행하세요.')
    return ['wsl', '--exec', 'env', 'PYTHONIOENCODING=utf-8', python, '-u',
            root_linux + '/dft_worker.py', '--config', _wsl_path(config_path),
            '--output', _wsl_path(directory)]


def run_calculation(config, output_root=None, backend='auto', progress=lambda message: None,
                    control=None, restart_from=None):
    """Return (summary, run directory). Nonconvergence is a result, not success."""
    if not isinstance(config, (DFTConfig, PeriodicConfig)):
        raise ValueError('DFTConfig 또는 PeriodicConfig가 필요합니다.')
    output_root = Path(output_root or ROOT / 'outputs' / ('periodic' if config.system == 'periodic' else 'molecular')).resolve()
    output_root.mkdir(parents=True, exist_ok=True)
    directory = Path(mkdtemp(prefix=datetime.now().strftime('%Y%m%d_%H%M%S_'), dir=output_root))
    save_dft_config(config, directory / 'settings.json')
    lines = []
    try:
        if (control is not None or restart_from is not None) and config.system != 'molecule':
            raise ValueError('중단/SCF 재시작은 현재 분자 엔진에서 지원합니다.')
        if restart_from is not None:
            source = Path(restart_from)
            if source.is_dir():
                source = source / 'checkpoint.npz'
            shutil.copyfile(source, directory / 'restart_input.npz')
            (directory / 'restart_source.json').write_text(json.dumps({'path': str(source.resolve())}), encoding='utf-8')
        if control is not None:
            control.bind(directory)
        command = worker_command(directory / 'settings.json', directory, backend)
        env = dict(os.environ, PYTHONIOENCODING='utf-8')
        with (directory / 'execution.log').open('w', encoding='utf-8') as log:
            with subprocess.Popen(command, cwd=ROOT, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                                  encoding='utf-8', errors='replace', bufsize=1, env=env,
                                  **process_options()) as process:
                for line in process.stdout:
                    line = line.rstrip()
                    log.write(line + '\n')
                    log.flush()
                    lines.append(line)
                    progress(line)
                code = process.wait()
        if code not in (0, 2, 3):
            raise RuntimeError('DFT 엔진 실행 실패.\n' + '\n'.join(lines[-14:]))
        summary = json.loads((directory / 'results.json').read_text(encoding='utf-8'))
        if (summary.get('status') != {0: 'converged', 2: 'not_converged', 3: 'cancelled'}[code]
                or summary.get('converged') is not (code == 0)):
            raise RuntimeError('DFT 엔진 종료 코드와 수렴 상태가 일치하지 않습니다.')
        return summary, directory
    except Exception as exc:
        data = {'schema_version': 1, 'status': 'error', 'converged': False,
                'created_at_utc': datetime.now(timezone.utc).isoformat(),
                'config': config.to_dict(), 'formula': config.formula, 'error': str(exc)}
        (directory / 'error.json').write_text(json.dumps(data, ensure_ascii=False, indent=2), encoding='utf-8')
        with (directory / 'execution.log').open('a', encoding='utf-8') as log:
            log.write('\n' + str(exc) + '\n')
        raise CalculationError(str(exc), directory) from exc


def list_dft_runs(root, limit=100):
    runs = []
    for directory in sorted(Path(root).glob('*'), reverse=True):
        if not directory.is_dir():
            continue
        path = directory / 'error.json' if (directory / 'error.json').exists() else directory / 'results.json'
        try:
            summary = json.loads(path.read_text(encoding='utf-8'))
            config_from_dict(summary['config'])
            if summary['status'] not in ('converged', 'not_converged', 'cancelled', 'error'):
                continue
            runs.append({'directory': directory, 'summary': summary})
        except (OSError, ValueError, KeyError, TypeError):
            continue
        if len(runs) >= limit:
            break
    return runs
