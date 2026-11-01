"""GUI와 명령행이 함께 쓰는 설정·결과 저장 계층. Tkinter에 의존하지 않는다."""
import csv
from dataclasses import asdict
from datetime import datetime, timezone
import json
from pathlib import Path
import platform

import numpy as np
import scipy

from solver import BOHR_TO_ANGSTROM, HARTREE_TO_EV, Config, Result, quality_warnings

VERSION = "0.9.3"
ORBITALS = ("s", "p", "d", "f", "g", "h")
PRESETS = {"빠른 확인": 2000, "권장": 6000, "정밀": 12000}
FIELD_NAMES = {
    "nuclear_charge": "핵전하 Z", "angular_momentum_l": "방위 양자수 l",
    "num_states": "계산할 상태 수", "r_max_bohr": "계산영역 반지름", "grid_points": "격자점 수",
}


def load_config(path: Path) -> Config:
    data = json.loads(Path(path).read_text(encoding="utf-8-sig"))
    if not isinstance(data, dict):
        raise ValueError("설정 파일에는 JSON 객체가 필요합니다.")
    unknown = set(data) - set(Config.__dataclass_fields__)
    if unknown:
        raise ValueError("알 수 없는 설정: " + ", ".join(sorted(unknown)))
    return Config(**data)


def save_config(config: Config, path: Path):
    Path(path).write_text(json.dumps(asdict(config), ensure_ascii=False, indent=2), encoding="utf-8")


def config_from_fields(fields: dict) -> Config:
    data = {}
    for key in Config.__dataclass_fields__:
        text = str(fields.get(key, "")).strip()
        try:
            data[key] = float(text) if key == "r_max_bohr" else int(text)
        except (ValueError, TypeError):
            kind = "숫자" if key == "r_max_bohr" else "정수"
            raise ValueError(f"{FIELD_NAMES[key]}에 {kind}를 입력해 주세요.") from None
    try:
        cfg = Config(**data)
    except ValueError as exc:
        message = str(exc)
        for key, label in FIELD_NAMES.items():
            message = message.replace(key, label)
        raise ValueError(message) from None
    if cfg.r_max_bohr / (cfg.grid_points + 1) * cfg.nuclear_charge < 1e-6:
        raise ValueError("격자 간격이 지나치게 작습니다. 계산영역을 늘리거나 권장 설정을 선택하세요.")
    return cfg


def recommended_config(z=1, ell=0, num_states=3, precision="권장") -> Config:
    max_n = ell + num_states
    radius = float(np.ceil(max(30, 6 * max_n ** 2) / 10) * 10 / z)
    return Config(z, ell, num_states, radius, PRESETS[precision])


def state_label(n, ell):
    return f"{int(n)}{ORBITALS[int(ell)]}"


def build_summary(result: Result, elapsed_seconds=None) -> dict:
    errors = np.abs(result.energies_hartree - result.exact_energies_hartree) * HARTREE_TO_EV
    return {
        "program": "QuantumStarter", "version": VERSION,
        "created_at_utc": datetime.now(timezone.utc).isoformat(),
        "model": "one-electron Coulomb radial Schrodinger equation; NOT DFT",
        "assumptions": ["fixed infinite-mass point nucleus", "nonrelativistic",
                        "isolated single electron", "Dirichlet radial boundaries"],
        "method": "second-order finite difference; tridiagonal eigensolver",
        "config": asdict(result.config), "dr_bohr": result.dr,
        "elapsed_seconds": elapsed_seconds,
        "units": {"energy": "hartree / eV", "distance": "bohr",
                  "reduced_radial_u": "bohr^(-1/2)", "radial_probability": "bohr^(-1)"},
        "constants": {"hartree_to_ev": HARTREE_TO_EV, "bohr_to_angstrom": BOHR_TO_ANGSTROM},
        "environment": {"python": platform.python_version(),
                        "numpy": np.__version__, "scipy": scipy.__version__},
        "warnings": quality_warnings(result),
        "states": [
            {"n": int(n), "l": result.config.angular_momentum_l,
             "label": state_label(n, result.config.angular_momentum_l),
             "energy_hartree": float(energy), "energy_ev": float(energy * HARTREE_TO_EV),
             "exact_energy_hartree": float(reference), "abs_error_ev": float(error),
             "radial_normalization": float(result.dr * np.sum(result.u[:, j] ** 2)),
             "mean_r_bohr": float(result.mean_r_bohr[j])}
            for j, (n, energy, reference, error) in enumerate(zip(
                result.principal_n, result.energies_hartree,
                result.exact_energies_hartree, errors))
        ],
    }


def save_result(result: Result, output_dir: Path, elapsed_seconds=None) -> dict:
    """원시 데이터와 재실행 설정 저장. 그래프는 UI/CLI의 주 스레드에서 저장한다."""
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    summary = build_summary(result, elapsed_seconds)
    save_config(result.config, output_dir / "settings.json")
    with (output_dir / "energies.csv").open("w", newline="", encoding="utf-8-sig") as stream:
        writer = csv.writer(stream)
        writer.writerow(["state", "n", "l", "energy_hartree", "energy_ev", "exact_energy_ev",
                         "abs_error_ev", "normalization", "mean_r_bohr"])
        for state in summary["states"]:
            writer.writerow([state["label"], state["n"], state["l"], state["energy_hartree"],
                             state["energy_ev"], state["exact_energy_hartree"] * HARTREE_TO_EV,
                             state["abs_error_ev"], state["radial_normalization"], state["mean_r_bohr"]])
    columns, headers = [result.r_bohr], ["r_bohr"]
    for j, n in enumerate(result.principal_n):
        columns.extend([result.u[:, j], result.radial_probability[:, j], result.u_exact[:, j]])
        suffix = f"n{n}_l{result.config.angular_momentum_l}"
        headers.extend([f"u_{suffix}", f"P_{suffix}", f"u_exact_{suffix}"])
    np.savetxt(output_dir / "radial.csv", np.column_stack(columns), delimiter=",",
               header=",".join(headers), comments="")
    temporary = output_dir / "results.json.tmp"
    temporary.write_text(json.dumps(summary, ensure_ascii=False, indent=2), encoding="utf-8")
    temporary.replace(output_dir / "results.json")
    return summary


def new_run_dir(output_root: Path) -> Path:
    from tempfile import mkdtemp
    output_root = Path(output_root)
    output_root.mkdir(parents=True, exist_ok=True)
    # Microseconds keep name order = run order; several fast runs can start within one second.
    prefix = datetime.now().strftime("%Y%m%d_%H%M%S_%f_")
    return Path(mkdtemp(prefix=prefix, dir=output_root))


def list_runs(output_root: Path, limit=100) -> list[dict]:
    runs = []
    for path in sorted(Path(output_root).glob("*/results.json"), reverse=True):
        try:
            data = json.loads(path.read_text(encoding="utf-8"))
            Config(**data["config"])
            if not isinstance(data["states"], list) or not data["states"]:
                continue
            if not all(isinstance(state, dict) and np.isfinite(float(state['energy_ev']))
                       for state in data['states']):
                continue
            runs.append({"directory": path.parent, "summary": data})
        except (OSError, ValueError, TypeError, KeyError):
            continue
        if len(runs) >= limit:
            break
    return runs
