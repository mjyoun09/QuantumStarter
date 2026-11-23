"""Harmonic dipole intensities and area-normalized line shapes."""
import numpy as np
from scipy.constants import Avogadro, atomic_mass, elementary_charge, epsilon_0, speed_of_light


# Isotropic integrated molar absorption for |dμ/dQ|² in e²/u, converted from m/mol.
KM_MOL_PER_E2_AMU = Avogadro * elementary_charge**2 / (12 * epsilon_0 * speed_of_light**2 * atomic_mass * 1000)
DEFAULT_FWHM_CM1 = 30.0


def mode_intensities(dipole_derivative, analysis):
    """Project a (3, 3N) Cartesian dipole derivative in e onto mass-normalized modes."""
    modes = analysis['mass_normalized_modes']
    try:
        tensor = np.asarray(dipole_derivative, dtype=float)
    except (TypeError, ValueError) as exc:
        raise ValueError('IR 쌍극자 미분을 실수 배열로 읽을 수 없습니다.') from exc
    if tensor.shape != (3, 3 * modes.shape[1]) or not np.isfinite(tensor).all():
        raise ValueError('IR 쌍극자 미분은 유한한 3 × 3N 배열이어야 합니다.')
    slopes = modes.reshape(len(modes), tensor.shape[1]) @ tensor.T
    strengths = KM_MOL_PER_E2_AMU * np.sum(slopes**2, axis=1)
    return [float(s) if f > 0 else None for f, s in zip(analysis['frequencies_cm1'], strengths)]


def validate_fwhm(fwhm):
    try:
        fwhm = float(fwhm)
    except (TypeError, ValueError):
        raise ValueError('가우스 선폭 FWHM(cm⁻¹)을 숫자로 입력하세요.') from None
    if not np.isfinite(fwhm) or fwhm <= 0:
        raise ValueError('가우스 선폭 FWHM(cm⁻¹)은 유한한 양수여야 합니다.')
    return fwhm


def broaden_spectrum(bands, fwhm):
    """Return the positive wavenumber grid and curves. Each band is (wavenumber, km/mol)."""
    fwhm = validate_fwhm(fwhm)
    sigma = fwhm / np.sqrt(8 * np.log(2))
    centers = [freq for series in bands for freq, area in series if freq > 0 and area is not None]
    if not centers:
        return np.array([]), np.zeros((len(bands), 0))
    # Sample around every line so narrow, separated peaks are resolved without a huge uniform grid.
    if sigma <= np.spacing(max(centers)) or not np.isfinite(max(centers) + 8 * sigma):
        raise ValueError('이 파수 범위에서 수치로 표현할 수 있는 선폭을 입력하세요.')
    lo, hi = max(0., min(centers) - 8 * sigma), max(centers) + 8 * sigma
    offsets = np.linspace(-8, 8, 161) * sigma
    grid = np.unique(np.concatenate([np.linspace(lo, hi, 1001), *(f + offsets for f in centers)]))
    grid = grid[grid >= 0]
    curves = np.zeros((len(bands), len(grid)))
    for curve, series in zip(curves, bands):
        for freq, area in series:
            if freq > 0 and area is not None:
                start, stop = np.searchsorted(grid, [freq - 8 * sigma, freq + 8 * sigma])
                z = (grid[start:stop] - freq) / sigma
                curve[start:stop] += (area / (sigma * np.sqrt(2 * np.pi))) * np.exp(-0.5 * z**2)
    return grid, curves
