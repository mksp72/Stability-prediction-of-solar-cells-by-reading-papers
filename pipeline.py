"""Minimal auditable data validation and auxiliary ISOS-L fitting."""
from dataclasses import dataclass
import numpy as np
import pandas as pd
from scipy.optimize import least_squares, brentq

REQUIRED_DEVICES = {'device_id', 'doi', 'pce0', 'voc0', 'jsc0', 'ff0'}
REQUIRED_CURVES = {'device_id', 'protocol', 'time_h', 'retention_pct'}


def validate(devices: pd.DataFrame, curves: pd.DataFrame):
    for name, frame, cols in [('devices', devices, REQUIRED_DEVICES), ('curves', curves, REQUIRED_CURVES)]:
        missing = cols - set(frame.columns)
        if missing:
            raise ValueError(f'{name}: missing columns {sorted(missing)}')
    if devices.device_id.isna().any() or devices.device_id.duplicated().any():
        raise ValueError('Device identifiers must be non-null and unique')
    if curves[list(REQUIRED_CURVES)].isna().any().any():
        raise ValueError('Curve required fields cannot be null')
    if not set(curves.protocol.unique()).issubset({'ISOS-D', 'ISOS-L'}):
        raise ValueError('Unexpected stability protocol')
    if not set(curves.device_id).issubset(set(devices.device_id)):
        raise ValueError('Orphan curve device IDs')
    for (_, _), part in curves.groupby(['device_id', 'protocol']):
        t = part.time_h.to_numpy(dtype=float)
        r = part.retention_pct.to_numpy(dtype=float)
        if (t < 0).any() or not np.isfinite(t).all() or not np.isfinite(r).all() or (r < 0).any():
            raise ValueError('Invalid curve time/retention')
        if len(np.unique(t)) != len(t):
            raise ValueError('Duplicate time points in one curve')
    return True


def observed_t80(time_h, retention_pct):
    """First interpolated crossing of 80% initial-reference retention, if observed."""
    order = np.argsort(time_h)
    t = np.asarray(time_h, float)[order]
    r = np.asarray(retention_pct, float)[order]
    if r[0] <= 80:
        return float(t[0])
    for i in range(1, len(t)):
        if r[i] <= 80 and r[i - 1] > 80:
            return float(t[i-1] + (80-r[i-1]) * (t[i]-t[i-1]) / (r[i]-r[i-1]))
    return None


def _retention(dt, amplitudes, rates):
    return np.exp(-np.outer(np.asarray(dt), rates)) @ amplitudes


@dataclass
class FitResult:
    n_components: int
    peak_time_h: float
    peak_retention_pct: float
    observed_t80_h: float | None
    fitted_t80_h: float | None
    fitted_t80_extrapolated: bool | None
    rmse_post_peak_pct: float
    amplitudes: list
    rates_per_h: list


def fit_post_peak(time_h, retention_pct, n_components=2, max_extrapolation_factor=10):
    """Double/triple-exponential post-peak fit; solve T80 without a time cap.

    T80 is against 80% of the ORIGINAL initial PCE/MPP reference.
    Model is anchored at the device-specific post-burn-in peak.
    """
    if n_components not in (2, 3):
        raise ValueError('Choose 2 or 3 exponential components')
    order = np.argsort(time_h)
    t = np.asarray(time_h, float)[order]
    r = np.asarray(retention_pct, float)[order]
    if len(t) < max(6, n_components * 3) or np.any(np.diff(t) <= 0):
        raise ValueError('Too few points or duplicate times')
    p = int(np.argmax(r))
    tp, rp = float(t[p]), float(r[p])
    dt = t[p:] - tp
    y = r[p:] / rp
    if len(dt) < max(5, n_components * 2) or dt[-1] <= 0:
        raise ValueError('Insufficient post-peak data')
    # Positive amplitudes summing to one, nonnegative rates.
    horizon = float(dt[-1])
    initial_rates = np.geomspace(0.001 / horizon, 0.2 / horizon, n_components)
    x0 = np.r_[np.zeros(n_components), np.log(initial_rates)]
    def unpack(x):
        z = x[:n_components]
        weights = np.exp(z - np.max(z))
        weights /= weights.sum()
        rates = np.exp(np.clip(x[n_components:], -35, 5))
        return weights, rates
    def residual(x):
        a, k = unpack(x)
        return _retention(dt, a, k) - y
    result = least_squares(residual, x0, max_nfev=10000)
    a, k = unpack(result.x)
    idx = np.argsort(k)
    a, k = a[idx], k[idx]
    rmse = float(np.sqrt(np.mean(residual(result.x)**2)) * rp)
    obs = observed_t80(t, r)
    threshold = 80 / rp
    predicted = None
    extrapolated = None
    # Reuse exactly the fitted post-peak 2- or 3-exponential function.
    # The old max_extrapolation_factor is kept in the signature for backward
    # compatibility, but is NOT used to suppress mathematically finite T80.
    from fit_lifetime_solver import solve_fitted_t80
    predicted, status = solve_fitted_t80(tp, rp, a, k,
                                         threshold=80., measured_last_h=float(t[-1]))
    extrapolated = None if predicted is None else status=='FITTED_EXTRAPOLATED'
    return FitResult(n_components, tp, rp, obs, predicted, extrapolated, rmse, a.tolist(), k.tolist())
