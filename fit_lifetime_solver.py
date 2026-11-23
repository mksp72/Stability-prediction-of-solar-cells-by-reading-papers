"""Solve fitted multi-exponential stability lifetimes with adaptive root finding.

No arbitrary multiple of measured duration is used as an upper bound.
T80 refers to 80% of the ORIGINAL PCE/MPP reference, not 80% of
post-burn-in peak. All time coordinates returned are absolute aging hours.
"""
from __future__ import annotations
import math
import numpy as np
from scipy.optimize import brentq


def solve_fitted_t80(peak_time_h, peak_retention_pct, amplitudes, rates_hinv,
                     threshold=80.0, measured_last_h=None):
    amps=np.asarray(amplitudes,float)
    rates=np.asarray(rates_hinv,float)
    peak=float(peak_retention_pct); start=float(peak_time_h)
    if not (len(amps)==len(rates) and len(amps)>=2 and np.all(np.isfinite(amps))
            and np.all(np.isfinite(rates)) and np.all(amps>=0) and np.all(rates>=0)
            and math.isfinite(peak) and math.isfinite(start) and peak>0
            and np.isclose(amps.sum(),1,atol=0.002)):
        return None,'INVALID_FIT_PARAMETERS'
    if peak <= threshold:
        # The crossing, if any, occurred during burn-in or before the peak.
        return None,'ALREADY_BELOW_THRESHOLD_AT_PEAK'
    floor=float(peak*amps[rates==0].sum())
    if floor >= threshold - 1e-10:
        return None,'NO_FINITE_T80_MODEL_ASYMPTOTE_AT_OR_ABOVE_80'

    def f(dt):
        return peak*float(np.dot(amps,np.exp(-rates*dt)))-threshold

    positives=rates[rates>0]
    bound=max(1.,float(measured_last_h or start)-start)
    if len(positives):
        bound=max(bound, math.log(max(peak,threshold)/threshold)/float(np.max(positives)))
    # Exponentially grow until crossing: no arbitrary 10x-measurement cutoff.
    # Max finite float prevents overflow; physical interpretability is recorded separately.
    for _ in range(1024):
        if f(bound)<=0:
            crossing=float(start+brentq(f,0.,bound,xtol=1e-10,rtol=1e-13))
            if not math.isfinite(crossing):return None,'NUMERIC_OVERFLOW'
            measured=(measured_last_h is not None and crossing<=float(measured_last_h)+1e-7)
            return crossing,'FITTED_WITHIN_MEASURED_INTERVAL' if measured else 'FITTED_EXTRAPOLATED'
        if bound>1e300:
            return None,'NUMERIC_TIME_OVERFLOW'
        bound*=2.
    return None,'NUMERIC_BRACKET_FAILURE'
