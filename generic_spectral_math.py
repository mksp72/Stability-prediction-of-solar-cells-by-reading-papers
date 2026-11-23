"""Generic red-edge interpolation with no article-specific data."""
import numpy as np

def spectral_50_bandgap(wavelength_nm, eqe_pct):
    x=np.asarray(wavelength_nm,dtype=float);y=np.asarray(eqe_pct,dtype=float)
    keep=np.isfinite(x)&np.isfinite(y)
    x,y=x[keep],y[keep]
    if len(x)<6:raise ValueError('Insufficient EQE points')
    order=np.argsort(x);x,y=x[order],y[order]
    high=np.nanpercentile(y,85)
    if high<=0:raise ValueError('No positive EQE plateau')
    target=0.5*high
    crossings=[]
    for i in range(1,len(x)):
        if x[i] > 500 and y[i-1]>=target>y[i] and y[i]!=y[i-1]:
            crossings.append(x[i-1]+(target-y[i-1])*(x[i]-x[i-1])/(y[i]-y[i-1]))
    if not crossings:raise ValueError('No red-edge half-height crossing')
    lam=float(max(crossings))
    return {'lambda50_nm':lam,'estimated_bandgap_ev':1239.841984/lam}
