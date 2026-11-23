from pathlib import Path
p=Path(__file__).with_name('general_curve_reader.py');s=p.read_text()
place='def locate_text_axes(page,words_override=None):'
func='''def _robust_ocr_y_ticks(ticks, bottom):
    """Fit raster y ticks despite occasional OCR junk (10 for 1.0, etc).

    Require at least 3 mutually consistent labels in one printed tick column,
    a negative slope, and a geometrically plausible intercept at the baseline.
    Never infer the y scale from unrelated chart text alone.
    """
    from itertools import combinations
    if len(ticks)<3:return None
    candidates=[]
    for tri in combinations(ticks,3):
        if max(q[0] for q in tri)-min(q[0] for q in tri)>10:continue
        values=np.array([v[2] for v in tri],float)
        yy=np.array([v[1] for v in tri],float)
        if np.ptp(yy)<32 or len(set(values))<3:continue
        slope,intercept=np.polyfit(yy,values,1)
        if slope>=0:continue
        residual=np.abs(slope*yy+intercept-values)
        if max(residual)>max(.25,.035*np.ptp(values)):continue
        at_bottom=slope*bottom+intercept
        max_v=max(values)
        if max_v<=0 or not(-.15*max_v<=at_bottom<=.25*max_v):continue
        aligned=[a for a in ticks if abs(a[0]-np.median([q[0] for q in tri]))<7]
        inliers=[a for a in aligned if abs(slope*a[1]+intercept-a[2])<=max(.13,.035*max_v)]
        if len(inliers)<3:continue
        score=(len(inliers), -abs(at_bottom), np.ptp([a[1] for a in inliers]),-max_v)
        candidates.append((score,{'ticks':sorted(inliers,key=lambda z:z[1]),'slope':float(slope),
                                  'intercept':float(intercept),'error':float(np.median(residual))}))
    return max(candidates,key=lambda z:z[0])[1] if candidates else None

'''
s=s.replace(place,func+place)
needle="yf=min(candidates_y,key=lambda x:x[0])[1] if candidates_y else _fit_axis(yt,'y',True)"
replacement="""yf=min(candidates_y,key=lambda x:x[0])[1] if candidates_y else _fit_axis(yt,'y',True)
        if words_override is not None and (not yf or (yf and max(a[2] for a in yf['ticks'])>120)):
            # Find correct tick values even when intermediate OCR glyphs are
            # malformed, but only when 3+ physical labels remain readable.
            robust=_robust_ocr_y_ticks(yt,bottom)
            if robust and (not yf or max(a[2] for a in robust['ticks'])<max(a[2] for a in yf['ticks'])):
                yf=robust"""
assert needle in s;s=s.replace(needle,replacement)
p.write_text(s)
