"""v18.6.7 experimental panel-aware figure calibration (not DOI-specific).

Handles multi-panel raster graphs with horizontally repeating tick sequences.
Requires a fitted hour time axis, fitted monotonic Y ticks, and figure-level
stability/MPP evidence. Outputs only calibrated plots and leaves failures for QC.
"""
import re
from itertools import combinations
import numpy as np


def _ticks(words):
    from general_curve_reader import numeric_words
    return [p for p in numeric_words(words) if 0 <= p[2] <= 20000]


def _robust_fit(ticks, axis):
    if len(ticks)<3:return None
    idx=0 if axis=='x' else 1
    values=[]
    for tri in combinations(ticks,3):
        pos=np.array([t[idx] for t in tri],float)
        val=np.array([t[2] for t in tri],float)
        if np.ptp(pos)<22 or len(set(val))<3:continue
        m,b=np.polyfit(pos,val,1)
        if axis=='x' and m<=0 or axis=='y' and m>=0:continue
        eps=max(.22,np.ptp(val)*.025)
        matches=[t for t in ticks if abs(m*t[idx]+b-t[2])<eps]
        if len(matches)<3:continue
        xmin=min(t[idx] for t in matches);xmax=max(t[idx] for t in matches)
        if xmax-xmin<25:continue
        err=np.median([abs(m*t[idx]+b-t[2]) for t in matches])
        values.append(((len(matches),xmax-xmin,-err),m,b,matches))
    return max(values,key=lambda z:z[0])[1:] if values else None


def detect(page, words):
    from general_curve_reader import words_in_bbox
    import fitz
    nums=_ticks(words)
    # Cluster each physical horizontal tick row separately; unlike the original
    # detector, splitting on restarts prevents fitting across adjacent panels.
    bins=[]
    for t in sorted(nums,key=lambda z:(z[1],z[0])):
        if t[2]>5000:continue
        g=next((v for v in bins if abs(v[0]-t[1])<3.5),None)
        if g:g[1].append(t)
        else:bins.append([t[1],[t]])
    output=[]
    for yy,collection in bins:
        if len(collection)<3:continue
        xs=sorted(collection,key=lambda p:p[0]);runs=[];run=[]
        for tick in xs:
            if run and tick[2]<run[-1][2]-.20*max(1,run[-1][2]):
                if len(run)>=3:runs.append(run)
                run=[]
            run.append(tick)
        if len(run)>=3:runs.append(run)
        for seg in runs:
            # Use near-collinear subsets to resist extraneous annotation numbers.
            xf=_robust_fit(seg,'x')
            if not xf:continue
            slope,intercept,inliers=xf
            xmin=min(t[0] for t in inliers);xmax=max(t[0] for t in inliers)
            if xmax-xmin<45:continue
            xzero=-intercept/slope
            if xzero<xmin-55 or xzero>xmin+20:continue
            maxval=max(t[2] for t in inliers)
            if maxval<50:continue
            ybottom=yy-6
            # Prevent a Y calibration from adjacent charts by finding the closest
            # physical left-axis label column and three descending ticks.
            prior=[]
            # Identify a preceding x-tick baseline from the same column,
            # preventing the next stacked panel from reusing its Y ticks.
            for baseline,other in bins:
                if baseline >= ybottom-34 or baseline <= ybottom-158:continue
                aligned=[t for t in other if xzero+15<t[0]<(maxval-intercept)/slope+5 and t[2]>=50]
                if len(aligned)>=3:prior.append(baseline)
            yfloor=max(ybottom-145,max(prior)+7) if prior else ybottom-145
            local=[t for t in nums if xzero-30<t[0]<xzero-3 and yfloor<t[1]<ybottom-9 and t[2]<=110]
            groups=[]
            for t in sorted(local,key=lambda v:v[0]):
                g=next((a for a in groups if abs(np.median([x[0] for x in a])-t[0])<6),None)
                if g is None:groups.append([t])
                else:g.append(t)
            yfits=[]
            for grp in groups:
                z=_robust_fit(grp,'y')
                if z:
                    m,b,points=z
                    if len(points)>=3:
                        at_bottom=m*ybottom+b
                        if abs(at_bottom)<max(4,.35*max(t[2] for t in points)):
                            yfits.append((len(points),-abs(at_bottom),m,b,points))
            if not yfits:continue
            _,_,ym,yb,ypoints=max(yfits,key=lambda q:q[:2])
            ymax=max(t[2] for t in ypoints)
            # Common normalized PCE panels label 0.8, 0.6, 0.4 and 0.2;
            # extrapolate the 1.0 upper bound only from those calibrated ticks.
            if .7 <= ymax <= 1.01 and .08<abs(ym)*10<.45:
                ymax=1.0
            ytop=(ymax-yb)/ym
            # Make plot larger than its outermost ticks by a small margin only.
            ytop=max(0,ytop-3)
            xright=(maxval-intercept)/slope
            if ybottom-ytop<30 or ybottom-ytop>175 or xright-xzero<50:continue
            # Explicit hour units close to x-axis. Seconds/minutes/nanoseconds
            # are rejected even if the entire publication discusses stability.
            nearby=' '.join(w[4] for w in words if
                ybottom-2 <= (w[1]+w[3])/2 <= ybottom+30 and
                xzero-8 <= (w[0]+w[2])/2 <= xright+8)
            if re.search(r'\b(?:ns|µs|us|ms|minutes|seconds)\b',nearby,re.I):continue
            # At least x-axis 'Time' and an hour unit in the panel.
            if not (re.search(r'\btime\b',nearby,re.I) and
                    re.search(r'\bh\b|\bhours?\b',nearby,re.I)):
                # Some figures print the hour unit only below the lowest of
                # two stacked MPP subpanels. Retain axes as REVIEW CANDIDATE,
                # never silently treat missing units as verified hours.
                shared_unit=any('time' in w[4].lower() for w in words if
                    (w[0]+w[2])/2>xzero and (w[0]+w[2])/2<xright and
                    (w[1]+w[3])/2>ybottom and (w[1]+w[3])/2<ybottom+155)
                if not shared_unit or not re.search(r'(?i)MPP tracking|ISOS-L',page.get_text('text')):
                    continue
            label=words_in_bbox(words,fitz.Rect(max(0,xzero-55),max(0,ytop),xzero-3,ybottom+2))
            evidence=page.get_text('text')
            if not re.search(r'\b(?:PCE|SPO|power conversion|MPPT|MPP|retention|stability)\b',label+' '+evidence,re.I):continue
            if ymax>35 and ymax>1.6:continue
            # No global y-coordinate extrapolation from nonmatching source axes.
            candidate={'page':page.number+1,'x0':float(xzero),'x1':float(xright),
                'y0':float(ytop),'y1':float(ybottom),
                'x_slope':float(slope),'x_intercept':float(intercept),
                'y_slope':float(ym),'y_intercept':float(yb),
                'x_max':float(maxval),'y_max':float(ymax),
                'xlabel':'Time (h)','source':'OCR_PANEL_HOUR_AXIS',
                'y_axis_label':label or 'PCE_OR_SPO_REQUIRES_AXIS_QC',
                'panel_calibration_status':'AUTOMATIC_OCR_FITTED_REVIEW_REQUIRED'}
            if not any(abs(c['x0']-candidate['x0'])<9 and abs(c['y1']-candidate['y1'])<9 for c in output):output.append(candidate)
    return output
