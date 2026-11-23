"""Layout-independent time-series plot reader for PDF text/vector/raster figures.

No PDF names, DOIs, or page coordinates are embedded. Extracts candidate plots only
where calibrated numeric axes and connected colored traces can be established.
Each curve retains its inferred color and evidence; unknown legend labels are not invented.
This is an experimental automated digitizer, not a guarantee for all figures.
"""
from __future__ import annotations
import re, json, math
from pathlib import Path
import cv2, fitz, numpy as np
from scipy.ndimage import binary_opening

TIME_LABEL=re.compile(r'^(?:time(?:\s*\(\s*h(?:ours?)?\s*\))?|t|duration|aging|ageing|hours?|hrs?)$',re.I)
NUM=re.compile(r'^[-+]?(?:\d{1,3}(?:,\d{3})+|\d+)(?:\.\d+)?$')
STABILITY=re.compile(r'stabilit|retention|aging|ageing|degrad|storage|ISOS|MPP|maximum power|PCE',re.I)
COLORS={
 'red':((0,13),(170,179)), 'blue':((96,135),), 'green':((33,88),),
 'cyan':((84,102),), 'magenta':((137,169),), 'violet':((117,155),), 'orange':((14,25),), 'gray':(), 'black':()
}

def numeric_words(words):
    out=[]
    for w in words:
        z=w[4].strip('()%,')
        if NUM.fullmatch(z):
            v=float(z.replace(',',''))
            if np.isfinite(v):out.append((float((w[0]+w[2])/2),float((w[1]+w[3])/2),v,w[4]))
    return out

def _fit_axis(ticks,axis,positive=False):
    # stable pair-based linear calibration, excluding duplicate/offset text.
    if len(ticks)<3:return None
    for ytol in ((9.,14.,20.) if axis=='y' else (2.,3.5,5.5)):
        bins=[]
        for t in ticks:
            near=next((b for b in bins if abs(b[0][1 if axis=='x' else 0]-t[1 if axis=='x' else 0])<ytol),None)
            if near is None:bins.append([t])
            else:near.append(t)
        for b in sorted(bins,key=len,reverse=True):
            if len(b)<3:continue
            ts=sorted(b,key=lambda x:x[0 if axis=='x' else 1]);coords=np.array([a[0 if axis=='x' else 1] for a in ts]); vals=np.array([a[2] for a in ts]);
            if len(np.unique(coords.round(1)))<3 or len(np.unique(vals))<3:continue
            slope,intercept=np.polyfit(coords,vals,1)
            if (axis=='x' and slope<=0) or (axis=='y' and slope>=0):continue
            predicted=slope*coords+intercept
            err=np.median(np.abs(predicted-vals)); span=max(vals)-min(vals)
            if span<=0 or err>span*.035:continue
            if positive and vals.min()<0:continue
            return {'ticks':ts,'slope':float(slope),'intercept':float(intercept),'error':float(err)}
    return None

def _robust_ocr_y_ticks(ticks, bottom):
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

def locate_text_axes(page,words_override=None):
    words=words_override if words_override is not None else page.get_text('words',sort=True)
    nums=numeric_words(words)
    candidates=[]
    for w in words:
        if not TIME_LABEL.fullmatch(w[4].strip(' :')):continue
        # Short JV/Mpp stabilization plots in seconds/minutes are not ISOS
        # aging curves measured in hours. Don't silently assign hours.
        if any(t[4].strip().lower() in ('(s)','(sec)','(seconds)','(min)','(mins)')
               and abs(t[1]-w[1])<8 and 0<=t[0]-w[2]<40 for t in words):continue
        cx=(w[0]+w[2])/2; ly=w[1]
        # x-axis numeric tick rows occur immediately above the 'Time' label
        xt=[(x,y,v,txt) for x,y,v,txt in nums if abs(x-cx)<340 and ly-38<y<ly-2 and v>=0]
        # Split adjacent plots sharing the same tick baseline when values restart
        # at zero (instead of treating the entire row as one malformed axis).
        groups=[]
        # OCR can shift labels in the *same* tick row by several pixels.
        for a in sorted(xt,key=lambda q:q[1]):
            near=next((g for g in groups if abs(np.median([t[1] for t in g])-a[1])<6.2),None)
            if near is None:groups.append([a])
            else:near.append(a)
        proposals=[]
        for group in groups:
            group=sorted(group,key=lambda a:a[0]); segments=[]; segment=[]
            for a in group:
                if segment and a[2]<segment[-1][2]-.1*max(1,segment[-1][2]):
                    segments.append(segment);segment=[]
                segment.append(a)
            if segment:segments.append(segment)
            for segment in segments:
                # OCR frequently mixes a neighboring JV-axis tick into the
                # time-tick row. Also test contiguous sub-runs of >=3 ticks.
                variants=[segment]
                if len(segment)>3:
                    for nsub in range(3,len(segment)):
                        for j in range(len(segment)-nsub+1):
                            variants.append(segment[j:j+nsub])
                for subset in variants:
                    xf=_fit_axis(subset,'x',True)
                    if xf and subset[0][0]-15<=cx<=subset[-1][0]+15:
                        proposals.append((abs(np.mean([a[0] for a in xf['ticks']])-cx)-3*len(xf['ticks']),-len(xf['ticks']),xf))
        xfit=sorted(proposals,key=lambda z:(z[0],z[1]))[0][2] if proposals else None
        if not xfit:continue
        ts=xfit['ticks']; x0=min(a[0] for a in ts);x1=max(a[0] for a in ts)
        # The x=0 tick center is also the y-axis x coordinate.
        sx,ix=xfit['slope'],xfit['intercept']
        xzero=-ix/sx
        if not x0-max(12,(1.4 if words_override is not None else 0.35)*(x1-x0))<=xzero<=x1+12:continue
        bottom=float(np.median([a[1] for a in ts]))-5.5
        # Numeric y-ticks are grouped by x at the left of the plot area.
        # Exclude other vertically stacked panels using previous axis labels.
        prev=[q[1] for q in words if TIME_LABEL.fullmatch(q[4].strip(' :')) and q[1]<ly-25 and abs((q[0]+q[2])/2-cx)<75]
        yfloor=max(0,max(prev)+12) if prev else max(0,bottom-(175 if words_override is not None else 500))
        yt=[a for a in nums if xzero-58<a[0]<xzero-3 and yfloor<a[1]<bottom+2]
        if words_override is not None:
            # If the axis explicitly contains fractional labels (e.g. 1.2,
            # 0.6, 0.0), OCR sometimes drops the decimal point from 1.0,
            # 0.8, 0.2. Correct only a single matching tick column.
            decimal_ticks=[a for a in yt if '.' in a[3] and 0<a[2]<=2]
            if len(decimal_ticks)>=2:
                major_x=float(np.median([a[0] for a in decimal_ticks]))
                revised=[]
                for a in yt:
                    if (abs(a[0]-major_x)<9 and '.' not in a[3]
                        and re.fullmatch(r'(?:0?[1-9]|1[0-2])',a[3]) and a[2]>2):
                        revised.append((a[0],a[1],a[2]/10.,a[3]+'[OCR_DECIMAL_REPAIRED]'))
                    else:revised.append(a)
                yt=revised
        if words_override is not None and len(yt)>=4:
            # Raster export can append a grid-line fragment to each y tick.
            # Example: 20, 154, 107, 54 printed as 20, 15, 10, 5.
            # Only repair when a low-valued upper tick and >=2 intermediate
            # oversized ticks yield a descending scale after truncation.
            descending=sorted(yt,key=lambda a:a[1]); top_tick=descending[0]
            if 0<top_tick[2]<=35:
                fixed=[];repairs=0
                for a in yt:
                    raw=a[3]
                    if (a[1]>top_tick[1]+8 and a[2]>top_tick[2]
                        and re.fullmatch(r'\d{2,3}',raw)
                        and 0<float(raw[:-1])<top_tick[2]):
                        fixed.append((a[0],a[1],float(raw[:-1]),raw+'[OCR_GRIDLINE_REPAIRED]'))
                        repairs+=1
                    else:fixed.append(a)
                if repairs>=2:yt=fixed
        if words_override is not None and len(yt)>=5:
            # Reject stray numeric labels from nearby plots. Proper y ticks
            # share approximately the same right-aligned text column.
            hist=sorted([a[0] for a in yt])
            groups_x=[]
            for xx in hist:
                grp=next((g for g in groups_x if abs(np.median(g)-xx)<6),None)
                if grp is None:groups_x.append([xx])
                else:grp.append(xx)
            major=max(groups_x,key=len)
            if len(major)>=3:
                cx_tick=float(np.median(major))
                yt=[a for a in yt if abs(a[0]-cx_tick)<7.5]
        yf=None
        # Split multiple vertical plots: y-tick values restart at the next chart.
        ordered=sorted(yt,key=lambda a:a[1]);segments=[];seg=[]
        for tick in ordered:
            if seg and (tick[1]-seg[-1][1]>46 or tick[2]>seg[-1][2]+1e-5):
                if len(seg)>=3:segments.append(seg)
                seg=[]
            seg.append(tick)
        if len(seg)>=3:segments.append(seg)
        candidates_y=[]
        for seg in segments:
            attempt=_fit_axis(seg,'y',True)
            if attempt:candidates_y.append((abs(bottom-max(a[1] for a in attempt['ticks'])),attempt))
        yf=min(candidates_y,key=lambda x:x[0])[1] if candidates_y else _fit_axis(yt,'y',True)
        if words_override is not None and (not yf or (yf and max(a[2] for a in yf['ticks'])>120)):
            # Find correct tick values even when intermediate OCR glyphs are
            # malformed, but only when 3+ physical labels remain readable.
            robust=_robust_ocr_y_ticks(yt,bottom)
            if robust and (not yf or max(a[2] for a in robust['ticks'])<max(a[2] for a in yf['ticks'])):
                yf=robust
        if not yf and words_override is not None:
            # Raster plots sometimes yield just two readable y ticks. Fit a
            # two-tick scale only if their extrapolated zero is near the
            # observed x-axis baseline; this rejects unrelated plot numbers.
            viable=[]
            for i,aa in enumerate(yt):
                for bb in yt[i+1:]:
                    if abs(aa[0]-bb[0])>6 or abs(aa[1]-bb[1])<17:continue
                    slope=(bb[2]-aa[2])/(bb[1]-aa[1])
                    if slope>=0 or abs(bb[2]-aa[2])<5:continue
                    off=aa[2]-slope*aa[1]
                    atbottom=slope*bottom+off
                    if not (-0.16*max(aa[2],bb[2],1)<atbottom<0.20*max(aa[2],bb[2],1)):continue
                    viable.append((abs(atbottom),{'ticks':[aa,bb],'slope':slope,'intercept':off,'error':0.0}))
            if viable:yf=min(viable,key=lambda z:z[0])[1]
        if not yf:continue
        yticks=yf['ticks'];ymin=min(a[1] for a in yticks); ytop=max(0,min(yticks,key=lambda a:a[1])[1]-1)
        # x max is from last tick; y upper based on first y tick.
        max_value=max(a[2] for a in ts)
        # A final OCR tick such as '1,000' or '500' can be shifted a few
        # pixels vertically and omitted from the best linear subset. Recover
        # only ticks that agree with its physical calibrated x position.
        for xx,yy,vv,_ in nums:
            if xx<=max(a[0] for a in ts)+4 or vv<=max_value:continue
            if abs(yy-np.median([a[1] for a in ts]))>8:continue
            if abs(sx*xx+ix-vv)<max(12,0.05*vv):max_value=float(vv)
        if max_value<10:continue
        xright=(max_value-ix)/sx
        ymaxv=max(a[2] for a in yticks)
        if ymaxv<=0:continue
        top=(ymaxv-yf['intercept'])/yf['slope']
        top=max(0,top-max(6,0.14*(bottom-top)))
        xright=min(xright,page.rect.width-5)
        if xright-xzero<50 or bottom-top<36 or top<0:continue
        candidates.append({'page':page.number+1,'x0':xzero,'x1':xright,'y0':top,'y1':bottom,
            'x_slope':sx,'x_intercept':ix,'y_slope':yf['slope'],'y_intercept':yf['intercept'],
            'x_max':max_value,'y_max':ymaxv,'xlabel':w[4], 'source':'OCR_IMAGE_AXIS' if words_override is not None else 'PDF_TEXT_AXIS'})
    # dedup labels centered on same chart
    unique=[]
    for c in sorted(candidates,key=lambda z:(z['page'],z['x0'],z['y0'])):
        if not any(abs(c['x0']-u['x0'])<10 and abs(c['y1']-u['y1'])<12 for u in unique):unique.append(c)
    return unique

def locate_numeric_time_axes(page, ocr_words):
    """Conservative raster-plot fallback when italic `t [h]` is unreadable.

    A candidate must have a monotonic approximately linear numeric x-tick
    row (at least four ticks; an origin near zero), contextual evidence of
    stability versus time and a valid PCE/retention-like y axis. No inferred
    units are exported without the nearby article mentioning time in hours.
    """
    caption=page.get_text('text')
    if not (re.search(r'stabil|degrad|tracking|retention',caption,re.I) and
            re.search(r'\b(?:h|hours?|MPP|maximum power|time)\b',caption,re.I)):
        return []
    nums=numeric_words(ocr_words)
    rows=[]
    for t in sorted(nums,key=lambda q:q[1]):
        group=next((g for g in rows if abs(np.median([a[1] for a in g])-t[1])<4),None)
        if group is None:rows.append([t])
        else:group.append(t)
    axes=[]
    for row in rows:
        if len(row)<4:continue
        ordered=sorted(row,key=lambda a:a[0]);segments=[];ss=[]
        for a in ordered:
            if ss and a[2]<ss[-1][2]:
                segments.append(ss);ss=[]
            ss.append(a)
        if ss:segments.append(ss)
        for seg in segments:
            if len(seg)<4:continue
            fit=_fit_axis(seg,'x',True)
            if not fit:continue
            ticks=fit['ticks']
            vmin=min(a[2] for a in seg);vmax=max(a[2] for a in ticks)
            if not (0<=vmin<=(.06*vmax) and 24<=vmax<=30000):continue
            if max(a[0] for a in ticks)-min(a[0] for a in ticks)<75:continue
            # Explicitly represent the inferred axis anchor for the existing
            # image/calibration function, while keeping evidence on source.
            cx=float(np.median([a[0] for a in ticks]));yt=float(np.median([a[1] for a in ticks]))+18
            virtual=(cx-9,yt-4,cx+9,yt+4,'Time',-1,-1,-1)
            for chart in locate_text_axes(page,ocr_words+[virtual]):
                if abs(chart['x0']-(-fit['intercept']/fit['slope']))>18:continue
                if abs(chart['y1']-(float(np.median([a[1] for a in ticks]))-5.5))>15:continue
                chart['source']='OCR_NUMERIC_TIME_AXIS'
                chart['time_axis_inferred_from_numeric_ticks']=True
                if not any(abs(chart['x0']-k['x0'])<12 and abs(chart['y1']-k['y1'])<12 for k in axes):
                    axes.append(chart)
    return axes

def _mask(hsv,color):
    h,s,v=cv2.split(hsv);mask=np.zeros(h.shape,dtype=np.uint8)
    if color=='gray':return ((s<46)&(v>=95)&(v<=229)).astype(np.uint8)
    if color=='black':return ((s<55)&(v>=25)&(v<95)).astype(np.uint8)
    for lo,hi in COLORS[color]:mask|=((h>=lo)&(h<=hi)&(s>(32 if color=='violet' else 72))&(v>57)).astype(np.uint8)
    return mask

def sample_colored_trace(image,plot,color,scale=3):
    hsv=cv2.cvtColor(image,cv2.COLOR_RGB2HSV)
    mask=_mask(hsv,color)
    x0,x1,y0,y1=[float(plot[k])*scale for k in ('x0','x1','y0','y1')]
    lo,hi=int(max(0,x0+4)),int(min(image.shape[1]-1,x1-4))
    ya,yb=int(max(0,y0+max(40 if color in ('black','gray') else 10,(y1-y0)*(.15 if color in ('black','gray') else .045)))),int(min(image.shape[0]-1,y1-max(4,(y1-y0)*(.085 if color in ('gray','black') else .035))))
    if hi<=lo or yb<=ya:return []
    # A graph line is a path through many x-columns, unlike legend symbols and text.
    sampled=[]
    for xi in range(lo,hi,4):
        ys=np.flatnonzero(np.max(mask[ya:yb,xi-1:xi+2],axis=1)>0)+ya
        if len(ys)<1:continue
        starts=np.r_[0,np.where(np.diff(ys)>3)[0]+1]
        ends=np.r_[starts[1:],len(ys)]
        candidates=[]
        for a,b in zip(starts,ends):
            if b-a<1:continue
            vv=ys[a:b];candidates.append(float(np.median(vv)))
        if candidates:sampled.append((float(xi),candidates))
    if len(sampled)<max(6,int((hi-lo)/40)):return []
    # Find smooth longest horizontal trajectory with dynamic programming. Max jump
    # across one step proportional to plot height and sample step.
    # Each segment may terminate at missing colors; short interruptions allowed.
    best=[]
    for start in range(min(14,len(sampled))):
        for ystart in sampled[start][1][:8]:
            path=[(sampled[start][0],ystart)];lastx,lasty=path[-1]
            for x,yy in sampled[start+1:]:
                dx=x-lastx
                if dx>max(32,(hi-lo)*.10):break
                nearest=min(yy,key=lambda y:abs(y-lasty))
                if abs(nearest-lasty)>max(7, min(22, dx*.45+3)):
                    continue
                path.append((x,nearest));lastx,lasty=x,nearest
            if len(path)>len(best):best=path
    if len(best)<max(6,int((hi-lo)/65)) or (best[-1][0]-best[0][0])<(hi-lo)*.16:return []
    # remove strong text contamination and compress into physical x/h and y units.
    rows=[]
    for xx,yy in best:
        x=plot['x_slope']*(xx/scale)+plot['x_intercept']
        y=plot['y_slope']*(yy/scale)+plot['y_intercept']
        if x<0 or x>plot['x_max']*1.04 or y < -.1*plot['y_max'] or y>plot['y_max']*1.15:continue
        rows.append((float(x),float(y)))
    if len(rows)<7:return []
    # 3-point rolling median rejects raster text artifacts
    out=[]
    for i,(x,y) in enumerate(rows):
        j=max(0,i-2);k=min(len(rows),i+3)
        z=float(np.median([q[1] for q in rows[j:k]]))
        if abs(y-z) <max(.9,.08*plot['y_max']):out.append((round(x,3),round(y,4)))
    return out

def ocr_page_words(page,zoom=3):
    """Last-resort OCR for a raster graph; optional Tesseract binary on PATH.

    Called once per possible figure page, never as part of per-pixel sampling.
    """
    try:
        import pytesseract
        from PIL import Image
        import os
        if os.name=='nt':
            exe=Path(r'C:\Program Files\Tesseract-OCR\tesseract.exe')
            if exe.is_file():pytesseract.pytesseract.tesseract_cmd=str(exe)
        pix=page.get_pixmap(matrix=fitz.Matrix(zoom,zoom),alpha=False)
        im=Image.frombytes('RGB',[pix.width,pix.height],pix.samples)
        data=pytesseract.image_to_data(im,config='--psm 11',output_type=pytesseract.Output.DICT,timeout=30)
    except (ImportError,RuntimeError,FileNotFoundError,OSError):return []
    words=[]
    for i,t in enumerate(data['text']):
        if not t.strip():continue
        try:confidence=float(data['conf'][i])
        except (TypeError,ValueError):continue
        if confidence<45:continue
        x,y=float(data['left'][i])/zoom,float(data['top'][i])/zoom
        w,h=float(data['width'][i])/zoom,float(data['height'][i])/zoom
        words.append((x,y,x+w,y+h,t, int(data['block_num'][i]),int(data['line_num'][i]),int(data['word_num'][i])))
    return words

def words_in_bbox(words,rect):
    a=[]
    for w in words:
        if rect.x0<=.5*(w[0]+w[2])<=rect.x1 and rect.y0<=.5*(w[1]+w[3])<=rect.y1:a.append(w[4])
    return ' '.join(a)

def is_stability_axis(page,plot,ocr_words=None):
    # Numeric time axes alone are insufficient: many figures show Voc, strain, etc.
    clip=fitz.Rect(max(0,plot['x0']-90),max(0,plot['y0']-8),plot['x0']-5,plot['y1']+8)
    label=page.get_text('text',clip=clip) if ocr_words is None else words_in_bbox(ocr_words,clip)
    keep=bool(re.search(r'\b(?:PCE|SPO|MPP|MPPT|retention|normalized|efficiency|power\s+output|power\s+conversion|P\s*norm)\b|η\s*MPPT|ηMPPT|Pnorm',label,re.I))
    reject=bool(re.search(r'FWHM|strain|V.?OC|J.?SC|current\s+density|fill\s+factor|\bFF\b|resistance|luminescence',label,re.I))
    if ocr_words is not None and not keep and not reject:
        # Scans frequently lose rotated PCE/SPO axis text. Accept only plots
        # with suitable numerical y ticks and an explicit stability/PCE caption.
        fulltext=page.get_text('text')
        if plot['y_max']<=110 and re.search(r'\b(?:stability|degradation|retention)\b',fulltext,re.I) and re.search(r'\b(?:PCEs?|power conversion efficiency|SPO|MPP)\b',fulltext,re.I):
            keep=True;label='PCE_OR_SPO_%_INFERRED_FROM_CAPTION'
    if ocr_words is not None and not keep and not reject and plot['y_max']<=1.6:
        # A *local* raster PCE/MPP criterion for figure-only manuscripts:
        # recognizable 0..1 normalized retention scale, hour-calibrated x
        # axis, and explicit MPP/retention evidence INSIDE THIS PANEL.
        # This does not infer stability from unrelated body text or other plots.
        region=fitz.Rect(max(0,plot['x0']),max(0,plot['y0']-8),
                         min(page.rect.width,plot['x1']+8),min(page.rect.height,plot['y1']+8))
        local=words_in_bbox(ocr_words,region)
        if (0.85<=plot['y_max']<=1.6 and plot['x_max']>=100
                and re.search(r'\bMPP\b|tracking|T\s*80|T\s*90|normalized',local,re.I)):
            keep=True;label='NORMALIZED_PCE_OCR_FIGURE_ONLY_PANEL'
    return keep and not reject,label.strip().replace('\n',' ')[:110]

def infer_legend(page,plot,image,scale=3,ocr_words=None):
    """Link a plotted color with the nearest same-color legend swatch + its text."""
    h,w=image.shape[:2];hsv=cv2.cvtColor(image,cv2.COLOR_RGB2HSV)
    labels={}
    d=page.get_text('dict')
    line_boxes=[]
    if ocr_words is None:
        for b in d['blocks']:
            for line in b.get('lines',[]):
                line_boxes.append((fitz.Rect(line['bbox']),' '.join(s['text'].strip() for s in line['spans']).strip()))
    else:
        groups={}
        for w0 in ocr_words:
            groups.setdefault((w0[5],w0[6]),[]).append(w0)
        for group in groups.values():
            group.sort(key=lambda v:v[0])
            text=' '.join(v[4] for v in group)
            rect=fitz.Rect(min(v[0] for v in group),min(v[1] for v in group),
                           max(v[2] for v in group),max(v[3] for v in group))
            line_boxes.append((rect,text))
    for rect,text in line_boxes:
            if not text or len(text)>75 or NUM.fullmatch(text) or re.search(r'time|hours|PCE|SPO|MPP|^\s*\(|^\s*>\s*\d|T\s*\d{2}\s*[=>]',text,re.I):continue
            cy=(rect.y0+rect.y1)/2
            if not (plot['x0']+10<rect.x0<plot['x1']+35 and plot['y0']-45<cy<plot['y1']+10):continue
            ax=max(0,int((rect.x0-35)*scale));bx=min(w,int((rect.x0-2)*scale))
            ay=max(0,int((cy-3)*scale));by=min(h,int((cy+3)*scale))
            if bx<=ax or by<=ay:continue
            for color in COLORS:
                n=int(_mask(hsv[ay:by,ax:bx],color).sum())
                if n>10 and (color not in labels or n>labels[color][0]):labels[color]=(n,text)
    if ocr_words is not None:
        # Raster legends often split each label into separate OCR words.
        # Match the pixel hue of the actual swatch *immediately preceding*
        # its category name, instead of using nearby arbitrary figure text.
        for word in ocr_words:
            name=word[4].strip('.,:')
            if not re.fullmatch(r'Control|Reference|Target(?:-\d+)?|Pristine|Theophylline|Passivated',name,re.I):continue
            if not (plot['x0']<word[0]<plot['x1']+10 and plot['y0']<word[1]<plot['y1']):continue
            xx=int(word[0]*scale);cy=int(.5*(word[1]+word[3])*scale)
            roi=hsv[max(0,cy-12):cy+12,max(0,xx-70):max(0,xx-3)]
            if roi.size==0:continue
            scores={color:int(_mask(roi,color).sum()) for color in ('red','orange','blue','cyan','green','magenta','black')}
            best=max(scores,key=scores.get)
            if scores[best]>16 and (best not in labels or scores[best]>labels[best][0]):
                labels[best]=(scores[best],name)
    out={color:v[1] for color,v in labels.items()}
    # Better read formula/additive labels from the legend itself if normal
    # line OCR sees only the swatch or an annotation above a trace.
    if ocr_words is not None and len([a for c,a in out.items() if c not in ('black','gray')])<2:
        for color,label in ocr_corner_legend(page,plot,image,scale).items():
            out[color]=label
    return out

def ocr_corner_legend(page,plot,image,scale=3):
    """Optional raster-legend OCR, with swatch color measured beside each line.

    Limited to the bottom-right plot region; original labels remain the
    preferred source. No publication-specific labels or colors are assumed.
    """
    try:
        import pytesseract
        from PIL import Image
        roi=fitz.Rect(plot['x0']+.45*(plot['x1']-plot['x0']),
                      plot['y0']+.40*(plot['y1']-plot['y0']),
                      plot['x1']-2,plot['y1']-3)
        if roi.width<75 or roi.height<25:return {}
        pix=page.get_pixmap(matrix=fitz.Matrix(4,4),clip=roi,alpha=False)
        im=Image.frombytes('RGB',(pix.width,pix.height),pix.samples)
        d=pytesseract.image_to_data(im,config='--psm 6',output_type=pytesseract.Output.DICT)
        groups={}
        for i,word in enumerate(d['text']):
            token=word.strip()
            if not token or float(d['conf'][i])<28:continue
            groups.setdefault((d['block_num'][i],d['par_num'][i],d['line_num'][i]),[]).append(i)
        hsv=cv2.cvtColor(image,cv2.COLOR_RGB2HSV)
        found={}
        for inds in groups.values():
            words=sorted(inds,key=lambda j:d['left'][j])
            meaningful=[j for j in words if re.search(r'[A-Za-z]{2,}|[A-Z][a-z]?[A-Z]?',d['text'][j])
                        and not re.fullmatch(r'[\W_]+',d['text'][j])]
            if not meaningful:continue
            line=' '.join(d['text'][j].strip() for j in words
                          if d['left'][j]>=d['left'][meaningful[0]] and
                          (j in meaningful or d['text'][j].strip() in ('+','-'))).strip()
            # Optical OCR often confuses a capital iodine with a lower-case l
            # in a chemical formula (PbI3 -> Pbl3); correct only the element
            # symbol context, never arbitrary scientific words.
            line=re.sub(r'Pb[l1](?=\d)', 'PbI',line)
            if not (3<len(line)<65) or re.search(r'\b(?:Time|hours|normalized|PCE|JSC|VOC)\b',line,re.I):continue
            x=roi.x0+(d['left'][meaningful[0]]/4)
            yc=roi.y0+((d['top'][meaningful[0]]+d['height'][meaningful[0]]*.5)/4)
            xa=max(0,int((x-31)*scale));xb=max(xa+1,int((x-3)*scale))
            ya=max(0,int((yc-3)*scale));yb=min(hsv.shape[0],int((yc+3)*scale))
            region=hsv[ya:yb,xa:xb]
            if region.size==0:continue
            scores={c:int(_mask(region,c).sum()) for c in COLORS if c not in ('black','gray')}
            color=max(scores,key=scores.get)
            if scores[color]>15:
                if color not in found or scores[color]>found[color][0]:found[color]=(scores[color],line)
        return {c:text for c,(score,text) in found.items()}
    except (ImportError,FileNotFoundError,OSError,RuntimeError,ValueError):return {}

def find_plots(pdf,max_pages=150):
    results=[]
    with fitz.open(str(pdf)) as doc:
        ocr_used=0
        for page in list(doc)[:max_pages]:
            # General detection works on figure axis labels, no DOI checking.
            plots=locate_text_axes(page)
            ocr_words=None
            # v18.6.4: Publisher preprints often append the figures as
            # image-only pages with only an ARTICLE IN PRESS text watermark.
            # These pages have no searchable caption; OCR the figures on those
            # pages instead of requiring a text caption on the SAME page.
            image_only_figure=(len(page.get_text('words'))<=12 and bool(page.get_images())
                and any((img['bbox'][2]-img['bbox'][0]>=.48*page.rect.width and img['bbox'][3]-img['bbox'][1]>=.30*page.rect.height)
                        for img in page.get_image_info()))
            if not plots and page.get_images() and (ocr_used<5 or image_only_figure):
                caption=page.get_text('text')
                # OCR only on pages with stability figure captions; avoid
                # multiple minutes of OCR on unrelated characterization pages.
                figcaption=re.search(r'(?:fig\.?|figure)\s*(?:S|\d)[^\n]{0,140}(?:stabil|retention|degrad|aging|operation)|(?:stabil|retention|degrad|aging|operation)[^\n]{0,140}(?:fig\.?|figure)\s*(?:S|\d)',caption,re.I)
                if (image_only_figure or figcaption or (re.search(r'\b(?:Fig\.?|Figure)\s*\d+',caption,re.I) and
                                  re.search(r'normalized power|stability|degradation|retention|MPP tracking|light stability|outdoor stability|damp heat',caption,re.I))):
                    if not image_only_figure:ocr_used+=1
                    ocr_words=ocr_page_words(page)
                    if ocr_words:
                        plots=locate_text_axes(page,ocr_words)
                        if not plots:plots=locate_numeric_time_axes(page,ocr_words)
            if not plots:continue
            zoom=3
            pix=page.get_pixmap(matrix=fitz.Matrix(zoom,zoom),alpha=False)
            im=np.frombuffer(pix.samples,np.uint8).reshape(pix.height,pix.width,pix.n).copy()[:,:,:3]
            for j,p in enumerate(plots):
                from figure_panel_context import panel_context
                p.update(panel_context(page,p,ocr_words if p['source'] in ('OCR_IMAGE_AXIS','OCR_NUMERIC_TIME_AXIS') else None))
                if p.get('panel_non_stability'):continue
                # Keep local OCR evidence from image-only publisher figure
                # pages, so protocol assignment does not depend on missing
                # searchable PDF captions.
                if ocr_words and p['source'] in ('OCR_IMAGE_AXIS','OCR_NUMERIC_TIME_AXIS'):
                    clip=fitz.Rect(max(0,p['x0']-8),max(0,p['y0']-15),
                        min(page.rect.width,p['x1']+18),min(page.rect.height,p['y1']+20))
                    p['ocr_panel_evidence']=words_in_bbox(ocr_words,clip)[:450]
                good,axis_label=is_stability_axis(page,p,ocr_words if p['source'] in ('OCR_IMAGE_AXIS','OCR_NUMERIC_TIME_AXIS') else None)
                if not good:continue
                p['y_axis_label']=axis_label
                legend=infer_legend(page,p,im,zoom,ocr_words if p['source'] in ('OCR_IMAGE_AXIS','OCR_NUMERIC_TIME_AXIS') else None)
                series=[]
                for color in COLORS:
                    points=sample_colored_trace(im,p,color,zoom)
                    if len(points)>8:
                        hsv=cv2.cvtColor(im,cv2.COLOR_RGB2HSV)
                        xa,xb=int(p['x0']*zoom),int(p['x1']*zoom)
                        ya,yb=int(p['y0']*zoom),int(p['y1']*zoom)
                        dens=float(_mask(hsv[ya:yb,xa:xb],color).mean())
                        if dens<.065:series.append({'color':color,'label':legend.get(color),'points':points,'mask_density':round(dens,4)})
                if not series:continue
                p['series']=series;p['figure_local_index']=j+1
                results.append(p)
    # A figure may label one category in panel (a) but put a numerical
    # annotation over the same-color trace in panel (b). Borrow a label only
    # from a sibling panel on the same PDF page, with exactly one consistent
    # and scientifically plausible color-to-name mapping.
    from collections import defaultdict
    by_page=defaultdict(list)
    for plot in results:by_page[plot['page']].append(plot)
    bad=re.compile(r'\bEQE\b|\d+(?:\.\d+)?\s*%|^.{1,2}$|T\s*80|T_?\d{2}',re.I)
    for candidates in by_page.values():
        if len(candidates)<2:continue
        for plot in candidates:
            for trace in plot['series']:
                current=trace.get('label')
                if current and not bad.search(current):continue
                alternatives={other['label'] for sibling in candidates if sibling is not plot
                              for other in sibling['series']
                              if other['color']==trace['color'] and other.get('label') and not bad.search(other['label'])}
                if len(alternatives)==1:
                    trace['label']=alternatives.pop()
                    trace['label_source']='MATCHED_SAME_COLOR_SIBLING_PANEL'
        for plot in candidates:
            unique=[]
            for trace in sorted(plot['series'],key=lambda x:len(x['points']),reverse=True):
                arr=np.asarray(trace['points'],float)
                duplicate=False
                for previous in unique:
                    prev=np.asarray(previous['points'],float)
                    if trace.get('label')!=previous.get('label'):continue
                    lo=max(arr[0,0],prev[0,0]);hi=min(arr[-1,0],prev[-1,0]);
                    if hi-lo<.6*min(arr[-1,0]-arr[0,0],prev[-1,0]-prev[0,0]):continue
                    xs=np.linspace(lo,hi,32)
                    difference=np.median(np.abs(np.interp(xs,arr[:,0],arr[:,1])-np.interp(xs,prev[:,0],prev[:,1])))
                    if difference<.018*plot['y_max']:
                        duplicate=True;break
                if not duplicate:unique.append(trace)
            plot['series']=unique
    return results

def write_overlays(pdf,plots,out_dir):
    out_dir=Path(out_dir);out_dir.mkdir(parents=True,exist_ok=True)
    colors={'red':(240,40,50),'blue':(20,50,230),'green':(20,200,50),'cyan':(10,210,210),'magenta':(210,30,210),'violet':(160,55,210),'orange':(240,140,40),'gray':(110,110,110),'black':(25,25,25)}
    with fitz.open(str(pdf)) as doc:
        for i,plot in enumerate(plots):
            page=doc[plot['page']-1];scale=3; pix=page.get_pixmap(matrix=fitz.Matrix(scale,scale),alpha=False)
            im=np.frombuffer(pix.samples,np.uint8).reshape(pix.height,pix.width,pix.n).copy()[:,:,:3]
            x0,x1,y0,y1=(round(plot[k]*scale) for k in ('x0','x1','y0','y1'))
            cv2.rectangle(im,(x0,y0),(x1,y1),(255,130,0),2)
            for series in plot['series']:
                rgb=colors[series['color']]
                for t,r in series['points'][::max(1,len(series['points'])//75)]:
                    x=int((t-plot['x_intercept'])/plot['x_slope']*scale)
                    y=int((r-plot['y_intercept'])/plot['y_slope']*scale)
                    cv2.circle(im,(x,y),3,rgb,-1)
            path=out_dir/f'{Path(pdf).stem}_p{plot["page"]}_{i+1}.png'
            cv2.imwrite(str(path),cv2.cvtColor(im,cv2.COLOR_RGB2BGR))
    return len(plots)
