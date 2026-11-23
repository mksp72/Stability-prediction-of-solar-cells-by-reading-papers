"""v18.6.8: label-preserving, marker-only PV aging reader for raster/vector PDF plots.

Experimental source-calibrated reader. Never mistakes connecting lines or error
bars for measurement points. Absolute PCE and retention remain separate arrays.
Unassigned/ambiguous devices are retained in audits, not promoted into ML.
"""
from __future__ import annotations
import re, json, csv
from pathlib import Path
from collections import defaultdict
import cv2, fitz, numpy as np
from general_curve_reader import ocr_page_words, numeric_words, _mask
from scipy.ndimage import maximum_filter1d
from axis_semantics import axis_kind, time_kind, retention, calendar_axis_candidates

COLOR_SET=('red','blue','orange','green','magenta','violet','cyan')
STAB=re.compile(r'\bstability\b|\bdegradation\b|\baging\b|\bageing\b|\bMPP\b|\bMPPT\b|\bretention\b|\bISOS[-– ]?[DLO]|\bstorage\b',re.I)
NUM=re.compile(r'^\d+(?:[,.]\d+)?$')

def x_axis_candidates(words):
    """RANSAC numeric tick rows, robust to missing 0 and OCR 1000->4000."""
    nums=[q for q in numeric_words(words) if 0<=q[2]<=30000]
    bins=[]
    for q in sorted(nums,key=lambda z:z[1]):
        group=next((g for g in bins if abs(np.median([t[1] for t in g])-q[1])<6.5),None)
        if group is None:bins.append([q])
        else:group.append(q)
    axes=[]
    for line in bins:
        if len(line)<3:continue
        # Fit each plotted panel independently, including adjacent plots with
        # identical x-tick baselines. A single malformed OCR token is an outlier.
        remaining=[q for q in line if q[2]>=0]
        for repetition in range(5):
            possibilities=[]
            for i,a in enumerate(remaining):
                for b in remaining[i+1:]:
                    dx=b[0]-a[0]
                    if dx<38 or dx>400 or b[2]<=a[2]+20:continue
                    m=(b[2]-a[2])/dx;intercept=a[2]-m*a[0]
                    zero=-intercept/m
                    if zero<15 or zero>800:continue
                    inliers=[q for q in remaining if zero-8<=q[0]<=min(zero+350,zero+30000/m)
                        and abs(m*q[0]+intercept-q[2])<max(8,.014*max(q[2],1))]
                    if len(inliers)<3:continue
                    xlo=min(q[0] for q in inliers);xhi=max(q[0] for q in inliers)
                    if xhi-xlo<45 or zero>xlo+15 or zero<xlo-60:continue
                    yy=np.median([q[1] for q in inliers])-5
                    context=' '.join(q[4] for q in words if zero-10<q[0]<xhi+20 and yy<q[1]<yy+32)
                    if not re.search(r'(?i)\btime\b|\bhours?\b',context):
                        # Some stacked panels print the time unit once under
                        # the bottom subplot. Accept numeric time ticks as a
                        # REVIEW candidate only; the y-axis still must pass.
                        if not re.search(r'(?i)\btime\b|\bhours?\b|aging|stability|MPP|MPPT',
                                ' '.join(w[4] for w in words)) or max(q[2] for q in inliers)<48:continue
                    if re.search(r'(?i)\b(?:ns|ms|seconds?|minutes?)\b',context):continue
                    # Penalize sweeping two distinct panels into one axis.
                    score=(len(inliers),-abs(zero-xlo),-abs(m),xhi-xlo)
                    possibilities.append((score,m,intercept,inliers,yy,zero,xhi))
            if not possibilities:break
            score,m,b,ins,yy,zero,xhi=max(possibilities,key=lambda a:a[0])
            axis=dict(x0=float(zero),x1=float(xhi),bottom=float(yy),x_m=float(m),x_b=float(b),x_max=float(m*xhi+b),tick_count=len(ins))
            if not any(abs(z['x0']-zero)<10 and abs(z['bottom']-yy)<8 for z in axes):axes.append(axis)
            remaining=[q for q in remaining if all(abs(q[0]-a[0])>3 or abs(q[1]-a[1])>3 for a in ins)]
    return sorted(axes,key=lambda a:(a['bottom'],a['x0']))

def y_axis(words,axis,page):
    # IMPORTANT: allow truncated y axes (e.g. 5..20 PCE) and 0..1 retention.
    nums=[q for q in numeric_words(words) if 0<=q[2]<=150]
    left,bottom=axis['x0'],axis['bottom']
    # Graph-axis OCR frequently fuses tick labels with small grid ticks:
    #  '—16', '~14', 'O12'. Repair only in the narrow LEFT-TICK band,
    #  never in running article text or legends, and retain numeric checks.
    for q in words:
        token=q[4].strip()
        if not (left-43<q[0]<left-5 and bottom-190<q[1]<bottom-5):continue
        clean=re.sub(r'^[\s~—–−_\|]+','',token).replace('O','0').replace('o','0')
        if re.fullmatch(r'0?\d{1,2}(?:[.,]\d+)?',clean):
            val=float(clean.replace(',','.'))
            if 0<=val<=110 and not any(abs(z[0]-q[0])<2 and abs(z[1]-q[1])<2 for z in nums):
                nums.append((q[0],q[1],val,q[3],q[4]))
    above=max([a['bottom'] for a in page if a['bottom'] < bottom-25 and abs(a['x0']-left)<25],default=bottom-170)
    lo=max(5,above+9,bottom-185)
    groups=[]
    for q in nums:
        if not (left-39<q[0]<left-4 and lo<q[1]<bottom-7):continue
        g=next((gg for gg in groups if abs(np.median([a[0] for a in gg])-q[0])<7),None)
        if g is None:groups.append([q])
        else:g.append(q)
    solutions=[]
    from itertools import combinations
    for group in groups:
        if len(group)<2:continue
        for k in (3,2):
            for chosen in combinations(group,k):
                yy=np.array([t[1] for t in chosen]); vv=np.array([t[2] for t in chosen]);
                if np.ptp(yy)<12 or len(set(vv))<k:continue
                m,b=np.polyfit(yy,vv,1)
                if m>=0:continue
                errs=abs(yy*m+b-vv)
                if max(errs)>max(.22,.032*np.ptp(vv)):continue
                inliers=[q for q in group if abs(m*q[1]+b-q[2])<max(.18,.03*max(vv))]
                # Allow absolute-PCE lower bound at 5 or higher.
                y_bottom=m*bottom+b
                if y_bottom < -2 or y_bottom > max(vv)+3:continue
                y_top=max(lo,min(q[1] for q in inliers)-35)
                if bottom-y_top<35 or bottom-y_top>185:continue
                score=(len(inliers),np.ptp([q[1] for q in inliers]),-max(errs))
                solutions.append((score,dict(y0=float(y_top),y1=float(bottom),y_m=float(m),y_b=float(b),y_min=float(m*bottom+b),y_max=float(m*y_top+b),y_tick_count=len(inliers))))
    return max(solutions,key=lambda v:v[0])[1] if solutions else None

def markers(image,axis,color,scale=4):
    """Detect *symbols* by 2-D shape correlation, never follow connecting paths.

    Candidate glyphs must match a filled diamond, open diamond, open triangle,
    filled circle, or filled square. Pure thin connector/error-bar segments do
    not meet the two-dimensional template threshold.
    """
    hsv=cv2.cvtColor(image,cv2.COLOR_RGB2HSV)
    mask=_mask(hsv,color).astype(np.float32)
    xa=max(0,int((axis['x0']+3)*scale));xb=min(image.shape[1],int((axis['x1']-3)*scale))
    ya=max(0,int((axis['y0']+3)*scale));yb=min(image.shape[0],int((axis['y1']-2)*scale))
    if xb-xa<50 or yb-ya<30:return []
    roi=mask[ya:yb,xa:xb]
    # Symbol sizes adapt to raster resolution (3--7 points in PDF space).
    candidates=[]
    for rad in (max(3,int(scale*.9)),max(4,int(scale*1.35))):
        sz=2*rad+1
        for shape in ('diamond','triangle','circle','square'):
            templ=np.zeros((sz,sz),np.uint8)
            if shape=='diamond':
                pts=np.array([[rad,0],[sz-1,rad],[rad,sz-1],[0,rad]],np.int32)
                cv2.fillConvexPoly(templ,pts,1)
            elif shape=='triangle':
                pts=np.array([[rad,0],[sz-1,sz-1],[0,sz-1]],np.int32)
                cv2.fillConvexPoly(templ,pts,1)
            elif shape=='circle':cv2.circle(templ,(rad,rad),rad,1,-1)
            else:templ[:,:]=1
            for open_shape in (False,True):
                if open_shape:
                    outline=np.zeros_like(templ)
                    if shape in ('diamond','triangle'):cv2.polylines(outline,[pts],True,1,2)
                    elif shape=='circle':cv2.circle(outline,(rad,rad),rad,1,2)
                    else:cv2.rectangle(outline,(0,0),(sz-1,sz-1),1,2)
                    templ=outline
                count=int(templ.sum())
                if count<10:continue
                matched=cv2.filter2D(roi,-1,templ.astype(np.float32),borderType=cv2.BORDER_CONSTANT)/count
                # Rule out thin lines (including error bars) using average
                # occupancy across orthogonal quadrants of the symbol window.
                local_max=cv2.dilate(matched,np.ones((max(5,rad*2),)*2,np.uint8))
                yy,xx=np.where((matched>=.66)&(matched==local_max))
                for cy,cx in zip(yy,xx):
                    if cx<rad or cy<rad or cx+rad>=roi.shape[1] or cy+rad>=roi.shape[0]:continue
                    patch=roi[cy-rad:cy+rad+1,cx-rad:cx+rad+1]
                    quads=[patch[:rad,:rad].sum(),patch[:rad,rad+1:].sum(),patch[rad+1:,:rad].sum(),patch[rad+1:,rad+1:].sum()]
                    if min(quads)<max(1,rad*.4):continue
                    xxp=(cx+xa)/scale;yyp=(cy+ya)/scale
                    t=axis['x_m']*xxp+axis['x_b'];v=axis['y_m']*yyp+axis['y_b']
                    if 0<=t<=axis['x_max']*1.02 and v>=0:
                        candidates.append((float(t),float(v),float(matched[cy,cx]),float(xxp),float(yyp)))
                if open_shape:break
    candidates.sort(key=lambda v:v[2],reverse=True)
    selected=[]
    for pt in candidates:
        # Nonmaximum suppression in physical plot position. Prevent error-bar
        # caps and duplicate template radii from becoming additional points.
        if any(abs(pt[3]-prev[3])<2.6 and abs(pt[4]-prev[4])<2.6 for prev in selected):continue
        selected.append(pt)
    selected.sort(key=lambda v:v[0])
    clusters=[]
    for pt in selected:
        group=next((g for g in clusters if abs(np.median([v[0] for v in g])-pt[0])<.013*axis['x_max']),None)
        if group is None:clusters.append([pt])
        else:group.append(pt)
    return [max(g,key=lambda v:v[2]) for g in clusters]

def infer_label(words,image,axis,color,scale):
    # A label must be visibly near the same-color legend swatch; never infer
    # device material from graph line color alone.
    hsv=cv2.cvtColor(image,cv2.COLOR_RGB2HSV)
    all_labels=[]
    for word in words:
        token=word[4].strip()
        if not re.search('[A-Za-z]',token) or len(token)<2 or token.lower() in ('time','pce','sop','spo','control:','unencapsulated'):continue
        x=(word[0]+word[2])/2; y=(word[1]+word[3])/2
        if not (axis['x0']+15<x<axis['x1']+10 and axis['y0']<y<axis['y1']):continue
        left=max(0,int((word[0]-23)*scale));right=max(0,int((word[0]-2)*scale));
        top=max(0,int((y-3)*scale));bottom=min(image.shape[0],int((y+3)*scale))
        if right<=left or bottom<=top:continue
        score=int(_mask(hsv[top:bottom,left:right],color).sum())
        if score>=max(8,scale*3):all_labels.append((score,token))
    if not all_labels:return None
    label=max(all_labels,key=lambda a:a[0])[1]
    if re.fullmatch(r'PCE|Voc|Jsc|FF|Normaliz(?:ed)?|SP[OC]|Time|Hours?',label,re.I):return None
    return label

def detect_pdf(pdf,source_name=None):
    """Audited candidate marker observations, never silently promoted."""
    records=[]; pdf=Path(pdf)
    with fitz.open(pdf) as doc:
        for page in doc:
            text=page.get_text('text')
            if not STAB.search(text):continue
            captions=re.findall(r'(?is)(?:Fig\.?|Figure)\s*(?:S\s*)?\d+\s*[.|:]?[^\n]{0,1000}',text)
            if not any(STAB.search(c) or (re.search(r'(?i)\btime\b|\bhours?\b',c) and re.search(r'(?i)\bPCE\b|\bSPO\b|\bpower\b|\bMPP\b',c)) for c in captions):
                # Some articles place a shared figure caption in a separate
                # end section. Preserve possible panels, but require local
                # time+PCE semantics and numeric axes below.
                if not re.search(r'(?i)time.{0,65}(?:PCE|MPP|retention|stability)',text,re.S):continue
            if not page.get_images() and len(page.get_drawings())<45:continue
            words=ocr_page_words(page,zoom=3)
            if not words:continue
            axes=x_axis_candidates(words)
            # Calendar tick labels (e.g. March 2024 or 25/11/2024) may
            # represent stability aging. Never mistake them for numeric hours.
            axes.extend(calendar_axis_candidates(words,text))
            if not axes:continue
            scale=4
            pix=page.get_pixmap(matrix=fitz.Matrix(scale,scale),alpha=False)
            im=np.frombuffer(pix.samples,np.uint8).reshape(pix.height,pix.width,pix.n)[:,:,:3].copy()
            for idx,ax in enumerate(axes):
                yy=y_axis(words,ax,axes)
                if not yy:continue
                ax.update(yy)
                # Must fit plausible metric: either absolute percentage PCE or normalized values.
                if not (0<ax['y_max']<=110):continue
                # Read the shared Y label from neighboring figure panels if
                # a subplot omits it. Numeric magnitude is only a fallback.
                nearby_y=' '.join(q[4] for q in words if
                    ax['x0']-65<q[0]<ax['x0']+8 and ax['y0']-12<q[1]<ax['y1']+8)
                kind,axis_conf=axis_kind([ax['y_min'],ax['y_max']],nearby_y,text[:800])
                if kind in ('UNKNOWN','NON_PCE_AXIS'):continue
                ax['axis_kind']=kind; ax['axis_confidence']=axis_conf
                for color in COLOR_SET:
                    data=markers(im,ax,color,scale=scale)
                    if len(data)<4:continue
                    # More than 35 apparent symbols on an aging plot usually
                    # indicates a continuous waveform/line, not discrete glyphs.
                    # Keep those for the existing continuous-signal digitizer.
                    if len(data)>35:continue
                    # Remove legend swatches; unlike experimental points, these
                    # sit immediately to the left of text within the same plot.
                    legend_boxes=[q for q in words if re.search(r'[A-Za-z]{2,}',q[4])
                        and ax['x0']+10<q[0]<ax['x1']+25 and ax['y0']<q[1]<ax['y1']]
                    data=[pt for pt in data if not any(0 < q[0]-pt[3] < 37 and abs((q[1]+q[3])/2-pt[4])<6
                        for q in legend_boxes)]
                    if len(data)<4:continue
                    label=infer_label(words,im,ax,color,scale)
                    vals=[{'time_h':round(max(0,t),3),'measured_y':round(y,5)} for t,y,_,xx,yy in data]
                    rawvals=np.array([v['measured_y'] for v in vals],dtype=float)
                    retained,normalization=retention(rawvals,kind)
                    retained=np.asarray(retained if retained is not None else np.full(len(rawvals),np.nan))
                    # This is not a device-specific burn-in or J-V initial PCE.
                    # Keep original absolute values so future verification can
                    # replace this provisional reference without data loss.
                    for item,rv in zip(vals,retained):item['provisional_retention_pct']=round(float(rv),4) if np.isfinite(rv) else None
                    panel_tokens=[]
                    for word in words:
                        token=word[4].strip(' :.')
                        mm=re.fullmatch(r'\(?([a-f])\)?',token,re.I)
                        if not mm:continue
                        if ax['x0']-57 < word[0] < ax['x0']+10 and ax['y0']-40 < word[1] < ax['y0']+17:
                            panel_tokens.append((abs(word[0]-ax['x0']+18)+abs(word[1]-ax['y0']+14),mm.group(1).lower()))
                    panel_letter=min(panel_tokens)[1] if panel_tokens else ''
                    records.append({'source_file':source_name or pdf.name,'source_page':page.number+1,
                      'subplot_letter':panel_letter,
                      'plot_number':idx+1,'color':color,'legend_candidate':label or '',
                      'point_count':len(vals),'raw_axis_lower':round(ax['y_min'],3),
                      'raw_axis_upper':round(ax['y_max'],3),'raw_measurement_type':{'NORMALIZED_FRACTION':'NORMALIZED_0_TO_1','ABSOLUTE_PCE_PERCENT':'ABSOLUTE_PCE_OR_POWER','RETENTION_PERCENT':'RETENTION_PERCENT'}.get(kind,'UNKNOWN'),
                      'provisional_normalization_basis':normalization, 'axis_confidence':axis_conf, 'x_time_unit':ax.get('x_unit','HOURS'), 'reference_date':ax.get('reference_date',''),
                      'status':'MARKER_ONLY_REQUIRES_DEVICE_LABEL_CONFIRMATION' if not label else 'MARKER_ONLY_AXIS_AND_LABEL_CANDIDATE_NOT_FITTED',
                      'axis_json':json.dumps({k:round(v,5) for k,v in ax.items() if isinstance(v,(int,float))}),
                      'points_json':json.dumps(vals),
                      'caption_context':text[max(0,text.lower().find('fig.')):][:600]})
    return records

def save(records,path):
    cols=('source_file','source_page','plot_number','color','legend_candidate','subplot_letter','point_count','raw_axis_lower','raw_axis_upper','raw_measurement_type','provisional_normalization_basis','status','axis_json','points_json','caption_context','axis_confidence','x_time_unit','reference_date')
    path=Path(path);path.parent.mkdir(parents=True,exist_ok=True)
    with path.open('w',newline='',encoding='utf-8-sig') as f:
        w=csv.DictWriter(f,fieldnames=cols,extrasaction='ignore');w.writeheader();w.writerows(records)
