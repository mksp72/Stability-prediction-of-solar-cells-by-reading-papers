"""Conservative, publication-independent raster EQE + integrated-Jsc recovery.

This additive fallback targets PDF pages with image-only EQE panels and no
searchable figure text. It calibrates wavelength tick locations, right-hand
integrated-current ticks and the actual red-edge from pixels. It never uses
DOI, page coordinates, a value lookup, or JV Jsc as a replacement.

Only *single explicitly identified* EQE devices are assigned. Both the red-edge
and cumulative-current plateau must be found on the SAME calibrated plot.
Nothing is filled when color, axes, device attribution or closure is ambiguous.
"""
from __future__ import annotations
import math,re
from pathlib import Path
import fitz, cv2, numpy as np
from general_curve_reader import ocr_page_words

HC=1239.841984
PALETTES={
    'red': [(0,12),(170,179)], 'orange':[(13,26)],
    'green':[(31,85)], 'cyan':[(86,102)], 'blue':[(103,135)],
    'magenta':[(136,169)]}


def _norm(s):return re.sub(r'[^a-z0-9]','',str(s).lower())
def _num(s):
    s=str(s).strip(' \t\n:;,.()[]<>~=')
    try:return float(s)
    except (ValueError,TypeError):return None

def _tick(word,axis):
    text=str(word[4]);val=_num(text)
    if axis=='y' and re.fullmatch(r'[^0-9]*[123]0[0-9][^0-9]*',text):return float(re.search(r'[123]0',text).group())
    if val is not None:return val
    if axis=='y':
        # OCR may conflate a red tick with a crossing curve/grid line.
        # Restricted to expected multiples of 10 and tested by 3+ point fit.
        s=re.sub(r'[^0-9]','',text)
        if s in ('200','100','300'):return float(s[:2])
        if s in ('10','20','30','0'):return float(s)
    return None

def _axes(words):
    """Find paired x-wavelength and right integrated-current axes by geometry."""
    axes=[]
    for label in words:
        if _norm(label[4])!='wavelength':continue
        cx=(label[0]+label[2])/2.; ly=label[1]
        xs=[]
        for w in words:
            val=_tick(w,'x');x=(w[0]+w[2])/2.;y=(w[1]+w[3])/2.
            if val is not None and 250<=val<=1600 and abs(x-cx)<110 and ly-31<y<ly-1:
                xs.append((x,y,val))
        if len(xs)<4:continue
        # Only the column's own shared tick row; avoids adjacent panel numbers.
        groups=[]
        for q in xs:
            group=next((g for g in groups if abs(np.median([a[1] for a in g])-q[1])<4),None)
            if group is None:groups.append([q])
            else:group.append(q)
        for group in groups:
            if len(group)<4:continue
            xx=np.array([q[0] for q in group]);vv=np.array([q[2] for q in group]);
            if np.ptp(xx)<35 or np.ptp(vv)<200:continue
            sx,ix=np.polyfit(xx,vv,1)
            if sx<=0 or np.median(abs(sx*xx+ix-vv))>4:continue
            wlmin=float(min(vv));wlmax=float(max(vv))
            if not 300<=wlmin<=700 or not 700<=wlmax<=1600:continue
            xleft=float((wlmin-ix)/sx);xright=float((wlmax-ix)/sx)
            ybottom=float(np.median([q[1] for q in group])-6)
            # Require three independently printed ticks on the *right* y-axis;
            # a left EQE % tick is not sufficient to infer the current scale.
            yt=[]
            for w in words:
                v=_tick(w,'y');x=(w[0]+w[2])/2.;y=(w[1]+w[3])/2.
                if v is None or v<0 or v>45:continue
                if xright-8<x<xright+25 and ybottom-125<y<ybottom-8:
                    yt.append((x,y,v))
            fits=[]
            from itertools import combinations
            for three in combinations(yt,3):
                y=np.array([q[1] for q in three]);val=np.array([q[2] for q in three]);
                if np.ptp(y)<25 or np.ptp(val)<10 or len(set(val))<3:continue
                sl,ic=np.polyfit(y,val,1)
                if sl>=-.07 or not -1.5<=sl*ybottom+ic<=3.0:continue
                err=max(abs(sl*y+ic-val))
                if err>1.5:continue
                fits.append((err,sl,ic,three))
            if not fits:continue
            _,sy,iy,ticks=min(fits,key=lambda t:t[0]);top=(30-iy)/sy
            if not 40<ybottom-top<180:continue
            if xright-xleft<60:continue
            axes.append({'x0':xleft,'x1':xright,'y0':top,'y1':ybottom,
                         'sx':float(sx),'ix':float(ix),'sy':float(sy),'iy':float(iy),
                         'wlmax':wlmax,'wlmin':wlmin,'y_ticks':[q[2] for q in ticks]})
    return axes

def _colored_masks(page,scale):
    pix=page.get_pixmap(matrix=fitz.Matrix(scale,scale),alpha=False)
    rgb=np.frombuffer(pix.samples,dtype=np.uint8).reshape(pix.height,pix.width,pix.n)[:,:,:3]
    hsv=cv2.cvtColor(rgb,cv2.COLOR_RGB2HSV)
    result={}
    for name,ranges in PALETTES.items():
        mask=np.zeros(hsv.shape[:2],dtype=np.uint8)
        for lo,hi in ranges:
            mask|=cv2.inRange(hsv,np.array([lo,58,75],np.uint8),np.array([hi,255,255],np.uint8))
        result[name]=mask
    return result

def _line_estimates(axis,mask,scale):
    """Integrated-current plateau and EQE λ50 from independent pixel windows."""
    x0,x1,y0,y1=(axis[k]*scale for k in ('x0','x1','y0','y1'))
    if not (x1-x0>180 and y1-y0>120):return None
    # Cumulative current flat region near far-right edge. Exclude border/axis;
    # an EQE spectrum here falls to zero so can be separated vertically.
    xlo=int(x0+.84*(x1-x0));xhi=int(x0+.95*(x1-x0))
    ya=max(0,int(y0+0.10*(y1-y0)));yb=min(mask.shape[0],int(y0+.60*(y1-y0)))
    if xhi>mask.shape[1] or xhi<=xlo or yb<=ya:return None
    # Multiple x-columns required, not isolated red letters or annotation.
    med=[]
    for x in range(xlo,xhi,3):
        ys=np.flatnonzero(np.max(mask[ya:yb,max(0,x-2):x+3],axis=1)>0)
        if len(ys)>0:
            h=np.bincount(ys//max(1,int(scale)),minlength=int((yb-ya)/scale)+2)
            if h.max()>=1:
                peak=np.argmax(h)*max(1,int(scale))+max(1,int(scale))/2
                med.append((x,ya+peak))
    if len(med)<8:return None
    yy=np.array([p[1] for p in med]);
    if np.percentile(yy,75)-np.percentile(yy,25)>5*scale:return None
    # Must have at least a visible continuous plateau over several wavelengths.
    m_y=float(np.median(yy));current=axis['sy']*m_y/scale+axis['iy']
    if not 5<current<35:return None
    # EQE wavelength at 50% of the measured high-energy plateau. Scan a
    # separate low-to-middle EQE band to avoid the integrated-current trace.
    # Estimate the EQE plateau at wavelengths 500--700 nm using bright pixels.
    pxy=[]
    for wl in np.arange(max(500,axis['wlmin']+100),min(700,axis['wlmax']-120),15):
        x=int((wl-axis['ix'])/axis['sx']*scale)
        if not 0<=x<mask.shape[1]:continue
        ya2=max(0,int(y0+0.04*(y1-y0)));yb2=min(mask.shape[0],int(y0+.36*(y1-y0)))
        sy=np.flatnonzero(mask[ya2:yb2,max(0,x-3):min(mask.shape[1],x+4)].max(axis=1)>0)
        if len(sy)>0:pxy.append(float(ya2+np.median(sy)))
    if len(pxy)<5:return None
    plateau_y=float(np.median(pxy));eqe_50_y=.5*(plateau_y+y1)
    # Seek the far-red descending EQE line at its 50% level. Require a
    # contiguous band of colored pixels, not an isolated legend or gridline.
    xlo=int(x0+.70*(x1-x0));xhi=int(x0+.96*(x1-x0));
    ymid=int(eqe_50_y);band=int(3.2*scale)
    patch=mask[max(0,ymid-band):min(mask.shape[0],ymid+band+1),xlo:xhi]
    if patch.size==0:return None
    count=np.count_nonzero(patch,axis=0)
    centers=np.where(np.convolve(count.astype(float),np.ones(max(3,int(scale*3))),'same')>=max(3,scale*2))[0]
    if len(centers)<2:return None
    # Ignore isolated extreme pixels; use densest red-edge cluster.
    hh=np.convolve(count.astype(float),np.ones(max(3,int(scale*4))),'same')
    peak=int(np.argmax(hh))
    use=centers[abs(centers-peak)<10*scale]
    if len(use)<2:return None
    edge_x=(xlo+float(np.median(use)))/scale
    lambda50=axis['sx']*edge_x+axis['ix']
    if not 550<lambda50<1150:return None
    return {'current_ma_cm2':round(float(current),3), 'lambda50_nm':round(float(lambda50),2),
            'bandgap_ev':round(float(HC/lambda50),4), 'samples':len(med)}

def recover(rows,main_pdf,supplements=(),max_raster_pages=6):
    """Mutate only missing matched-category fields. Return recovered-field count."""
    if not rows:return 0
    count=0
    for path in [main_pdf,*supplements]:
        if Path(path).suffix.lower()!='.pdf':continue
        try:doc=fitz.open(str(path))
        except Exception:continue
        with doc:
            scanned=0
            for pno,page in enumerate(doc):
                # Only raster-only figures: native text recoveries already ran.
                if len(page.get_text('text').strip())>130:continue
                if scanned>=max_raster_pages:break
                scanned+=1
                words=ocr_page_words(page,zoom=2.5)
                if not any(_norm(w[4])=='wavelength' for w in words):continue
                axes=_axes(words)
                if not axes:continue
                masks=_colored_masks(page,3)
                for axis in axes:
                    # Require the single device's complete legend label INSIDE
                    # or immediately adjacent to this wavelength plot.
                    hits=[]
                    for row in rows:
                        category=str(row.get('device_category') or '')
                        name=_norm(category)
                        if len(name)<5:continue
                        found=[w for w in words if _norm(w[4])==name and
                               axis['x0']-12<w[0]<axis['x1']+10 and
                               axis['y0']-12<w[1]<axis['y1']+15]
                        if found:hits.append(row)
                    if len(hits)!=1:continue
                    row=hits[0]
                    old_current=row.get('jsc_eqe_integrated');old_gap=row.get('bandgap_ev')
                    if old_current is not None and old_gap is not None:continue
                    candidates=[]
                    for color,mask in masks.items():
                        val=_line_estimates(axis,mask,3)
                        if not val:continue
                        # JV-to-EQE closure is a *cross-check*, never a source.
                        jv=row.get('jsc0')
                        if jv is not None:
                            try:
                                discrepancy=abs(float(jv)-val['current_ma_cm2'])/float(jv)
                                if discrepancy>.13:continue
                            except (TypeError,ValueError,ZeroDivisionError):pass
                        candidates.append((color,val))
                    # Competing colored traces with different estimates are
                    # ambiguous unless a single color has a valid calibration.
                    if len(candidates)!=1:continue
                    color,val=candidates[0]
                    evidence=f'{Path(path).name} p.{pno+1}; raster EQE panel, {color}; calibrated wavelength and integrated-current axes'
                    if old_current is None:
                        row['jsc_eqe_integrated']=val['current_ma_cm2']
                        row['jsc_eqe_source']=evidence
                        row['jsc_eqe_extraction_method']='raster_calibrated_right_axis_plateau'
                        jv=row.get('jsc0')
                        if jv is not None:
                            jv=float(jv);delta=100*(jv-val['current_ma_cm2'])/jv
                            for col,v in [('jsc_difference_signed_pct',delta),('jsc_discrepancy_signed_pct',delta),
                                          ('jsc_difference_absolute_pct',abs(delta)),('jsc_discrepancy_absolute_pct',abs(delta)),
                                          ('jsc_eqe_to_jv_ratio',val['current_ma_cm2']/jv)]:row[col]=round(v,5)
                        count+=1
                    if old_gap is None:
                        row['bandgap_ev']=val['bandgap_ev']
                        row['bandgap_lambda50_nm']=val['lambda50_nm']
                        row['bandgap_source']=evidence+'; EQE red-edge 50%-of-plateau'
                        count+=1
    return count
