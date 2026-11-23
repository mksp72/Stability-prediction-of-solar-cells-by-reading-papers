"""Conservative, figure-only raster J-V digitization fallback (v18.6.4).

Only runs on image-only PDF pages with OCR-calibrated V and J tick scales,
color-matched device labels, and published same-category J-V PCE to cross-check
maximum power. Never uses a DOI, page map or hard-coded measured value.

We do NOT substitute inferred Jsc from an EQE/MPP annotation for J-V Jsc.
Recorded FF is computed from source J-V PCE and digitized Jsc/Voc and is
marked *derived*, not falsely called a measured FF.
"""
from __future__ import annotations
from pathlib import Path
import json,re,math
import fitz,numpy as np,cv2
from general_curve_reader import ocr_page_words,_fit_axis,_mask

COLOR_CLASSES=('red','violet','magenta','blue','cyan','green','orange')

def _key(s):return re.sub(r'[^a-z0-9]','',str(s).lower())

def _numeric(text,voltage=False):
    t=text.strip('()[]{}.,:')
    if voltage:
        if re.fullmatch(r'0[2468]',t):t='0.'+t[1]
        elif t=='12':t='1.2'
    else:
        # A thin graph gridline can fuse to the end of a tick in OCR.
        if re.fullmatch(r'(?:[123]0|[123]5)\d',t):t=t[:-1]
    try:return float(t)
    except (TypeError,ValueError):return None

def _candidate_axes(words):
    """Locate Voltage x-axis by OCR, then align its own current-density ticks."""
    found=[]
    for label in words:
        if not re.fullmatch(r'Voltage|Potential',label[4],re.I):continue
        cx=.5*(label[0]+label[2]);yt=label[1]
        xt=[]
        for v in words:
            x,y,xr,yr,txt=v[:5];n=_numeric(txt,voltage=True)
            if n is not None and 0<=n<=2.5 and abs(.5*(x+xr)-cx)<105 and yt-27<y<yt-2:
                xt.append((.5*(x+xr),.5*(y+yr),n,txt))
        if len(xt)<4:continue
        # Discard adjacent plotted panels by monotonic tick-value resets.
        ordered=sorted(xt,key=lambda v:v[0]);segments=[];segment=[]
        for q in ordered:
            if segment and q[2]<segment[-1][2]:segments.append(segment);segment=[]
            segment.append(q)
        if segment:segments.append(segment)
        for seg in segments:
            xf=_fit_axis(seg,'x',positive=True)
            if not xf:continue
            vals=[x[2] for x in xf['ticks']]
            if min(vals)>0.1 or max(vals)<.8 or max(vals)>2.4:continue
            sx,ix=xf['slope'],xf['intercept'];xzero=-ix/sx
            if not(xf['ticks'][0][0]-20<xzero<xf['ticks'][-1][0]+5):continue
            bottom=float(np.median([a[1] for a in xf['ticks']]))-6
            yticks=[]
            for w in words:
                x,y,xr,yr,txt=w[:5];v=_numeric(txt)
                if v is None or not 0<=v<=45:continue
                xx=.5*(x+xr);yy=.5*(y+yr)
                if xzero-18<xx<xzero-3 and bottom-140<yy<bottom-15:
                    yticks.append((xx,yy,v,txt))
            yf=_fit_axis(yticks,'y',positive=True)
            if not yf:continue
            sy,iy=yf['slope'],yf['intercept']
            if not(-.9<sy<-.04):continue
            if abs(sy*bottom+iy)>max(3,.13*max(q[2] for q in yf['ticks'])):continue
            # The last labeled tick can be vertically offset a few pixels
            # and omitted by the best-line subset. Include it only when its
            # value agrees with the calibrated horizontal coordinate.
            vmax=max(vals)
            for xx,yy,vv,_ in seg:
                if vv>vmax and abs(sx*xx+ix-vv)<max(.035,.05*vv) and abs(yy-bottom-6)<7:
                    vmax=vv
            # A trace may run slightly ABOVE the highest printed J tick.
            # Include that physical headroom, rather than cropping Jsc.
            ytop=(max(q[2] for q in yf['ticks'])+min(5,.20*max(q[2] for q in yf['ticks']))-iy)/sy
            xmin=xzero;xmax=(vmax-ix)/sx
            if xmax-xmin<50 or bottom-ytop<35:continue
            found.append({'x_slope':sx,'x_intercept':ix,'y_slope':sy,'y_intercept':iy,
                'xmin':xmin,'xmax':xmax,'ytop':ytop,'ybottom':bottom,'vmax':vmax,
                'jmax':max(q[2] for q in yf['ticks']),'tick_error_x':xf['error'],
                'tick_error_y':yf['error']})
    # Cluster same plot so 'Voltage (V)' cannot duplicate it.
    unique=[]
    for a in found:
        if not any(abs(a['xmin']-b['xmin'])<8 and abs(a['ybottom']-b['ybottom'])<8 for b in unique):unique.append(a)
    return unique

def _labels(words,hsv,axes,rows,scale):
    linked={}
    for word in words:
        raw=word[4].strip(' :;,.')
        if len(raw)<3:continue
        for row in rows:
            name=str(row.get('device_category') or '')
            if _key(raw)!=_key(name):continue
            x0,y0,x1,y1=word[:4]
            if not(axes['xmin']+2<x0<axes['xmax'] and axes['ytop']<y0<axes['ybottom']):continue
            a,b,c,d=[int(v*scale) for v in (x0,y0,x1,y1)]
            roi=hsv[max(0,b):max(b+1,d),max(0,a):max(a+1,c)]
            if not roi.size:continue
            scores=sorted(((int(_mask(roi,color).sum()),color) for color in COLOR_CLASSES),reverse=True)
            if not scores or scores[0][0]<max(20,roi.shape[0]*roi.shape[1]*.06):continue
            if len(scores)>1 and scores[0][0]<scores[1][0]*1.30:continue
            color=scores[0][1]
            if name not in linked or scores[0][0]>linked[name][0]:linked[name]=(scores[0][0],color)
    return {name:color for name,(score,color) in linked.items()}

def _infer_jv(hsv,axes,color,pce,scale):
    mask=_mask(hsv,color)
    sx,ix,sy,iy=[axes[v] for v in ('x_slope','x_intercept','y_slope','y_intercept')]
    ytop=max(0,int(axes['ytop']*scale));ybottom=min(mask.shape[0],int(axes['ybottom']*scale))
    if ybottom-ytop<80:return None
    def js_at(v):
        x=int((v-ix)/sx*scale)
        if x-4<0 or x+5>mask.shape[1]:return np.array([],float)
        ys=np.flatnonzero(np.max(mask[ytop:ybottom,x-4:x+5],axis=1)>0)
        if len(ys)<2:return np.array([],float)
        out=sy*(ytop+ys)/scale+iy
        return out[np.isfinite(out)]
    plate=[];curve=[]
    # The flat part nearest the knee avoids OCR color-legend annotations in
    # the left half, but is still well before the open-circuit rollover.
    for v in np.arange(max(.72,.66*axes['vmax']),.985*axes['vmax'],.005):
        js=js_at(v)
        high=js[(js>.7*axes['jmax']) & (js<axes['jmax']*1.25)]
        if len(high)>=2:
            current=float(np.median(high))
            if v<.8*axes['vmax']:plate.append(current)
            curve.append((v,current))
    if len(plate)<10 or len(curve)<18:return None
    jsc=float(np.median(plate))
    if not 5<=jsc<=42:return None
    tails=[]
    for v in np.arange(.72*axes['vmax'],min(axes['vmax']*1.04,2.5),.003):
        js=js_at(v)
        if np.any((js>.3)&(js<2.5)):tails.append(v)
    if len(tails)<7:return None
    voc=float(np.percentile(tails,94))
    if not .45<=voc<=2.4:return None
    ff=100*pce/(voc*jsc)
    pmax=max(v*j for v,j in curve)
    if not(50<=ff<=94) or abs(pmax-pce)>max(1.4,.08*pce):return None
    return {'voc0':round(voc,3),'jsc0':round(jsc,3),'ff0':round(ff,3),
       'jv_graph_pmax_pct':round(pmax,3),'jv_graph_pce_closure_error_pctpts':round(pmax-pce,3),
       'jv_graph_flat_samples':len(plate),'jv_graph_tail_samples':len(tails)}

def recover(rows,main_pdf,supplements=()):
    updates=0
    for src in [main_pdf,*supplements]:
        if Path(src).suffix.lower()!='.pdf':continue
        try:doc=fitz.open(str(src))
        except Exception:continue
        with doc:
            for page in doc:
                if len(page.get_text('words'))>12 or not page.get_images():continue
                # The J-V raster fallback is used only where the PDF itself
                # has no selectable source data. Never overwrite source tables.
                words=ocr_page_words(page,zoom=3)
                if not any(re.fullmatch(r'Voltage|Potential',w[4],re.I) for w in words):continue
                axeslist=_candidate_axes(words)
                if not axeslist:continue
                zoom=8
                pix=page.get_pixmap(matrix=fitz.Matrix(zoom,zoom),alpha=False)
                image=np.frombuffer(pix.samples,np.uint8).reshape(pix.height,pix.width,pix.n).copy()[:,:,:3]
                hsv=cv2.cvtColor(image,cv2.COLOR_RGB2HSV)
                for axes in axeslist:
                    linked=_labels(words,hsv,axes,rows,zoom)
                    for row in rows:
                        label=str(row.get('device_category') or '')
                        if label not in linked:continue
                        if all(row.get(k) not in (None,'') for k in ('voc0','jsc0','ff0')):continue
                        pce=row.get('pce_jv_pct')
                        try:pce=float(pce)
                        except (ValueError,TypeError):continue
                        metrics=_infer_jv(hsv,axes,linked[label],pce,zoom)
                        if not metrics:continue
                        source=f'{Path(src).name} p.{page.number+1}, OCR calibrated raster J-V axis'
                        for field in ('voc0','jsc0','ff0'):
                            if row.get(field) not in (None,''):continue
                            row[field]=metrics[field]
                            row[field+'_source']=(source+', inferred from the plotted curve' if field!='ff0' else
                                source+', FF derived using same-category J-V PCE and digitized Voc/Jsc')
                            updates+=1
                        row['raster_jv_recovery_status']='GRAPH_DIGITIZED_VALIDATED_BY_PUBLISHED_PCE'
                        row['raster_jv_recovery_confidence_json']=json.dumps(
                            {'source':source,'legend_color':linked[label],**metrics},ensure_ascii=False)
    return updates
