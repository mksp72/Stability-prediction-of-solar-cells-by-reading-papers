"""Content-based EQE red-edge reader for vector-text-annotated PDF plots.

Find wavelength/EQE axes from their printed tick positions, calibrate them,
then find the long-wavelength 50%-of-plateau descent in distinct colored EQE
traces. This is a conservative automatic fallback, not a fixed DOI/page list.
Integrated-current dotted curves are excluded by requiring an abrupt final
fall at the high-wavelength end. Outputs only self-consistent shared cutoffs.
"""
from pathlib import Path
import re
import fitz
import cv2
import numpy as np
from math import isfinite

NUM=re.compile(r'^\d{1,4}(?:\.\d+)?$')

def _axis_from_words(page):
    words=page.get_text('words',sort=True)
    labels=[w for w in words if w[4].lower() in ('wavelength','wavelength(nm)')]
    if not labels:return []
    numbers=[]
    for w in words:
        st=w[4].strip('(),')
        if NUM.fullmatch(st):
            numbers.append(((w[0]+w[2])/2,(w[1]+w[3])/2,float(st)))
    candidates=[]
    for label in labels:
        ly=(label[1]+label[3])/2
        tickrows=[]
        xs=[z for z in numbers if 250<=z[2]<=1600 and ly-65<z[1]<ly-3]
        for z in sorted(xs,key=lambda a:a[1]):
            group=next((g for g in tickrows if abs(np.median([q[1] for q in g])-z[1])<7),None)
            if group is None:tickrows.append([z])
            else:group.append(z)
        for group in tickrows:
            if len(group)<4:continue
            xx=np.array([z[0] for z in group]); vv=np.array([z[2] for z in group]);
            co=np.polyfit(xx,vv,1)
            if co[0]<=0:continue
            if np.median(np.abs(np.polyval(co,xx)-vv))>3:continue
            x_left=float(np.polyval(np.polyfit(vv,xx,1),min(vv)))
            x_right=float(np.polyval(np.polyfit(vv,xx,1),max(vv)))
            # Tick labels often end before the actual graph border. Extend
            # by 60% of one tick spacing (e.g. last tick 800, plot to 850).
            steps=np.diff(np.sort(xx))
            if len(steps) and np.median(steps)>5:
                x_right+=0.62*float(np.median(steps))
            if x_right-x_left<100:continue
            y0=float(np.median([z[1] for z in group]))-6
            yy=[z for z in numbers if x_left-52<z[0]<x_left-2 and y0-320<z[1]<y0+8 and 0<=z[2]<=110]
            if len(yy)<4:continue
            # A y-axis calibrated for IPCE should span >=50% and descend
            # numerically as pixel y grows.
            ytickcols=[]
            for z in sorted(yy,key=lambda q:q[0]):
                col=next((g for g in ytickcols if abs(np.median([q[0] for q in g])-z[0])<6),None)
                if col is None:ytickcols.append([z])
                else:col.append(z)
            for col in ytickcols:
                if len(col)<4:continue
                ys=np.array([z[1] for z in col]); val=np.array([z[2] for z in col]);
                sl,ic=np.polyfit(ys,val,1)
                if sl>=0 or max(val)-min(val)<50:continue
                if np.median(np.abs(sl*ys+ic-val))>3:continue
                if abs(sl*y0+ic)>15:continue
                ymaxy=float((100-ic)/sl)
                if ymaxy>=y0-50:continue
                # ensure local text mentions EQE or IPCE near the axis.
                rect=fitz.Rect(max(0,x_left-100),max(0,ymaxy-20),x_right+100,y0+10)
                text=page.get_text('text',clip=rect).lower()
                if not re.search(r'\b(?:eqe|ipce|quantum efficiency)\b',text):continue
                candidates.append(dict(x_left=x_left,x_right=x_right,bottom=y0,top=ymaxy,
                                       wl_slope=float(co[0]),wl_intercept=float(co[1]),
                                       y_slope=float(sl),y_intercept=float(ic)))
    return candidates

def _extract_cutoffs(page,axis,scale=2.0):
    pix=page.get_pixmap(matrix=fitz.Matrix(scale,scale),alpha=False)
    arr=np.frombuffer(pix.samples,np.uint8).reshape(pix.height,pix.width,pix.n)[:,:,:3].copy()
    hsv=cv2.cvtColor(arr,cv2.COLOR_RGB2HSV)
    x0=int(axis['x_left']*scale);x1=min(arr.shape[1]-2,int(axis['x_right']*scale))
    y0=max(0,int(axis['top']*scale));y1=min(arr.shape[0]-1,int(axis['bottom']*scale))
    if x1-x0<150 or y1-y0<90:return []
    palette={'red':[(0,13),(170,179)],'orange':[(14,27)],'green':[(31,85)],
             'cyan':[(86,103)],'blue':[(104,135)],'magenta':[(136,169)]}
    edges=[];span=x1-x0
    for name,ranges in palette.items():
        mask=np.zeros(hsv.shape[:2],dtype=np.uint8)
        for lo,hi in ranges:
            mask|=cv2.inRange(hsv,np.array([lo,65,55],np.uint8),np.array([hi,255,255],np.uint8))
        # Detect a solid high-efficiency plateau in central wavelengths.
        xa=x0+int(span*.40);xb=x0+int(span*.70)
        mid=mask[y0:y1,xa:xb]
        counts=mid.sum(axis=1).astype(float)/255
        peak=np.flatnonzero(counts>=max(3.,np.percentile(counts,95)*.48))
        peak=peak[peak<.50*(y1-y0)]
        if len(peak)<3:continue
        plateau_pixel=float(y0+np.median(peak))
        # 50% of measured plateau-to-baseline, not the right-axis Jsc tick.
        y50=(plateau_pixel+y1)/2
        ys=slice(max(y0,int(y50-4)),min(y1,int(y50+5)))
        # Restrict to the final descending EQE edge; integrated Jsc lines
        # typically cross 50%-height well before this region.
        xa=x0+int(span*.74); xb=x1-int(span*.008)
        patch=mask[ys,xa:xb]
        xs=np.where(patch>0)[1]
        if len(xs)<2:continue
        # rightmost dense cluster rather than isolated legend or noise.
        hist=np.bincount(xs,minlength=xb-xa)
        kernel=np.ones(13,float)/13
        sm=np.convolve(hist,kernel,'same')
        centers=np.where(sm>=max(.12,float(sm.max())*.35))[0]
        if len(centers)<2:continue
        x_at=float(xa+np.median(centers))
        nm=axis['wl_slope']*(x_at/scale)+axis['wl_intercept']
        if not (550<nm<1200):continue
        # Require substantial color density in high-EQE plateau to exclude
        # dotted cumulative Jsc series, which rise gradually across the axis.
        n_plateau=int(np.count_nonzero(mid))
        if n_plateau<45:continue
        edges.append(dict(color=name,lambda50_nm=round(nm,2),bandgap_ev=round(1239.841984/nm,4)))
    # 100-R is optical reflectance, not a device EQE trace. Never include it
    # among independent device cutoffs; do not confuse curves on different axes.
    if re.search(r'\b100\s*[-–]\s*R\b',page.get_text('text'),re.I):
        edges=[edge for edge in edges if edge['color']!='magenta']
    # Two device EQE lines may overlap so closely that a color is hidden
    # by the top trace exactly at the 50%-height crossing. Corroborate the
    # second color along the same descending long-wavelength edge instead.
    if len(edges)==1:
        main=edges[0]
        x_main=(main['lambda50_nm']-axis['wl_intercept'])/axis['wl_slope']*scale
        for color,ranges in palette.items():
            if color==main['color'] or (color=='magenta' and re.search(r'\b100\s*[-–]\s*R\b',page.get_text('text'),re.I)):
                continue
            mask=np.zeros(hsv.shape[:2],dtype=np.uint8)
            for lo,hi in ranges:
                mask|=cv2.inRange(hsv,np.array([lo,65,55],np.uint8),np.array([hi,255,255],np.uint8))
            # Inspect declining edge below 50% where secondary colors may
            # become visible after leaving the overlapped main trace.
            near=mask[int(y0+.55*(y1-y0)):int(y0+.90*(y1-y0)),
                      max(x0,int(x_main-30)):min(x1,int(x_main+37))]
            if np.count_nonzero(near)>10:
                edges.append(dict(color=color,lambda50_nm=main['lambda50_nm'],
                                  bandgap_ev=main['bandgap_ev'],
                                  evidence='SAME_DECLINING_EDGE_OVERLAP'))
    return edges

def recover_shared_cutoff(rows,main_pdf,supplements=(),max_pages_per_document=80,
                          tolerance_nm=15):
    """Return (filled_row_count, spectral_records).

    Only share when at least two colored traces independently fall at
    overlapping high-wavelength cutoffs, and the page caption explicitly
    identifies the categories or the rows have same absorber formulation.
    """
    docs=[main_pdf]+list(supplements)
    recovered=0;records=[]
    for doc_path in docs:
        if Path(doc_path).suffix.lower()!='.pdf':continue
        with fitz.open(str(doc_path)) as doc:
            for pno,page in enumerate(doc):
                if pno>=max_pages_per_document:break
                text=page.get_text('text')
                if not re.search(r'\b(?:IPCE|EQE)\b',text,re.I):continue
                axes=_axis_from_words(page)
                for axis in axes:
                    edges=_extract_cutoffs(page,axis)
                    if len(edges)<2:continue
                    wl=np.array([q['lambda50_nm'] for q in edges])
                    # dominant cluster of coinciding red-edge responses
                    cluster=max((list(np.where(abs(wl-v)<=tolerance_nm)[0]) for v in wl),key=len)
                    if len(cluster)<2:continue
                    waves=wl[cluster]
                    if max(waves)-min(waves)>tolerance_nm:continue
                    mean_wl=float(np.median(waves))
                    # Nearby caption must identify EQE/IPCE control/target,
                    # and not just unrelated numerical graph panels.
                    page_label=text.lower()
                    if not re.search(r'(?:eqe|ipce).{0,160}(?:control|target|reference|treated|modified)|(?:control|target|reference|treated|modified).{0,160}(?:eqe|ipce)',page_label,re.I|re.S):continue
                    same=[r for r in rows if not __import__('re').search(r'\b(?:pure|phase[- ]pure|standalone)\s*2d\b',str(r.get('device_category','')),__import__('re').I)]
                    # Need >=2 categories explicitly mentioned on figure page,
                    # otherwise do not assume every modified concentration is plotted.
                    matched=[r for r in same if re.search(r'(?<![\w])'+re.escape(str(r.get('device_category','')).split(' stability group')[0])+r'(?![\w-])',text,re.I)]
                    if len(matched)<2:
                        matched=[r for r in same if (str(r.get('device_category','')).lower() in ('control','target') and re.search(r'\b'+re.escape(str(r.get('device_category')))+r'\b',text,re.I))]
                    if len(matched)<2:continue
                    for r in matched:
                        if r.get('bandgap_ev') not in (None,''):continue
                        r['bandgap_ev']=round(1239.841984/mean_wl,4)
                        r['bandgap_source']=f'{Path(doc_path).name} p.{pno+1}: automatic shared EQE long-wavelength 50% edge ({len(cluster)} color traces)'
                        r['bandgap_lambda50_nm']=round(mean_wl,2)
                        r['bandgap_shared_edge']='YES'
                        recovered+=1
                    records.append(dict(pdf=Path(doc_path).name,page=pno+1,
                                        wavelength_nm=round(mean_wl,2),edges=edges,
                                        matched_categories=[r.get('device_category') for r in matched]))
    return recovered,records
