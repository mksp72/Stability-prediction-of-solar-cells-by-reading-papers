"""Detect raster traces when graph axes are selectable PDF text.

Graph panels are found by numerical tick geometry, not page location.  This
reader is for *continuous measured signals* only; marker/errorbar plots are
handled independently by marker_stability_recovery.
"""
from pathlib import Path
from collections import defaultdict
from itertools import combinations
import re, fitz, numpy as np, cv2
from general_curve_reader import numeric_words, _mask

NUMERIC_LABEL=re.compile(r'(?i)normalized|retention|PCE|SPO|power|efficiency')

def _split_line(groups,tol=3):
 rows=[]
 for q in sorted(groups,key=lambda a:(a[1],a[0])):
  near=next((v for v in rows if abs(v[0]-q[1])<tol),None)
  if near:near[1].append(q)
  else:rows.append([q[1],[q]])
 return rows

def _fit(nums,direction):
 if len(nums)<3:return None
 best=None
 for a,b in combinations(nums,2):
  i=0 if direction=='x' else 1
  dx=b[i]-a[i];dy=b[2]-a[2]
  if abs(dx)<12 or abs(dy)<.15:continue
  slope=dy/dx;off=a[2]-slope*a[i]
  if direction=='x' and slope<=0 or direction=='y' and slope>=0:continue
  match=[t for t in nums if abs(slope*t[i]+off-t[2])<max(.12,.02*max(t[2],1))]
  if len(match)<3:continue
  err=np.median([abs(slope*t[i]+off-t[2]) for t in match])
  score=(len(match),np.ptp([q[i] for q in match]),-err)
  if best is None or score>best[0]:best=(score,slope,off,match)
 return best[1:] if best else None

def _dense_continuous(rgb,plot,color,scale=3):
 """Only sample densely rendered *measured* normalized MPP trajectories.

 This does not digitize straight connector lines between discrete markers.
 Candidate traces require >70% horizontal color coverage; error-bar plots
 belong to the marker-only reader.
 """
 hsv=cv2.cvtColor(rgb,cv2.COLOR_RGB2HSV)
 m=_mask(hsv,color)
 x0=max(0,int((plot['x0']+3)*scale));x1=min(m.shape[1],int((plot['x1']-3)*scale))
 y0=max(0,int((plot['y0']-6)*scale));y1=min(m.shape[0],int((plot['y1']-5)*scale))
 if x1-x0<80 or y1-y0<30:return []
 pts=[];present=0;total=0
 for xx in range(x0,x1,7):
  total+=1
  ys=np.where(m[y0:y1,max(0,xx-2):xx+3].any(axis=1))[0]
  if not len(ys):continue
  # The selected normalized signal is the locally dominant horizontal
  # colored band, not text or a single error-bar endpoint.
  counts=np.bincount(ys//max(2,scale),minlength=max(1,(y1-y0)//scale+1))
  peak=int(np.argmax(counts))
  yy=(peak+.5)*max(2,scale)+y0
  if counts[peak]<1:continue
  xxpdf=xx/scale;yypdf=yy/scale
  t=plot['x_slope']*xxpdf+plot['x_intercept'];v=plot['y_slope']*yypdf+plot['y_intercept']
  if -0.02<=v<=1.25 and -0.01*plot['x_max']<=t<=1.02*plot['x_max']:
   pts.append((float(t),float(v)));present+=1
 if present/max(total,1)<.70:return []
 # Series near the normalized top should not be rejected as out of axis bounds.
 arr=np.array(pts)
 if len(arr)<20:return []
 # Sparse text contamination causes large spikes. Remove after local median.
 filtered=[]
 for i,(t,v) in enumerate(arr):
  ref=np.median(arr[max(0,i-4):i+5,1])
  if abs(v-ref)<.09:filtered.append((round(t,3),round(v,5)))
 return filtered if len(filtered)>=20 else []


def detect(pdf):
 found=[]
 with fitz.open(pdf) as doc:
  for page in doc:
   txt=page.get_text('text')
   if not re.search('(?i)MPP|tracking|aging|stability|retention',txt):continue
   words=page.get_text('words')
   nums=numeric_words(words)
   xlines=[]
   for yy,row in _split_line([n for n in nums if 0<=n[2]<=100000]):
    if len(row)<3:continue
    for segment in [sorted(row,key=lambda a:a[0])]:
     # Multiple x axes in one horizontal row are split at tick-value resets.
     parts=[];part=[]
     for q in segment:
      if part and q[2]<part[-1][2]:parts.append(part);part=[]
      part.append(q)
     parts.append(part)
     for part in parts:
      if len(part)<3 or max(q[2] for q in part)<75:continue
      fit=_fit(part,'x')
      if fit is None:continue
      sx,bx,inl=fit
      xmax=max(q[2] for q in inl)
      xmin=-bx/sx
      if xmin<5 or xmin>page.rect.width-60:continue
      right=(xmax-bx)/sx
      if right-xmin<85 or right>page.rect.width+15:continue
      bottom=np.median([q[1] for q in inl])-5
      # Separate horizontally neighboring charts and reject non-hour plots.
      near=' '.join(w[4] for w in words if xmin-20<.5*(w[0]+w[2])<right+20 and bottom<=w[1]<bottom+35)
      if re.search(r'(?i)\bns\b|\bms\b|\bmin\b|seconds',near):continue
      # Some vertical multi-panels share one Time (hours) label printed only
      # underneath their bottom panel. Require hour units elsewhere in the
      # SAME source figure and an explicitly described hour-based aging test.
      if not re.search(r'(?i)\btime\b|\bhours?\b|\(h\)',near):
       shared=' '.join(w[4] for w in words if xmin-30<w[0]<right+20 and bottom<w[1]<bottom+140)
       if not re.search(r'(?i)time.*hours?',shared):continue
      if not re.search(r'(?i)\bhours?\b|\(h\)',near):
       shared=' '.join(w[4] for w in words if xmin-30<w[0]<right+20 and bottom<w[1]<bottom+140)
       if not re.search(r'(?i)hours?',shared):continue
      # For each x baseline, only left-axis printed y ticks above it.
      leftnums=[q for q in nums if xmin-25<q[0]<xmin-3 and bottom-180<q[1]<bottom-9 and 0<=q[2]<=120]
      ycolumns=[]
      for xpos,grp in _split_line([(q[1],q[0],q[2],q[3]) for q in leftnums],tol=4):
       # Here swapped x/y to cluster by source x coordinate.
       ycolumns.append(grp)
      candidates=[]
      for group in ycolumns:
       zz=sorted([(q[1],q[0],q[2],q[3]) for q in group],key=lambda q:q[1])
       # Two vertically stacked plots can share identical left tick columns.
       # Split whenever the printed Y values restart at a larger number.
       batches=[];batch=[]
       for q in zz:
        if batch and q[2]>batch[-1][2]+0.25:
         batches.append(batch);batch=[]
        batch.append(q)
       batches.append(batch)
       for zz in batches:
        fit_y=_fit(zz,'y')
        if not fit_y:continue
        sy,by,matched=fit_y
        if len(matched)<3:continue
        # Tick group must include last axis tick near baseline.
        if bottom-max(v[1] for v in matched)>40:continue
        valbot=sy*bottom+by
        ytop=min(t[1] for t in matched)-2
        ymax=sy*ytop+by
        if ymax>120 or ymax<.6 or valbot<-.8*ymax:continue
        if bottom-ytop<26:continue
        candidates.append((len(matched),-abs(valbot),sy,by,ytop,ymax))
      if not candidates:continue
      n,_,sy,by,top,ymax=max(candidates)
      label=' '.join(w[4] for w in words if xmin-70<w[0]<xmin-5 and top-5<w[1]<bottom)
      if not NUMERIC_LABEL.search(label):
       continue
      normalized=bool(re.search(r'(?i)normalized',label)) or ymax<=1.6
      if not normalized:continue  # prevent connector lines on discrete-PCE/error-bar charts
      plot=dict(page=page.number+1,x0=xmin,x1=right,y0=top,y1=bottom,
        x_slope=sx,x_intercept=bx,y_slope=sy,y_intercept=by,x_max=xmax,y_max=ymax,
        y_axis_label=label,source='PDF_NATIVE_TEXT_TICKS')
      pix=page.get_pixmap(matrix=fitz.Matrix(3,3),alpha=False)
      rgb=np.frombuffer(pix.samples,np.uint8).reshape(pix.height,pix.width,pix.n)[:,:,:3].copy()
      from general_curve_reader import COLORS
      for color in ['red','orange','blue','cyan','magenta','violet','green']:
       track=_dense_continuous(rgb,plot,color,scale=3)
       if len(track)<12:continue
       if track[-1][0]-track[0][0]<.65*xmax:continue
       arr=np.array(track,float)
       if not np.isfinite(arr).all() or np.any(arr[:,1]<0):continue
       # Thin connecting lines are not measurements on marker-style plots.
       # This reader accepts only curves with no distinct marker glyphs.
       found.append(dict(source_file=Path(pdf).name,source_page=page.number+1,
         axis='NORMALIZED_0_TO_1',color=color,number_of_points=len(track),
         t_end_h=round(float(arr[-1,0]),2),retention_final_pct=round(float(100*arr[-1,1]),3),
         points=track,plot_bounds=[round(float(plot[k]),2) for k in ('x0','y0','x1','y1')],
         status='CONTINUOUS_RASTER_TRACE_CANDIDATE_CONFIRM_DEVICE_IDENTITY'))
 return found
