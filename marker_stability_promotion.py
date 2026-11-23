"""Qualified marker-center extraction for stability figures with printed error bars.

No connecting lines or error-bar heights are sampled as measured values.
Qualified series carry explicit REVIEW status until their source overlay passes QC.
"""
from pathlib import Path
import json,re,difflib,statistics
from collections import defaultdict
import numpy as np, fitz
from marker_stability_recovery import detect_pdf

PANEL=re.compile(r'(?i)^\(?([a-f])\)?$')
MATERIAL=re.compile(r'(?i)(?:control\s+)?([A-Za-z][A-Za-z0-9₀-₉₆_./()\-]{2,40})\s+devices?')

def _key(s):return re.sub(r'[^a-z0-9]','',str(s).lower()).replace('6','6')

def _name_choices(text):
 # Printed stability figure captions frequently introduce the tested
 # formulations as "control XYZ devices ... and ABC devices".
 # Materials are taken from THIS figure caption, never a DOI dictionary.
 names=[]
 for m in MATERIAL.finditer(text):
  name=m.group(1)
  if name.lower() in ('perovskite','solar','encapsulated','unencapsulated','optimized','treated','different','best','the','these','those'):continue
  if name not in names:names.append(name)
 return names

def _figure_caption(text):
 match=re.search(r'(?i)Fig(?:ure)?\.?\s*\d+[.:]',text)
 if match:return text[match.start():match.start()+2100]
 return text

def _panel_letter(ws,record):
 a=json.loads(record['axis_json']);x=a['x0']; y=a['y0']
 options=[]
 for w in ws:
  t=w[4].strip(' :.')
  m=PANEL.fullmatch(t)
  if not m:continue
  if x-58<w[0]<x+10 and y-35<w[1]<y+18:
   options.append((abs(w[0]-x+18)+abs(w[1]-y+14),m.group(1)))
 return min(options)[1] if options else None

def _clean_points(record):
 source=json.loads(record['points_json'])
 pts=np.array([[p['time_h'],p['measured_y']] for p in source],float)
 if len(pts)<8:return None,0
 pts=pts[np.argsort(pts[:,0])]
 # Duplicate template matches at one timestamp: retain one center. Error bars
 # are not drawn into an interpolated polyline.
 keep=np.r_[True,np.diff(pts[:,0])>1e-7];pts=pts[keep]
 if len(pts)<8:return None,0
 valid=np.ones(len(pts),dtype=bool)
 # Reject isolated false positives only when both neighboring centers jointly
 # support a very different value. Never enforce monotonic degradation.
 for i in range(1,len(pts)-1):
  prev,nxt=pts[i-1],pts[i+1]
  if nxt[0]<=prev[0] or nxt[0]-prev[0]>.24*(pts[-1,0]-pts[0,0]):continue
  expected=prev[1]+(nxt[1]-prev[1])*(pts[i,0]-prev[0])/(nxt[0]-prev[0])
  if abs(pts[i,1]-expected)>max(1.15,.065*np.median(pts[:,1])):
   valid[i]=False
 clean=pts[valid]
 return clean,int((~valid).sum())

def obtain(pdf,raw_records=None):
 records=raw_records if raw_records is not None else detect_pdf(pdf)
 if not records:return records,[]
 verified=[];p=Path(pdf)
 with fitz.open(p) as doc:
  for page_index in set(r['source_page']-1 for r in records):
   page=doc[page_index]
   caption=_figure_caption(page.get_text())
   names=_name_choices(caption)
   if len(names)!=2:continue
   plotrecords=[r for r in records if r['source_page']-1==page_index]
   # For this page each figure might have panels (a), (b), (c): do not
   # allow a weak violet/blue duplicate to become a third device.
   candidates=[]
   for r in plotrecords:
    if r['raw_measurement_type']!='ABSOLUTE_PCE_OR_POWER':continue
    if r['color'] not in ('red','blue','orange','green'):continue
    pts,dropped=_clean_points(r)
    if pts is None or len(pts)<8 or np.ptp(pts[:,0])<150:continue
    label=r['legend_candidate']
    if label:
     best=max(names,key=lambda q:difflib.SequenceMatcher(None,_key(label),_key(q)).ratio())
     score=difflib.SequenceMatcher(None,_key(label),_key(best)).ratio()
     if score<.57:label=None
     else:label=best
    candidates.append((r,pts,label,dropped))
   # Infer a stable color-to-device mapping only from printed, source-matched
   # legends on THIS figure page. Many scientific papers repeat the same red
   # and blue traces across four subpanels but print a legend only once.
   grouped=defaultdict(list)
   for r,pts,label,dropped in candidates:grouped[r['plot_number']].append((r,pts,label,dropped))
   color_evidence=defaultdict(set)
   for group in grouped.values():
    for r,pts,label,dropped in group:
     if label:color_evidence[r['color']].add(label)
   color_map={col:next(iter(nms)) for col,nms in color_evidence.items() if len(nms)==1}
   # If exactly two source-caption device identities recur and one color is
   # anchored unambiguously, the other repeated color is the complement.
   paired=[g for g in grouped.values() if len(g)==2 and set(x[0]['color'] for x in g)=={'red','blue'}]
   if paired and len(names)==2:
    if 'red' in color_map and 'blue' not in color_map:
     color_map['blue']=next(n for n in names if n!=color_map['red'])
    elif 'blue' in color_map and 'red' not in color_map:
     color_map['red']=next(n for n in names if n!=color_map['blue'])
   if len(set(color_map.values()))!=len(color_map):color_map={}
   for plot_id,group in grouped.items():
    group=[(r,pts,label or color_map.get(r['color']),dropped) for r,pts,label,dropped in group]
    if len(group)!=2 or len(set(x[0]['color'] for x in group))!=2:continue
    resolved=[x[2] for x in group if x[2]]
    if len(resolved)==1:
     complement=next((n for n in names if n!=resolved[0]),None)
     group=[(r,pts,label or complement,dropped) for r,pts,label,dropped in group]
    if not all(x[2] for x in group) or len(set(x[2] for x in group))!=2:continue
    for r,pts,label,dropped in group:
     r['matched_device']=label
     r['panel_letter']=r.get('subplot_letter') or None
     r['rejected_outlier_marker_positions']=dropped
     r['selected_time_h']=[round(float(x),3) for x in pts[:,0]]
     r['selected_absolute_pce_pct']=[round(float(y),4) for y in pts[:,1]]
     r['status']='MARKER_CENTERS_SOURCE_LABEL_MATCHED_VISUAL_QC_REQUIRED'
     verified.append(r)
 return records,verified

def to_rows(records,doi,title,main_file):
 from pipeline import observed_t80
 studies=defaultdict(list)
 for r in records:
  label=r['matched_device']
  t=np.array(r['selected_time_h']);yy=np.array(r['selected_absolute_pce_pct'])
  if len(t)<8 or yy[0]<=0:continue
  rt=100*yy/yy[0]
  protocol={ 'a':'LIGHT_STRESS_NON_MPP', 'b':'DARK_THERMAL_STORAGE',
    'c':'DAMP_HEAT_STORAGE','d':'OUTDOOR_ISOS_O_1',
    'e':'OPERATIONAL_MPP'}.get(r.get('panel_letter'),'UNRESOLVED_STABILITY_MODE')
  trace={'pdf':r['source_file'],'page':r['source_page'],
   'figure':f'Figure panel {r.get("panel_letter") or r["plot_number"]}',
   'protocol':protocol,'measurement_mode':'INTERMITTENT_JV' if protocol!='OPERATIONAL_MPP' else 'MPP',
   'conditions':'SOURCE FIGURE CAPTION - confirm aging protocol for this panel',
   'metric':'ABSOLUTE_PCE','n_points':len(t),'t_start_h':float(t[0]),'t_end_h':float(t[-1]),
   'time_h':t.tolist(),'raw_plot_y':yy.tolist(),'retention_pct':rt.round(4).tolist(),
   'retention_last_pct':float(rt[-1]),
   'normalization_source':'RELATIVE_TO_FIRST_DIGITIZED_ABSOLUTE_PCE',
   'observed_t80_h':observed_t80(t,rt),
   'digitization_quality_flags':['MARKER_ONLY','AUTO_VISUAL_QC_REQUIRED',
         f'REJECTED_{r["rejected_outlier_marker_positions"]}_POTENTIAL_ARTIFACTS']}
  # No numerical fits on provisional glyph matches: false glyph centers
  # can otherwise cause both misleading T80 and excessive optimizer time.
  trace['fits']=[]
  trace['fit_status']='WAITING_FOR_MARKER_OVERLAY_REVIEW'
  studies[label].append(trace)
 rows=[]
 for label,traces in studies.items():
  primary=max(traces,key=lambda x:x['t_end_h']);raw_json=json.dumps(traces,ensure_ascii=False)
  r=dict(doi=doi,article_title=title,pdf_file=main_file,device_category=label,
    category_id='marker_'+_key(label),extraction_method='RASTER_MARKER_CENTERS_PRELIMINARY',
    stability_experiment_count=len(traces),primary_stability_figure=primary['figure'],
    primary_stability_protocol=primary['protocol'],primary_stability_conditions=primary['conditions'],
    final_retention_pct=primary['retention_last_pct'],final_aging_time_h=primary['t_end_h'],
    stability_endpoint=str((primary['retention_last_pct'],primary['t_end_h'])),
    stability_data_json=raw_json,observed_t80_h=primary['observed_t80_h'],
    t80_status='OBSERVED' if primary['observed_t80_h'] is not None else 'NOT_REACHED_DURING_MEASUREMENT',
    retention_reference=primary['normalization_source'],stability_measurement='ABSOLUTE_PCE',
    fit_traces_successful=0,
    digitization_quality_status='PRELIMINARY_EXCLUDE_FROM_ML_UNTIL_VISUAL_REVIEW',
    ml_training_ready=False,
    digitization_quality_flags_json=json.dumps(primary['digitization_quality_flags']),
    category_legend_source='SOURCE_FIGURE_CAPTION_AND_COLOR_LEGEND',
    pce0=None,voc0=None,jsc0=None,ff0=None)
  r.update(primary.get('fit_summary',{}));rows.append(r)
 return rows
