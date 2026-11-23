"""Conservatively promote explicit single-device normalized continuous MPP panels.

This is not for marker/error-bar plots. Only source captions naming a unique
single formulation qualify. Never infer control/optimized pairing from colors.
"""
from pathlib import Path
import re, json
import fitz
import numpy as np
from vector_time_plot import detect

NAMED=re.compile(r'(?i)performance\s+of\s+(.{3,85}?)\s+(?:PSCs?|solar\s+cells?)\b')

def source_label(text):
 # Prefer the caption of this FIGURE, not a label mentioned elsewhere in
 # related-work paragraphs or a paper title.
 for m in re.finditer(r'(?i)Fig(?:ure)?\.?\s*\d+\s*[|:.]',text):
  fragment=text[m.end():m.end()+180]
  z=NAMED.search(fragment)
  if not z:continue
  label=' '.join(z.group(1).split()).strip(' .:;')
  if not label or len(label)>65 or re.search(r'(?i)control\s+and|versus|vs\.|different|various|multiple|several',label):continue
  return label
 return None

def obtain(main,supplements=()):
 candidates=[]; output=[]
 for p in [Path(main),*[Path(s) for s in supplements]]:
  if p.suffix.lower()!='.pdf':continue
  candidates.extend(detect(p))
 # De-duplicate overlapping HSV masks, e.g. blue and cyan detecting identical pixels.
 sorted_candidates=sorted(candidates,key=lambda a:a['number_of_points'],reverse=True)
 unique=[]
 for c in sorted_candidates:
  if any(c['source_file']==old['source_file'] and c['source_page']==old['source_page'] and
     abs(c['t_end_h']-old['t_end_h'])<6 and abs(c['retention_final_pct']-old['retention_final_pct'])<1.5
     and np.median(np.abs(np.array(c['points'])[:min(len(c['points']),len(old['points'])),1]-
     np.array(old['points'])[:min(len(c['points']),len(old['points'])),1]))<.03
     for old in unique):continue
  unique.append(c)
 for c in unique:
  pdf=next((p for p in [Path(main),*[Path(s) for s in supplements]] if p.name==c['source_file']),None)
  if pdf is None:continue
  with fitz.open(pdf) as doc:
   page=doc[c['source_page']-1]
   cap=source_label(page.get_text('text'))
   # Single-formulation caption is necessary, not sufficient: local panel must
   # explicitly describe continuous operating PCE/MPP (not a simulated line).
   rect=c['plot_bounds']; box=fitz.Rect(rect[0],rect[1],rect[2],rect[3])
   local=page.get_text('text',clip=box)
   local=' '.join(local.split())
   if not re.search(r'(?i)MPP|tracking|sun\s+illumination',local):
    c['status']='REVIEW_NO_LOCAL_OPERATING_EVIDENCE';continue
   if not cap:
    c['status']='REVIEW_CATEGORY_NOT_UNAMBIGUOUS_IN_CAPTION';continue
   c['device_label']=cap;c['local_panel_text']=local[:350]
   c['status']='SOURCE_CALIBRATED_CONTINUOUS_MPP_SINGLE_CAPTION_DEVICE'
   output.append(c)
 return unique,output

def to_rows(candidates,doi,title,main_file):
 """Create ML records only for caption-verified normalized MPP series."""
 if not candidates:return []
 from perovskite_onefile_final import fit_every_trace,choose_decay_fit
 from pipeline import observed_t80
 from collections import defaultdict
 grouped=defaultdict(list)
 for c in candidates:grouped[c['device_label']].append(c)
 results=[]
 for label,series in grouped.items():
  studies=[]
  for c in series:
   raw=np.asarray(c['points'],float)
   xx=raw[:,0]; yy=raw[:,1]*100
   if len(xx)<25 or np.ptp(xx)<100 or np.any(yy<0) or np.any(yy>125):continue
   order=np.argsort(xx); xx=xx[order];yy=yy[order]
   # Avoid duplicate timestamps introduced by image sampling.
   keep=np.r_[True,np.diff(xx)>0];xx=xx[keep];yy=yy[keep]
   if len(xx)<20:continue
   isos='ISOS-L3' if re.search(r'(?i)65\s*°?\s*C|65C',c.get('local_panel_text','')) else 'ISOS-L_OR_MPP'
   initial_match=re.search(r'(?i)initial\s+PCE\s*[:=]\s*(\d+(?:\.\d+)?)\s*%',c.get('local_panel_text',''))
   aging_initial_pce=float(initial_match.group(1)) if initial_match else None
   study=dict(aging_initial_pce_pct=aging_initial_pce,pdf=c['source_file'],page=c['source_page'],figure='Figure - normalized MPP',
     protocol=isos,measurement_mode='CONTINUOUS_MPP',
     conditions=c.get('local_panel_text','')[:250],metric='NORMALIZED_PCE',
     n_points=len(xx),t_start_h=float(xx[0]),t_end_h=float(xx[-1]),
     time_h=[round(float(x),3) for x in xx],
     retention_pct=[round(float(x),4) for x in yy],
     raw_plot_y=[round(float(x),5) for x in yy/100],
     normalization_source='PUBLISHED_NORMALIZED_0_TO_1',
     retention_last_pct=round(float(yy[-1]),3),
     observed_t80_h=observed_t80(xx,yy),
     digitization_quality_flags=['RASTER_AUTOMATIC_REVIEW_RECOMMENDED',
       'CONTINUOUS_MPP_NO_ERROR_BAR_MARKERS','TEXT_IDENTIFIED_CATEGORY'])
   study['fits']=fit_every_trace(xx,yy)
   study['fit_summary']=choose_decay_fit(study)
   studies.append(study)
  if not studies:continue
  primary=max(studies,key=lambda s:s['t_end_h'])
  result=dict(doi=doi,article_title=title,pdf_file=str(main_file),device_category=label,
   category_id='continuous_mpp_'+re.sub('[^a-z0-9]+','_',label.lower()).strip('_'),
   extraction_method='NATIVE_AXIS_CONTINUOUS_MPP_RASTER',
   stability_experiment_count=len(studies),primary_stability_figure=primary['figure'],
   primary_stability_protocol=primary['protocol'],
   primary_stability_conditions=primary['conditions'],
   final_retention_pct=primary['retention_last_pct'],
   final_aging_time_h=primary['t_end_h'],
   stability_endpoint=str((primary['retention_last_pct'],primary['t_end_h'])),
   stability_data_json=json.dumps(studies,ensure_ascii=False),
   observed_t80_h=primary['observed_t80_h'],
   t80_status='OBSERVED' if primary['observed_t80_h'] is not None else 'NOT_REACHED_DURING_MEASUREMENT',
   retention_reference='PUBLISHED_NORMALIZED_0_TO_1',
   primary_stability_metric='PCE_NORMALIZED',
   pce0=None,voc0=None,jsc0=None,ff0=None,
   initial_aging_pce_pct=primary.get('aging_initial_pce_pct'),
   category_legend_source='FIGURE_CAPTION_SINGLE_DEVICE',
   digitization_quality_status='RASTER_AUTO_NEEDS_VISUAL_QC',
   ml_training_ready=False,
   digitization_quality_flags_json=json.dumps(primary['digitization_quality_flags']),
   fit_traces_successful=sum('fit_summary' in s for s in studies))
  result.update(primary.get('fit_summary',{}))
  results.append(result)
 return results
