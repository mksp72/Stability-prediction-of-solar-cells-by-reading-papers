"""Transparent source-grounded long-format measurement/coverage exports.

Does not promote preliminary marker data to validated ML training data.
"""
from __future__ import annotations
import csv,json,re
from collections import defaultdict
from pathlib import Path

POINT_FIELDS=['source_pdf','device_category','category_quality','protocol','figure','page','trace_index','point_index','time_h','raw_y','raw_y_type','retention_pct','source_digitization','usable_for_ml']
PV_FIELDS=['source_pdf','device_category','pce_pct','voc_v','jsc_ma_cm2','ff_pct','source_digitization','quality_status']
COVER_FIELDS=['source_pdf','detected_stability_figure_pages','text_endpoint_mentions','raw_marker_candidates','calibrated_continuous_candidates','exported_device_categories','exported_stability_traces','exported_points','ml_ready_categories','coverage_status','unresolved_reason']

def _write(path,fields,items):
 path=Path(path);path.parent.mkdir(parents=True,exist_ok=True)
 with path.open('w',encoding='utf-8-sig',newline='') as f:
  w=csv.DictWriter(f,fieldnames=fields,extrasaction='ignore');w.writeheader();w.writerows(items)

def _safe_json(s):
 if isinstance(s,list):return s
 if isinstance(s,dict):return [s]
 try:
  if isinstance(s,str) and s.startswith('zlib64:'):
   import zlib,base64
   s=zlib.decompress(base64.b64decode(s[7:])).decode('utf-8')
  data=json.loads(s or '[]')
  return data if isinstance(data,list) else [data]
 except (ValueError,TypeError):return []

def produce(dest,all_rows,all_figures,all_reported_endpoints,markers,continuous,source_files):
 dest=Path(dest);output=dest.parent
 point_rows=[];pv_rows=[];coverage=[]
 perfile=defaultdict(list)
 for r in all_rows:perfile[Path(str(r.get('pdf_file',''))).name].append(r)
 for r in all_rows:
  file=Path(str(r.get('pdf_file',''))).name;device=str(r.get('device_category',''))
  quality=str(r.get('digitization_quality_status',r.get('extraction_method','')))
  pv_rows.append(dict(source_pdf=file,device_category=device,pce_pct=r.get('pce0'),voc_v=r.get('voc0'),
   jsc_ma_cm2=r.get('jsc0'),ff_pct=r.get('ff0'),source_digitization=r.get('extraction_method',''),quality_status=quality))
  for i,trace in enumerate(_safe_json(r.get('stability_data_json'))):
   if not isinstance(trace,dict):continue
   time=trace.get('time_h',[]);ret=trace.get('retention_pct',[]);raw=trace.get('raw_plot_y',[])
   if not isinstance(time,list):continue
   for j,t in enumerate(time):
    point_rows.append(dict(source_pdf=file,device_category=device,category_quality=quality,protocol=trace.get('protocol',''),
     figure=trace.get('figure',''),page=trace.get('page',''),trace_index=i+1,point_index=j+1,time_h=t,
     raw_y=raw[j] if isinstance(raw,list) and j<len(raw) else '',raw_y_type=trace.get('metric',''),
     retention_pct=ret[j] if isinstance(ret,list) and j<len(ret) else '',
     source_digitization=r.get('extraction_method',''),usable_for_ml=bool(r.get('ml_training_ready'))))
 all_files=list(dict.fromkeys(str(Path(x).name) for x in source_files))
 for file in all_files:
  rr=perfile[file];figs=[f for f in all_figures if f.get('source_file')==file];ends=[e for e in all_reported_endpoints if e.get('source_file')==file]
  mm=[m for m in markers if m.get('source_file')==file];cc=[c for c in continuous if c.get('source_file')==file]
  tr=sum(len(_safe_json(r.get('stability_data_json'))) for r in rr)
  pts=sum(x['source_pdf']==file for x in point_rows)
  ml=sum(bool(r.get('ml_training_ready')) for r in rr)
  missing=[]
  if figs and not tr:missing.append('STABILITY_FIGURES_PRESENT_BUT_NO_EXPORTED_TRACES')
  if mm and len(mm)>tr:missing.append('SOME_MARKER_SERIES_NOT_DEVICE_ASSIGNED')
  if tr and not ml:missing.append('EXPORTED_TRACES_ARE_PRELIMINARY_NOT_ML_READY')
  if not figs and not tr:missing.append('NO_STABILITY_EVIDENCE_RECOVERED')
  coverage.append(dict(source_pdf=file,detected_stability_figure_pages=len(figs),text_endpoint_mentions=len(ends),
   raw_marker_candidates=len(mm),calibrated_continuous_candidates=len(cc),exported_device_categories=len(rr),
   exported_stability_traces=tr,exported_points=pts,ml_ready_categories=ml,
   coverage_status='PARTIAL_REVIEW_REQUIRED' if missing else 'NO_UNRESOLVED_COVERAGE_FLAGS',unresolved_reason='; '.join(missing)))
 _write(output/'extracted_stability_points.csv',POINT_FIELDS,point_rows)
 _write(output/'extracted_pv_metrics.csv',PV_FIELDS,pv_rows)
 _write(output/'extraction_coverage.csv',COVER_FIELDS,coverage)
 return coverage
