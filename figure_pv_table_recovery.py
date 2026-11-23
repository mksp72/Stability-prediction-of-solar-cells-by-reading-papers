"""v18.6.7 generic same-row J-V tables printed inside multipanel figures.

Does not require literal Table S1 text. Uses real PDF word coordinates,
recognizes four figure header columns and requires physical same-row values,
plausible units, and PCE closure. Independent of DOI/material/contact names.
"""
from pathlib import Path
import csv,math,re
import fitz

METRICS=('voc0','jsc0','ff0','pce0')

def _metric(s):
 s=re.sub(r'[^a-z0-9η]+','',str(s).lower())
 if s in ('voc','opencircuitvoltage'):return 'voc0'
 if s in ('jsc','shortcircuitcurrentdensity'):return 'jsc0'
 if s in ('ff','fillfactor'):return 'ff0'
 if s in ('pce','η','eta','efficiency'):return 'pce0'
 return None

def _num(s):
 x=str(s).replace('−','-').replace(',','.').replace('%','').strip()
 if not re.fullmatch(r'\d{1,3}(?:\.\d{1,5})?',x):return None
 v=float(x)
 return v if math.isfinite(v) else None

def _key(s):
 from hysteresis_recovery import canonical_key
 return re.sub(r'(?:treated|treatment)$','',canonical_key(s))


def extract(path):
 found=[];path=Path(path)
 if path.suffix.lower()!='.pdf' or not path.is_file():return found
 try:
  with fitz.open(path) as doc:
   for pi,page in enumerate(doc):
    words=page.get_text('words')
    if len(words)<20:continue
    lines=[]
    for w in sorted(words,key=lambda w:w[1]):
     cy=(w[1]+w[3])/2
     grp=next((g for g in lines if abs(g[0]-cy)<2.8),None)
     if grp is None:lines.append([cy,[w]])
     else:grp[1].append(w)
    lines.sort(key=lambda g:g[0]);lines=[(yy,sorted(ws,key=lambda w:w[0])) for yy,ws in lines]
    for hi,(hy,hw) in enumerate(lines):
     anchors={}
     for w in hw:
      field=_metric(w[4])
      if field and field not in anchors:anchors[field]=(w[0]+w[2])/2
     if set(anchors)!=set(METRICS):continue
     xvals=list(anchors.values());cols=sorted(anchors,key=anchors.get)
     if max(xvals)-min(xvals)<35:continue
     distances=[abs(x-y) for x,y in zip(sorted(xvals)[:-1],sorted(xvals)[1:])]
     if min(distances)<6:continue
     for yy,ws in lines[hi+1:hi+18]:
      if yy-hy>100:break
      nums={};labels=[]
      for w in ws:
       xx=(w[0]+w[2])/2
       value=_num(w[4])
       nearest=min(anchors,key=lambda f:abs(anchors[f]-xx))
       if value is not None and abs(xx-anchors[nearest])<max(8,min(distances)*.55):
        if nearest in nums:nums={};break
        nums[nearest]=value
       elif min(xvals)-34<xx<min(xvals)-6 and not str(w[4]).strip().isdigit():labels.append(w[4])
      if set(nums)!=set(METRICS):continue
      label=' '.join(labels).strip()
      if not label or len(label)>90:continue
      v,j,ff,pce=(nums[k] for k in METRICS)
      if ff<=1:ff*=100
      if not (.3<=v<=2.8 and 1<j<45 and 25<=ff<=95 and 2<pce<=40):continue
      closure=v*j*ff/100
      if abs(closure-pce)>max(.75,.055*pce):continue
      found.append({'file':path.name,'page':pi+1,'device_label':label,
       'voc0':v,'jsc0':j,'ff0':ff,'pce0':pce,'jv_closure_error':round(closure-pce,4),
       'status':'SOURCE_ROW_CONFIRMED_IDENTITY_NOT_YET_MATCHED'})
 except Exception:return found
 # Eliminate duplicate recognition from the same page, retaining the exact set.
 unique={}
 for r in found:unique[(r['page'],_key(r['device_label']),r['pce0'],r['voc0'],r['jsc0'],r['ff0'])]=r
 return list(unique.values())


def recover(rows,main_pdf,supplements=()):
 candidates=[c for f in [main_pdf,*supplements] for c in extract(f)]
 changed=0
 for row in rows:
  name=_key(row.get('device_category',''))
  matches=[x for x in candidates if _key(x['device_label'])==name]
  if not matches:continue
  # Different champion/scan/area sets are NOT interchangeable.
  sets={tuple(m[f] for f in METRICS) for m in matches}
  if len(sets)!=1:
   row['figure_table_jv_status']='MULTIPLE_SETS_REQUIRE_SCAN_AREA_REVIEW';continue
  m=matches[0]
  # A complete four-metric printed figure-table row, with exact formulation
  # identity and PCE closure, is internally consistent evidence. Do not build
  # a hybrid row by keeping unrelated values from broad narrative recovery.
  conflicts=[]
  for field in METRICS:
   before=row.get(field)
   try:
    if before not in (None,'') and abs(float(before)-float(m[field]))>1e-6:
     conflicts.append(f'{field}: previous={before} figure_table={m[field]}')
   except (ValueError,TypeError):pass
   if before!=m[field]:changed+=1
   row[field]=m[field]
   row[field+'_source']=f'Complete same-row native figure J-V table, {m["file"]} p.{m["page"]}; {m["device_label"]}'
  row['pce_jv_pct']=m['pce0']
  if conflicts:row['pv_source_conflicts']=' | '.join(conflicts)
  row['figure_table_jv_status']='SOURCE_ROW_EXACT_LABEL_MATCH'
 return changed,candidates


def save_candidates(candidates,path):
 if not candidates:return
 cols=['file','page','device_label','voc0','jsc0','ff0','pce0','jv_closure_error','status']
 with Path(path).open('w',encoding='utf-8-sig',newline='') as f:
  w=csv.DictWriter(f,cols);w.writeheader();w.writerows(candidates)
