"""Additive, publication-independent evidence scan for PDFs, including end-of-file figures.

Scans every page of the main article and matched supplements. No DOI whitelist.
Evidence/marker candidates are NEVER silently treated as validated fitted curves.
"""
from pathlib import Path
import re, csv, json
import fitz

STAB = re.compile(r'\bstabilit\w*|\bdegrad\w*|\baging\b|\bageing\b|\bretention\b|\bMPP\b|\bMPPT\b|\bISOS[-\s]?[DLO]\b|\bdamp.heat\b|\bdark.storage\b',re.I)
CAP = re.compile(r'(?:Fig(?:ure)?\.?\s*[S]?\s*\d+[a-z]?\s*[.:|]?[^\n]{0,700})',re.I)
REPORT = re.compile(r'(?i)(?:\b(?:retain(?:s|ed|ing)?|remaining|preserv(?:es|ed)|maintain(?:s|ed)?)[^\n.]{0,90}?)?\b(\d{1,3}(?:\.\d+)?)\s*%[^\n.]{0,130}?\b(?:after|for|over|at|within)\s*(\d+(?:\.\d+)?)\s*(h|hours?|days?)\b')
T80 = re.compile(r'(?i)\bT\s*[_ ]?(80|90|95)\b[^\n.]{0,30}?(\d+(?:\.\d+)?)\s*(h|hours?|days?)\b')

def _window(text, match, n=175):
 return ' '.join(text[max(0,match.start()-n):min(len(text),match.end()+n)].split())

def scan(main, supplements=()):
 records=[]; figures=[]
 for path in [Path(main),*(Path(p) for p in supplements)]:
  if path.suffix.lower()!='.pdf' or not path.exists():continue
  with fitz.open(path) as doc:
   for page in doc:
    t=page.get_text('text',sort=True)
    figcaps=[m.group(0) for m in CAP.finditer(t) if STAB.search(m.group(0))]
    # Scan all pages, not just graph-bearing pages in the body: figures and
    # captions may be at the very end of appendices or supplementary PDFs.
    if figcaps or STAB.search(t) and (page.get_images() or len(page.get_drawings())>30):
     figures.append(dict(source_file=path.name,page=page.number+1,
                         figure_captions=' | '.join(figcaps)[:1100],
                         embedded_images=len(page.get_images()),vector_drawings=len(page.get_drawings()),
                         evidence_state='FOUND_NOT_NECESSARILY_DIGITIZED'))
    for regex,typ in ((REPORT,'PERCENT_AT_DURATION'),(T80,'REPORTED_T_LIFETIME')):
     for m in regex.finditer(t):
      context=_window(t,m)
      if not STAB.search(context):continue
      if typ=='PERCENT_AT_DURATION':
       value=float(m.group(1));duration=float(m.group(2));unit=m.group(3)
       if not 0<value<=130:continue
       metric='REPORTED_PERCENT_NEEDS_REFERENCE_CHECK'
      else:
       value=float(m.group(1));duration=float(m.group(2));unit=m.group(3)
       metric='T'+str(int(value));value=None
      if unit.lower().startswith('day'):duration*=24
      records.append(dict(source_file=path.name,source_page=page.number+1,
                          record_kind=typ,measurement=metric,reported_percent=value,
                          reported_hours=duration,source_context=context,
                          status='PUBLISHED_TEXT_EVIDENCE_ONLY_DEVICE_AND_PROTOCOL_REVIEW'))
 return figures,records

def save(records,path,fields):
 path=Path(path);path.parent.mkdir(parents=True,exist_ok=True)
 with path.open('w',newline='',encoding='utf-8-sig') as f:
  writer=csv.DictWriter(f,fieldnames=fields,extrasaction='ignore');writer.writeheader();writer.writerows(records)
