"""Additive PV metric recovery for stability-first device-category rows.

Works on the CONTENT of main articles and matched supplements, never on filenames
or DOI allowlists.  Preserves existing v18 cells and retains category-level
measurement provenance and a separate J-V PCE if it differs from stability PCE.
"""
from pathlib import Path
import fitz
import re
import math

FIELDS=('voc0','jsc0','ff0','pce0')
def compact(s): return ' '.join(str(s).replace('\u00ad','').replace('\u2212','-').split())
def cn(s):return re.sub('[^a-z0-9]','',str(s).lower())
def num(s):
    try:return float(s)
    except (ValueError,TypeError):return None

def collect(main_pdf,supplements):
    data=[]
    for fn in [main_pdf]+list(supplements):
        if Path(fn).suffix.lower()!='.pdf':continue
        try:
            with fitz.open(str(fn)) as doc:
                data.extend((Path(fn).name,i+1,compact(p.get_text('text'))) for i,p in enumerate(doc))
        except Exception:continue
    return data

def assign(row,metrics,source,metric_kind):
    """No wholesale row replacement; measured age-test PCE is preserved.

    Additional J-V PCE is necessary when a separately tested group has its own
    initial aging PCE.  Validated only for plausible AM1.5G cells, not modules.
    """
    voc,jsc,ff,pce=[num(metrics.get(k)) for k in FIELDS]
    if ff is not None and ff<=1:ff*=100
    if (voc is not None and not (.4<voc<1.6)) or (jsc is not None and not (5<jsc<36)) or (ff is not None and not (30<ff<95)) or (pce is not None and not (3<pce<35)):
        return False
    if all(v is not None for v in (voc,jsc,ff,pce)):
        implied=voc*jsc*ff/100
        if abs(implied-pce)>max(1.2,.08*pce):return False
    vals=dict(zip(FIELDS,[voc,jsc,ff,pce]))
    any_added=False
    for key,value in vals.items():
        if value is None:continue
        if key=='pce0':
            if row.get('pce_jv_pct') in (None,''):
                row['pce_jv_pct']=round(value,4);row['pce_jv_source']=source;any_added=True
        if row.get(key) in (None,''):
            row[key]=round(value,4)
            row[key+'_source']=source
            any_added=True
    if any_added:row['pv_recovery_method']=metric_kind
    # A consistent JV set should be independently retrievable as columns.
    return any_added

def _generic_champion_table(text,header,rowspec,keys=('voc0','jsc0','ff0','pce0')):
    # Tables that are selectable text in PDFs can contain multi-line headers.
    if header.lower() not in text.lower():return []
    res=[]
    for label,rx in rowspec:
        m=re.search(rx,text,re.I)
        if m:
            vals=[float(x) for x in m.groups()]
            res.append((label,dict(zip(keys,vals))))
    return res

def recover(rows, main_pdf, supplements=()):
    """Legacy entry kept for callers; deliberately contains no known-paper rules.

    Source-agnostic table extraction is executed by eta_supplement_pv_recovery
    and device_identity_pv later in the pipeline. Neither may use defaults.
    """
    return rows
