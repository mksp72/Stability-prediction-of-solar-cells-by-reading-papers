"""Generic device-linked complete photovoltaic sets in publication prose.

A reported PCE/Voc/Jsc/FF set is linked only when the source associates
all four numbers with the same named device, possibly AFTER the values.
Avoid guessing a scan direction; do not confuse etaMPPT with J-V PCE.
"""
from __future__ import annotations
import re
from pathlib import Path
import fitz
from hysteresis_recovery import canonical_key, canonical_scan_label
from general_pv_measurement_recovery import metrics_of, valid_set

PCE_CLAUSE = re.compile(r'\bPCE\s+(?:can\s+be\s+)?(?:further\s+)?(?:improved|increased|enhanced|raised|reached|achieved|obtained|reported)\s+to\s+\d',re.I)
NAME_AFTER = re.compile(r'\b(?:by\s+using|using|for|of|with)\s+(?:the\s+)?(?P<name>[\w\d₀-₉/+().\s-]{3,110}?)\s+(?:structure|architecture|device|devices|cells|solar\s+cells)\b',re.I)
NAME_BEFORE = re.compile(r'(?P<name>[\w\d₀-₉/+().\s-]{3,110}?)\s+(?:structure|architecture|device|devices|cells)\s+(?:exhibit|show|achieve|yield|deliver|reach|has|have)',re.I)

def _chem_tokens(s):
    s=str(s).lower().replace('₆','6').replace('₀','0').replace('₂','2')
    return re.findall(r'[a-z]+\d*[a-z]*\d*',s)

def _source_identity_match(device, source):
    """Exact canonical key OR strong, non-generic material ingredient evidence."""
    d=canonical_key(device);s=canonical_key(source)
    if not d or not s:return False
    if d==s:return True
    # Multi-component device labels have ingredient tokens; require every
    # discriminating ingredient, so plain C60 cannot match Cl6SubPc/C60.
    ingredients=[p for p in re.split(r'[/+,]',str(device).lower().replace('₆','6').replace('₂','2'))
                 if len(re.sub(r'[^a-z0-9]','',p))>=3]
    if len(ingredients)<2:return False
    return all(canonical_key(p) in s for p in ingredients)

def _figure_scan(doc,reference_text):
    """Source-caption scan identity for prose values cited to a J-V panel."""
    ref=re.search(r'(?i)\bFig(?:ure)?\.?\s*(\d+)\s*([a-z])\b',reference_text)
    if not ref:return None
    figure,letter=ref.group(1),ref.group(2).lower()
    for page in doc:
        t=' '.join(page.get_text('text').split())
        m=re.search(r'(?i)\bFig(?:ure)?\.?\s*'+re.escape(figure)+r'\s*[.|:]\s',t)
        if not m:continue
        cap=t[m.end():m.end()+1450]
        pm=re.search(r'\('+re.escape(letter)+r'\)\s',cap,re.I)
        if not pm:continue
        rest=cap[pm.end():]
        nxt=re.search(r'\([a-z](?:[–-][a-z])?\)\s',rest,re.I)
        if nxt:rest=rest[:nxt.start()]
        rev=bool(re.search(r'(?i)\breverse\s+scan\b|\bbackward\s+scan\b',rest))
        fw=bool(re.search(r'(?i)\bforward\s+scan\b',rest))
        if rev!=fw:return 'reverse' if rev else 'forward'
    return None

def extract(path):
    path=Path(path);out=[]
    if not path.is_file() or path.suffix.lower()!='.pdf':return out
    with fitz.open(path) as doc:
        for pi,page in enumerate(doc):
            t=' '.join(page.get_text('text').split())
            if not re.search(r'\bPCE\b',t,re.I):continue
            # Search near the numerical PCE block, including a trailing device
            # identifier; avoid crossing independent PV measurement blocks.
            for match in PCE_CLAUSE.finditer(t):
                frag=t[match.start():min(len(t),match.start()+420)]
                boundary=re.search(r'\b(?:In contrast|On the other hand|By comparison|Compared with)\b',frag,re.I)
                if boundary:frag=frag[:boundary.start()]
                vals=metrics_of(frag)
                pce_value=re.search(r'\bPCE\s+(?:can\s+be\s+)?(?:further\s+)?(?:improved|increased|enhanced|raised|reached|achieved|obtained|reported)\s+to\s+(\d{1,2}(?:\.\d+)?)\s*%',frag,re.I)
                if pce_value:vals['pce0']=float(pce_value.group(1))
                if not all(k in vals for k in ('pce0','voc0','jsc0','ff0')):continue
                if not valid_set(vals,allow_partial=False):continue
                names=list(NAME_AFTER.finditer(frag))
                if not names:continue
                for nm in names:
                    name=nm.group('name').strip(' .,:;')
                    # Crop conjunctions and article text from name, retaining
                    # the terminal chemical architecture clause.
                    name=re.split(r'\b(?:as shown|which|respectively)\b',name,re.I)[0].strip()
                    if not re.search(r'[A-Za-z].*\d|\d.*[A-Za-z]',name):continue
                    direction=_figure_scan(doc,frag)
                    record=dict(device_label=name,scan_direction=direction,values=vals,
                        source_file=path.name,source_page=pi+1,
                        source_excerpt=frag[:360],method='explicit_complete_prose')
                    if record not in out:out.append(record)
    return out

def recover(rows,main_pdf,supplements=()):
    candidates=[c for path in [main_pdf,*supplements] for c in extract(path)]
    for row in rows:
        match=[c for c in candidates if _source_identity_match(row.get('device_category',''),c['device_label'])]
        row.setdefault('narrative_jv_candidate_count',len(match))
        if not match:continue
        distinct={tuple(round(c['values'][k],5) for k in ('pce0','voc0','jsc0','ff0')) for c in match}
        if len(distinct)!=1:
            row['narrative_jv_status']='AMBIGUOUS_MULTIPLE_COMPLETE_SETS';continue
        c=match[0]
        # Narrative without scan direction is permitted only if the user has
        # not supplied forward/reverse measurements for this identity.
        direction=c.get('scan_direction') or canonical_scan_label(c['device_label'])[1]
        if row.get('jv_selected_scan_direction') in ('forward','reverse'):
            row['narrative_jv_status']='SCAN_QUALIFIED_MEASUREMENT_TAKES_PRECEDENCE'
            continue
        for k,v in c['values'].items():
            if k not in ('pce0','voc0','jsc0','ff0'):continue
            if row.get(k) not in ('',None):
                try:
                    if abs(float(row[k])-float(v))>0.00001:
                        row[k+'_previous_value']=row[k]
                        row[k+'_previous_source']=row.get(k+'_source')
                except (TypeError,ValueError):pass
            row[k]=float(v)
            row[k+'_source']=f'Complete source-aligned prose, scan {direction or "not specified"}; {c["source_file"]} p.{c["source_page"]}'
        row['narrative_jv_status']='SOURCE_LINKED_COMPLETE_SET_SCAN_'+(direction or 'UNSPECIFIED').upper()
        if direction:
            row['jv_selected_scan_direction']=direction
            row['jv_scan_selection_status']='NARRATIVE_JV_FIGURE_CAPTION_SCAN_'+direction.upper()
            row['jv_selected_measurement_source']=c['source_file']+' p.'+str(c['source_page'])
        row['narrative_jv_evidence']=c['source_excerpt']
    return candidates
