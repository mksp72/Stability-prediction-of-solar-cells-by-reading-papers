"""HI extraction for perovskite papers. No paper identifiers or fixed numeric values.

Signed HI convention: 100 * (PCE_reverse - PCE_forward)/PCE_reverse.
Absolute HI is its magnitude. Only same-category, same-condition F/R pairs
from source-aligned J-V table rows may generate a calculated HI.
"""
from __future__ import annotations
import re, math
from pathlib import Path
from collections import defaultdict
from typing import Optional

FORWARD={'f','fw','fwd','forward','forwards','fs','forwardscan','forwardscanning'}
REVERSE={'r','re','rw','rev','reverse','reversed','rs','b','bw','bwd','backward','backwards','back','backwardscan','reversescan'}

def canonical_scan_label(label):
    """Return (device formulation, scan direction); strip only separate qualifiers.

    Never strip 'F' within FAPbI3, 'B' from BCP or trailing R from a material.
    """
    s=str(label or '').strip()
    direction=None
    # Qualifiers in brackets or after a delimiter are supported, but do not
    # erase the remainder of a chemical name.
    tokens=re.findall(r'[A-Za-z0-9]+|[^A-Za-z0-9]+',s)
    lexical=[(i,v) for i,v in enumerate(tokens) if re.fullmatch(r'[A-Za-z]+',v)]
    if lexical:
        i,word=lexical[-1]
        w=word.lower()
        if w in FORWARD|REVERSE and (len(lexical)>1):
            direction='forward' if w in FORWARD else 'reverse'
            # discard this tail and a trailing 'scan' token if present
            start=sum(len(x) for x in tokens[:i])
            s=s[:start].rstrip(' :-_/([,')
        elif len(lexical)>1 and lexical[-1][1].lower() in ('scan','scanning'):
            j,w2=lexical[-2];w2=w2.lower()
            if w2 in FORWARD|REVERSE:
                direction='forward' if w2 in FORWARD else 'reverse'
                s=s[:sum(len(x) for x in tokens[:j])].rstrip(' :-_/([,')
    return s, direction

def canonical_key(label):
    s,_=canonical_scan_label(label)
    s=s.casefold().replace('₂','2').replace('₀','0')
    s=re.sub(r'(?i)solution[ -]?processed','solution',s)
    s=re.sub(r'(?i)(?:co[ -]?)?evaporated','evaporation',s)
    s=re.sub(r'(?i)\bsno\s*2\b','sno2',s)
    s=re.sub(r'(?i)\b(devices?|solar cells?|pscs?)\b','',s)
    key=re.sub(r'[^a-z0-9]+','',s)
    return re.sub(r'(?:treated|treatment)$','',key)

def calculate_hi(reverse_pce,forward_pce):
    r=float(reverse_pce);f=float(forward_pce)
    if not (math.isfinite(r) and math.isfinite(f) and r>0 and f>0):
        raise ValueError('Paired positive, finite PCE values required')
    signed=100*(r-f)/r
    return round(signed,6),round(abs(signed),6)

def _read_table_candidates(path):
    import fitz
    from eta_supplement_pv_recovery import _from_pymupdf_tables,_from_spatial,_from_lines
    extracted=[]
    if path.suffix.lower()!='.pdf' or not path.is_file():return extracted
    with fitz.open(path) as doc:
        for i in range(len(doc)):
            # Only page-native (no costly OCR here). Each result is one aligned JV row.
            text=doc[i].get_text('text')
            if not re.search(r'(?:\bPCE\b|efficien|\bVoc\b|\bJsc\b|fill factor|\bFF\b|\u03b7)',text,re.I):continue
            for method in (_from_pymupdf_tables,_from_spatial,_from_lines):
                try: records=method(doc,i)
                except Exception:records=[]
                for name,values,page,evidence in records:
                    base,scan=canonical_scan_label(name)
                    p=values.get('pce0')
                    if scan and p is not None and 0<float(p)<45:
                        extracted.append(dict(device_label=base,scan_direction=scan,pce_pct=float(p),
                            source_file=path.name,source_page=page,evidence=evidence,
                            methodology=method.__name__))
    return extracted

def _reported_hi_candidates(paths):
    """Locate explicitly published HI with a percent sign; do not guess its convention."""
    import fitz
    found=[]
    rx=re.compile(r'(?i)(?:\b(?:hysteresis\s+(?:index|factor)|HI|HF)\b)\s*(?:\([^)]{0,40}\))?\s*(?:=|:|of|was|is)\s*([-+]?\d+(?:\.\d+)?)\s*(%)?')
    for path in paths:
        path=Path(path)
        if path.suffix.lower()!='.pdf' or not path.is_file():continue
        try:
            with fitz.open(path) as doc:
                for page_index,page in enumerate(doc):
                    raw=page.get_text('text')
                    for match in rx.finditer(raw):
                        try: val=float(match.group(1))
                        except ValueError:continue
                        if not math.isfinite(val):continue
                        context=raw[max(0,match.start()-160):min(len(raw),match.end()+130)]
                        found.append(dict(value=val,unit='percent' if match.group(2) else 'unspecified',
                            context=context,source=f'{path.name} p.{page_index+1}'))
        except Exception:pass
    return found

def recover(rows,main_pdf,supplements=(),extra_candidates=()):
    """Fill ML HI features only for exact unambiguous paired scan measurements.

    Never synthesize a missing scan from champion PCE or from etaMPPT.
    Does not modify J-V metrics or the stability-fit pipeline.
    Returns full per-publication scan candidates for audit export.
    """
    candidates=[]
    for source in [main_pdf,*supplements]:
        try:candidates.extend(_read_table_candidates(Path(source)))
        except Exception:pass
    for c in extra_candidates:
        label=c.get('device_label','');base,scan=canonical_scan_label(label)
        if scan and c.get('pce0') not in ('',None):
            try:val=float(c['pce0'])
            except (ValueError,TypeError):continue
            candidates.append(dict(device_label=base,scan_direction=scan,pce_pct=val,
                source_file=c.get('file',''),source_page=c.get('page',''),
                methodology='figure_aligned_table',evidence=c.get('status','')))
    # Deduplicate three independent text extraction routes for the same source row.
    distinct={}
    for c in candidates:
        key=(c['source_file'],str(c['source_page']),canonical_key(c['device_label']),
             c['scan_direction'],round(float(c['pce_pct']),4))
        distinct.setdefault(key,c)
    candidates=list(distinct.values())
    reported=_reported_hi_candidates([main_pdf,*supplements])
    for row in rows:
        row.setdefault('hi_reported_pct',None)
        row.setdefault('hi_signed_pct',None)
        row.setdefault('hi_absolute_pct',None)
        row.setdefault('hi_reverse_pce_pct',None)
        row.setdefault('hi_forward_pce_pct',None)
        row.setdefault('hi_reported_source',None)
        row.setdefault('hi_reported_raw',None)
        row.setdefault('hi_reported_unit',None)
        row.setdefault('hi_source',None)
        row.setdefault('hi_definition','100*(PCE_reverse-PCE_forward)/PCE_reverse')
        row.setdefault('hi_status','NO_PAIRED_FORWARD_REVERSE_JV_PCE')
        row.setdefault('scan_direction','category_unqualified')
        row['device_category_canonical']=canonical_scan_label(row.get('device_category',''))[0]
        key=canonical_key(row['device_category_canonical'])
        # Explicitly reported HI remains separate from the calculated standard.
        # Require one unambiguous source-context match to this formulation.
        matched_reported=[z for z in reported if key and key in canonical_key(z['context'])]
        if len(matched_reported)==1:
            z=matched_reported[0]
            row['hi_reported_raw']=z['value'];row['hi_reported_unit']=z['unit']
            row['hi_reported_source']=z['source']
            if z['unit']=='percent':row['hi_reported_pct']=z['value']
        matching=[c for c in candidates if canonical_key(c['device_label'])==key]
        # Do not pair across differing pages/source files, or mix multiple tests.
        groups=defaultdict(lambda: defaultdict(set))
        for c in matching:
            groups[(c['source_file'],str(c['source_page']))][c['scan_direction']].add(round(c['pce_pct'],5))
        pairs=[]
        for (file,page),scans in groups.items():
            if len(scans['forward'])==1 and len(scans['reverse'])==1:
                pairs.append((file,page,next(iter(scans['reverse'])),next(iter(scans['forward']))))
        if len(pairs)!=1:
            if len(pairs)>1:row['hi_status']='MULTIPLE_PAIRED_MEASUREMENTS_SELECT_CONDITION'
            elif matching:row['hi_status']='MISSING_OR_AMBIGUOUS_SCAN_PAIR'
            continue
        file,page,r,f=pairs[0]
        signed,absol=calculate_hi(r,f)
        row.update(hi_signed_pct=signed,hi_absolute_pct=absol,
                   hi_reverse_pce_pct=r,hi_forward_pce_pct=f,
                   hi_status='CALCULATED_FROM_SAME_SOURCE_JV_PAIR',
                   hi_source=f'{file} p.{page}: forward/reverse paired JV table')
    return candidates
