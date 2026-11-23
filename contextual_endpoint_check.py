"""v18.6.4: optional category-specific cross-check against original article prose.

This is a diagnostic only. No digitized points, rates, AIC/RMSE models, T80,
or original published claims are replaced. A named category, percentage, and
ageing time must be in the same short passage of the same publication.
"""
from pathlib import Path
import re, math, fitz, numpy as np


def _clean(s):
    s=s.replace('−','-').replace('–','-')
    s=re.sub(r'(?<=[A-Za-z])-\s+(?=[a-z])','-',s)
    return ' '.join(s.split())


def _extract_claim(text,label,all_labels=()):
    if len(label)<3 or not re.search(r'[A-Za-z]',label):return []
    # Prevent matching a base material FACs within FACs@PNIC.
    match=re.compile(r'(?<![A-Za-z0-9@])'+re.escape(label)+r'(?![A-Za-z0-9@])',re.I)
    claims=[]
    for m in match.finditer(text):
        passage=text[m.end():m.end()+220]
        # Stop before *another named formulation*, including cases where the
        # next category is embedded in 'outperforming FACs@PIC (85%)'.
        for other in all_labels:
            if other.lower()==label.lower():continue
            nm=re.search(r'(?<![A-Za-z0-9@])'+re.escape(other)+r'(?![A-Za-z0-9@])',passage,re.I)
            if nm:passage=passage[:nm.start()]
        # Never cross into a second comparative clause or another device.
        nextnamed=re.search(r'\b(?:whereas|by comparison|while the|however, the|in contrast, the)\b',passage,re.I)
        if nextnamed:passage=passage[:nextnamed.start()]
        r=re.search(r'\b(?:retaining|maintains?|retains?|degrades?|fall(?:s|ing)?|drops?|remain(?:s|ing)?)\b'
                    r'.{0,65}?(\d{2,3}(?:\.\d+)?)\s*%'
                    r'.{0,65}?\b(?:after|within|over|for)\s*'
                    r'([\d,]{2,7}(?:\.\d+)?)\s*(?:hours?|hrs?|h)\b',passage,re.I)
        if r:
            retention=float(r.group(1));hours=float(r.group(2).replace(',',''))
            if 0<=retention<=125 and 10<=hours<=100000:
                claims.append((retention,hours,passage[:r.end()]))
    # If one candidate is clearly identified, return it; if competing claims
    # for the same device are present, do NOT assume they are the same test.
    unique={(round(v,2),round(t,2)) for v,t,_ in claims}
    return claims if len(unique)==1 else []


def verify(rows,main_pdf,supplements=()):
    try:
        with fitz.open(str(main_pdf)) as d:
            pages=[(i+1,_clean(p.get_text('text'))) for i,p in enumerate(d)]
    except Exception:return 0
    changed=0
    for row in rows:
        if not row.get('stability_data_json'):continue
        category=str(row.get('device_category') or '')
        names=[str(r.get('device_category') or '') for r in rows]
        claims=[(pg,*claim) for pg,text in pages for claim in _extract_claim(text,category,names)]
        if len({(round(q[1],2),round(q[2],2)) for q in claims})!=1:continue
        page,retention,hours,evidence=claims[0]
        # Do not compare to an unrelated dark/thermal trace.
        import json,base64,zlib
        payload=row.get('stability_data_json')
        if isinstance(payload,str) and payload.startswith('zlib64:'):
            payload=zlib.decompress(base64.b64decode(payload[7:])).decode('utf8')
        try:studies=json.loads(payload) if isinstance(payload,str) else payload
        except (TypeError,ValueError):continue
        preferred=next((s for s in studies if s.get('figure')==row.get('primary_stability_figure')),None)
        if not preferred:continue
        proto=str(preferred.get('protocol') or '').upper()
        if 'DARK' in proto or 'ISOS-D' in proto:continue
        try:
            tx=np.asarray(preferred['time_h'],float);ry=np.asarray(preferred['retention_pct'],float)
            ok=np.isfinite(tx)&np.isfinite(ry);tx=tx[ok];ry=ry[ok]
            if len(tx)<4 or tx[-1]<=0:continue
            if tx[0]<=hours<=tx[-1]:interp=float(np.interp(hours,tx,ry));type_='IN_CURVE'
            elif tx[-1]<hours<=tx[-1]*1.05:
                interp=float(ry[-1]);type_='NEAR_LAST_POINT_NOT_EXACT_INTERPOLATION'
            else:continue
        except (TypeError,ValueError,IndexError):continue
        err=round(interp-retention,3)
        row['source_endpoint_claim_pct']=retention
        row['source_endpoint_claim_h']=hours
        row['source_endpoint_digitized_pct']=round(interp,3)
        row['source_endpoint_abs_difference_pctpts']=round(abs(err),3)
        row['source_endpoint_interpolation_type']=type_
        row['source_endpoint_claim_page']=page
        row['source_endpoint_claim_evidence']=f'{Path(main_pdf).name} p.{page}: '+evidence[:180]
        row['source_endpoint_claim_validation']=('AGREES_WITHIN_3PCTPTS' if abs(err)<=3
                                                else 'DISCREPANCY_REQUIRES_REVIEW')
        changed+=1
    return changed
