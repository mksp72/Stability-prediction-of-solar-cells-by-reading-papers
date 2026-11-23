"""Filename-independent discovery of probable stability figures in unfamiliar PDFs.

This does NOT digitize curves; it only identifies evidence locations, so no
synthetic ML rows are created. A future general curve decoder can consume the
candidate page/figure list. No preconfigured page numbers or DOI allowlists.
"""
from __future__ import annotations
import re
from pathlib import Path

AGE=re.compile(r'\b(stabilit(?:y|ies)|degrad(?:ation|ed|ing)|aging|ageing|retention|maximum.power.point|MPP|ISOS[- ]?[DLTO]|operational|T80|T90|shelf.life|storage.life)\b',re.I)
FIG=re.compile(r'\b(?:Fig(?:ure)?\.?\s*(?:S?\d+[A-Za-z]?(?:[-–][A-Za-z])?))',re.I)
NUMTIME=re.compile(r'\b(?:\d{2,6}\s*(?:h|hours?|hrs?|days?)|hours?\s+under\s+illumination)\b',re.I)

def stability_figures(pdf:Path,max_results=20):
    import fitz
    found=[]
    with fitz.open(str(pdf)) as document:
        for p,page in enumerate(document):
            text=page.get_text(sort=True)
            if not AGE.search(text):continue
            chunks=re.split(r'(?=\b(?:Fig(?:ure)?\.?\s*S?\d{1,3}[A-Za-z]?(?:[:.\s]|\))))',text,flags=re.I)
            for chunk in chunks:
                if not AGE.search(chunk):continue
                if not FIG.search(chunk) and not (AGE.search(chunk) and NUMTIME.search(chunk)):continue
                labels=list(dict.fromkeys(x.group(0) for x in FIG.finditer(chunk)))[:5]
                snippet=re.sub(r'\s+',' ',chunk[:360]).strip()
                score=int(bool(labels))*2+int(bool(NUMTIME.search(chunk)))+int(bool(re.search(r'(?:ISOS|MPP|retention|stability)',chunk,re.I)))
                found.append({'page':p+1,'figure_labels':', '.join(labels),
                              'candidate_score':score,'caption_excerpt':snippet[:280]})
    # Keep highest-evidence candidate per (page, labels), stable output order.
    unique={}
    for c in found:
        key=(c['page'],c['figure_labels'])
        if key not in unique or c['candidate_score']>unique[key]['candidate_score']:unique[key]=c
    return sorted(unique.values(),key=lambda x:(-x['candidate_score'],x['page']))[:max_results]
