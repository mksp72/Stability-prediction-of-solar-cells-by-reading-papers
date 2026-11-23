"""Reusable same-red-edge EQE band-gap completion for a matched publication.

Priority: actual calibrated EQE spectra -> own 50% red edge.  If curves
share their long-wavelength cutoff (within tolerance), the same Eg may be
used for their common absorber.  Never inherit a 3D gap into pure 2D,
or a layer-gap/model parameter into a different absorber.  Existing measured
values are never overwritten. No DOI or sample-specific numeric lookups.
"""
from __future__ import annotations
import re
from collections import defaultdict
from pathlib import Path

HC=1239.841984

def number(x):
    try:
        from math import isfinite
        v=float(x);return v if isfinite(v) else None
    except (ValueError,TypeError):return None

def norm(x):return re.sub(r'[^a-z0-9]+','',str(x).lower())

def _pure_2d(row):
    label=str(row.get('device_category') or '').lower()
    composition=str(row.get('absorber_formula') or '').lower()
    return bool(re.search(r'\b(?:pure|phase[- ]pure|standalone)\s*(?:2d|quasi[- ]?2d)\b',label))

def _absorber_compatible(a,b):
    if _pure_2d(a)!=_pure_2d(b):return False
    af=norm(a.get('absorber_formula') or '')
    bf=norm(b.get('absorber_formula') or '')
    if af and bf and af!=bf:return False
    aa=str(a.get('architecture') or '').lower();bb=str(b.get('architecture') or '').lower()
    if aa and bb and aa not in ('unknown','none') and bb not in ('unknown','none') and aa!=bb:return False
    return True

def _aliases(row):
    name=str(row.get('device_category') or '')
    choices=[name]
    if re.search(r'^control\b',name,re.I) and not _pure_2d(row):choices+=['control','reference','ref']
    name=re.sub(r'\s+(?:stability group|champion|average)$','',name,flags=re.I)
    if len(norm(name))>=3: choices.append(name)
    return list(dict.fromkeys(c for c in choices if len(norm(c))>=3))

def _label_contains(group,row):
    text=str(group).lower()
    if not text:return False
    if _pure_2d(row) and not re.search(r'\b(?:pure|phase[- ]pure|standalone)\s*(?:2d|quasi[- ]?2d)\b',text,re.I):
        return False
    for alias in _aliases(row):
        if norm(alias)==norm(text):return True
        # Delimited list: "Control / Target-1 / Target-2", but never
        # match generic "target" inside Target-2 in an unrelated group.
        if re.search(r'(?<![a-z0-9])'+re.escape(alias.lower())+r'(?![a-z0-9-])',text):return True
    return False

def _numeric_spectral_group(csv_path):
    """Return actual red-edge results from numerical EQE trace rows, when present."""
    import csv
    from generic_spectral_math import spectral_50_bandgap
    groups=defaultdict(lambda: {'x':[],'y':[]})
    if not Path(csv_path).is_file():return []
    with open(csv_path,encoding='utf-8-sig',newline='') as f:
        for row in csv.DictReader(f):
            x=number(row.get('wavelength_nm'))
            y=number(row.get('ipce_pct'))
            if x is None or y is None:continue
            key=(str(row.get('doi') or ''),str(row.get('sample_group') or ''),str(row.get('source_page') or ''))
            groups[key]['x'].append(x);groups[key]['y'].append(y)
    records=[]
    for (doi,label,page),xy in groups.items():
        # Removing duplicates is essential for monotonic wavelength interpolation.
        xy_pairs=sorted(set(zip(xy['x'],xy['y'])))
        filtered={}
        for x,y in xy_pairs:filtered.setdefault(x,[]).append(y)
        try:
            from numpy import median
            out=spectral_50_bandgap(list(filtered),[float(median(v)) for v in filtered.values()])
            records.append({'doi':doi,'label':label,'page':page,'lambda50_nm':out['lambda50_nm'],
                            'bandgap_ev':out['estimated_bandgap_ev'],'source':'digitized EQE spectrum'})
        except (ValueError,TypeError,ImportError):pass
    return records

def recover_eqe_bandgaps(rows, main_pdf, supplementary_paths=(), spectral_root=None,
                         cutoff_overlap_tolerance_nm=8.0):
    """Fill only genuinely supported missing gap fields, preserving source."""
    if not rows:return 0
    records=[]
    if spectral_root:
        root=Path(spectral_root)
        records=_numeric_spectral_group(root/'spectral_curve_points.csv')
        # Existing EQE edge extraction may store the calibrated 50% result in
        # bandgap_results.csv rather than a full trace (particularly Science).
        file=root/'bandgap_results.csv'
        if file.is_file():
            import csv
            for b in csv.DictReader(file.open(encoding='utf-8-sig',newline='')):
                lam=number(b.get('lambda50_nm'))
                ev=number(b.get('estimated_bandgap_ev'))
                if lam and ev and 0.9<ev<3.3:
                    records.append({'doi':b.get('doi') or '','label':b.get('sample_group') or '',
                                    'page':b.get('source_page') or '',
                                    'lambda50_nm':lam,'bandgap_ev':ev,
                                    'source':'calibrated EQE 50% cutoff'})
    # Main-paper DOI is the authoritative scope, so a shared label in another
    # PDF in a batch can NEVER contribute a cutoff to this publication.
    wanted=norm(rows[0].get('doi') or '')
    if wanted:
        records=[rec for rec in records if norm(rec.get('doi'))==wanted]
    recovered=0
    # Detect overlapping edges across actual calibrated traces before sharing.
    # A single trace explicitly labeled for multiple formulations is also
    # strong shared-cutoff evidence; cannot share to unrelated 2D controls.
    for i,r in enumerate(rows):
        if number(r.get('bandgap_ev')) is not None:continue
        compatible=[rec for rec in records if _label_contains(rec['label'],r)]
        if not compatible:
            continue
        # Distinct curves assigned to the same category with different cutoffs
        # do NOT justify averaging or a false single gap.
        wave=[v['lambda50_nm'] for v in compatible]
        if max(wave)-min(wave)>cutoff_overlap_tolerance_nm:continue
        import statistics
        lam=statistics.median(wave)
        r['bandgap_ev']=round(HC/lam,4)
        r['bandgap_source']='EQE red-edge half-height; lambda50='+str(round(lam,2))+' nm; '+', '.join(sorted({v['source'] for v in compatible}))
        r['bandgap_lambda50_nm']=round(lam,3)
        r['bandgap_shared_edge']='YES' if any('/' in str(v['label']) for v in compatible) else 'NO'
        recovered+=1

    # For new chart layouts, allow SAME observed edge to be propagated only
    # where both categories have numerical EQE-cutoff evidence with matching
    # wavelengths. A single paper-wide Eg number alone is not proof of overlap.
    for i,a in enumerate(rows):
        if number(a.get('bandgap_ev')) is not None:continue
        compatible=[]
        for j,b in enumerate(rows):
            if i==j or number(b.get('bandgap_ev')) is None or not _absorber_compatible(a,b):continue
            if number(b.get('bandgap_lambda50_nm')) is None:continue
            sources=[r for r in records if _label_contains(r['label'],a)]
            if sources and all(abs(r['lambda50_nm']-number(b['bandgap_lambda50_nm']))<=cutoff_overlap_tolerance_nm for r in sources):
                compatible.append((b,sources))
        if len(compatible)==1:
            source,spectra=compatible[0]
            a['bandgap_ev']=number(source['bandgap_ev'])
            a['bandgap_source']='Overlapping calibrated EQE long-wavelength cutoff with '+str(source.get('device_category'))
            a['bandgap_lambda50_nm']=number(source.get('bandgap_lambda50_nm'))
            a['bandgap_shared_edge']='YES'
            recovered+=1
    return recovered
