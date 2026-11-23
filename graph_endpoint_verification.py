"""Generic, source-backed checks for automatic stability digitization (v18.6.2).

Does NOT alter measured points, fitting models, PV values or T80. Extracts
reported time/retention anchors from the SAME publication, matches by category,
architecture and source-page context, interpolates the digitized retention, and
adds diagnostics used in primary-curve selection. No per-DOI configurations.
"""
from __future__ import annotations
import json
import math
import re
from pathlib import Path

SPACE = re.compile(r'\s+')
NUM = r'(?:\d{1,3}(?:,\d{3})+|\d{1,7})(?:\.\d+)?'
TIME = rf'(?P<time>{NUM})\s*(?:h(?:ours?|rs?)?\b)'
PCT = rf'(?P<pct>\d{{1,3}}(?:\.\d+)?)\s*%'
# Do not interpret arbitrary efficiency or JV results as stability anchors.
STABILITY_V = r'(?:retain(?:ed|s|ing)?|maintain(?:ed|s|ing)?|remain(?:ed|s|ing)?|drop(?:ped|s|ping)?\s+to|declin(?:ed|es|ing)\s+to|degrad(?:ed|es|ing)\s+to|surviv(?:ed|es|ing)?\s+with)'
ANCHOR_RX = re.compile(
    rf'(?P<prefix>(?:{STABILITY_V})\s+(?:[^.;]{{0,70}}?))?{PCT}'
    rf'(?P<between>[^.;]{{0,105}}?)\b(?:after|at|for|over)\s+(?:around|approximately|nearly|about|~)?\s*{TIME}',
    re.I,
)
# Source statements like "T80 = 640 h" or "T80,Champ ≈ 1,080 h".
T80_RX = re.compile(rf'\bT\s*80(?:\s*[,\- ]\s*[A-Za-z]{{2,12}})?\s*(?:=|≈|~|:|of|after|\s+is)\s*(?:around|approximately|nearly|about|~)?\s*(?P<time>{NUM})\s*h\b', re.I)
RESPECTIVE_T80_RX = re.compile(
    rf'\bT\s*80\b[^.;]{{0,40}}?\bafter\s*[~≈]?\s*(?P<t1>{NUM})\s+and\s*[~≈]?\s*(?P<t2>{NUM})\s*h\s+for\s+the\s+(?P<label1>[A-Za-z0-9]+)\s*[-–]?\s+and\s+(?P<label2>[A-Za-z][A-Za-z\-\s]{{2,38}}?)(?:\s*\((?P<alias2>[A-Z0-9]{{2,8}})\))?[-–]?based\s+devices,?\s+respectively',re.I)
FIG_RX = re.compile(r'\b(?:Fig(?:ure)?\.?\s*|panel\s*)(\d+)([a-z])?\b',re.I)


def normalize_text(text):
    return SPACE.sub(' ',str(text or '').replace('\u2009',' ').replace('\u202f',' ')
                      .replace('\u00a0',' ').replace('ﬁ','fi').replace('ﬃ','ffi').replace('₈','8').replace('₀','0'))


def _aliases(label):
    s=normalize_text(label).lower()
    tokens=re.findall(r'[a-z]+(?:-[a-z0-9]+)?|\d+[a-z]+',s)
    banned={'with','without','based','device','devices','evaporation','solution','sample','cells','cell','treated','and','the'}
    out={x for x in tokens if len(x)>2 and x not in banned}
    if 'control' in out or 'pristine' in out or 'reference' in out or 'ref' in out:
        out.update(('control','pristine','reference'))
    # Don't let 'target' match every different 'Target-1' category.
    if 'target' in out:out.remove('target')
    return out


def _local_matches(text, label, other_labels=()):
    words=_aliases(label)
    if not words:return False
    found={x for x in words if re.search(r'(?<![a-z0-9])'+re.escape(x)+r'(?![a-z0-9])',text,re.I)}
    if not found:return False
    # Specific distinctions such as BCP versus SnO2 or TPTA versus control.
    competitors=set()
    for other in other_labels:
        other_words=_aliases(other)
        if other_words==words and normalize_text(other).casefold()!=normalize_text(label).casefold():
            # A claim merely saying 'BCP' cannot distinguish solution BCP
            # from evaporated BCP: do not claim to validate either one.
            return False
        competitors|=other_words-words
    for term in competitors:
        if re.search(r'(?<![a-z0-9])'+re.escape(term)+r'(?![a-z0-9])',text,re.I):
            return False
    return True


def _page_claims(text,page,source):
    clean=normalize_text(text)
    claims=[]
    # Sentence windows preserve sample names without searching entire article.
    for match in ANCHOR_RX.finditer(clean):
        prefix=match.group('prefix') or ''
        if not re.search(STABILITY_V, prefix, re.I):
            continue
        # The stability verb must directly govern this percentage, not an
        # unrelated efficiency elsewhere in the paragraph.
        if re.search(r'(?:\b(?:PCE|MPPT|efficiency)\s*(?:of|=|≈|:)?\s*)$', prefix[-24:], re.I):
            continue
        t=float(match.group('time').replace(',',''))
        pct=float(match.group('pct'))
        if not (1<=t<=200000 and 0<=pct<=120):continue
        w=clean[max(0,match.start()-160):min(len(clean),match.end()+75)]
        immediate=clean[max(0,match.start()-75):match.end()]
        immediate=re.split(r'\b(?:whereas|while|compared\s+with|but)\b',immediate,flags=re.I)[-1]
        # Exclude efficiency statements and unrelated statistics.
        if re.search(r'\b(?:certified|champion|initial PCE|J[ -]?V)\b',w[:max(0,match.start()-max(0,match.start()-160))],re.I) and 'retain' not in w.lower():continue
        claims.append({'type':('NEAR_T80_RETENTION' if 80.0<=pct<81.0 else 'RETENTION'),'time_h':t,'retention_pct':pct,'page':page,
                       'source':source,'context':w,'local_context':immediate,'page_start':match.start()})
    # Two devices and two explicitly ordered lifetimes in a 'respectively'
    # statement: map the pairs instead of assigning both T80s to all curves.
    for match in RESPECTIVE_T80_RX.finditer(clean):
        for idx,(time,lab) in enumerate(((match.group('t1'),match.group('label1')),
                                        (match.group('t2'),match.group('alias2') or match.group('label2')))):
            claim_time=float(time.replace(',',''))
            label=normalize_text(lab)
            claims.append({'type':'T80','time_h':claim_time,'retention_pct':80.,'page':page,
                 'source':source,'context':clean[max(0,match.start()-100):match.end()+50],
                 'local_context':label+' device T80 '+time+' h',
                 'page_start':match.start(),'respective_pair':True})
    for match in T80_RX.finditer(clean):
        t=float(match.group('time').replace(',',''))
        if not 1<=t<=200000:continue
        w=clean[max(0,match.start()-115):min(len(clean),match.end()+95)]
        claims.append({'type':'T80','time_h':t,'retention_pct':80.,'page':page,
                       'source':source,'context':w,'local_context':w,'page_start':match.start()})
    # De-duplicate repeated PDF text lines / same source statements.
    seen=set();result=[]
    for c in claims:
        key=(c['type'],c['time_h'],c['retention_pct'],c['page'])
        if key not in seen:result.append(c);seen.add(key)
    return result


def publication_claims(main_pdf, supplement_files=()):
    try:import fitz
    except ImportError:return []
    all_claims=[]
    for path in [Path(main_pdf)]+[Path(f) for f in supplement_files if str(f).lower().endswith('.pdf')]:
        try:
            with fitz.open(str(path)) as doc:
                for page_num,page in enumerate(doc):
                    all_claims.extend(_page_claims(page.get_text('text'),page_num+1,path.name))
        except (OSError,ValueError,RuntimeError):continue
    return all_claims


def _architecture_matches(study,context):
    kind=normalize_text(study.get('plot_architecture') or '').lower()
    s=context.lower()
    if any(x in kind for x in ('tandem','2t','2-terminal')):
        if 'single-junction' in s or 'single junction' in s:return False
    if 'single' in kind and ('tandem' in s or 'ptsc' in s or '2-terminal' in s):return False
    return True


def _claim_score(claim,study,category,other_labels):
    ctx=claim.get('local_context') or claim['context'];low=ctx.lower()
    full=claim['context'].lower()
    if not _architecture_matches(study,ctx):return -1
    if not _local_matches(ctx,category,other_labels):return -1
    # Restrict matching to exact document and figure neighborhood; distant
    # abstract prose is lower confidence than same-page Figure results.
    src_score=2 if claim['source']==study.get('pdf') else 0
    page_gap=abs(claim['page']-int(study.get('page') or 0))
    page_score=4 if page_gap==0 else 2 if page_gap<=2 else 0
    if page_gap>4:return -1
    if re.search(r'\b(?:dark|shelf storage|long.term storage|glovebox)\b',full) and str(study.get('protocol','')).startswith('ISOS-L'):return -1
    if re.search(r'continuous.*(?:mpp|operation)|illumination|mpp tracking|isos-l',full) and str(study.get('protocol','')).startswith('ISOS-D'):return -1
    # Avoid turning a generic result about 'TPTA' into the wrong figure's
    # endpoint when multiple panels show different experiments.
    return src_score+page_score+int(bool(claim['type']=='RETENTION'))


def verify_study(study,category,all_claims,other_labels=()):
    """Attach independent measured/source consistency checks to a study.

    A result is a cross-check, not an edited trace. A published checkpoint
    inside the curve's time span is compared at that exact time (interpolated).
    """
    points=study.get('retention_pct') or []
    times=study.get('time_h') or []
    if len(points)<2 or len(points)!=len(times):return {'status':'INSUFFICIENT_TRACE'}
    import numpy as np
    xs=np.asarray(times,float);ys=np.asarray(points,float)
    if not np.all(np.isfinite(xs)) or not np.all(np.isfinite(ys)) or not np.all(np.diff(xs)>0):
        return {'status':'INVALID_POINTS_OR_TIME_ORDER'}
    candidates=[]
    for c in all_claims:
        score=_claim_score(c,study,category,other_labels)
        if score<4:continue
        if c['time_h']<xs[0] or c['time_h']>xs[-1]*1.03:
            # A source checkpoint beyond the extracted graph's reach is a
            # coverage warning, not a measured/fit mismatch.
            candidates.append((score,{**c,'status':'OUTSIDE_DIGITIZED_WINDOW'}));continue
        yi=float(np.interp(c['time_h'],xs,ys))
        difference=round(yi-c['retention_pct'],3)
        # One-sided user-defined near-T80 endpoint allowance: [80, 81).
        # Never call 79.x a near-T80 observation. This check is not the
        # mathematical T80 root of the fitted model.
        near_80=(80.0<=yi<81.0)
        if c['type']=='T80':
            status='PASS' if near_80 else 'BELOW_80_NOT_NEAR_T80' if yi<80 else 'T80_RETENTION_MISMATCH'
        elif c['type']=='NEAR_T80_RETENTION':
            status='PASS' if near_80 else 'BELOW_80_NOT_NEAR_T80' if yi<80 else 'T80_RETENTION_MISMATCH'
        else:
            status='PASS' if abs(difference)<=7.5 else 'MISMATCH'
        candidates.append((score,{**c,'digitized_retention_pct':round(yi,3),
             'difference_pctpts':difference,'near_t80_one_sided_tolerance':near_80,
             'status':status}))
    valid=[(score,c) for score,c in candidates if c['status'] in ('PASS','MISMATCH','BELOW_80_NOT_NEAR_T80','T80_RETENTION_MISMATCH')]
    if valid:
        # Higher context score first; within same context confidence, prefer
        # closest reported checkpoint to curve's measured endpoint.
        score,c=max(valid,key=lambda z:(z[0],z[1]['time_h']/max(xs[-1],1)))
        return {'status':c['status'],'claim_type':c['type'],'source':c['source'],
                'page':c['page'],'reported_time_h':c['time_h'],
                'reported_retention_pct':c['retention_pct'],
                'digitized_retention_pct':c['digitized_retention_pct'],
                'difference_pctpts':c['difference_pctpts'],
                'evidence_context':c['context'][:270],'context_score':score,
                'supporting_claims':len(valid)}
    outside=[(sc,c) for sc,c in candidates if c['status']=='OUTSIDE_DIGITIZED_WINDOW']
    if outside:
        _,c=max(outside,key=lambda z:z[0])
        return {'status':'REPORTED_CHECKPOINT_OUTSIDE_DIGITIZED_WINDOW',
                'reported_time_h':c['time_h'],'source':c['source'],
                'page':c['page'],'evidence_context':c['context'][:270]}
    return {'status':'NO_UNAMBIGUOUS_MATCHING_PUBLISHED_CHECKPOINT'}




def _selected_fit_t80(study):
    summary=study.get('fit_summary') or {}
    value=summary.get('selected_fit_t80_h')
    if value is None:
        n=summary.get('selected_fit_components')
        value=summary.get(f'fit{n}_t80_h') if n in (2,3) else None
    try:
        t=float(value)
        return t if math.isfinite(t) and t>=0 else None
    except (ValueError,TypeError):return None


def _published_t80_candidates(study,category,claims,other_labels):
    """Only source-backed category-matched explicit T80 values or near-80 endpoints.

    A reported 80.7% at 1000 h is classified as a NEAR threshold
    observation, not misrepresented as an exact published T80 crossing.
    """
    found=[]
    for claim in claims:
        if claim['type'] not in ('T80','NEAR_T80_RETENTION'):continue
        score=_claim_score(claim,study,category,other_labels)
        if score<4:continue
        found.append((score,claim))
    return sorted(found,key=lambda x:(-x[0],abs(x[1]['page']-int(study.get('page') or 0))))


def compare_published_t80(study,category,claims,other_labels=()):
    """Reported-vs-selected-fit T80 diagnostic without refitting or replacing T80.

    Exactly published T80 and approximate 80%-81% endpoint are distinct.
    Refuses to choose between incompatible highest-confidence citations.
    """
    candidates=_published_t80_candidates(study,category,claims,other_labels)
    if not candidates:return {'status':'NO_CATEGORY_MATCHED_SOURCE_T80'}
    top=candidates[0][0]
    top_claims=[c for score,c in candidates if score==top]
    # Multiple incompatible T80s often arise from different stress conditions.
    # Preserve ambiguity instead of attaching the wrong value to a curve.
    if len({round(c['time_h'],3) for c in top_claims})>1:
        return {'status':'AMBIGUOUS_REPORTED_T80','candidate_times_h':sorted({c['time_h'] for c in top_claims}),
                'candidate_sources':[{'source':c['source'],'page':c['page']} for c in top_claims]}
    c=top_claims[0]
    is_exact=c['type']=='T80'
    reported=float(c['time_h'])
    fitted=_selected_fit_t80(study)
    result={'status':'VERIFIED_COMPARISON' if fitted is not None else 'NO_FINITE_SELECTED_FIT_T80',
        'evidence_type':'EXPLICIT_PUBLISHED_T80' if is_exact else 'NEAR_80_RETENTION_ENDPOINT',
        'reported_time_h':reported,'reported_retention_pct':float(c['retention_pct']),
        'source':c['source'],'page':c['page'],'evidence_context':c['context'][:300],
        'selected_fit_t80_h':fitted,'selected_fit_components':(study.get('fit_summary') or {}).get('selected_fit_components'),
        't80_tolerance_rule':'80.0 <= reported/digitized retention < 81.0 (one-sided; <80 excluded)'}
    if fitted is not None and reported>0:
        err=fitted-reported
        result['signed_error_h']=round(err,4)
        result['absolute_error_h']=round(abs(err),4)
        result['relative_absolute_error_pct']=round(100*abs(err)/reported,4)
        result['fit_discrepancy_flag_20pct']=bool(abs(err)/reported>0.2)
        result['comparison_interpretation']=('FIT_VS_PUBLISHED_T80' if is_exact
             else 'FIT_VS_NEAR_THRESHOLD_ENDPOINT_NOT_EXACT_T80')
    # Check actual curve at the published time separately from fitted T80.
    xs=study.get('time_h') or [];ys=study.get('retention_pct') or []
    if len(xs)>=2 and len(xs)==len(ys) and float(xs[0])<=reported<=float(xs[-1]):
        import numpy as np
        digitized=float(np.interp(reported,np.asarray(xs,float),np.asarray(ys,float)))
        result['digitized_retention_at_reported_time_pct']=round(digitized,4)
        result['digitized_near_t80_one_sided']=bool(80.0<=digitized<81.0)
        result['digitized_below_80_not_near_threshold']=bool(digitized<80.0)
    else:result['digitized_checkpoint_status']='OUTSIDE_DIGITIZED_TIME_WINDOW'
    return result

def annotate_publication(studies_by_category,main_pdf,supplement_files=()):
    claims=publication_claims(main_pdf,supplement_files)
    labels=[obj['label'] for obj in studies_by_category.values()]
    counts={}
    for obj in studies_by_category.values():
        other=[x for x in labels if x!=obj['label']]
        for s in obj['studies']:
            q=verify_study(s,obj['label'],claims,other)
            s['published_endpoint_verification']=q
            s['published_t80_comparison']=compare_published_t80(s,obj['label'],claims,other)
            stat=q['status'];counts[stat]=counts.get(stat,0)+1
            if stat in ('MISMATCH','BELOW_80_NOT_NEAR_T80','T80_RETENTION_MISMATCH'):
                s.setdefault('digitization_quality_flags',[]).append('PUBLISHED_CHECKPOINT_MISMATCH')
            if stat=='REPORTED_CHECKPOINT_OUTSIDE_DIGITIZED_WINDOW':
                s.setdefault('digitization_quality_flags',[]).append('PUBLISHED_CHECKPOINT_OUTSIDE_TRACE')
            if stat=='INVALID_POINTS_OR_TIME_ORDER':
                s.setdefault('digitization_quality_flags',[]).append('INVALID_TIME_MONOTONICITY')
    return counts


def choose_primary_verified(studies):
    """Prefer explicit source agreement within same protocol/identity quality.

    Changes ONLY choice among already-digitized traces, never point arrays or
    fitted parameters. Otherwise delegates to the original v18.6.1 selector.
    """
    from curve_integrity_guard import choose_primary, _quality_score
    if not studies:raise ValueError('No studies')
    # Protocol preference remains first. Verified graph checkpoints are used
    # only to adjudicate otherwise comparable traces within the same category.
    def rank(s):
        islight=int(str(s.get('protocol','')).startswith('ISOS-L'))
        q=s.get('published_endpoint_verification') or {}
        ver={'PASS':2,'NO_UNAMBIGUOUS_MATCHING_PUBLISHED_CHECKPOINT':1,
             'REPORTED_CHECKPOINT_OUTSIDE_DIGITIZED_WINDOW':1,
             'MISMATCH':0,'BELOW_80_NOT_NEAR_T80':0,'T80_RETENTION_MISMATCH':0}.get(q.get('status'),1)
        quality=_quality_score(s)
        # Pale/borrowed graph masks cannot win solely from a coincidental
        # published endpoint match.
        return (islight,quality[:5],ver,quality[5],float(s.get('t_end_h') or 0))
    return max(studies,key=rank)


def annotate_calibrated_rows(rows, main_pdf, supplementary_files=()):
    """Add the same source-T80 diagnostics to calibrated/v15 rows.

    Generic rows are already verified per trace and not rescanned. Nothing is
    re-fitted, renormalized, or changed in the stability JSON.
    """
    pending=[r for r in rows if not r.get('published_t80_comparison_json')]
    if not pending:return 0
    claims=publication_claims(main_pdf,supplementary_files)
    labels=[str(r.get('device_category') or '') for r in rows]
    count=0
    for row in pending:
        try:
            payload=row.get('stability_data_json') or '[]'
            if isinstance(payload,str) and payload.startswith('zlib64:'):
                import base64,zlib
                payload=zlib.decompress(base64.b64decode(payload[7:])).decode('utf-8')
            studies=json.loads(payload) if isinstance(payload,str) else payload
            if not isinstance(studies,list) or not studies:continue
        except (TypeError,ValueError):continue
        target=str(row.get('primary_stability_figure') or '')
        chosen=next((s for s in studies if str(s.get('figure'))==target),studies[0])
        # Do not mutate the original stored curve/fitting JSON even indirectly.
        s=dict(chosen)
        s.setdefault('pdf',Path(main_pdf).name)
        if not s.get('fit_summary'):
            s['fit_summary']={k:row.get(k) for k in ('selected_fit_t80_h','selected_fit_components')}
        label=str(row.get('device_category') or '')
        comparison=compare_published_t80(s,label,claims,[v for v in labels if v!=label])
        row['published_t80_comparison_json']=json.dumps(comparison,ensure_ascii=False)
        row['published_t80_comparison_status']=comparison.get('status')
        row['reported_t80_source_type']=comparison.get('evidence_type')
        if comparison.get('evidence_type')=='EXPLICIT_PUBLISHED_T80':
            row.setdefault('reported_t80_from_source_h',comparison.get('reported_time_h'))
            if row.get('published_t80_h') is None:row['published_t80_h']=comparison.get('reported_time_h')
        elif comparison.get('evidence_type')=='NEAR_80_RETENTION_ENDPOINT':
            row['reported_near_t80_endpoint_h']=comparison.get('reported_time_h')
        row['selected_fit_vs_source_t80_abs_error_h']=comparison.get('absolute_error_h')
        row['selected_fit_vs_source_t80_rel_error_pct']=comparison.get('relative_absolute_error_pct')
        row['selected_fit_vs_source_t80_discrepancy_flag']=comparison.get('fit_discrepancy_flag_20pct')
        row['digitized_near_t80_one_sided']=comparison.get('digitized_near_t80_one_sided')
        ret=row.get('final_retention_pct')
        try:
            ret=float(ret)
            row['final_retention_near_t80_one_sided']=bool(80.0<=ret<81.0)
        except (ValueError,TypeError):pass
        if comparison.get('status')=='VERIFIED_COMPARISON':count+=1
    return count
