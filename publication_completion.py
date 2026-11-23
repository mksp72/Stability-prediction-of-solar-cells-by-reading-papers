"""v18.3 additive publication-level completion.

Called ONCE after all extractors have finished ONE matched main-paper/SI bundle.
No changes to stability digitizers, fitting equations, or previously captured values.
Retains all source/normalization evidence. Unknown is never interpreted as zero.
"""
from __future__ import annotations
import json, math, re
from pathlib import Path
import fitz
import numpy as np

_NUM = r'(\d(?:\.\d{1,3})?)'
_GAP_RX = [
    re.compile(r'(?i)(?:optical\s+|electronic\s+)?band[ -]?gap\s+(?:of|is|was|at|to|~|≈|about|around|measured\s+to\s+be|determined\s+to\s+be)\s*'+_NUM+r'\s*eV'),
    re.compile(r'(?i)(\d(?:\.\d{1,3})?)\s*[-– ]*eV[-– ]*band[ -]?gap'),
    re.compile(r'(?i)(?:band[ -]?gap|gap\s+energy)\s*(?:\([^)]*\))?\s*(?::|=)\s*'+_NUM+r'\s*eV'),
]
_FORMULA_RX = re.compile(
    r'(?<![A-Za-z0-9])((?:Rb|Cs|FA|MA)[A-Za-z0-9().\[\]+-]{0,70}?Pb(?:\([^\s)]{2,25}\)|\[[^\s]]{2,25}\]|I\d*(?:\.\d+)?(?:Br\d*(?:\.\d+)?)?)(?:\d+)?)(?=\s|[,;:.\)]|$)', re.I)


def _number(v):
    try:
        x=float(v)
        return x if math.isfinite(x) else None
    except (ValueError,TypeError):return None


def _missing(v):return v is None or v=='' or (isinstance(v,float) and not math.isfinite(v))


def _set(row,key,val,origin=None):
    if _missing(row.get(key)) and not _missing(val):
        row[key]=val
        if origin and _missing(row.get(key+'_source')):row[key+'_source']=origin
        return True
    return False


def _study_list(row):
    payload=row.get('stability_data_json')
    if not payload:return []
    if isinstance(payload,list):return payload
    if isinstance(payload,str) and payload.startswith('zlib64:'):
        import base64,zlib
        payload=zlib.decompress(base64.b64decode(payload[7:])).decode('utf8')
    try:
        data=json.loads(payload)
        return data if isinstance(data,list) else []
    except (ValueError,TypeError):return []


def _protocol(study):return str(study.get('protocol') or '')


def _is_dark(s):return bool(re.search(r'ISOS[-_ ]?D|DARK|STORAGE|DAMP.?HEAT',_protocol(s),re.I))


def _is_light(s):return bool(re.search(r'ISOS[-_ ]?L|LIGHT|ILLUM|MPP',_protocol(s),re.I)) and not _is_dark(s)


def _primary(studies):
    return max(studies,key=lambda s:(int(_is_light(s)), int(str(s.get('measurement_mode','')).upper()=='MPP'),float(s.get('t_end_h') or 0))) if studies else None


def _interp(study, hour):
    try:
        x=np.asarray(study.get('time_h',[]),float); y=np.asarray(study.get('retention_pct',[]),float)
        ok=np.isfinite(x)&np.isfinite(y);x=x[ok];y=y[ok]
        if len(x)<2:return None
        idx=np.argsort(x);x=x[idx];y=y[idx]
        return round(float(np.interp(hour,x,y)),4) if x[0]<=hour<=x[-1] else None
    except (TypeError,ValueError):return None


def _t80(t,y):
    for i in range(1,len(t)):
        if y[i-1]>80>=y[i]:
            return round(float(t[i-1]+(80-y[i-1])*(t[i]-t[i-1])/(y[i]-y[i-1])),3)
    return None


def _complete_stability(row):
    studies=_study_list(row)
    if not studies:return
    best=_primary(studies)
    _set(row,'stability_experiment_count',len(studies))
    _set(row,'stability_figures','; '.join(dict.fromkeys(str(s.get('figure')) for s in studies if s.get('figure'))))
    if not best:return
    fields={
      'primary_stability_figure':'figure','primary_stability_protocol':'protocol',
      'primary_stability_conditions':'conditions','primary_measurement_mode':'measurement_mode',
      'primary_stability_metric':'metric','retention_last_pct':'retention_last_pct',
      'stability_last_time_h':'t_end_h','final_retention_pct':'retention_last_pct',
      'final_aging_time_h':'t_end_h','observed_t80_h':'observed_t80_h',
      'retention_reference':'normalization_source'
    }
    for dest,src in fields.items():_set(row,dest,best.get(src))
    for hr in (100,500,1000):_set(row,f'retention_{hr}h_pct',_interp(best,hr))
    times=best.get('time_h',[]);vals=best.get('retention_pct',[])
    if len(times)==len(vals) and times:
        try:
            idx=int(np.nanargmax(np.asarray(vals,dtype=float)))
            _set(row,'burn_in_peak_time_h',round(float(times[idx]),4))
            _set(row,'burn_in_peak_retention_pct',round(float(vals[idx]),4))
            _set(row,'observed_t80_h',_t80(times,vals))
        except (ValueError,TypeError):pass
    if _missing(row.get('t80_status')):
        row['t80_status']='OBSERVED' if not _missing(row.get('observed_t80_h')) else 'NOT_REACHED_DURING_MEASUREMENT'
    if _missing(row.get('stability_endpoint')) and not _missing(row.get('retention_last_pct')):
        row['stability_endpoint']=str((row.get('retention_last_pct'),row.get('stability_last_time_h')))
    dark=max((s for s in studies if _is_dark(s)),key=lambda s:float(s.get('t_end_h') or 0),default=None)
    if dark:
        for d,k in [('dark_stability_protocol','protocol'),('dark_retention_last_pct','retention_last_pct'),('dark_duration_h','t_end_h')]:_set(row,d,dark.get(k))
    # Copy fitted parameters out of existing per-trace JSON, without refitting.
    fs=best.get('fits') or []
    for fit in fs:
        n=fit.get('n_components')
        if n not in (2,3) or fit.get('status')!='ok':continue
        _set(row,f'fit{n}_rmse_pct',fit.get('rmse_post_peak_pct'))
        for i in range(1,n+1):_set(row,f'fit{n}_rate{i}_hinv',fit.get(f'rate_{i}_per_h'))
    summary=best.get('fit_summary') or {}
    for key,value in summary.items():_set(row,key,value)
    _set(row,'fit_traces_successful',sum(any(f.get('status')=='ok' for f in s.get('fits',[])) for s in studies))
    _set(row,'data_origin',best.get('figure_source_type') or 'DIGITIZED_PUBLISHED_STABILITY')
    if _missing(row.get('primary_stability_metric')):_set(row,'primary_stability_metric',row.get('stability_measurement'))
    # Reported T80 is independent of fitted T80 and only used for validation.
    _set(row,'published_t80_h',best.get('published_t80_h'))
    valid=[]
    for s in studies:
        t80=_number(s.get('published_t80_h'))
        if t80 is not None:
            valid.append({'figure':s.get('figure'),'reported_t80_h':t80,'fit2_t80_h':(s.get('fit_summary') or {}).get('fit2_t80_h'),'fit3_t80_h':(s.get('fit_summary') or {}).get('fit3_t80_h')})
    if valid:
        _set(row,'reported_t80_validation_json',json.dumps(valid,ensure_ascii=False))
        _set(row,'reported_t80_validation_count',len(valid))
        for v in valid:
            if v.get('figure')==best.get('figure'):
                errors={n:abs(v['reported_t80_h']-v[f'fit{n}_t80_h']) for n in (2,3) if _number(v.get(f'fit{n}_t80_h')) is not None}
                if errors:_set(row,'closest_model_to_published_t80',min(errors,key=errors.get))


def _journal_from_doi(doi):
    """Journal is taken from source cover and metadata, never guessed from DOI."""
    return None


def _year_from_cover(cover):
    patterns=[r'(?i)Published\s+(?:online\s*)?(?:on\s*)?:?\s*(?:\d{1,2}\s+)?[A-Za-z]+\s+(20\d{2})',
              r'(?i)(?:Published|Copyright|©|\(c\))\s*(?:by\s+)?(?:\S+\s+){0,4}?(20\d{2})',
              r'(?i)Science\s+\d+\s*,\s*\d+\s*\((20\d{2})\)',
              r'(?i)©\s*(20\d{2})',
              r'(?i)Document\s+status\s+and\s+date\s*:\s*Published\s*:\s*\d{1,2}/\d{1,2}/(20\d{2})',
              r'(?i)\bNat(?:ure)?\s+Energy\s+\d+\s*,\s*\d+[-–]\d+\s*\((20\d{2})\)',
              r'(?i)\((20\d{2})\)\.\s*[A-Z][a-z].{0,70}\bperovskite',
              r'(?i)\b(?:publisher.s PDF|version of record).{0,100}?Published.{0,30}?(20\d{2})']
    for rx in patterns:
        m=re.search(rx,cover[:13000]);
        if m:return int(m.group(1))
    return None


def _evidence(main_pdf,supplements):
    docs=[]
    for f in [main_pdf,*supplements]:
        f=Path(f)
        if f.suffix.lower()!='.pdf':continue
        try:
            with fitz.open(str(f)) as doc:
                pages=[(i+1,' '.join(doc[i].get_text('text').split())) for i in range(len(doc))]
                docs.append((f.name,pages))
        except Exception:pass
    return docs


def _formula_candidates(docs):
    candidates=[]
    for file,pages in docs:
        for p,txt in pages:
            for m in _FORMULA_RX.finditer(txt):
                f=re.sub(r'\s+','',m.group(1)).strip('().,;:')
                if len(f)<8 or len(f)>95 or not re.search(r'(?:I\d|Br\d|\[[^]]+\]|\([^)]*[IBr][^)]*\))',f):continue
                around=txt[max(0,m.start()-170):m.end()+110]
                score=1
                if re.search(r'(?i)(?:absorber|composition|fabricat|precursor|perovskite\s+film|PSC|device)',around):score+=2
                if re.search(r'(?i)(?:the\s+composition\s+(?:is|of)|light\s+absorber|baseline\s+perovskite|perovskite\s+recipe)',around):score+=3
                if re.search(r'(?i)\b(?:DFT|simulation|simulated|modelling|modelled)\b',around):score-=3
                if p>12:score-=2
                if score>=3:candidates.append((score,f,file,p))
    # Also recognize author-reported mixed-halide 3D recipes written with
    # symbolic mole fractions, e.g. (FAPbI3)x(MAPbBr3)1-x.
    symbolic=re.compile(r'\(FA\s*PbI3\)\s*x\s*\(MA\s*PbBr3\)\s*1\s*[-−]\s*x',re.I)
    for file,pages in docs:
        for p,txt in pages[:12]:
            for m in symbolic.finditer(txt):
                candidates.append((7,re.sub(r'\s+','',m.group(0)),file,p))
    # The actual absorber composition is usually declared in first pages;
    # candidates in references or theoretical surrogate sections are penalized.
    return sorted(candidates,key=lambda c:(-c[0],c[3]))


def _category_formula(row,candidates):
    if _number(row.get('absorber_formula')) is not None or not _missing(row.get('absorber_formula')):return
    if not candidates:return
    cat=str(row.get('device_category','')).lower()
    candidates=[c for c in candidates if not re.search(r'(?i)\b(?:PbI2|PbBr2)\b',c[1])]
    # Pure-2D controls cannot inherit common 3D composition.
    if 'pp-2d' in cat or 'pure-2d' in cat:
        specialized=[c for c in candidates if 'BA2' in c[1]]
        if not specialized:return
        candidates=specialized
    # Choose only when one well-supported formulation dominates evidence.
    scores={}
    for score,f,file,page in candidates:
        key=re.sub(r'[^a-z0-9.]','',f.lower())
        scores[key]=scores.get(key,0)+max(1,score)
    order=sorted(scores,key=scores.get,reverse=True)
    if not order:return
    # No chemistry guessing when two different strong formulas coexist.
    if len(order)>1 and scores[order[0]] < 1.6*scores[order[1]]:return
    for score,f,file,page in candidates:
        if re.sub(r'[^a-z0-9.]','',f.lower())==order[0]:
            _set(row,'absorber_formula',f,f'{file} p.{page}, formulation from article text')
            return


def _category_gap(row,docs):
    if _number(row.get('bandgap_ev')) is not None:return
    cat=str(row.get('device_category','')).lower()
    if 'pp-2d' in cat or 'pure-2d' in cat:return
    matches=[]
    for file,pages in docs:
        for page,txt in pages:
            # Experimental optical/EQE gap, not a modelled HTL/ETL gap.
            for m in re.finditer(r'(?i)\b(?:optical\s+|EQE[- ]?derived\s+)?band[ -]?gap\b.{0,115}?(?<![0-9])([1-2]\.\d{1,3})\s*eV',txt):
                v=_number(m.group(1))
                if v is None or not 1.1<=v<=2.8:continue
                span=txt[max(0,m.start()-70):m.end()+35]
                if re.search(r'(?i)\b(?:DFT|computed|theoretical|simulat|HTL|ETL|transport\s+layer|barrier)\b',span):continue
                matches.append((v,file,page,span))
            for rx in _GAP_RX:
                for m in rx.finditer(txt):
                    v=_number(m.group(1))
                    if v is None or not 1.1<=v<=2.8:continue
                    span=txt[max(0,m.start()-70):m.end()+35]
                    if re.search(r'(?i)\b(?:DFT|computed|theoretical|simulat|HTL|ETL|transport\s+layer|barrier)\b',span):continue
                    matches.append((v,file,page,span))
    if not matches:return
    vals={round(x[0],2) for x in matches}
    # One dominant gap in the publication => shared 3D absorber.
    from collections import Counter
    counts=Counter(round(x[0],2) for x in matches)
    v,count=counts.most_common(1)[0]
    if len(vals)>1 and count<2*counts.most_common(2)[1][1]:return
    for value,file,page,_ in matches:
        if round(value,2)==v:
            _set(row,'bandgap_ev',round(value,4))
            _set(row,'bandgap_source',f'{file} p.{page}, experimental optical/EQE gap in text')
            break


def finalize_publication(rows,main_pdf,supplements=(),metadata=None):
    """Complete each category once, before writing or advancing to next PDF."""
    main_pdf=Path(main_pdf);docs=_evidence(main_pdf,supplements)
    cover=docs[0][1][0][1] if docs and docs[0][1] else ''
    supplementary_covers='\n'.join(pages[0][1] for name,pages in docs[1:] if pages)
    journal=_journal_from_doi((metadata or {}).get('doi') or (rows[0].get('doi') if rows else ''))
    year=_year_from_cover(cover) or _year_from_cover(supplementary_covers)
    candidates=_formula_candidates(docs)
    for row in rows:
        _complete_stability(row)
        _set(row,'supplementary_pdf','; '.join(Path(s).name for s in supplements if Path(s).suffix.lower()=='.pdf'))
        _set(row,'supplementary_evidence','Main article and supplementary documents processed together; matching details in supplement_match_evidence' if supplements else None)
        # Correct a clearly wrong inherited journal from metadata parser.
        if journal and row.get('journal_name')!=journal:
            row['journal_name']=journal
            row['journal_source']='DOI journal-family identification'
        if year and row.get('article_year')!=year:
            row['article_year']=year
            row['article_year_source']='Main article cover / publication line'
        _category_formula(row,candidates)
        _category_gap(row,docs)
        f=str(row.get('absorber_formula') or '').lower()
        if f:
            if _missing(row.get('mixed_cation')):
                row['mixed_cation']=bool(sum(tok in f for tok in ('cs','fa','ma','rb'))>1)
            if _missing(row.get('mixed_halide')):
                row['mixed_halide']='br' in f and 'i' in f
        # Preserve stability-test initial PCE separately from representative J-V PCE.
        # FF validation MUST use J-V PCE when available; otherwise the values
        # may come from different devices within the same category (e.g. champion
        # JV vs initial stability group) and produce an artificial discrepancy.
        p_initial,v,j=(_number(row.get(k)) for k in ('pce0','voc0','jsc0'))
        p_jv=_number(row.get('pce_jv_pct'))
        p=p_jv if p_jv is not None else p_initial
        basis=('JV_PCE_VOC_JSC_100mWcm2' if p_jv is not None
               else 'CATEGORY_PCE_VOC_JSC_100mWcm2_NOT_SAME_CELL')
        if _number(row.get('ff0')) is None and None not in (p,v,j) and v*j>0:
            ff=100*p/(v*j)
            if 5<=ff<=100:_set(row,'ff0',round(ff,4),
                'Calculated using '+basis)
        if None not in (p,v,j) and v*j>0:
            _set(row,'ff_calculated_pct',round(100*p/(v*j),5))
            _set(row,'ff_calculation_basis',basis)
            _set(row,'ff_calculation_status',
                 'JV_METRIC_CROSSCHECK' if p_jv is not None else 'CATEGORY_APPROXIMATION')
            ff=_number(row.get('ff0'))
            if ff is not None:_set(row,'ff_reported_vs_calculated_diff_pctpts',round(ff-100*p/(v*j),4))
        jv,eqe=(_number(row.get(k)) for k in ('jsc0','jsc_eqe_integrated'))
        if jv is not None and eqe is not None and jv>0:
            for k,v in [('jsc_difference_signed_pct',100*(jv-eqe)/jv),('jsc_difference_absolute_pct',100*abs(jv-eqe)/jv),('jsc_eqe_to_jv_ratio',eqe/jv)]:
                _set(row,k,round(v,5))
        # Recover a category-specific literature page number from evidence strings.
        if _missing(row.get('metric_page')):
            for src_key in ('pce0_source','voc0_source','jsc0_source','metric_source'):
                source=str(row.get(src_key) or '')
                mm=re.search(r'\bp\.?\s*(\d+)\b',source,re.I)
                if mm:_set(row,'metric_page',int(mm.group(1)));break
        if _missing(row.get('jv_figure_assignment')) and not _missing(row.get('metric_source')):
            _set(row,'jv_figure_assignment',str(row['metric_source']))
        if _missing(row.get('metric_source')):
            _set(row,'metric_source',row.get('pce0_source') or row.get('voc0_source') or row.get('jsc0_source'))
        if _missing(row.get('jsc_eqe_integrated_source')):
            _set(row,'jsc_eqe_integrated_source',row.get('jsc_eqe_source') or row.get('jsc_eqe_integrated_source'))
        if _missing(row.get('data_origin')):_set(row,'data_origin','PUBLISHED_STABILITY_DIGITIZED')
    return rows
