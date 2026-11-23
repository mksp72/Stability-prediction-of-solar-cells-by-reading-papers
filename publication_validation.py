"""v18.5 additive, evidence-preserving per-publication validation.

No changes to curve points, decay fits, or T80. This module is run after
all v18.4 extraction/enrichment for a publication, before its rows are saved.
"""
from __future__ import annotations
import re, math
from collections import defaultdict
from pathlib import Path
from curve_integrity_guard import suspicious_legend

BAD_LABEL = re.compile(r'^[\W_]+$|^(?:\d+(?:[.,]\d+)?\s*(?:h|hour|hours|nm|%|v|ma|cm2)?|time|retention|pce|jsc|voc|ff|axis|legend|figure|fig\.?|nature|science|experimental|simulated|fitted|fit|sample|device|unknown|unassigned|not assigned)$', re.I)
MAIN_CHEMICAL_KEYS = ('control','reference','pristine','target','treated','passivat','modified','3d','2d')


def invalid_category(label):
    s = re.sub(r'\s+', ' ', str(label or '')).strip()
    if len(s) < 2 or len(s)>110 or BAD_LABEL.fullmatch(s) or suspicious_legend(s):return True
    if s.lower().startswith(('unassigned_', 'unknown_', 'axis_')):return True
    if s in ('|','-','—','/','I','l'):return True
    # Obvious isolated tick axis instead of a chemical or control label.
    if re.fullmatch(r'[\d.,\s%°]+',s):return True
    if re.search(r'\b(?:time\s*\(?h\)?|wavelength\s*\(?nm\)?|normalized\s*pce|relative\s*humidity)\b',s,re.I):return True
    return False


def review_article(document):
    """Conservatively skip *explicitly identified* review articles, never by 'review' in citations.

    Requires 'review' (or 'perspective', 'mini-review') as a front-matter
    article-type label; 'reviewed' in body or 'literature review' alone is not enough.
    """
    title=str(getattr(document,'title','') or '')
    header=str(getattr(document,'header','') or '')
    top='\n'.join(header.splitlines()[:32])
    title_review=bool(re.match(r'^\s*(?:a\s+)?(?:systematic\s+)?(?:review|mini[ -]?review|perspective)\s*[:—-]',title,re.I) or
                      re.search(r'\b(?:a review|review article|systematic review|meta-analysis)\b',title,re.I))
    tag_review=bool(re.search(r'(?im)^\s*(?:article\s*type\s*[:\-]\s*)?(?:review\s+article|reviews?|mini[ -]?review|perspectives?)\s*$',top))
    if not (title_review or tag_review):return False
    # A review may contain an original experiment; require more explicit evidence
    # before skipping it. A review with 'we fabricated/measured' goes to extractor.
    claims_original=bool(re.search(r'\bwe\s+(?:fabricated|deposited|prepared|synthesi[sz]ed|measured|recorded|constructed)\b',top,re.I))
    return not claims_original


def _number(v):
    try:
        if v is None or str(v).strip() in ('','nan','None'):return None
        x=float(v)
        return x if math.isfinite(x) else None
    except (ValueError,TypeError):return None


def _canonical_alias(label):
    s=re.sub(r'\s+', ' ',str(label or '')).strip()
    if re.fullmatch(r'(?i)(ctrl|ctrl\.|ref|ref\.|reference|pristine)',s):return 'Control'
    # Strip sampling-statistic suffixes ONLY; keep distinct Target-1/2 and
    # different passivants or concentrations as distinct formulations.
    s=re.sub(r'\s*\((?:champion|champ\.?|average|avg\.?)\)\s*$','',s,flags=re.I).strip()
    return s


def _representative_pv_check(row):
    """Audit same-scan J-V closure; never equate stability-initial PCE to champion J-V PCE."""
    p0=_number(row.get('pce0'));pj=_number(row.get('pce_jv_pct'))
    v=_number(row.get('voc0'));j=_number(row.get('jsc0'));ff=_number(row.get('ff0'))
    row['jv_pce_from_voc_jsc_ff_pct']=round(v*j*ff/100,5) if None not in (v,j,ff) else None
    comparison=pj
    # If both PV and stability PCE are the same provenance use PCE0 only when
    # the source explicitly states a J-V table/curve, not an aging starting PCE.
    if comparison is None and p0 is not None:
        origin=' '.join(str(row.get(k) or '') for k in ('metric_source','pce0_source','pce_source')).lower()
        if re.search(r'j[–-]v|current.voltage|photovoltaic\s+table|table\s+s?\d+',origin) and not re.search(r'initial\s+(?:stability|aging)|stability\s+initial',origin):
            comparison=p0
    if comparison is None:
        row['jv_closure_status']='NO_MATCHED_JV_PCE_FOR_COMPARISON' if row['jv_pce_from_voc_jsc_ff_pct'] is not None else 'INCOMPLETE_JV_COMPONENTS'
        row['jv_closure_error_pctpts']=None
        return
    if row['jv_pce_from_voc_jsc_ff_pct'] is None:
        row['jv_closure_status']='INCOMPLETE_JV_COMPONENTS';row['jv_closure_error_pctpts']=None;return
    gap=abs(row['jv_pce_from_voc_jsc_ff_pct']-comparison)
    row['jv_closure_error_pctpts']=round(gap,5)
    row['jv_closure_status']='CONSISTENT_WITH_ROUNDING' if gap<=0.6 else 'POSSIBLE_MISMATCH_DIFFERENT_JV_DEVICES'


def _protocol_check(row):
    proto=str(row.get('primary_stability_protocol') or '')
    mode=str(row.get('primary_measurement_mode') or '')
    conditions=str(row.get('primary_stability_conditions') or '')
    field=' '.join((proto,conditions)).lower()
    dark=bool(re.search(r'isos[– -]?d\b|\bdark\b|dark.storage|shelf',field))
    light=bool(re.search(r'isos[– -]?l\b|light.?soak|illumina|one.?sun|\bmpp\b|maximum.power',field))
    if dark and light:
        row['protocol_validation_status']='CONFLICT_DARK_AND_LIGHT_EVIDENCE'
    elif dark and mode.upper() in ('MPP','MPPT','MAXIMUM POWER POINT'):
        row['protocol_validation_status']='CONFLICT_DARK_AGING_WITH_MPP_LABEL'
        # Do not change measured curve's mode or protocol: source-based evidence
        # needed. Consumers should treat this feature as ambiguous.
    else:
        row['protocol_validation_status']='NO_CONTRADICTION_DETECTED' if proto else 'PROTOCOL_UNSPECIFIED'


def validate_publication(rows, article=None):
    """Return filtered rows and diagnostic counts. Retains all original fitting data."""
    surviving=[];stats={'invalid_category_removed':0,'category_alias_normalized':0,
                         'jv_mismatch_flags':0,'protocol_conflict_flags':0,'duplicate_records_removed':0}
    seen=set()
    for row in rows:
        original=str(row.get('device_category') or '')
        if invalid_category(original):
            stats['invalid_category_removed']+=1
            continue
        new=_canonical_alias(original)
        if new!=original:
            row['device_category_original']=original
            row['device_category']=new
            stats['category_alias_normalized']+=1
        # Keep separate concentrations, absorber stacks, and numbered targets.
        # Only de-duplicate exact repetitions of identical same-paper category
        # and identical stable curve payload; don't merge non-identical studies.
        identity=(str(row.get('doi') or ''),str(row.get('category_id') or new),
                  str(row.get('stability_data_json') or ''))
        if identity in seen:
            stats['duplicate_records_removed']+=1
            continue
        seen.add(identity)
        _representative_pv_check(row)
        _protocol_check(row)
        if row['jv_closure_status']=='POSSIBLE_MISMATCH_DIFFERENT_JV_DEVICES':stats['jv_mismatch_flags']+=1
        if row['protocol_validation_status'].startswith('CONFLICT_'):stats['protocol_conflict_flags']+=1
        surviving.append(row)
    return surviving,stats
