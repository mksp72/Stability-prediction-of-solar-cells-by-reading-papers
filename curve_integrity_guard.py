"""Publication-local curve identity checks for the v18.6 general PDF reader.

Only removes demonstrably redundant/suspect automatic traces. This is a
validation-layer change: numerical traces, existing fits, AIC/RMSE, and the
root solver are never modified. Excluded evidence is retained as diagnostics.
No DOI, figure, filename, condition or fixed curve values are configured.
"""
from __future__ import annotations
import math
import re
from collections import defaultdict

PALE = {'gray'}
# Failed OCR fragments that are not chemical names or sample formulations.
BAD_FRAGMENTS = re.compile(
    r'^(?:[a-z]{1,3}\s*[,;:]\s*(?:of|for|and)\b|(?:of|for|and)\s+|[+,:;/|]\s*)', re.I
)


def suspicious_legend(label: str) -> bool:
    s = ' '.join(str(label or '').split())
    if not s or BAD_FRAGMENTS.match(s):
        return True
    if len(s) > 90 or re.match(r'(?i)^\s*(?:normalized|retention|current|stability|time|power)\b', s):
        return True
    return False


def _quality_score(study):
    """Order by independent identification and trace fidelity; not duration."""
    f = set(study.get('digitization_quality_flags') or [])
    color = str(study.get('color') or '').lower()
    legend = str(study.get('legend_evidence') or '')
    # A named, directly matched non-gray curve is stronger than a pale-gray
    # overlapping mask. Never use the longest span as a proxy for accuracy.
    return (0 if suspicious_legend(legend) else 1,
            0 if 'CATEGORY_LABEL_BORROWED_FROM_SIBLING_PANEL' in f else 1,
            0 if color in PALE else 1,
            0 if 'HOUR_UNIT_INFERRED_FROM_CAPTION_REQUIRES_QC' in f else 1,
            0 if 'PALE_GRAY_TRACE_REQUIRES_VISUAL_QC' in f else 1,
            len(study.get('time_h') or []))


def _span_overlap(a, b):
    a0, a1 = float(a['t_start_h']),float(a['t_end_h'])
    b0, b1 = float(b['t_start_h']),float(b['t_end_h'])
    shared = max(0., min(a1,b1)-max(a0,b0))
    short = min(max(0,a1-a0),max(0,b1-b0))
    return shared / short if short else 0


def prune_duplicate_color_traces(studies):
    """Dedupe competing colors assigned to ONE label in the SAME graph.

    Only suppress lower-confidence traces when there is a stronger trace of
    the exact same category and figure covering >=50% of the shorter span.
    Differently labelled samples, distinct plots, and replicate stability
    experiments are preserved.
    """
    grouped=defaultdict(list)
    for idx, study in enumerate(studies):
        graph=(study.get('pdf'),study.get('page'),study.get('figure'))
        grouped[graph].append((idx,study))
    dropped={}
    for graph, members in grouped.items():
        for i,study in members:
            for j,other in members:
                if i==j or study.get('color')==other.get('color'):
                    continue
                if _span_overlap(study,other)<0.50:
                    continue
                if study.get('protocol')!=other.get('protocol'):
                    continue
                # Only suppress a pale-gray candidate if another color is
                # independently labelled as this category on the SAME plot.
                if study.get('color') not in PALE or other.get('color') in PALE:
                    continue
                if _quality_score(other)>_quality_score(study):
                    dropped[i]={'status':'REJECTED_DUPLICATE_PALE_COLOR',
                                'reason':f'Gray mask competes with stronger {other.get("color")} curve on the same graph',
                                'figure':study.get('figure'),'source':study.get('pdf'),
                                'candidate_color':study.get('color'),'preferred_color':other.get('color')}
                    break
    return [s for i,s in enumerate(studies) if i not in dropped],list(dropped.values())


def choose_primary(studies):
    """Operational preference, reliable graph identity, then measurement span."""
    assert studies
    def rank(s):
        is_light=str(s.get('protocol','')).startswith('ISOS-L')
        # Quality before duration prevents a long axis contamination from
        # overriding a shorter but reliably identified degradation trace.
        return (int(is_light),*_quality_score(s),float(s.get('t_end_h') or 0))
    return max(studies,key=rank)


def trajectory_warnings(study):
    f=set(study.get('digitization_quality_flags') or [])
    points=study.get('retention_pct') or []
    start=float(study.get('t_start_h') or 0)
    end=float(study.get('t_end_h') or 0)
    if end>0 and start/end > 0.10:
        f.add('CURVE_START_MISSING_GT_10PCT_MEASUREMENT_WINDOW')
    if points:
        try:
            vals=[float(p) for p in points]
            if not all(math.isfinite(v) for v in vals):f.add('NONFINITE_RETENTION')
            if min(vals)<-1 or max(vals)>140:f.add('IMPLAUSIBLE_RETENTION_RANGE')
            if len(vals)>=3:
                sorted_d=sorted(abs(vals[i]-vals[i-1]) for i in range(1,len(vals)))
                if sorted_d[-1]>25 and sorted_d[-1]>5*sorted_d[len(sorted_d)//2]:
                    f.add('SUSPECTED_GRAPH_ANNOTATION_JUMP')
        except (TypeError,ValueError):
            f.add('UNREADABLE_RETENTION_VALUES')
    return sorted(f)


def review_only_label(label,studies):
    """Reject certainly corrupted identities, or ambiguous borrowed gray only."""
    if suspicious_legend(label):return 'MALFORMED_DEVICE_LABEL'
    if not studies:return 'NO_NUMERICAL_TRACES'
    if all(s.get('color') in PALE and
           'CATEGORY_LABEL_BORROWED_FROM_SIBLING_PANEL' in
           set(s.get('digitization_quality_flags') or []) for s in studies):
        return 'ONLY_BORROWED_PALE_GRAY_LABEL'
    return None
