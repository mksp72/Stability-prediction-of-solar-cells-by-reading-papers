"""Figure-panel metadata inference. No DOI, figure-number, or document allowlists.

The returned data are evidence for validation rather than a license to invent a
sample label. Ambiguous legends stay unassigned.
"""
from __future__ import annotations
import re

PANEL_TOKEN = re.compile(r'^\(?([a-z])\)$', re.I)

def identify_subplot(plot, words):
    x0,y0,y1=plot['x0'],plot['y0'],plot['y1']
    choices=[]
    for w in words:
        token=w[4].strip()
        match=PANEL_TOKEN.fullmatch(token)
        if not match:continue
        x,y=(w[0]+w[2])*.5,(w[1]+w[3])*.5
        if x0-105 < x < x0-6 and y0-105 < y < y0+20:
            distance=(abs((x0-40)-x)+.45*abs((y0-4)-y))
            choices.append((distance,match.group(1).lower()))
    return min(choices)[1] if choices else None


def _caption_block(text):
    # Prefer a real figure caption, not mentions of a figure in the body.
    matches=list(re.finditer(r'(?m)^(?:Fig\.?|Figure)\s*\d+\s*[|.:]',text,re.I))
    return text[matches[-1].start():matches[-1].start()+1700] if matches else ''


def panel_context(page, plot, ocr_words=None):
    words=ocr_words if ocr_words is not None else page.get_text('words')
    panel=identify_subplot(plot,words)
    caption=_caption_block(page.get_text('text'))
    segment=''
    if panel and caption:
        # Caption syntax commonly uses either '(b)' or 'b Normalized ...'.
        pattern=re.compile(r'(?<!\w)(?:\(([a-z])\)|([a-z])\s+(?=[A-Z]))\s*',re.I)
        markers=[(m.start(),m.end(),(m.group(1) or m.group(2)).lower())
                 for m in pattern.finditer(caption)]
        for i,(a,b,letter) in enumerate(markers):
            if letter==panel:
                boundary=markers[i+1][0] if i+1<len(markers) else min(len(caption),b+360)
                segment=caption[b:boundary].strip()
                break
    text=(segment or caption[:450])
    # Optical and JV axes must not be accepted just because a caption nearby
    # mentions a stability experiment in another panel.
    panel_is_non_stability=bool(segment and re.search(
        r'\b(?:V\s*OC|V_OC|J\s*SC|J_SC|current density|J[–-]?V curves?|histogram|distribution|XRD|spectra)\b',
        segment,re.I) and not re.search(r'normalized power|retention|MPP tracking|stability curve',segment,re.I))
    protocol=None
    if segment:
        if re.search(r'ISOS[- ]?D|dark|stored|storage|glovebox',segment,re.I):
            protocol='ISOS-D_OR_DARK_AGING'
        elif re.search(r'ISOS[- ]?L|MPP|illumination|light.soak|continuous operation',segment,re.I):
            protocol='ISOS-L_OR_LIGHT_AGING'
    sample_architecture=None
    if re.search(r'\b(?:single[ -]?junction|single[ -]?cell|NBG single[ -]?junction)\b',caption,re.I):
        sample_architecture='single_junction'
    elif re.search(r'\b(?:tandem|PTSCs|2\s*T\s+PTSC)\b',caption,re.I):
        sample_architecture='tandem'
    return dict(subplot_label=panel,panel_caption=segment[:420] if segment else None,
                panel_protocol=protocol,panel_non_stability=panel_is_non_stability,
                plot_architecture=sample_architecture)
