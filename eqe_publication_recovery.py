"""Publication-agnostic EQE/IPCE integrated Jsc recovery.

Read the actual main PDF and all matched PDF supplements, selecting only
values accompanied by EQE/IPCE + *integrated* current evidence.  Assign
only to explicit device categories, never to the closest J-V current or
some unrelated stability category.  This is additive to the v18.1 pipeline:
no DOI, filename, figure number, or known device current is embedded here.

Sources: explicit one- and two-category numerical prose; printed endpoint
annotation on a single-category IPCE figure (conservative); paper/figure
provenance. Curves without numeric endpoint labels need calibrated digitization
and are deliberately left unfilled by this text/annotation stage.
"""
from __future__ import annotations
from pathlib import Path
import math, re, unicodedata
import fitz

FLOAT = r'(?<![\w.])(?:\d{1,2}\.\d{1,3}|\d{1,2})(?![\w.])'
VALUE = re.compile(FLOAT)
CONTEXT = re.compile(r'(?i)(?:integrat(?:ed|ion|ing)|cumulative).{0,90}(?:J\s*[_-]?\s*SC|photocurrent|current(?:\s*density)?)|(?:J\s*[_-]?\s*SC|photocurrent).{0,90}integrat(?:ed|ion|ing)')
EQE = re.compile(r'(?i)\b(?:EQE|IPCE|external\s+quantum\s+efficiency|incident\s+photon)\b')
UNIT = re.compile(r'(?i)(?:mA\s*(?:cm|/cm)|mA\s*cm|mAcm|mA\s*×\s*cm)')
RESPECTIVELY = re.compile(r'(?i)respectively')

def _clean(s):
    s=unicodedata.normalize('NFKC',s).replace('\u00ad','').replace('\u2212','-').replace('−','-')
    s=re.sub(r'(?<=\w)-\s+(?=\w)', '-', s)  # wrapped chemical names
    return re.sub(r'\s+',' ',s).strip()

def _num(v):
    try:
        val=float(v)
        return val if math.isfinite(val) and 4 <= val <= 32 else None
    except (ValueError,TypeError):return None

def _key(s):return re.sub('[^a-z0-9]','',str(s).lower())

def _aliases(row):
    label=str(row.get('device_category') or '')
    cat=str(row.get('category_id') or '')
    choices=[label,re.sub(r'\s+(?:stability\s+group|champion|average)$','',label,flags=re.I)]
    # Generic formulation categories.  A treatment code such as 'PSP' stays
    # distinct from the control, and numerically different targets do not merge.
    canon=_key(label)
    if canon in ('control','reference','ref','pristine') or cat=='control':
        choices += ['control','reference','ref.','pristine']
    if canon in ('target','treated','modified'):
        choices += ['target','treated','modified']
    # For 'Target-1/2', never infer identity from the word 'target' alone.
    return sorted({c.strip() for c in choices if len(_key(c))>=3},key=len,reverse=True)

def _mentions(text, alias):
    # Literal tokens not arbitrary prefix substrings ('target' != 'target-1').
    escaped=re.escape(alias.strip().rstrip('.'))
    return bool(re.search(r'(?<![\w])'+escaped+r'(?![\w\-\d])',text,re.I))

def _matching(rows,phrase):
    hits=[]
    for i,r in enumerate(rows):
        if any(_mentions(phrase,a) for a in _aliases(r)):
            hits.append(i)
    return hits

def _assign(rows,idx,val,source,kind,counts):
    v=_num(val)
    if v is None:return False
    row=rows[idx]
    present=_num(row.get('jsc_eqe_integrated'))
    if present is not None:
        # Recovered value must not silently overwrite existing literature data.
        if abs(present-v)>1.25:
            row.setdefault('eqe_recovery_conflict',f'existing {present} vs {v} ({source})')
        return False
    row['jsc_eqe_integrated']=round(v,4)
    row['jsc_eqe_source']=source
    row['jsc_eqe_extraction_method']=kind
    jv=_num(row.get('jsc0'))
    if jv is not None:
        sign=round(100*(jv-v)/jv,5)
        absolute=round(abs(sign),5)
        row['jsc_difference_signed_pct']=sign
        row['jsc_difference_absolute_pct']=absolute
        row['jsc_discrepancy_signed_pct']=sign
        row['jsc_discrepancy_absolute_pct']=absolute
        row['jsc_eqe_to_jv_ratio']=round(v/jv,6)
    counts[0]+=1
    return True

def _sentence_regions(s):
    # Paragraphs have references and superscripts; retain a broad overlapping
    # range around integrated-current anchors to capture multi-line 'respectively'.
    for m in CONTEXT.finditer(s):
        start=max(0,m.start()-75);end=min(len(s),m.end()+420)
        yield s[start:end],m.start()-start

def _units_near(s, span):
    return bool(UNIT.search(s[max(0,span[0]-20):min(len(s),span[1]+45)]))

def _parse_pair(region,rows,source,counts):
    """Two currents and two named categories, explicitly linked by 'respectively'."""
    if not RESPECTIVELY.search(region) or not EQE.search(region):return
    # Keep association within one local clause; look for (num ... and num ...
    # for name and name, respectively).  Numeric current candidates require
    # unit near each number or the second, to avoid axis ticks and PCE pairs.
    candidates=[]
    rx=re.compile(r'(?P<v1>\d{1,2}\.\d{1,3})\s*(?:mA\s*(?:/|×|\s)?\s*cm\s*[-−]?[23])?\s*(?:,|and|&)\s*(?P<v2>\d{1,2}\.\d{1,3})\s*(?:mA\s*(?:/|×|\s)?\s*cm\s*[-−]?[23])?\s+for\s+(?P<names>.{4,120}?)\s*,?\s*respectively',re.I)
    for m in rx.finditer(region):
        names=m.group('names')
        # 'CsFAMA and CsFAMA-F PSCs' or 'the control and the target devices'
        sep=re.split(r'\s+(?:and|&)\s+',names)
        if len(sep)!=2:continue
        l,r=sep
        l=re.sub(r'(?i)^\s*(?:the|both|devices?\s+with)\s+','',l).strip()
        r=re.sub(r'(?i)\s+(?:PSCs?|cells?|devices?|solar\s+cells?)\s*$','',r).strip()
        r=re.sub(r'(?i)^\s*(?:the|both)\s+','',r).strip()
        a=_matching(rows,l);b=_matching(rows,r)
        if len(a)==1 and len(b)==1 and a[0]!=b[0]:
            candidates.append((a[0],m.group('v1')));candidates.append((b[0],m.group('v2')))
    for idx,val in candidates:_assign(rows,idx,val,source,'explicit_paired_prose',counts)

def _parse_single(region,rows,source,counts):
    if not (EQE.search(region) or re.search(r'(?i)integrat(?:ed|ion|ing)\s+J\s*[_-]?\s*SC',region)):return
    # Pattern: "integrated Jsc of 24.88 mA cm-2 for the target device
    # and 23.92 mA cm-2 for the control device".  Each value is tied
    # to an explicit following category, including later currents in clause.
    anchors=list(CONTEXT.finditer(region))
    for a in anchors:
        tail=region[a.start():min(len(region),a.end()+245)]
        pairs=list(re.finditer(
            r'(?P<value>\d{1,2}\.\d{1,3})\s*'
            r'(?:mA\s*(?:/|×|\s)?\s*cm\s*[-−]?[23])'
            r'\s*(?:for|from)\s+(?:the\s+)?'
            r'(?P<label>[a-zA-Z][\w/-]{2,35})(?:\s+(?:devices?|cells?|PSCs?))?',
            tail,re.I))
        for mm in pairs:
            v=_num(mm.group('value'))
            if v is None:continue
            label=mm.group('label')
            hits=_matching(rows,label)
            if not hits and _key(label) in ('target','treated','modified'):
                # Safe only when there is exactly one non-reference category
                # in this paper's stability set.  Never merge Target-1/2.
                noncontrol=[i for i,r in enumerate(rows)
                            if not any(_key(r.get('device_category'))==a for a in
                                ('control','reference','ref','pristine'))]
                if len(noncontrol)==1:hits=noncontrol
            if len(hits)==1:
                _assign(rows,hits[0],v,source,'explicit_labeled_prose',counts)
        # Form: "the control device showed an integrated Jsc of 23.56 mA cm-2"
        after=tail
        mm=re.search(r'integrated\s+(?:J\s*SC|photocurrent|current)\s+of\s+'
                     r'(\d{1,2}\.\d{1,3})\s*mA',after,re.I)
        if mm:
            before=region[max(0,a.start()-125):a.start()]
            labels=re.findall(r'(?i)([\w/-]{3,35})\s+(?:device|cell|PSC)\s+(?:showed|has|yielded|exhibited|gave)',before)
            for label in labels:
                hits=_matching(rows,label)
                if len(hits)==1:
                    _assign(rows,hits[0],mm.group(1),source,'explicit_labeled_prose',counts)


def _figure_annotation(page,rows,source,counts):
    """Exact printed numeric integrated-Jsc endpoint in labeled EQE/IPCE panel.

    This is NOT estimating a value from the axis maximum. The figure must print
    the endpoint number, and the same figure caption must clearly name exactly
    one target formulation. Multiple categories => no guessing.
    """
    text=page.get_text('text')
    if not EQE.search(text) or not CONTEXT.search(text):return
    # Captions uniquely naming one series by 'IPCE plots for the X solar cells'.
    caps=re.findall(r'(?i)(?:IPCE|EQE)\s+(?:plots?|curves?|spectra)\s+for\s+(?:the\s+)?(.{2,55}?)\s+(?:solar\s+cells?|PSCs?|devices?)',text)
    if not caps:return
    indices=set()
    for cap in caps:indices.update(_matching(rows,cap))
    if len(indices)!=1:return
    idx=next(iter(indices))
    words=page.get_text('words')
    # Require a non-axis number with two decimal places in the integrated
    # current plot near its heading.  Reject dense numeric panels and legend.
    nums=[]
    for w in words:
        raw=w[4].strip('()%,')
        if not re.fullmatch(r'(?:[12]\d|\d)\.\d{2}',raw):continue
        v=_num(raw)
        if v is None:continue
        rect=fitz.Rect(w[:4]);nums.append((v,rect))
    if not nums:return
    labels=page.search_for('Integrated Jsc')+page.search_for('Integrated JSC')
    if not labels:return
    scored=[]
    for v,rect in nums:
        dist=min(abs(rect.x0-l.x0)+abs(rect.y0-l.y0) for l in labels)
        if dist<190:scored.append((dist,v))
    # Avoid interpreting unrelated champion JV numbers from nearby panels.
    if len(scored)==1:
        _assign(rows,idx,scored[0][1],source,'printed_plot_endpoint',counts)


def _docx_text(filepath):
    """Read Word supplementary paragraphs and table cells with python-docx."""
    from docx import Document
    doc=Document(str(filepath))
    for ix,p in enumerate(doc.paragraphs,1):
        if p.text.strip():yield p.text,f'paragraph {ix}'
    for ix,table in enumerate(doc.tables,1):
        for jx,row in enumerate(table.rows,1):
            yield ' | '.join(c.text for c in row.cells),f'table {ix} row {jx}'

def _text_candidates(rows, raw, source, counts):
    text=_clean(raw)
    for region,_ in _sentence_regions(text):
        if not (EQE.search(region) or re.search(r'(?i)integrat(?:ed|ion|ing)\s+J\s*[_-]?\s*SC',region)):
            continue
        _parse_pair(region,rows,source,counts)
        _parse_single(region,rows,source,counts)

def recover(rows, main_pdf, supplements=()):
    """Mutate only missing EQE currents; preserve existing values and evidence."""
    counter=[0]
    for filepath in [main_pdf,*supplements]:
        suffix=Path(filepath).suffix.lower()
        if suffix=='.docx':
            try:
                for part,where in _docx_text(filepath):
                    _text_candidates(rows,part,f'{Path(filepath).name} {where}',counter)
            except Exception as exc:
                print(f'  EQE_WORD_READ_SKIPPED: {Path(filepath).name}: {exc}',flush=True)
            continue
        if suffix!='.pdf':continue
        try:
            doc=fitz.open(str(filepath))
        except Exception:continue
        with doc:
            for page in doc:
                source=f'{Path(filepath).name} p.{page.number+1}'
                _text_candidates(rows,page.get_text('text'),source,counter)
                _figure_annotation(page,rows,source,counter)
    return counter[0]
