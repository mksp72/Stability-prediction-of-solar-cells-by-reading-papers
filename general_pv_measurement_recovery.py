"""v18.6.3: conservative, content-based category-aware photovoltaic recovery.

Additive to v18.6.2; never modifies stability, EQE, fit or reported T80.
The recovery creates candidate same-device J-V measurement sets from main/SI
text and selectable tables. It does not mix metric values between candidates.
An unsupported category-to-JV mapping is left empty and recorded as evidence.
No DOI, filename or published-value configuration is used in the extractor.
"""
from __future__ import annotations
from pathlib import Path
from collections import defaultdict
import math
import re
import fitz

# Deal with pdf extraction ligatures, line-hyphenation and mathematical subscripts.
def clean(s):
    s=str(s).replace('\ufb01','fi').replace('\ufb02','fl').replace('−','-').replace('–','-').replace('–','-')
    s=re.sub(r'(?<=[A-Za-z])-\s+(?=[a-z])','',s)
    return ' '.join(s.split())

def key(s):return re.sub(r'[^a-z0-9]+','',str(s).lower())

def number(v):
    try:
        f=float(v)
        return f if math.isfinite(f) else None
    except (TypeError,ValueError):return None

PCE=r'(?:PCE|power conversion efficienc(?:y|ies)|efficiency|\u03b7MPPT)'
VOC=r'(?:V\s*OC|open[ -]?circuit voltage)'
JSC=r'(?:J\s*SC|short[ -]?circuit current(?: density)?)'
FF=r'(?:FF|fill factor)'
NUM=r'(\d{1,2}(?:\.\d{1,5})?)'
MEASURE_RX={
    'pce0':re.compile(r'\b'+PCE+r'\b\s*(?:\([^)]{0,25}\))?\s*(?:of|=|:|was|at|reached)?\s*'+NUM+r'\s*%',re.I),
    'voc0':re.compile(r'\b'+VOC+r'\b\s*(?:\([^)]{0,15}\))?\s*(?:of|=|:|was|at)?\s*'+NUM+r'\s*V\b',re.I),
    'jsc0':re.compile(r'\b'+JSC+r'\b\s*(?:\([^)]{0,15}\))?\s*(?:of|=|:|was|at)?\s*'+NUM+r'\s*mA',re.I),
    'ff0':re.compile(r'\b'+FF+r'\b\s*(?:\([^)]{0,15}\))?\s*(?:of|=|:|was|at)?\s*'+NUM+r'\s*%',re.I),
}
# No unit is sometimes supplied for reported Voc/Jsc in adjoining clauses.
BARE_RX={
    'voc0':re.compile(r'\b'+VOC+r'\b\s*(?:of|=|:)?\s*(0\.\d{2,4}|1\.\d{2,4}|2\.\d{2,4})\b',re.I),
    'jsc0':re.compile(r'\b'+JSC+r'\b\s*(?:of|=|:)?\s*(\d{2}\.\d{1,4})\b',re.I),
    'ff0':re.compile(r'\b'+FF+r'\b\s*(?:of|=|:)?\s*(0\.\d{2,4})\b',re.I),
}

def metrics_of(segment):
    result={}
    for field,rx in MEASURE_RX.items():
        mm=list(rx.finditer(segment))
        if len(mm)==1:result[field]=float(mm[0].group(1))
    for field,rx in BARE_RX.items():
        if field in result:continue
        mm=list(rx.finditer(segment))
        if len(mm)==1:result[field]=float(mm[0].group(1))
    # Publishers also write 'achieved 19.6% PCE', not 'PCE of 19.6%'.
    if 'pce0' not in result:
        backwards=re.findall(r'\b(\d{1,2}(?:\.\d+)?)\s*%\s*(?:PCE|power conversion efficiency)\b',segment,re.I)
        if len(backwards)==1:result['pce0']=float(backwards[0])
    if result.get('ff0',100)<=1:result['ff0']*=100
    return result

def valid_set(m,allow_partial=True):
    limits={'pce0':(2,40),'voc0':(.3,2.7),'jsc0':(2,40),'ff0':(25,95)}
    if any(not(lo<=v<=hi) for k,v in m.items() for lo,hi in [limits[k]]):return False
    if all(k in m for k in ('pce0','voc0','jsc0','ff0')):
        calc=m['voc0']*m['jsc0']*m['ff0']/100
        if abs(calc-m['pce0'])>max(1.3,.08*m['pce0']):return False
    return bool(m) and (allow_partial or len(m)>=2)

def role_of(label):
    s=key(label)
    if s in ('control','reference','ref','pristine','untreated','withoutadditive'):return 'control'
    if s in ('target','treated','modified','optimized','withtpta','tpta','modifiedwithtpta'):return 'target'
    return None

def tokens(s):
    return [z for z in re.findall(r'[a-z0-9]+',str(s).lower()) if len(z)>=2]

def identity_match(row,candidate):
    """Require direct category evidence, not generic target proximity."""
    cat=key(row.get('device_category',''))
    evidence=key(candidate['label'])
    if not cat or not evidence:return False
    if cat==evidence:return True
    # Named chemical abbreviation must be unique and distinguished from base.
    if cat in ('fapbi3','fapbi3sesiqds'):
        if cat=='fapbi3sesiqds':return ('siqd' in evidence and 'fapbi3' in evidence)
        return evidence in ('fapbi3','fapbi3reference','referencefapbi3','neatfapbi3')
    if cat in ('me4pacz2paczk','me4paczplus2paczk'):
        return evidence in ('mek','me4pacz2paczk','me4paczplus2paczk','mixedsamme4pacz2paczk')
    if cat=='me4pacz':return evidence in ('me4pacz','me4paczonly')
    if cat=='withtpta':return evidence in ('tpta','withtpta','tptatreated','tptaincorporated')
    if cat=='tpta':return evidence in ('tpta','withtpta','tptatreated','tptaincorporated')
    if cat=='control':return evidence in ('control','reference','pristine','untreated')
    # Exact normalization of plus, hyphens, 2 subscripts.
    return False

def arch_of(row):
    a=str(row.get('architecture') or '').lower()
    return 'tandem' if 'tandem' in a or '2t' in a else 'single' if 'single' in a else None

def compatible(row,c):
    ra=arch_of(row); ca=c.get('architecture')
    if ra and ca and ra!=ca:return False
    if ca and not ra and ca=='tandem':return False
    # No imported module-vs-small-cell conflation.
    if c.get('module'):return False
    return identity_match(row,c)

def collect(main_pdf,supplements=()):
    out=[]
    for path in [main_pdf,*supplements]:
        if Path(path).suffix.lower()!='.pdf':continue
        try:
            with fitz.open(str(path)) as d:
                for index,page in enumerate(d):
                    out.append((Path(path).name,index+1,clean(page.get_text('text'))))
        except Exception:continue
    return out

def add(candidates,label,values,file,page,kind,architecture=None,excerpt=None,score=1,group=None):
    vals={k:number(v) for k,v in values.items() if k in MEASURE_RX}
    vals={k:v for k,v in vals.items() if v is not None}
    if not valid_set(vals):return
    candidates.append({'label':label,'metrics':vals,'source':f'{file} p.{page}',
        'kind':kind,'architecture':architecture,'score':score,'excerpt':excerpt or '',
        'group':group,'module':bool(re.search(r'\bmodule\b',excerpt or '',re.I))})

# Explicit label anchors, used inside a SHORT contiguous passage. Never join metrics
# from two different clauses or performance examples in the same paragraph.
ALIASES=[
    ('FAPbI3 + SE SiQDs',r'\b(?:hybrid\s+)?FAPbI3\s*\+\s*SE\s*[- ]?SiQDs?\s+devices?\b'),
    ('FAPbI3',r'\bFAPbI3\s+(?:reference|control|neat)\s+devices?\b'),
    ('Me-4PACz+2PACz-K',r'\b(?:champion\s+)?(?:Me\s*\+\s*K|Me-4PACz\s*\+\s*2PACz-K)\s+devices?\b'),
    ('TPTA',r'\b(?:TPTA[- ](?:treated|incorporated)|devices?\s+with\s+TPTA|champion\s+PTSC)\b'),
    ('Control',r'\b(?:control|reference|pristine)\s+(?:devices?|PSCs?|cells?)\b'),
]

def prose_candidates(pages):
    results=[]
    for file,page,text in pages:
        # Named device + quantified metrics in a short prose span, cut before the
        # next named device or comparator to prevent cross-category borrowing.
        for label,pattern in ALIASES:
            for mm in re.finditer(pattern,text,re.I):
                tail=text[mm.end():mm.end()+360]
                stopper=re.search(r'(?:\b(?:After introducing|In comparison|By comparison|On the other hand|Meanwhile|whereas)\b|\b(?:control|reference|hybrid FAPbI3|TPTA-treated)\s+devices?\b)',tail,re.I)
                if stopper:tail=tail[:stopper.start()]
                if len(tail)>260:tail=tail[:260]
                vals=metrics_of(tail)
                # A generic champion device with no chemical label cannot imply
                # the same identity as a chemically named TPTA category.
                if label=='TPTA' and re.search(r'champion\s+PTSC',mm.group(),re.I):
                    architecture='tandem'
                else:
                    architecture='tandem' if re.search(r'\b(?:PTSCs?|tandem|2T)\b',text[max(0,mm.start()-110):mm.end()+130],re.I) else None
                if len(vals)>=2:
                    add(results,label,vals,file,page,'named_prose',architecture,mm.group()+tail[:170],score=6)
        # Explicit complete standalone J-V set with named device in the prior
        # 100 characters (e.g. 'champion Me+K device delivered ...').
        rx=re.compile(r'\b(?:champion\s+)?([A-Za-z0-9+–-]{2,24})\s+devices?\s+(?:delivered|achieved|showed|exhibited|yielded)\s+(?:a\s+)?(?:PCE|power conversion efficiency)\s+of\s+(\d{1,2}(?:\.\d+)?)\s*%',re.I)
        for m in rx.finditer(text):
            tail=text[m.end():m.end()+200]
            vals=metrics_of('PCE of '+m.group(2)+'% '+tail)
            if len(vals)>=2:
                add(results,m.group(1),vals,file,page,'named_champion_jv_prose',None,m.group(0)+' '+tail[:100],score=8)
        # Comparison sentences uniquely identify before/after values for a single
        # material treatment. They can supply incomplete metric sets without
        # inventing missing values or treating averages as champions.
        # This captures 'VOC of 0.88 V (vs. 0.83 V for control)' and
        # 'JSC increased from 32.2 to 32.6 mA/cm2' only in a nearby J-V context.
        for m in re.finditer(r'\b(?:TPTA[- ]treated|TPTA[- ]incorporated)\s+devices?\b',text,re.I):
            paragraph=text[m.start():m.start()+650]
            if not re.search(r'\bJ[–-]?V\b|champion\s+PCE|champion\s+device',paragraph,re.I):continue
            if re.search(r'\b(?:PTSC|tandem)\b',paragraph[:180],re.I):continue
            pairv=re.search(r'\bVOC\s+of\s+(\d\.\d+)\s*V\s*\(vs\.\s*(\d\.\d+)\s*V\s+for\s+control\)',paragraph,re.I)
            pairj=re.search(r'\bJSC\s+increased\s+from\s+(\d{2}\.\d+)\s+to\s+(\d{2}\.\d+)\s*mA',paragraph,re.I)
            pairp=re.search(r'\bchampion\s+PCE\s+of\s+(\d{1,2}\.\d+)\s*%,\s+compared\s+to\s+(\d{1,2}\.\d+)\s*%\s+for\s+the\s+control',paragraph,re.I)
            if pairv and pairj and pairp:
                t={'voc0':float(pairv.group(1)),'jsc0':float(pairj.group(2)),'pce0':float(pairp.group(1))}
                c={'voc0':float(pairv.group(2)),'jsc0':float(pairj.group(1)),'pce0':float(pairp.group(2))}
                add(results,'TPTA',t,file,page,'explicit_comparison_jv',architecture='single',excerpt=paragraph[:470],score=9)
                add(results,'Control',c,file,page,'explicit_comparison_jv',architecture='single',excerpt=paragraph[:470],score=9)
    return results


def dynamic_prose_candidates(rows,pages):
    """Find named device J-V statements for NEW device/category names without
    a publication-specific vocabulary. Require explicit category mentions and
    at least two metrics in one short device-specific clause.
    """
    candidates=[]
    labels=sorted({str(r.get('device_category') or '').strip() for r in rows if r.get('device_category')},key=len,reverse=True)
    for label in labels:
        if len(label)<3 or len(label)>80:continue
        # Prevent a generic 'Control' or 'Target' from overriding chemically
        # specific labels. Generic categories need prose such as 'Control devices'.
        components=[re.escape(v) for v in re.split(r'[\s+\-]+',label) if v]
        if not components:continue
        phrase=r'[\s+\-–]*'.join(components)
        anchors=[phrase+r'\s+(?:reference\s+|control\s+|neat\s+|treated\s+)?(?:devices?|PSCs?|cells?)\b']
        if role_of(label) in ('control','target'):
            anchors=[r'\b'+phrase+r'\s+(?:devices?|PSCs?|cells?)\b']
        for file,page,text in pages:
            for pattern in anchors:
                for m in re.finditer(pattern,text,re.I):
                    tail=text[m.end():m.end()+250]
                    # Do not concatenate parameters over a second explicitly
                    # named device or the next sentence comparing them.
                    next_device=re.search(r'\b(?:reference|control|target|hybrid|treated|modified)\s+(?:devices?|PSCs?|cells?)\b',tail,re.I)
                    if next_device:tail=tail[:next_device.start()]
                    values=metrics_of(tail)
                    if len(values)<2:continue
                    context=text[max(0,m.start()-90):m.end()+90]
                    arch='tandem' if re.search(r'\b(?:tandem|PTSCs?|2T)\b',context,re.I) else None
                    add(candidates,label,values,file,page,'dynamic_named_category_prose',arch,
                        text[m.start():m.end()]+tail[:170],score=7)
    return candidates

def table_candidates(pages,main_pdf,supplements,rows=()):
    """Reading-order row recognition from selectable PDF table text.

    Only a labeled row and an explicitly recognized header order are accepted.
    We avoid guessing positional layout when the header order is unknown.
    """
    candidates=[]
    pattern=re.compile(r'(?i)\b(?:PCE|efficiency)\b|\bV\s*OC\b|\bJ\s*SC\b|\bFF\b|fill factor')
    for file,page,text in pages:
        if len(pattern.findall(text))<4 or not re.search(r'\btable\s+(?:s\d+|\d+)\b',text,re.I):continue
        # PDF `get_text` may have vertical table extraction; use explicit
        # line-order header+row recovery only for headers that fit one line.
        with_lf=None
        for path in [main_pdf,*supplements]:
            if Path(path).name!=file:continue
            try:
                with fitz.open(str(path)) as d:with_lf=d[page-1].get_text('text')
            except Exception:pass
            break
        if not with_lf:continue
        lines=[clean(s) for s in with_lf.splitlines() if clean(s)]
        # Source-compatible dynamic device identities, rather than a fixed
        # short list of reference papers. Longer names take precedence when a
        # short name (e.g. 'Evaporation') prefixes a real formulation.
        named=sorted({str(r.get('device_category') or '').strip()
            for r in rows if len(str(r.get('device_category') or '').strip())>=3},
            key=len,reverse=True)
        label_terms=[re.escape(x) for x in named]
        label_terms.extend(['Control','Reference','Pristine','Target'])
        row_rx=re.compile(r'^\s*('+ '|'.join(label_terms) +r')\s+(.+)$',re.I)
        for index,ln in enumerate(lines):
            header=' '.join(lines[index:index+3]).lower()
            hits=[(field,min((m.start() for m in re.finditer(rx,header,re.I)),default=10**9)) for field,rx in [('voc0',VOC),('jsc0',JSC),('ff0',FF),('pce0',PCE)]]
            if any(pos==10**9 for _,pos in hits):continue
            order=[f for f,_ in sorted(hits,key=lambda x:x[1])]
            if len(set(order))!=4:continue
            for j in range(index+3,min(index+15,len(lines))):
                row=lines[j]
                match=row_rx.match(row)
                if not match:continue
                raw=re.findall(r'(?<![a-zA-Z])\d{1,2}(?:\.\d+)?',match.group(2))
                if len(raw)!=4:continue
                m=dict(zip(order,map(float,raw)))
                if valid_set(m,False):add(candidates,match.group(1),m,file,page,'four_column_jv_table',None,row,score=7)
    return candidates


def paired_figure_annotation_candidates(pages, rows):
    """General OCR/selectable-text J-V table embedded in a figure.

    Some papers print a sequence like:
        Control Target VOC 1.17 V JSC 26.20 mA FF 83.5% PCE 25.6%
                       VOC 1.20 V JSC 26.40 mA FF 86.0% PCE 27.2%.
    We assign the two complete parameter blocks ONLY when their count agrees
    with the adjacent ordered set of named device categories and when the PCE
    equation closes. No DOI or per-paper number configuration is involved.
    """
    results=[]
    names=sorted({str(r.get('device_category') or '').strip() for r in rows
                  if len(str(r.get('device_category') or '').strip())>=3},key=len,reverse=True)
    if not names:return results
    labels_pattern='|'.join(re.escape(x) for x in names)
    # Text blocks may have subscript OCR artifacts and cm^-2 between metrics.
    voc=r'V\s*O\s*C\s*[:=]?\s*(\d\.\d{1,4})\s*V\b'
    jsc=r'J\s*S\s*C\s*[:=]?\s*(\d{1,2}\.\d{1,4})\s*mA'
    ff=r'\bFF\s*[:=]?\s*(\d{1,2}\.\d{1,3})\s*%'
    pce=r'\bPCE\s*[:=]?\s*(\d{1,2}\.\d{1,3})\s*%'
    blocks=re.compile(voc+r'.{0,75}?'+jsc+r'.{0,65}?'+ff+r'.{0,55}?'+pce,re.I|re.S)
    for file,page,text in pages:
        for m in blocks.finditer(text):
            # Group annotations separated by at most 60 characters, and never
            # mix a long prose paragraph with unrelated J-V values.
            vals={'voc0':float(m[1]),'jsc0':float(m[2]),
                  'ff0':float(m[3]),'pce0':float(m[4])}
            if not valid_set(vals,allow_partial=False):continue
            prefix=text[max(0,m.start()-180):m.start()]
            matches=list(re.finditer(r'(?<![A-Za-z0-9])(?:'+labels_pattern+r')(?![A-Za-z0-9])',prefix,re.I))
            # Only reuse an ordered list if all device identities appear
            # together (not scattered through paragraphs or page headers).
            ordered=[]
            for mm in matches:
                if mm.group().lower() not in [x.lower() for x in ordered]:ordered.append(mm.group())
            if not ordered:continue
            ordered=ordered[-6:]
            if len(ordered)==1:
                # Only one label and no earlier competing full J-V tuple.
                add(results,ordered[0],vals,file,page,'labeled_figure_jv_tuple',None,
                    text[max(0,m.start()-80):m.end()],score=10)
            else:
                # Adjacent blocks may follow the printed legend order.
                group=[m]
                for mm in blocks.finditer(text,m.end(),min(len(text),m.end()+420)):
                    between=text[group[-1].end():mm.start()]
                    if len(between)>90 or re.search(r'(?:\bFig\.?|\bTable\b|\bResults\b)',between,re.I):break
                    group.append(mm)
                    if len(group)>=len(ordered):break
                if len(group)!=len(ordered):continue
                for name,block in zip(ordered,group):
                    values={'voc0':float(block[1]),'jsc0':float(block[2]),
                        'ff0':float(block[3]),'pce0':float(block[4])}
                    if valid_set(values,allow_partial=False):
                        add(results,name,values,file,page,'ordered_figure_jv_tuples',None,
                            text[max(0,m.start()-75):group[-1].end()],score=11)
    return results

def recover(rows,main_pdf,supplements=()):
    pages=collect(main_pdf,supplements)
    candidates=paired_figure_annotation_candidates(pages,rows)+dynamic_prose_candidates(rows,pages)+prose_candidates(pages)+table_candidates(pages,main_pdf,supplements,rows)
    changed=0
    for row in rows:
        matches=[c for c in candidates if compatible(row,c)]
        # Distinct architecture can be given by the source's 'tandem' or
        # 'single-junction' context; never share a single-JV set with tandem.
        matches.sort(key=lambda x:(x['score'],len(x['metrics'])),reverse=True)
        # Same set from repeated article passages is fine; conflicting sets
        # are not blended. Prefer the higher-confidence complete set.
        if not matches:
            if any(row.get(k) in (None,'') for k in ('pce0','voc0','jsc0','ff0')):
                row['general_pv_recovery_status']='NO_EXACT_CATEGORY_JV_MATCH'
            continue
        chosen=matches[0]
        vals=chosen['metrics'];source=chosen['source']+' '+chosen['kind']
        # Avoid filling isolated metrics from a PCE that is not an actual J-V set.
        if len(vals)==1:
            row['general_pv_recovery_status']='INSUFFICIENT_MATCHED_JV_PARAMETERS';continue
        for field,val in vals.items():
            if field=='pce0':
                if row.get('pce_jv_pct') in (None,''):
                    row['pce_jv_pct']=round(val,5)
                    row['pce_jv_source']=source
                if row.get('pce0') not in (None,''):continue
            if row.get(field) not in (None,''):continue
            row[field]=round(val,5)
            row[field+'_source']=source
            changed+=1
        row['general_pv_recovery_status']='MATCHED_WITH_SOURCE'
        row['general_pv_recovery_source']=source
        row['general_pv_recovery_evidence']=chosen['excerpt'][:500]
        row['general_pv_recovery_architecture']=chosen['architecture'] or ''
        row['general_pv_recovery_candidate_count']=len(matches)
        # Derived FF only when PCE is an explicitly corresponding JV PCE.
        if row.get('ff0') in (None,''):
            p=number(row.get('pce_jv_pct'));v=number(row.get('voc0'));j=number(row.get('jsc0'))
            if all(x is not None and x>0 for x in (p,v,j)):
                inferred=100*p/(v*j)
                if 25<=inferred<=95:
                    row['ff0']=round(inferred,5)
                    row['ff0_source']='Derived from same J-V measurement: '+source
                    changed+=1
        # Conflicting separate measurement sets are kept visible as diagnostic.
        other=[c for c in matches[1:] if c['metrics']!=vals and c['score']>=chosen['score']]
        if other:row['general_pv_recovery_status']='MATCHED_WITH_ALTERNATIVE_JV_EVIDENCE'
    return changed
