"""v18.6.6: additive, DOI-independent supplementary J-V table recovery.

Handles η/ηMPP(T)/efficiency headings and PDFs whose cells are positioned
independently; matches FULL device+contact+scan identities. Never borrows
champion values from a different contact or tandem architecture.

PCE calculation uses same-row Voc[V]*Jsc[mA/cm²]*FF[%]/100 under 100 mW/cm²;
an η_MPPT output is deliberately not assumed identical to J-V PCE.
"""
from __future__ import annotations
import math, re
from pathlib import Path
import fitz

_NUM=re.compile(r'^[-+]?\d{1,3}(?:[.,]\d{1,5})?%?$')
METRICS=('voc0','jsc0','ff0','pce0','mppt_pce')

def clean(s):
    s=str(s).replace('\u00ad','').replace('−','-').replace('–','-').replace('₂','2').replace('₀','0').replace('ₛ','s').replace('ᶜ','c')
    s=s.replace('η','eta').replace('Η','eta').replace('𝜂','eta').replace('𝛈','eta').replace('𝜼','eta')
    s=s.replace('𝑉','V').replace('𝑱','J').replace('𝒇','f')
    return re.sub(r'\s+',' ',s).strip()

def namekey(s):
    s=clean(s).lower().replace('solution-processed','solution').replace('co-evaporated','evaporation').replace('evaporated','evaporation')
    s=re.sub(r'\b(sno\s*2)\b','sno2',s)
    s=re.sub(r'\b(solar cell|solar cells|psc|pscs|device|devices)\b',' ',s)
    return re.sub('[^a-z0-9]','',s)

def number(v):
    v=clean(v).replace('%','').replace(',','.')
    if not _NUM.fullmatch(v):return None
    try:
        x=float(v)
        return x if math.isfinite(x) else None
    except ValueError:return None

def field_for_header(s):
    s=clean(s).lower()
    s=re.sub(r'\s+','',s)
    if re.search(r'(?:voc|v_oc|opencircuitvoltage)',s):return 'voc0'
    if re.search(r'(?:jsc|j_sc|shortcircuitcurrent)',s):return 'jsc0'
    if re.search(r'(?:fillfactor|^ff(?:\b|\(|\[|$))',s):return 'ff0'
    if re.search(r'(?:etamppt|etamp(p|t)|etastabilized|stabilizedpower)',s):return 'mppt_pce'
    if re.search(r'(?:^eta(?:\b|\(|\[|$)|etajv|pce|efficiency)',s):return 'pce0'
    return None

def classify_header(words):
    """Header tokens can be split over several words, e.g. η MPPT (%)."""
    terms=[clean(w).lower() for w in words]
    result={}
    for i,w in enumerate(terms):
        probes=[w, ' '.join(terms[i:i+2]), ' '.join(terms[i:i+3])]
        for probe in probes:
            f=field_for_header(probe)
            if f and f not in result:result[f]=i
    return result

def _identity_ok(row_label, source_label):
    """Match contact/process identity, ignoring standalone J-V scan qualifiers.

    F/R/B/Fw/Bw/etc are measurement directions, never material identities.
    """
    from hysteresis_recovery import canonical_key
    r=canonical_key(row_label);s=canonical_key(source_label)
    if not r or not s:return False
    return r==s

def _derive(values):
    m=dict(values)
    for k,lo,hi in [('voc0',.25,2.8),('jsc0',2,45),('ff0',.25,100),('pce0',2,40),('mppt_pce',2,40)]:
        if k in m and not(lo<=m[k]<=hi):return None
    if m.get('ff0',100)<1.05:m['ff0']*=100
    if all(k in m for k in ('voc0','jsc0','ff0')):
        calculated=m['voc0']*m['jsc0']*m['ff0']/100
        if not 1<calculated<40:return None
        if 'pce0' in m:
            # if η is MPPT, not J-V, it must be in distinct column
            if abs(m['pce0']-calculated)>max(.9,.055*m['pce0']):return None
        else:
            m['pce0']=round(calculated,5)
            m['pce_calculated']=True
    else:
        if 'pce0' not in m:return None
    return m

def _from_lines(doc,page_index):
    """Read selectable plain text, including split-header and wrapped row labels."""
    page=doc[page_index]; lines=[clean(s) for s in page.get_text('text').splitlines() if clean(s)]
    out=[]
    for i,ln in enumerate(lines):
        if not (re.search(r'\b(?:table|photovoltaic|performance|j.v)\b', ' '.join(lines[max(0,i-9):i+1]),re.I) or i<15):continue
        words=ln.split()
        if len(words)<3:continue
        nums=[]
        for j,w in enumerate(words):
            n=number(w)
            if n is not None:nums.append((j,n))
        if len(nums)<3:continue
        label=' '.join(words[:nums[0][0]])
        if not label or len(label)>72:continue
        # Look backwards for a compact four-column header. Tables have many
        # formatting styles; never infer column order when no header exists.
        for k in range(max(0,i-9),i):
            head=' '.join(lines[k:i]); hdr=classify_header(head.split())
            present=[x for x in METRICS if x in hdr]
            if len(set(present))<3:continue
            # Only interpret a simple row when number of numeric columns
            # equals the number of explicitly named header metrics.
            ordered=[f for f,_ in sorted(((f,idx) for f,idx in hdr.items()),key=lambda x:x[1])]
            if len(nums)!=len(ordered) or len(set(ordered))!=len(ordered):continue
            result=_derive(dict(zip(ordered,(v for _,v in nums))))
            if result:
                out.append((label,result,page_index+1,'line table: '+ln[:180]))
            break
    return out

def _from_spatial(doc,page_index, words_override=None):
    """Fit text columns using word x-coordinates, grouping nearby baselines."""
    page=doc[page_index]; words=words_override if words_override is not None else page.get_text('words')
    if not words:return []
    ys=[]
    for w in words:
        cy=(w[1]+w[3])/2
        group=next((g for g in ys if abs(g[0]-cy)<3.3),None)
        if group:group[1].append(w);group[0]=(group[0]*(len(group[1])-1)+cy)/len(group[1])
        else:ys.append([cy,[w]])
    lines=[sorted(g[1],key=lambda w:w[0]) for g in sorted(ys,key=lambda g:g[0])]
    out=[]
    for hi,header_words in enumerate(lines):
        header=' '.join(str(w[4]) for w in header_words)
        if not re.search(r'\b(?:Voc|Jsc|J\s*sc|V\s*oc|fill factor|FF)\b',header,re.I):continue
        # Detect *separate* metric column anchors rather than assuming order.
        columns={}
        for j,w in enumerate(header_words):
            token=clean(w[4]);f=field_for_header(token)
            if token.lower() in ('eta','η') and j+1<len(header_words) and clean(header_words[j+1][4]).lower().startswith('mppt'):
                f='mppt_pce'
            if f is None and token.lower() in ('v','j','eta','η','open','short','fill','power') and j+1<len(header_words):
                f=field_for_header(token+' '+clean(header_words[j+1][4]))
            if f and f not in columns:columns[f]=.5*(w[0]+w[2])
        if not all(k in columns for k in ('voc0','jsc0','ff0')):continue
        xpos=sorted(columns.values()); minx=min(xpos)
        for line in lines[hi+1:hi+19]:
            numbers=[];label=[]
            for w in line:
                x=.5*(w[0]+w[2]);n=number(w[4])
                if n is not None and x>=minx-27:numbers.append((x,n))
                elif x<minx-8:label.append(str(w[4]))
            if len(numbers)<3 or not label:continue
            label=' '.join(label)
            if len(label)>80:continue
            values={}
            for x,v in numbers:
                f=min(columns,key=lambda field:abs(columns[field]-x))
                if abs(columns[f]-x)<max(38,.45*(max(xpos)-min(xpos))):
                    if f in values:values={};break
                    values[f]=v
            if not all(k in values for k in ('voc0','jsc0','ff0')):continue
            result=_derive(values)
            if result:out.append((label,result,page_index+1,'spatial table: '+label+' '+str(values)))
    return out

def _from_pymupdf_tables(doc,page_index):
    """Read ruled and borderless detected tables with multi-line metric headers.

    A row can have separate material and scan columns, or a multi-cell material
    label.  Preserve all non-numeric label cells and reject value mixing.
    """
    page=doc[page_index]
    try:tables=page.find_tables().tables
    except Exception:return []
    output=[]
    for table in tables:
        try:grid=table.extract()
        except Exception:continue
        if not grid:continue
        grid=[[clean(x or '') for x in line] for line in grid]
        width=max(map(len,grid))
        grid=[r+['']*(width-len(r)) for r in grid]
        for hi in range(min(len(grid),7)):
            for nh in (1,2,3):
                if hi+nh>len(grid):continue
                header=[' '.join(grid[j][i] for j in range(hi,hi+nh)).strip()
                        for i in range(width)]
                names=[field_for_header(h) for h in header]
                if not {'voc0','jsc0','ff0'}.issubset(set(names)):continue
                if len([f for f in names if f in METRICS])!=len(set(f for f in names if f in METRICS)):
                    continue
                for line in grid[hi+nh:]:
                    vals={};labels=[]
                    for val,field in zip(line,names):
                        n=number(val)
                        if field and n is not None:vals[field]=n
                        elif not field and val and (n is None):labels.append(val)
                    if not labels or not {'voc0','jsc0','ff0'}.issubset(vals):continue
                    # Combine independent scan column with material label.
                    label=' '.join(labels)
                    if len(label)>150:continue
                    result=_derive(vals)
                    if result:output.append((label,result,page_index+1,
                        'PyMuPDF multiheader table: '+str(line)[:180]))
                break
            else:continue
            break
    return output

def recover(rows,main_pdf,supplements=()):
    source_files=[Path(p) for p in [main_pdf,*supplements] if str(p).lower().endswith('.pdf')]
    candidates=[]
    for path in source_files:
        if not path.exists():continue
        try:
            with fitz.open(path) as doc:
                scanned_ocr=0
                for i in range(len(doc)):
                    native=_from_pymupdf_tables(doc,i)+_from_spatial(doc,i)+_from_lines(doc,i)
                    for label,values,page,evidence in native:
                        candidates.append((label,values,f'{path.name} p.{page}: {evidence}'))
                    # Optional OCR for genuinely rasterized SI table/figure panels.
                    # Limited to pages with negligible native text, and maximum
                    # three pages per PDF to avoid lengthy large-folder runs.
                    if not native and scanned_ocr<3 and len(doc[i].get_text('words'))<12:
                        try:
                            from general_curve_reader import ocr_page_words
                            ocr_words=ocr_page_words(doc[i],zoom=2.5)
                        except Exception:ocr_words=[]
                        scanned_ocr+=1
                        for label,values,page,evidence in _from_spatial(doc,i,ocr_words):
                            candidates.append((label,values,f'{path.name} p.{page}: OCR {evidence}'))
        except Exception:continue
    changed=0
    for row in rows:
        label=row.get('device_category') or ''
        matches=[c for c in candidates if _identity_ok(label,c[0])]
        if not matches:
            if not candidates:
                row['eta_supplement_jv_status']='NO_MACHINE_READABLE_JV_TABLE_FOUND'
            else:
                row['eta_supplement_jv_status']='NO_EXACT_DEVICE_CONTACT_SCAN_MATCH'
                row['eta_supplement_jv_candidate_count']=len(candidates)
            continue
        # Differing J-V sets of same identity/scan are not averaged or merged.
        distinct={(tuple((k,round(v,4)) for k,v in sorted(m.items()) if isinstance(v,(int,float)) and k!='pce_calculated')) for _,m,_ in matches}
        if len(distinct)!=1:
            row['eta_supplement_jv_status']='CONFLICTING_SAME_CATEGORY_SETS_REVIEW_SCAN_AREA';row['eta_supplement_jv_candidate_count']=len(matches);continue
        src_label,vals,source=matches[0]
        for field in ('voc0','jsc0','ff0'):
            if field in vals and row.get(field) in (None,''):
                row[field]=round(vals[field],5);row[field+'_source']='Measured, same-category '+source;changed+=1
        if 'pce0' in vals:
            if row.get('pce_jv_pct') in (None,''):
                row['pce_jv_pct']=round(vals['pce0'],5);row['pce_jv_source']=('Calculated from same-row Voc*Jsc*FF/100; ' if vals.get('pce_calculated') else 'Reported J-V efficiency; ')+source
            if row.get('pce0') in (None,''):
                row['pce0']=round(vals['pce0'],5);row['pce0_source']='Representative J-V PCE, not confirmed aging-initial; '+source;changed+=1
        if 'mppt_pce' in vals:
            row.setdefault('mppt_eta_pct',vals['mppt_pce'])
        row['eta_supplement_jv_status']='SOURCE_MATCHED_'+('DERIVED_PCE' if vals.get('pce_calculated') else 'REPORTED_PCE')
        row['eta_supplement_jv_evidence']=source
    return changed
