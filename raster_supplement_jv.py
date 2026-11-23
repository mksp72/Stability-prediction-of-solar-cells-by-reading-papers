"""Source-verified OCR recovery for vector-outlined/raster J-V supplement tables.

Activated only on pages containing a table caption but no selectable tabular
numbers. Uses printed SCAN, material, and FOUR physically aligned PV metrics.
Supports multiple adjacent tables and active-area headings. No DOI-specific
numbers, filenames, column coordinates, or device-name whitelist.
Tesseract is an optional system executable. If unavailable no data are invented.
"""
from __future__ import annotations
from pathlib import Path
import math, os, re, shutil, subprocess, tempfile
import fitz
from hysteresis_recovery import canonical_scan_label, canonical_key
from eta_supplement_pv_recovery import _derive

NUMBER=re.compile(r'^\d{1,3}(?:[.,]\d{1,4})?$')

def _tesseract_path():
    exe=shutil.which('tesseract')
    if exe:return exe
    # Windows users commonly install Tesseract without adding it to PATH.
    if os.name=='nt':
        for root in (os.environ.get('ProgramFiles',''),
                     os.environ.get('ProgramFiles(x86)',''),
                     r'C:\Program Files'):
            if not root:continue
            candidate=Path(root)/'Tesseract-OCR'/'tesseract.exe'
            if candidate.is_file():return str(candidate)
    return None

def _ocr_words(page):
    engine=_tesseract_path()
    if engine is None:return [],0
    with tempfile.TemporaryDirectory(prefix='jv_ocr_') as td:
        raster=Path(td)/'page.png'
        # Only the first 75% is rendered: J-V tables may be above any notes.
        clip=fitz.Rect(0,0,page.rect.width,page.rect.height*.75)
        page.get_pixmap(matrix=fitz.Matrix(3,3),clip=clip,alpha=False).save(raster)
        try:
            # Tesseract emits UTF-8 TSV. Do not set text=True here: on Windows
            # subprocess otherwise decodes pipe bytes using the process locale
            # (often cp1252), and a UnicodeDecodeError in its _readerthread
            # leaves stdout=None and can break the entire PDF batch.
            proc=subprocess.run([engine,str(raster),'stdout','--psm','6','tsv'],
                stdout=subprocess.PIPE,stderr=subprocess.PIPE,
                timeout=32,check=False)
        except (OSError,subprocess.SubprocessError):
            return [],0
        if proc.returncode != 0:
            return [],0
        # Invalid bytes from native OCR engines should not crash extraction.
        # Keeping replacement characters local to the affected OCR token
        # prevents a malformed character from damaging other numeric rows.
        raw=proc.stdout if proc.stdout is not None else b''
        if isinstance(raw,bytes):
            tsv=raw.decode('utf-8',errors='replace')
        else:
            tsv=str(raw)
        words=[]
        for record in tsv.splitlines()[1:]:
            f=record.split('\t',11)
            if len(f)!=12 or not f[11].strip():continue
            try:
                x,y,w,h=float(f[6]),float(f[7]),float(f[8]),float(f[9])
            except ValueError:continue
            words.append({'text':f[11].strip(),'x':x+w/2,'y':y+h/2,'w':w,'h':h})
        return words,clip.width*3

def _row(words,y,dy=17):
    return sorted((w for w in words if abs(w['y']-y)<dy),key=lambda w:w['x'])

def _is_num(s):return NUMBER.fullmatch(s.replace(',','.')) is not None

def _num(s):
    try:return float(s.replace(',','.'))
    except ValueError:return None

def _area(words,cx,first_y):
    # Explicit area heading above this table, rather than guessed from device.
    near=[w for w in words if first_y-240<w['y']<first_y-80 and abs(w['x']-cx)<320]
    for n in near:
        if not re.fullmatch(r'\d+(?:[.,]\d+)?',n['text']):continue
        unit=next((w for w in near if w['x']>n['x'] and w['x']-n['x']<100
                   and abs(w['y']-n['y'])<14 and re.match(r'(?i)cm',w['text'])),None)
        if unit:
            try:return float(n['text'].replace(',','.'))
            except ValueError:pass
    return None

def _normalize_material(raw):
    text=raw.strip().replace('₂','2')
    # OCR may render a subscript 2 as a comma on a source whose material
    # labels elsewhere visibly include SnO2. This repair is source-audited.
    text=re.sub(r'(?i)\bsno\s*[,.;]$', 'SnO2',text)
    text=re.sub(r'(?i)\bsno\s*[,.;]\b', 'SnO2',text)
    return text

def _extract_page(page,source_file,pi):
    words,width=_ocr_words(page)
    if not words:return []
    scans=sorted((w for w in words if w['text'].lower().strip(' |]}') in
        ('reverse','forward','backward','rev','fw','bw','rw','re')
        and sum(_is_num(z['text']) for z in _row(words,w['y']) if z['x']>w['x']+70)>=4),
        key=lambda w:(w['x'],w['y']))
    if len(scans)<4:return []
    # Multiple table panels have separate scan columns. Infer their positions
    # from the data, never from fixed page coordinates.
    centers=[]
    for w in scans:
        cx=w['x'];g=next((z for z in centers if abs(z[0]-cx)<60),None)
        if g:g[1].append(w);g[0]=sum(t['x'] for t in g[1])/len(g[1])
        else:centers.append([cx,[w]])
    centers.sort(key=lambda z:z[0]);out=[]
    for ci,(scan_x,group) in enumerate(centers):
        left=0 if ci==0 else (centers[ci-1][0]+scan_x)/2
        right=width if ci==len(centers)-1 else (scan_x+centers[ci+1][0])/2+95
        area=_area([w for w in words if left<w['x']<right],(left+right)/2,min(t['y'] for t in group))
        sc=sorted(group,key=lambda w:w['y'])
        for i,w in enumerate(sc):
            if w['text'].lower().strip(' |]}') not in ('reverse','backward','rev','rw','re','bw'):continue
            # Table categories can span adjacent Reverse and Forward rows:
            # 'Evaporation' appears by Reverse; 'BCP' appears by Forward.
            next_y=next((q['y'] for q in sc if q['y']>w['y']+17),w['y']+40)
            row_words=_row(words,w['y'])
            label_words=[q for q in row_words if left<q['x']<scan_x-25
                         and re.search('[A-Za-z]',q['text'])]
            below=[q for q in words if left<q['x']<scan_x-25
                   and w['y']+18<q['y']<min(w['y']+60,next_y+16)
                   and re.search('[A-Za-z]',q['text'])]
            part_a=' '.join(q['text'].strip(' |]}') for q in sorted(label_words,key=lambda z:z['x']))
            part_b=' '.join(q['text'].strip(' |]}') for q in sorted(below,key=lambda z:(z['y'],z['x'])))
            # 'Steady-State' and the following device heading are not parts
            # of this device name; prefer closest next source row.
            part_b=re.split(r'\b(?:Steady|State|Reverse|Forward)\b',part_b,flags=re.I)[0].strip()
            if not part_a or not part_b:continue
            label=_normalize_material(part_a+' '+part_b)
            if any(x.lower() in label.lower() for x in ('devices','direction')):continue
            # One row has four numeric PV metrics, in the printed column order.
            numbers=[q for q in row_words if scan_x+70<q['x']<right and _is_num(q['text'])]
            if len(numbers)!=4:continue
            metrics=list(map(lambda q:_num(q['text']),numbers))
            if any(x is None for x in metrics):continue
            voc,jsc,ff,pce=metrics
            decimal_repaired=False
            if 100<ff<1000 and float(ff).is_integer():
                ff=ff/10;decimal_repaired=True
            if 100<pce<1000 and float(pce).is_integer():
                pce=pce/10;decimal_repaired=True
            vals=_derive({'voc0':voc,'jsc0':jsc,'ff0':ff,'pce0':pce})
            if not vals:continue
            out.append(dict(device_label=label,scan_direction='reverse',
                area_cm2=area,values=vals,source_file=source_file,
                source_page=pi+1,method='OCR_SAME_ROW_TABLE',
                evidence='OCR four aligned numeric columns with PV closure'+
                  ('; FF decimal corrected' if decimal_repaired else '')))
            # Recover forward scan from same block for HI, with exactly the
            # same category identity and active area.
            fwd=next((q for q in sc if q['text'].lower().strip(' |]}') in ('forward','fw','fwd')
                and 16<q['y']-w['y']<80),None)
            if fwd:
                valswords=[q for q in _row(words,fwd['y']) if scan_x+70<q['x']<right and _is_num(q['text'])]
                if len(valswords)==4:
                    vv=list(map(lambda q:_num(q['text']),valswords))
                    if all(v is not None for v in vv):
                        v,j,f,p=vv
                        if 100<f<1000 and float(f).is_integer():f/=10
                        if 100<p<1000 and float(p).is_integer():p/=10
                        fvals=_derive({'voc0':v,'jsc0':j,'ff0':f,'pce0':p})
                        if fvals:out.append(dict(device_label=label,scan_direction='forward',
                            area_cm2=area,values=fvals,source_file=source_file,
                            source_page=pi+1,method='OCR_SAME_ROW_TABLE',
                            evidence='OCR aligned forward J-V set with PV closure'))
    # Deduplicate same page/area/direction/identity/values across OCR rows.
    dedup={}
    for d in out:
        key=(d['source_page'],canonical_key(d['device_label']),d['scan_direction'],
             d['area_cm2'],tuple(round(d['values'][k],3) for k in ('voc0','jsc0','ff0','pce0')))
        dedup.setdefault(key,d)
    return list(dedup.values())

def extract(main_pdf,supplements=()):
    result=[]
    for path in map(Path,[main_pdf,*supplements]):
        if not path.is_file() or path.suffix.lower()!='.pdf':continue
        with fitz.open(path) as doc:
            inspected=0
            for pi,page in enumerate(doc):
                t=page.get_text('text')
                if not re.search(r'(?i)\b(?:Supplementary\s+)?Table\s+[S]?\d+',t):continue
                if inspected>=8:break
                inspected+=1
                # Native selectable table rows take priority, OCR only when
                # invisible/vector-outlined cells prevent PDF text parsing.
                from eta_supplement_pv_recovery import _from_spatial,_from_pymupdf_tables
                if _from_spatial(doc,pi) or _from_pymupdf_tables(doc,pi):continue
                result+=_extract_page(page,path.name,pi)
    return result
