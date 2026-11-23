"""Non-destructive evidence-based enrichment of AUTO curve categories.

Never overwrite successfully extracted v15 metric values; use same-named
published category mentions and preserve the provenance in each dataset row.
"""
import re
from pathlib import Path
import fitz

NUM=r'(\d{1,3}(?:\.\d+)?)'

def _norm_label(t):
    return re.sub(r'[^a-z0-9]+','',str(t).lower())

def _parse_percentage_list(segment):
    # This is deliberately narrow: an explicit enumerated ordered list
    vals=re.findall(r'(\d{1,2}(?:\.\d+)?)\s*(?:±\s*\d+(?:\.\d+)?\s*)?%', segment)
    return [float(v) for v in vals]

def extract_caption_initial_pces(text):
    clean=' '.join(text.replace('\u2009',' ').replace('\u2002',' ').replace('\xa0',' ').split())
    out=[]
    for m in re.finditer(r'(?:The\s+)?initial\s+PCEs?\s+of\s+(.{1,130}?)\s+devices?\s+(?:is|are|were)\s+(.{1,160}?)(?:respectively\.|\.(?=\s+[a-zA-Z]))',clean,re.I):
        label_phrase=m.group(1).strip(); values=_parse_percentage_list(m.group(2))
        # Require the explicit ordered 'respectively' or exactly one category.
        if 'respectively' not in m.group(0).lower() and len(values)>1:continue
        labels=[x.strip() for x in re.split(r',\s*|\s+and\s+',label_phrase,flags=re.I)]
        labels=[re.sub(r'^(?:the|both|of)\s+','',v,flags=re.I).strip() for v in labels]
        if len(labels)==len(values) and len(labels)<=5:
            out.extend(zip(labels,values))
    return out

def _replace_missing(row,key,value,source):
    if value is None or row.get(key) not in (None,''):return False
    row[key]=value
    row[key+'_source']=source
    return True

def enrich_auto_rows(rows,main_pdf,supplements=()):
    pages=[]
    for filepath in [main_pdf,*[p for p in supplements if Path(p).suffix.lower()=='.pdf']]:
        try:
            with fitz.open(str(filepath)) as d:
                for i,p in enumerate(d):pages.append((Path(filepath).name,i+1,p.get_text('text')))
        except Exception:continue
    text='\n'.join(t for _,_,t in pages)
    # Only an explicitly named and ordered caption PV list may fill a PCE.
    for file,page,content in pages:
        for label,pce in extract_caption_initial_pces(content):
            for row in rows:
                if _norm_label(row['device_category'])==_norm_label(label):
                    _replace_missing(row,'pce0',pce,f'{file} p.{page}, initial PCE ordered category list')
    for file,page,content in pages:
        c=' '.join(content.replace('\u2009',' ').replace('\u2002',' ').split())
        # Two named, separately measured stability groups without the word
        # 'respectively' (e.g. damp-heat control and target).
        m=re.search(r'initial\s+PCE\s+of\s+control\s+and\s+target\s+devices?\s+is\s+(\d{1,2}(?:\.\d+)?)\s*(?:±\s*\d+(?:\.\d+)?)?%\s+and\s+(\d{1,2}(?:\.\d+)?)\s*(?:±\s*\d+(?:\.\d+)?)?%',c,re.I)
        if m:
            for row in rows:
                if _norm_label(row['device_category'])=='target':
                    _replace_missing(row,'pce0',float(m.group(2)),f'{file} p.{page}, target stability caption initial PCE')
    # Other information is applicable to the common absorber/architecture, not
    # an unrelated device JV champion, and is labeled as document-derived.
    main_text=' '.join(text[:130000].split())
    inverter=bool(re.search(r'\b(?:inverted|p[–-]i[–-]n)\b',main_text[:12000],re.I))
    normal=bool(re.search(r'\b(?:regular|n[–-]i[–-]p)\b',main_text[:12000],re.I))
    for row in rows:
        if not row.get('architecture'):
            if inverter and not normal:row['architecture']='p-i-n';row['architecture_source']='Main article text'
            if normal and not inverter:row['architecture']='n-i-p';row['architecture_source']='Main article text'
    # Do not infer a category PCE merely from nearby sentences: other
    # formulations (e.g., caffeine versus theophylline) may intervene.
    # Derive FF only if SAME category PCE, Voc and JV Jsc are available.
    for r in rows:
        p,v,j=(r.get(k) for k in ('pce0','voc0','jsc0'))
        try:
            p,v,j=float(p),float(v),float(j)
            if not r.get('ff0') and p>0 and v>0 and j>0:
                ff=100*p/(v*j)
                if 20<=ff<=100:_replace_missing(r,'ff0',round(ff,4),'Calculated from matched PCE, Voc and JV Jsc under 100 mW/cm2')
        except (ValueError,TypeError):pass
        try:
            jv=float(r['jsc0']);eqe=float(r['jsc_eqe_integrated'])
            if jv>0:
                r['jsc_difference_signed_pct']=round(100*(jv-eqe)/jv,5)
                r['jsc_difference_absolute_pct']=round(abs(100*(jv-eqe)/jv),5)
                r['jsc_eqe_to_jv_ratio']=round(eqe/jv,6)
        except (ValueError,TypeError,KeyError):pass
    recover_explicit_jv_prose(rows,main_pdf)
    recover_explicit_integrated_jsc(rows,main_pdf,supplements)
    renormalize_measured_curves(rows)
    return rows

# Published prose frequently contains an entire J-V parameter block for
# Control and Target in the same sentence rather than a machine-readable table.
def recover_explicit_jv_prose(rows,main_pdf):
    try:
        with fitz.open(str(main_pdf)) as d:
            pages=[(i+1,' '.join(p.get_text().split())) for i,p in enumerate(d)]
    except Exception:return rows
    jv_re=re.compile(r'(?:\b(?:the|a)\s+)?(control|reference|target|treated|modified)\s+devices?\s+(?:showed|yielded|exhibited|achieved|had)\s+(.{0,240}?)(?:[,;.](?=\s*(?:whereas|while|the\s+target|the\s+control|the\s+enhancement|in\s+comparison))|\.\s)',re.I)
    def get(measure,frag):
        rx={
           'voc0':r'(?:\bVOC\b|open[ -]?circuit voltage)\s*(?:\([^)]{1,12}\))?\s*(?:of|=|:)\s*(\d(?:\.\d+)?)',
           'jsc0':r'(?:\bJSC\b|short[ -]?circuit current(?:\s+density)?)\s*(?:\([^)]{1,12}\))?\s*(?:of|=|:)\s*(\d{1,2}(?:\.\d+)?)',
           'ff0':r'(?:\bFF\b|fill factor)\s*(?:\([^)]{1,12}\))?\s*(?:of|=|:)\s*(\d{1,2}(?:\.\d+)?)\s*%',
           'pce0':r'\bPCE\s*(?:of|=|:)\s*(\d{1,2}(?:\.\d+)?)\s*%'
        }[measure]
        m=re.search(rx,frag,re.I)
        return float(m.group(1)) if m else None
    for page,txt in pages:
        for m in jv_re.finditer(txt):
            name=m.group(1).lower(); chunk=m.group(2)
            actual={k:get(k,chunk) for k in ('pce0','voc0','jsc0','ff0')}
            if sum(v is not None for v in actual.values())<2:continue
            for r in rows:
                key=_norm_label(r['device_category'])
                role = (name in ('control','reference') and key in ('reference','control','ref')) or (name in ('target','modified','treated') and key in ('target','treated','modified'))
                if not role and name in ('target','modified','treated') and key=='theophylline':
                    role=bool(re.search(r'theophylline',txt[max(0,m.start()-350):m.end()+500],re.I))
                if not role:continue
                src=f'{Path(main_pdf).name} p.{page}, named J-V prose: {name} device'
                for k,v in actual.items():_replace_missing(r,k,v,src)
    # Explicit compact J-V block for a named treated/control PSC. Source text
    # example: 'target PSCs delivered a PCE of 23.8% ... VOC of 1.17 V,
    # JSC of 24.1 mA cm-2, and a fill factor (FF) of 0.842'.
    for page,txt in pages:
        for m in re.finditer(r'\b(control|reference|target)\s+PSCs?\s+(?:delivered|achieved|showed|exhibited)\s+(?:a\s+)?PCE\s+of\s+(\d{1,2}(?:\.\d+)?)%',txt,re.I):
            label=m.group(1).lower(); fragment=txt[m.start():m.start()+460]
            end=re.search(r'[.;]\s+(?:The|Figure|In|Our)\s',fragment)
            if end:fragment=fragment[:end.start()]
            voc=re.search(r'(?:VOC|open[ -]?circuit voltage)\)?\s+of\s+(\d(?:\.\d+)?)\s*V',fragment,re.I)
            jsc=re.search(r'(?:JSC|short[ -]?circuit current density)\)?\s+of\s+(\d{1,2}(?:\.\d+)?)',fragment,re.I)
            ff=re.search(r'(?:fill factor\s*\(FF\)|\bFF\b)\s+of\s+(0?\.\d+|\d{1,2}(?:\.\d+)?)',fragment,re.I)
            for row in rows:
                cat=_norm_label(row['device_category'])
                if not ((label in ('control','reference') and cat in ('reference','control','ref','pristine')) or
                        (label=='target' and cat in ('target','treated','modified'))):continue
                source=f'{Path(main_pdf).name} p.{page}, explicitly named PSC J-V prose'
                _replace_missing(row,'pce0',float(m.group(2)),source)
                if voc:_replace_missing(row,'voc0',float(voc.group(1)),source)
                if jsc:_replace_missing(row,'jsc0',float(jsc.group(1)),source)
                if ff:
                    value=float(ff.group(1));_replace_missing(row,'ff0',value*100 if value<=1 else value,source)
    for row in rows:
        if row.get('pce0') in (None,''):
            try:
                v,j,ff=(float(row[k]) for k in ('voc0','jsc0','ff0'))
                value=v*j*ff/100
                if 0<value<=35:_replace_missing(row,'pce0',round(value,4),'Calculated from category-matched Voc, Jsc, FF at 100 mW/cm2')
            except (TypeError,ValueError,KeyError):pass
    return rows


def renormalize_measured_curves(rows):
    """Use the author's category-specific initial PCE, not the first pixel,
    when the graph reports absolute PCE. Preserve all digitized raw values.
    """
    import json,numpy as np
    from perovskite_onefile_final import fit_every_trace,choose_decay_fit
    for r in rows:
        src=str(r.get('pce0_source') or '')
        if not ('initial PCE ordered category list' in src or 'caption MPP initial PCE' in src or 'stability caption initial PCE' in src):continue
        try:pce=float(r['pce0'])
        except (ValueError,TypeError,KeyError):continue
        if pce<=0:continue
        studies=json.loads(r.get('stability_data_json') or '[]')
        changed=False
        for study in studies:
            if study.get('normalization_source')!='RELATIVE_TO_FIRST_DIGITIZED_POINT':continue
            raw=np.asarray(study.get('raw_plot_y',[]),float)
            time=np.asarray(study.get('time_h',[]),float)
            if len(raw)!=len(time) or len(raw)<8:continue
            values=raw/pce*100
            if np.any(~np.isfinite(values)) or np.max(values)>200:continue
            study['retention_pct']=[round(float(v),4) for v in values]
            study['retention_last_pct']=round(float(values[-1]),4)
            study['normalization_source']='AUTHOR_REPORTED_CATEGORY_INITIAL_PCE'
            cross=None
            for j in range(1,len(values)):
                if values[j-1]>80>=values[j]:
                    cross=float(time[j-1]+(80-values[j-1])*(time[j]-time[j-1])/(values[j]-values[j-1]));break
            study['observed_t80_h']=cross
            study['fits']=fit_every_trace(time,values)
            study['fit_summary']=choose_decay_fit(study)
            changed=True
        if changed:
            primary=max(studies,key=lambda x:(x.get('protocol','').startswith('ISOS-L'),x.get('t_end_h',0)))
            r['stability_data_json']=json.dumps(studies,ensure_ascii=False)
            r['final_retention_pct']=primary['retention_last_pct']
            r['stability_endpoint']=str((primary['retention_last_pct'],primary['t_end_h']))
            r['observed_t80_h']=primary.get('observed_t80_h')
            r['t80_status']='OBSERVED' if primary.get('observed_t80_h') is not None else 'NOT_REACHED_DURING_MEASUREMENT'
            r['retention_reference']=primary['normalization_source']
            r.update(primary.get('fit_summary',{}))


def recover_explicit_integrated_jsc(rows, main_pdf, supplements=()):
    """Match published integrated-EQE Jsc sentences to stability formulations.

    Examples supported:
       integrated Jsc of 24.88 ... for target ... 23.92 ... control
       integrated Jsc of 25.17 and 25.74 ... control and target respectively
       integrated Jsc of 24.42 ... target ... whereas control ... 23.56
    No transfer from a generic Target to concentration-specific Target-1/2.
    """
    pattern=r'(?:integrated|EQE[- ]?integrated)\s+J\s*(?:SC|sc)'
    files=[main_pdf]+[x for x in supplements if Path(x).suffix.lower()=='.pdf']
    for file in files:
        try:
            with fitz.open(str(file)) as doc:
                pages=[(i+1,' '.join(p.get_text().split())) for i,p in enumerate(doc)]
        except Exception:
            continue
        for page, text in pages:
            for m in re.finditer(pattern,text,re.I):
                # Restrict matching to a short scientific statement rather than
                # harvesting unconnected numbers elsewhere on the page.
                span=text[m.start():m.start()+340]
                value_rx=r'(\d{1,2}(?:\.\d+)?)\s*mA\s*(?:cm|×\s*cm|/\s*cm)'
                mentions=[(float(mm.group(1)),mm.start()) for mm in re.finditer(value_rx,span,re.I)]
                labels=list(re.finditer(r'\b(control|reference|target)\b',span,re.I))
                assignments={}
                # Two ordered values followed by 'control and target ... respectively'.
                ordered=re.search(r'\b(?:control|reference)\s+and\s+(?:the\s+)?target\b.{0,30}\brespectively\b',span,re.I)
                if ordered and len(mentions)>=2 and mentions[1][1]<ordered.start():
                    assignments={'control':mentions[0][0],'target':mentions[1][0]}
                else:
                    # Each explicit 'for the target/control' must attach its
                    # immediately preceding current value in that phrase.
                    for v,pos in mentions:
                        phrase=span[pos:pos+115]
                        lm=re.search(r'\b(?:for|from)\s+(?:the\s+)?(control|reference|target)\b',phrase,re.I)
                        if lm:assignments[lm.group(1).lower()]=v
                    # Example: control device showed an integrated Jsc of ...
                    for label in ('control','reference','target'):
                        q=re.search(r'\b'+label+r'\s+device\s+showed\s+an?\s+integrated\s+J\s*(?:SC|sc)\s+of\s+(\d{1,2}(?:\.\d+)?)',span,re.I)
                        if q:assignments[label]=float(q.group(1))
                if not assignments:
                    continue
                for row in rows:
                    category=_norm_label(row.get('device_category',''))
                    role='control' if category in ('control','reference','pristine','ref') else 'target' if category in ('target','treated','modified') else None
                    if category=='theophylline' and ('theophylline' in text.lower()):role='target'
                    if role is None:continue
                    v=assignments.get(role, assignments.get('reference') if role=='control' else None)
                    if v is None or not (0<v<40):continue
                    _replace_missing(row,'jsc_eqe_integrated',v,f'{Path(file).name} p.{page}, explicitly paired integrated Jsc from text')
    for row in rows:
        try:
            jv=float(row['jsc0']); eqe=float(row['jsc_eqe_integrated'])
            if jv<=0:continue
            row['jsc_difference_signed_pct']=round(100*(jv-eqe)/jv,5)
            row['jsc_difference_absolute_pct']=round(100*abs(jv-eqe)/jv,5)
            row['jsc_eqe_to_jv_ratio']=round(eqe/jv,6)
        except (KeyError,TypeError,ValueError):pass
    return rows
