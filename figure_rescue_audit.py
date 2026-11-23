"""v18.6.7 audit of stability panels hidden by compound raster/PDF layouts.

Does NOT turn uncalibrated candidates into model-ready data. Writes deterministic
source/page evidence and panel-level OCR calibration for manual confirmation.
This prevents false 'no stability' statements when the graph reader misses plots.
"""
from pathlib import Path
import re, csv, json
import fitz

KEYWORDS=re.compile(r'(?i)\b(?:stability|degradation|ageing|aging|retention|MPP\s*tracking|MPPT|damp\s*heat|outdoor\s*test|light\s*stability|ISOS[-– ]?[LD]|shelf\s*life)\b')
FIG_RE=re.compile(r'(?i)\b(?:Fig\.?|Figure)\s*(?:S\s*)?\d{1,3}\b')


def scan_paper(main_pdf,supplements=(),preview_dir=None):
    """Return potential figure panels; no DOI-specific shortcuts or fake curves."""
    from general_curve_reader import ocr_page_words
    from panel_axis_fallback import detect
    records=[]
    for path in [main_pdf,*supplements]:
        path=Path(path)
        if path.suffix.lower()!='.pdf' or not path.is_file():continue
        try:
            with fitz.open(path) as doc:
                for i,page in enumerate(doc):
                    full=page.get_text('text')
                    if not KEYWORDS.search(full) or not FIG_RE.search(full):continue
                    # Inspect figure captions, not unrelated discussion text:
                    # otherwise full journal pages trigger expensive OCR.
                    anchors=list(re.finditer(r'(?im)^\s*(?:fig\.?|figure)\s*(?:S\s*)?\d{1,3}\s*[.|:]',full))
                    captions=[full[m.start():min(len(full),m.start()+1150)] for m in anchors]
                    if not any(KEYWORDS.search(c) for c in captions):continue
                    # A caption can be hundreds of characters after the Fig.
                    # token. Inspect the ENTIRE caption/page, not 140 characters.
                    paragraphs=[s.strip() for s in re.split(r'\n\s*\n',full) if s.strip()]
                    evidence=[c for c in captions if KEYWORDS.search(c)]
                    if not evidence and not re.search(r'(?i)\b(?:PCE|SPO|MPPT|normalized\s+PCE)\b',full):continue
                    if not page.get_images() and len(page.get_drawings())<25:continue
                    rec={'source_file':path.name,'source_page':i+1,
                        'record_type':'STABILITY_FIGURE_EVIDENCE',
                        'calibrated_panels':0,'status':'UNRESOLVED_AXIS_OR_LEGEND',
                        'caption_evidence':' '.join((evidence or [full[:600]])[:2])[:800],
                        'preview_file':'','panel_json':''}
                    if preview_dir:
                        folder=Path(preview_dir);folder.mkdir(parents=True,exist_ok=True)
                        target=folder/(path.stem[:70]+'_p'+str(i+1)+'.png')
                        if not target.exists():
                            pix=page.get_pixmap(matrix=fitz.Matrix(1.4,1.4),alpha=False)
                            pix.save(target)
                        rec['preview_file']=str(target)
                    try:
                        ocr=ocr_page_words(page,zoom=2.5)
                        panels=detect(page,ocr) if ocr else []
                        rec['calibrated_panels']=len(panels)
                        rec['panel_json']=json.dumps([{'x0':round(p['x0'],2),'x1':round(p['x1'],2),
                            'y0':round(p['y0'],2),'y1':round(p['y1'],2),
                            'max_time_h':p['x_max'],'y_max':p['y_max'],
                            'quality':'AXES_CALIBRATED_ONLY_TRACES_REQUIRE_LABEL_QC'} for p in panels])
                        if panels:rec['status']='AXES_CALIBRATED_LABELS_UNVERIFIED'
                    except Exception as ex:
                        rec['status']='PANEL_OCR_ERROR';rec['panel_json']=str(ex)[:250]
                    records.append(rec)
        except Exception as e:
            records.append({'source_file':path.name,'source_page':'','record_type':'READ_ERROR',
                            'calibrated_panels':0,'status':str(e)[:200],
                            'caption_evidence':'','preview_file':'','panel_json':''})
    return records


def save(records,path):
    if not records:return
    cols=('source_file','source_page','record_type','calibrated_panels','status',
          'caption_evidence','preview_file','panel_json')
    path=Path(path);path.parent.mkdir(parents=True,exist_ok=True)
    with path.open('w',newline='',encoding='utf-8-sig') as f:
        writer=csv.DictWriter(f,fieldnames=cols,extrasaction='ignore')
        writer.writeheader();writer.writerows(records)
