"""Evidence-first fallback from automatically calibrated PDF plots.

Only creates records for figures with a time axis, a stability y axis, and a
traced curve. The plot source is not restricted by DOI or filename.
"""
from pathlib import Path
import re, json
import numpy as np
import fitz
from general_curve_reader import find_plots
from pipeline import observed_t80
from curve_integrity_guard import (suspicious_legend, prune_duplicate_color_traces,
                                   choose_primary, trajectory_warnings, review_only_label)
from graph_endpoint_verification import annotate_publication, choose_primary_verified

REJECT_LEGENDS=re.compile(r'\bEQE\b|\d+(?:\.\d+)?\s*%|^.{1,2}$|^(?:OC|SC|MPP|open.?circuit|short.?circuit|-)$|confidence|error.?bar|95%|guide|fit|simulation|fitted|trend|nature|science|\bref\.?\s*\d+\b|damp.?heat|\b(?:SPO|PCE|Time|RH)\b|^\s*[0-9]|^\s*[a-z]\s*$',re.I)

def _clean_label(text,color):
    if not text:return 'UNASSIGNED_'+color.upper()+'_TRACE'
    x=re.sub(r'\s+',' ',text).strip()
    # Champion/average are statistics within one material formulation, not
    # separate device categories. Preserve the original legend in study data.
    x=re.sub(r'\s*\((?:average|avg\.?|champ\.?|champion)\)\s*$','',x,flags=re.I)
    if x.lower() in ('ctrl','ref','ref.','reference'):return 'Control'
    return x

def _retention(raw,label,ymx):
    arr=np.asarray(raw,float)
    if ('retention' in label.lower() or 'normalized' in label.lower()
        or re.search(r'P\s*norm|Norm\.|Norm\s*PCE',label,re.I)):
        if ymx<=1.5:return arr*100,'PLOTTED_NORMALIZED_0_TO_1'
        return arr,'PLOTTED_RETENTION_PERCENT'
    if ymx<=1.6 and np.median(arr)<1.3:return arr*100,'PLOTTED_NORMALIZED_0_TO_1'
    # Absolute PCE and SPO graph: use earliest digitized value as a transparent
    # proxy rather than fabricating an unplotted initial measurement.
    baseline=arr[0]
    if baseline==0:return None,None
    return arr/baseline*100,'RELATIVE_TO_FIRST_DIGITIZED_POINT'

def _protocol(text):
    if re.search(r'ISOS[- ]?D|dark|shelf|storage',text,re.I):return 'ISOS-D_OR_DARK_AGING'
    if re.search(r'ISOS[- ]?L|MPP|light|sun|illumination',text,re.I):return 'ISOS-L_OR_LIGHT_AGING'
    return 'STABILITY_PROTOCOL_UNSPECIFIED'

def auto_extract_publication(main_pdf,supplement_files=(),doi='',article_title='',overlay_root=None):
    """Digitize recognized stability plots across the matched publication bundle.

    Output: (category_rows, graph_diagnostics). If a legend is not identified,
    unassigned-colored trace still appears as a distinct experimental category,
    explicitly labelled UNASSIGNED. No fabricated PV metrics.
    """
    studies_by_category={};diagnostics=[]
    files=[Path(main_pdf)]+[Path(p) for p in supplement_files if Path(p).suffix.lower()=='.pdf']
    for source in files:
        try:plots=find_plots(source)
        except Exception as e:
            diagnostics.append({'file':source.name,'record_type':'AUTO_GRAPH_ERROR','reason':str(e)})
            continue
        if not plots:
            # Make the crucial external OCR binary dependency explicit when
            # a preprint stores all plots as raster-only figure pages.
            with fitz.open(str(source)) as probe:
                has_raster_figures=any(len(pg.get_text('words'))<=12 and bool(pg.get_images())
                    for pg in probe)
            if has_raster_figures:
                try:
                    import pytesseract
                    from general_curve_reader import ocr_page_words
                    import os
                    if os.name=='nt':
                        exe=Path(r'C:\Program Files\Tesseract-OCR\tesseract.exe')
                        if exe.exists():pytesseract.pytesseract.tesseract_cmd=str(exe)
                    pytesseract.get_tesseract_version()
                except Exception as ex:
                    diagnostics.append({'file':source.name,
                        'record_type':'RASTER_ONLY_FIGURE_OCR_UNAVAILABLE',
                        'reason':'Install Tesseract OCR executable, not only pip pytesseract: '+str(ex)[:180]})
        if overlay_root and plots:
            from general_curve_reader import write_overlays
            write_overlays(source,plots,overlay_root)
        with fitz.open(str(source)) as document:
            page_context={}
            for plot in plots:
                p=document[plot['page']-1]
                caption=p.get_text('text',clip=fitz.Rect(0,min(p.rect.height,plot['y1']+12),p.rect.width,min(p.rect.height,plot['y1']+270)))
                page_context[(plot['page'],plot['figure_local_index'])]=caption[:2500]
        for plot in plots:
            page,plot_idx=plot['page'],plot['figure_local_index']
            caption=page_context.get((page,plot_idx),'')
            # Page may contain just an embedded image with no PDF-text
            # caption. OCR evidence is specific to this calibrated panel.
            if not caption.strip() and plot.get('ocr_panel_evidence'):
                caption=plot['ocr_panel_evidence']
            with fitz.open(str(source)) as localdoc:
                localpage=localdoc[page-1]
                local_text=localpage.get_text('text',clip=fitz.Rect(max(0,plot['x0']-7),plot['y0'],plot['x1']+7,plot['y1'])).strip()
            protocol_text=local_text if re.search(r'damp.?heat|ISOS|\bMPP\b|maximum power|one.?sun|storage|dark',local_text,re.I) else caption
            for trace in plot['series']:
                name=_clean_label(trace.get('label'),trace['color'])
                # Never treat an unassigned colored curve as a named device category.
                # Keep the graph in discovery evidence instead of corrupting ML rows.
                if name.startswith('UNASSIGNED_') or suspicious_legend(name) or REJECT_LEGENDS.search(name) or (trace['color']=='gray' and not trace.get('label')) or (trace['color']=='black' and (name.casefold() not in ('reference','control','pristine') or np.std(np.asarray(trace['points'])[:,1])<0.04*plot['y_max'])):continue
                ar=np.asarray(trace['points'],float)
                tt=ar[:,0];raw=ar[:,1]
                if len(tt)<8 or tt[-1]-tt[0]<5:continue
                # Reject a clipped fragment of an absolute-PCE decay as a full
                # device trace. This otherwise fabricates initial retention and T80.
                if (plot['y_max']<=35 and tt[0]>0.05*plot['x_max']
                    and (tt[-1]-tt[0])<0.40*plot['x_max']
                    and raw[0]<0.50*plot['y_max']):
                    diagnostics.append({'file':source.name,'record_type':'CLIPPED_TRACE_REJECTED',
                         'page':page,'reason':f'{name}: plotted initial section not detected'})
                    continue
                ret,norm=_retention(raw,plot['y_axis_label'],plot['y_max'])
                if ret is None or np.any(~np.isfinite(ret)) or np.any(ret<0) or np.any(ret>200):continue
                # Evidence is retained at plot/legend resolution; independent graphs
                # with the same device category become separate experiments.
                key=(plot.get('plot_architecture') or 'unspecified')+'::'+re.sub(r'[^a-z0-9]','',name.lower())
                trace_flags=[]
                if trace['color']=='gray':trace_flags.append('PALE_GRAY_TRACE_REQUIRES_VISUAL_QC')
                if plot.get('time_axis_inferred_from_numeric_ticks'):
                    trace_flags.append('HOUR_UNIT_INFERRED_FROM_CAPTION_REQUIRES_QC')
                if trace.get('label_source')=='MATCHED_SAME_COLOR_SIBLING_PANEL':
                    trace_flags.append('CATEGORY_LABEL_BORROWED_FROM_SIBLING_PANEL')
                study={'figure':f'AUTO p{page}, chart {plot_idx}','page':page,
                    'pdf':source.name,'protocol':(plot.get('panel_protocol') or ('DAMP_HEAT' if re.search(r'damp.?heat|85\s*°?\s*C\s*[/,]\s*85\s*%',protocol_text[:400]+' '+plot['y_axis_label'],re.I) else _protocol(plot['y_axis_label']+' '+protocol_text))),
                    'measurement_mode':'MPP' if re.search(r'maximum.power.point|\bMPP\b',caption,re.I) else 'UNSPECIFIED',
                    'conditions':' '.join(caption.split())[:250] if caption else 'Not machine-extracted',
                    'metric':plot['y_axis_label'],'n_points':len(tt),
                    't_start_h':round(float(tt[0]),3),'t_end_h':round(float(tt[-1]),3),
                    'retention_last_pct':round(float(ret[-1]),3),
                    'time_h':[round(float(x),3) for x in tt],
                    'retention_pct':[round(float(x),3) for x in ret],
                    'raw_plot_y':[round(float(x),4) for x in raw],
                    'normalization_source':norm,
                    'observed_t80_h':observed_t80(tt,ret),
                    'color':trace['color'],'legend_evidence':trace.get('label'),
                    'axis_calibration':plot.get('source') or 'PDF_NUMERIC_TICKS',
                    'figure_source_type':('AUTOMATIC_RASTER_IMAGE_FIGURE'
                        if plot.get('source') in ('OCR_IMAGE_AXIS','OCR_NUMERIC_TIME_AXIS')
                        else 'AUTOMATIC_UNSEEN_LAYOUT'),
                    'subplot_label':plot.get('subplot_label'),
                    'panel_caption':plot.get('panel_caption'),
                    'plot_architecture':plot.get('plot_architecture'),
                    'axis_inferred':bool(plot.get('time_axis_inferred_from_numeric_ticks')),
                    'digitization_quality_flags':trace_flags}
                try:
                    from perovskite_onefile_final import fit_every_trace,choose_decay_fit
                    study['fits']=fit_every_trace(tt,ret)
                    study['fit_summary']=choose_decay_fit(study)
                except Exception as e:study['fit_error']=str(e)
                study['digitization_quality_flags']=trajectory_warnings(study)
                studies_by_category.setdefault(key,{'label':name,'studies':[],
                    'architecture':plot.get('plot_architecture')})['studies'].append(study)
                diagnostics.append({'file':source.name,'record_type':'AUTO_DIGITIZED_TRACE',
                    'page':page,'figure_labels':f'chart {plot_idx}',
                    'reason':f'{name}: {len(tt)} points; {norm}; source={source.name}'})
    # v18.6.2: independently verify digitized points against published,
    # category-matched time/retention checkpoints. No points or fits are edited.
    verification_counts=annotate_publication(studies_by_category,main_pdf,supplement_files)
    diagnostics.append({'file':Path(main_pdf).name,
                        'record_type':'SOURCE_ENDPOINT_VERIFICATION_SUMMARY',
                        'reason':json.dumps(verification_counts,sort_keys=True)})
    for obj in studies_by_category.values():
        for study in obj['studies']:
            quality=study.get('published_endpoint_verification') or {}
            if quality.get('status') in ('MISMATCH','REPORTED_CHECKPOINT_OUTSIDE_DIGITIZED_WINDOW','BELOW_80_NOT_NEAR_T80','T80_RETENTION_MISMATCH'):
                diagnostics.append({'file':study.get('pdf'),
                    'record_type':'SOURCE_ENDPOINT_VERIFICATION_'+quality['status'],
                    'page':study.get('page'), 'figure_labels':study.get('figure'),
                    'reason':json.dumps({'device_category':obj['label'],**quality},ensure_ascii=False)})
    rows=[]
    for key,obj in studies_by_category.items():
        studies=obj['studies']
        studies,removed=prune_duplicate_color_traces(studies)
        for event in removed:
            diagnostics.append({'file':event['source'],'record_type':event['status'],
                                'reason':event['reason'],'figure_labels':event['figure']})
        # Evidence stays attached to each raw trace. The selected model and fit
        # parameters are never altered by the paper-reported T80 comparison.
        status=review_only_label(obj['label'],studies)
        if status:
            diagnostics.append({'file':Path(main_pdf).name,'record_type':'CATEGORY_HELD_FOR_REVIEW',
                                'reason':f'{obj["label"]}: {status}'})
            continue
        studies.sort(key=lambda s:(s['pdf'],s['page'],s['figure']))
        preferred=choose_primary_verified(studies)
        result={'doi':doi,'article_title':article_title,'pdf_file':Path(main_pdf).name,
            'device_category':obj['label'],'category_id':key,
            'architecture':obj.get('architecture'),
            'extraction_method':'GENERAL_AUTO_DIGITIZED','stability_experiment_count':len(studies),
            'primary_stability_figure':preferred['figure'],
            'primary_stability_protocol':preferred['protocol'],
            'primary_stability_conditions':preferred['conditions'],
            'final_retention_pct':preferred['retention_last_pct'],
            'final_aging_time_h':preferred['t_end_h'],
            'stability_endpoint':str((preferred['retention_last_pct'],preferred['t_end_h'])),
            'stability_data_json':json.dumps(studies,ensure_ascii=False),
            'observed_t80_h':preferred['observed_t80_h'],
            't80_status':'OBSERVED' if preferred['observed_t80_h'] is not None else 'NOT_REACHED_DURING_MEASUREMENT',
            'published_t80_h':((preferred.get('published_t80_comparison') or {}).get('reported_time_h')
              if (preferred.get('published_t80_comparison') or {}).get('evidence_type')=='EXPLICIT_PUBLISHED_T80' else None),
            'retention_reference':preferred['normalization_source'],
            'stability_measurement':preferred['metric'],
            'published_endpoint_verification_json':json.dumps(preferred.get('published_endpoint_verification',{}),ensure_ascii=False),
            'published_endpoint_verification_status':(preferred.get('published_endpoint_verification') or {}).get('status'),
            'published_t80_comparison_json':json.dumps(preferred.get('published_t80_comparison',{}),ensure_ascii=False),
            'published_t80_comparison_status':(preferred.get('published_t80_comparison') or {}).get('status'),
            'reported_t80_source_type':(preferred.get('published_t80_comparison') or {}).get('evidence_type'),
            'reported_t80_from_source_h':((preferred.get('published_t80_comparison') or {}).get('reported_time_h')
              if (preferred.get('published_t80_comparison') or {}).get('evidence_type')=='EXPLICIT_PUBLISHED_T80' else None),
            'reported_near_t80_endpoint_h':((preferred.get('published_t80_comparison') or {}).get('reported_time_h')
              if (preferred.get('published_t80_comparison') or {}).get('evidence_type')=='NEAR_80_RETENTION_ENDPOINT' else None),
            'selected_fit_vs_source_t80_abs_error_h':(preferred.get('published_t80_comparison') or {}).get('absolute_error_h'),
            'selected_fit_vs_source_t80_rel_error_pct':(preferred.get('published_t80_comparison') or {}).get('relative_absolute_error_pct'),
            'selected_fit_vs_source_t80_discrepancy_flag':(preferred.get('published_t80_comparison') or {}).get('fit_discrepancy_flag_20pct'),
            'digitized_near_t80_one_sided':(preferred.get('published_t80_comparison') or {}).get('digitized_near_t80_one_sided'),
            'final_retention_near_t80_one_sided':(80.0<=float(preferred['retention_last_pct'])<81.0),
            'final_retention_near_t80_note':('NEAR_T80_ENDPOINT_NOT_PROVEN_CROSSING'
               if 80.0<=float(preferred['retention_last_pct'])<81.0 else
               'BELOW_80_NOT_WITHIN_ONE_SIDED_TOLERANCE' if float(preferred['retention_last_pct'])<80.0
               else 'ABOVE_ONE_SIDED_NEAR_80_WINDOW'),
            'fit_traces_successful':sum('fit_summary' in s for s in studies),
            'digitization_quality_flags_json':json.dumps(sorted(set(
                flag for st in studies for flag in st.get('digitization_quality_flags',[])))),
            'digitization_quality_status':('VISUAL_REVIEW_RECOMMENDED' if any(
                st.get('digitization_quality_flags') for st in studies) else 'AUTO_NO_SPECIFIC_WARNING'),
            'category_legend_source':'PLOT_LEGEND' if not obj['label'].startswith('UNASSIGNED') else 'UNRESOLVED_LEGEND',
            'pce0':None,'voc0':None,'jsc0':None,'ff0':None,'jsc_eqe_integrated':None,
            'jsc_difference_signed_pct':None,'jsc_difference_absolute_pct':None,
            'jsc_eqe_to_jv_ratio':None}
        if 'fit_summary' in preferred:result.update(preferred['fit_summary'])
        rows.append(result)
    return rows,diagnostics
