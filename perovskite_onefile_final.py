"""Stability-first ML dataset: one row per paper x stability-tested formulation.

Run in VS Code (folder picker) or: python perovskite_onefile_final.py --pdf-dir PATH
Uses only file-derived measurement evidence and source-linked categories.
No DOI-specific measurement defaults or calibrated figure-coordinate catalogs.
The stability measurements are linked at published device-CATEGORY level, not cell serial number.
"""
from __future__ import annotations
import argparse, json, math, tempfile, sys, contextlib, io
from pathlib import Path
import pandas as pd
import numpy as np
from openpyxl import Workbook
from openpyxl.styles import Font, PatternFill, Alignment
from openpyxl.utils import get_column_letter


def blank(v):
    if v is None: return None
    try:
        if pd.isna(v):return None
    except (TypeError,ValueError):pass
    if isinstance(v,(float,np.floating)) and not math.isfinite(float(v)):return None
    if isinstance(v,(np.integer,np.floating)):return v.item()
    if isinstance(v,bool):return v
    if isinstance(v,str) and not v.strip():return None
    return v


def category(doi,label):
    from hysteresis_recovery import canonical_key
    return canonical_key(str(label))




def choose_decay_fit(study):
    """Compare post-peak fits by AIC, with RMSE as tie-breaker.

    k=2*n-1 because amplitudes sum to one. Compute RSS from fitted
    retention percentages and same post-peak points for both candidates.
    A missing T80 does not disqualify an otherwise valid decay fit.
    """
    from scipy.optimize import brentq
    t=np.asarray(study.get('time_h',[]),float)
    y=np.asarray(study.get('retention_pct',[]),float)
    if not len(t) or len(t)!=len(y):return {}
    imax=int(np.argmax(y));tp=float(t[imax]);rp=float(y[imax]); dt=t[imax:]-tp; target=y[imax:]
    nobs=len(dt)
    if nobs < 8 or dt[-1] <= 0 or not np.isfinite(rp) or rp<=0:return {}
    out={};possible=[]
    for f in study.get('fits',[]):
        n=blank(f.get('n_components'))
        if n not in (2,3) or f.get('status')!='ok':continue
        amps=np.array([float(f.get(f'amplitude_{j}')) for j in range(1,n+1)])
        rates=np.array([float(f.get(f'rate_{j}_per_h')) for j in range(1,n+1)])
        if not np.isfinite(amps).all() or not np.isfinite(rates).all() or not np.isclose(amps.sum(),1,atol=.002):continue
        model=lambda x:float(rp*np.sum(amps*np.exp(-rates*x)))
        pred=np.array([model(x) for x in dt]);rss=float(np.sum((target-pred)**2))
        k=2*n-1
        rmse=float(np.sqrt(rss/nobs))
        aic=float(nobs*np.log(max(rss/nobs,1e-24))+2*k)
        aicc=float(aic+2*k*(k+1)/(nobs-k-1)) if nobs>k+1 else None
        # Solve the SELECTED mathematical curve at 80% of the ORIGINAL
        # retention reference, expanding the root bracket as necessary.
        # No 10x-measurement-duration cutoff or unreported T80 is invented.
        from fit_lifetime_solver import solve_fitted_t80
        t80,status=solve_fitted_t80(tp,rp,amps,rates,
                                    threshold=80.0,measured_last_h=float(t[-1]))
        out[f'fit{n}_rmse_pct']=rmse
        out[f'fit{n}_aic']=aic
        out[f'fit{n}_aicc']=aicc
        out[f'fit{n}_t80_h']=t80
        out[f'fit{n}_t80_status']=status
        ref=blank(study.get('published_t80_h'))
        out[f'fit{n}_published_t80_error_h']=abs(float(t80)-float(ref)) if t80 is not None and ref is not None else None
        out[f'fit{n}_published_t80_error_pct']=100*abs(float(t80)-float(ref))/float(ref) if t80 is not None and ref is not None and float(ref)>0 else None
        out[f'fit{n}_amplitudes_json']=json.dumps(amps.tolist())
        possible.append((aic,rmse,n))
    if possible:
        _,_,best=min(possible)
        out['selected_fit_components']=best
        out['selected_fit_aic']=out[f'fit{best}_aic']
        out['selected_fit_rmse_pct']=out[f'fit{best}_rmse_pct']
        out['selected_fit_t80_h']=out[f'fit{best}_t80_h']
        out['selected_fit_t80_status']=out[f'fit{best}_t80_status']
        # Numerical existence is different from predictive trustworthiness.
        # Extremely long extrapolations remain numeric (as requested),
        # but must not be confused with demonstrated experimental lifetime.
        if out['selected_fit_t80_h'] is not None:
            measured_end=float(t[-1]); measured_span=max(measured_end-tp,1e-9)
            ratio=(out['selected_fit_t80_h']-tp)/measured_span
            out['selected_fit_t80_to_measurement_span_ratio']=float(ratio)
            out['selected_fit_t80_interpretation']=(
                'MEASURED_RANGE_FITTED' if out['selected_fit_t80_h']<=measured_end
                else 'EXTREME_EXTRAPOLATION_UNRELIABLE' if ratio>10
                else 'EXTRAPOLATED_MODEL_ESTIMATE')
        out['selected_fit_criterion']='MIN_AIC_THEN_RMSE'
        out['selected_fit_published_t80_error_h']=out.get(f'fit{best}_published_t80_error_h')
        out['selected_fit_published_t80_error_pct']=out.get(f'fit{best}_published_t80_error_pct')
        ref=blank(study.get('published_t80_h'))
        if ref is not None:
            errs={n:out.get(f'fit{n}_published_t80_error_h') for n in (2,3)}
            eligible={n:v for n,v in errs.items() if v is not None}
            out['closest_model_to_published_t80']=min(eligible,key=eligible.get) if eligible else None
            out['published_t80_h']=ref
    return out


def fit_every_trace(time_h,retention_pct):
    """Refit each digitized category trace with 2 and 3 exponentials, T80 or not.

    Rates>=0, positive amplitudes normalized to one. Fits are after own peak.
    AIC uses equal-support post-peak residuals. No reported T80 information in fitting.
    """
    from pipeline import fit_post_peak
    res=[]
    for n in (2,3):
        try:
            f=fit_post_peak(time_h,retention_pct,n_components=n)
            rec={'n_components':n,'status':'ok','rmse_post_peak_pct':f.rmse_post_peak_pct,
                 'peak_time_h':f.peak_time_h,'peak_retention_pct':f.peak_retention_pct}
            for j,(a,k) in enumerate(zip(f.amplitudes,f.rates_per_h),1):
                rec[f'amplitude_{j}']=a;rec[f'rate_{j}_per_h']=k
            res.append(rec)
        except (ValueError,OverflowError,FloatingPointError) as exc:
            res.append({'n_components':n,'status':'not_fit','reason':str(exc)})
    return res


def _excel_clean(value):
    """PDF text can contain ASCII control codes invalid in XLSX XML."""
    import re
    if isinstance(value,str):
        value=re.sub(r'[\x00-\x08\x0b\x0c\x0e-\x1f\ufffe\uffff]', '', value)
        # Excel cells have a 32,767-character hard limit. Keep full audit in CSV.
        if len(value)>32760:value=value[:32725]+' [TRUNCATED: SEE SOURCE CSV]'
        if value.startswith(('=','+','-','@')) and not value.lstrip('+-').replace('.','',1).isdigit():
            value="'"+value
    return value


def _append_audit_sheet(wb, name, path):
    """Never crash at the worksheet boundary on invalid publisher glyphs."""
    from csv import DictReader
    if not Path(path).exists():return
    with Path(path).open(encoding='utf-8-sig',newline='') as f:
        reader=DictReader(f)
        if not reader.fieldnames:return
        ws=wb.create_sheet(name[:31]);ws.append([_excel_clean(x) for x in reader.fieldnames])
        for rec in reader:
            ws.append([_excel_clean(rec.get(k,'')) for k in reader.fieldnames])
        ws.freeze_panes='A2';ws.auto_filter.ref=ws.dimensions
        for c in ws[1]:
            c.fill=PatternFill('solid',fgColor='15324F')
            c.font=Font(color='FFFFFF',bold=True)
        for col in ws.columns:
            ws.column_dimensions[col[0].column_letter].width=27 if col[0].column<5 else 42


def write_excel(rows, destination, diagnostics=None):
    """Single readable XLSX using openpyxl (Windows pip-installable).

    One row per stability-tested category; no nonstandard spreadsheet dependency.
    """
    if not rows:
        raise RuntimeError('No stability-tested device-category records found')
    headers = list(dict.fromkeys(k for r in rows for k in r))
    wb = Workbook()
    ws = wb.active
    ws.title = 'ML_dataset'
    ws.append(headers)
    long_columns = {'pdf_file': 37, 'doi': 33, 'device_category': 30,
                    'stability_data_json': 45, 'reported_t80_validation_json': 40,
                    'absorber_formula': 37, 'primary_stability_conditions': 42}
    numeric_columns = {'pce0', 'voc0', 'jsc0', 'ff0', 'jsc_eqe_integrated',
                       'fit2_t80_h', 'fit3_t80_h', 'selected_fit_t80_h',
                       'fit2_published_t80_error_h', 'fit3_published_t80_error_h',
                       'selected_fit_aic', 'selected_fit_rmse_pct'}
    for row in rows:
        values = [_excel_clean(blank(row.get(h))) for h in headers]
        for i,value in enumerate(values):
            if isinstance(value,str) and len(value)>32000:
                if headers[i]=='stability_data_json':
                    import zlib,base64
                    values[i]='zlib64:'+base64.b64encode(zlib.compress(value.encode('utf8'),9)).decode('ascii')
                else:
                    # Metadata must not be silently truncated.
                    raise ValueError('Excel metadata cell too long: '+headers[i])
        ws.append(values)
    ws.freeze_panes = 'A2'
    ws.auto_filter.ref = ws.dimensions
    ws.row_dimensions[1].height = 29
    for cell in ws[1]:
        cell.fill = PatternFill(fill_type='solid', fgColor='15324F')
        cell.font = Font(bold=True, color='FFFFFF')
        cell.alignment = Alignment(horizontal='center', vertical='center')
    for col_idx, name in enumerate(headers, 1):
        ws.column_dimensions[get_column_letter(col_idx)].width = long_columns.get(name, 17)
        if name in numeric_columns:
            for row_idx in range(2, ws.max_row + 1):
                ws.cell(row_idx, col_idx).number_format = '0.000'
    for row_idx in range(2, ws.max_row + 1):
        ws.row_dimensions[row_idx].height = 23
    # Keep exactly one worksheet as requested; unmatched-file diagnostics
    # are reported in the console and preserved as source fields per row.
    destination = Path(destination)
    destination.parent.mkdir(parents=True, exist_ok=True)
    for tab, filename in (
        ('PV_table_candidates','figure_jv_table_candidates.csv'),
        ('Stability_marker_review','marker_only_stability_candidates.csv'),
        ('Continuous_MPP','continuous_mpp_candidates.csv'),
        ('Stability_figures','stability_all_figure_pages.csv'),
        ('Reported_endpoints','stability_reported_endpoints.csv'),
        ('Stability_figure_review','stability_figures_needing_review.csv'),
        ('Stability_points','extracted_stability_points.csv'),
        ('PV_measured_values','extracted_pv_metrics.csv'),
        ('Extraction_coverage','extraction_coverage.csv'),
        ('HI_scan_sources','hysteresis_scan_candidates.csv')):
        _append_audit_sheet(wb,tab,destination.parent/filename)
    wb.save(str(destination))
    return len(rows)


def is_supplementary_pdf(path):
    """Avoid treating SI as an independent experimental paper."""
    import re
    name=path.stem.lower()
    return bool(re.search(r'(?:_sm|_si|_supp|_supporting|_supplementary)$',name))


def find_supplement(main_pdf, all_pdfs):
    """Only pair a clearly matching DOI manuscript token; never pair arbitrary SI."""
    import re
    stem=main_pdf.stem.lower()
    tokens=re.findall(r'[a-z]{3}[0-9]{3,}',stem)
    if not tokens:return None
    matches=[q for q in all_pdfs if q!=main_pdf and is_supplementary_pdf(q)
             and any(token in q.stem.lower() for token in tokens)]
    return matches[0] if len(matches)==1 else None


def derive_missing_ff(rows):
    """Compute FF% if, and only if, PCE%, Voc[V], Jsc[mA/cm2] exist.

    AM1.5G Pin=100 mW/cm2; FF_fraction=PCE_percent/(Voc*Jsc).
    Keep *reported* FF instead of overwriting it. Annotate formula and source.
    Check physical plausibility; don't fabricate missing Voc/Jsc from unrelated curves.
    """
    for row in rows:
        p=blank(row.get('pce_jv_pct'));v=blank(row.get('voc0'));j=blank(row.get('jsc0'))
        if p is None:
            source=str(row.get('metric_source') or '').lower()
            if ('j-v' in source or 'jv' in source or 'photovoltaic table' in source) and not ('stability' in source or 'aging' in source):
                p=blank(row.get('pce0'))
        row['ff_calculation_basis']='AM1.5G 100 mW/cm2'
        row['ff0_source']='Reported / digitized from matching category' if blank(row.get('ff0')) is not None else None
        if None in (p,v,j) or not all(isinstance(x,(int,float)) for x in (p,v,j)) or v<=0 or j<=0:
            row['ff_calculated_pct']=None
            row['ff_calculation_status']='INSUFFICIENT_MATCHED_JV_INPUTS' if blank(row.get('ff0')) is None else 'REPORTED_FF_AVAILABLE'
            continue
        ff=float(100*p/(v*j))
        row['ff_calculated_pct']=round(ff,4) if 0<=ff<=100 else None
        if not (0<=ff<=100):
            row['ff_calculation_status']='PV_INPUTS_INCONSISTENT_FF_GT_100'
        elif blank(row.get('ff0')) is None:
            row['ff0']=round(ff,4)
            row['ff0_source']='Calculated: 100*matched_JV_PCE/(voc0*jsc0); original stability initial PCE never substituted'
            row['ff_calculation_status']='CALCULATED_FROM_AVAILABLE_PV'
        else:
            row['ff_calculation_status']='REPORTED_FF_RETAINED; CALCULATED_COMPARISON_AVAILABLE'
        # Category-level literature measurements can come from different representative cells.
        # Keep the discrepancy between FF calculated from PCE0 and reported JV-FF visible.
        if blank(row.get('ff0')) is not None and row['ff_calculated_pct'] is not None:
            row['ff_reported_vs_calculated_diff_pctpts']=round(float(row['ff0'])-row['ff_calculated_pct'],4)


def main():
    parser=argparse.ArgumentParser(description='One XLSX with one row per publication-device category')
    parser.add_argument('--pdf-dir');parser.add_argument('--output')
    args=parser.parse_args()
    if args.pdf_dir:folder=Path(args.pdf_dir)
    else:
        import tkinter as tk
        from tkinter import filedialog
        root=tk.Tk();root.withdraw();root.attributes('-topmost',True)
        chosen=filedialog.askdirectory(title='Choose folder with perovskite research PDFs');root.destroy()
        if not chosen:return
        folder=Path(chosen)
    dest=Path(args.output) if args.output else folder.parent/'Perovskite_ML_dataset.xlsx'
    from tempfile import TemporaryDirectory
    from publication_matcher import discover
    documents,associations,unmatched_supplements,metadata_errors=discover(folder)
    # Word figures need real page geometry: convert each DOCX to an ephemeral
    # PDF when a supported converter is available. Never alter originals.
    from tempfile import mkdtemp
    _word_bridge_dir=Path(mkdtemp(prefix='psc_word_figure_bridge_'))
    try:
        from word_pdf_bridge import replace_paths
        associations,word_bridge_warnings=replace_paths(documents,associations,_word_bridge_dir)
    except Exception as exc:
        word_bridge_warnings=[('all DOCX',str(exc))]
    for file,reason in word_bridge_warnings:
        print(f'  WORD_FIGURE_CONVERSION_UNAVAILABLE: {file}: {reason}',flush=True)

    # DOCX converted through the optional bridge are included in PDF processing.
    pdfs=[doc for doc in documents if doc.path.suffix.lower()=='.pdf']
    word_mains=[doc.path.name for doc in documents if doc.path.suffix.lower()=='.docx']
    if metadata_errors:
        for name,reason in metadata_errors:print(f'  DOCUMENT_READ_ERROR: {name}: {reason}',flush=True)
    if unmatched_supplements:
        for name,reason,choices in unmatched_supplements:
            print(f'  SUPPLEMENT_NOT_PAIRED: {name}: {reason}; candidates={choices}',flush=True)
    if word_mains:
        print(f'  WORD_MAIN_ARTICLES_METADATA_ONLY (PDF needed for graph digitization): {word_mains}',flush=True)
    if not pdfs:raise ValueError('No main research PDF identified; only supplementary or Word files were found')
    all_rows=[];skipped=[];diagnostics=[];rescue_records=[];figure_pv_candidates=[]
    for file,reason in word_bridge_warnings:diagnostics.append(dict(file=file,record_type='WORD_FIGURE_CONVERSION_UNAVAILABLE',reason=reason))
    all_marker_candidates=[];all_figures=[];all_reported_endpoints=[];all_continuous_candidates=[];all_hi_candidates=[]
    from figure_discovery import stability_figures
    for name,reason,choices in unmatched_supplements:
        diagnostics.append(dict(file=name,record_type='UNMATCHED_SUPPLEMENT',reason=reason,
                                match_evidence=json.dumps(choices)))
    for name,reason in metadata_errors:
        diagnostics.append(dict(file=name,record_type='UNREADABLE_DOCUMENT',reason=reason))
    for name in word_mains:
        diagnostics.append(dict(file=name,record_type='WORD_MAIN_TEXT_ONLY',reason='Current graph digitizer requires PDF main article'))
    from publication_validation import review_article,validate_publication
    for index,publication in enumerate(pdfs,1):
        pdf=publication.path
        if review_article(publication):
            print(f'[{index}/{len(pdfs)}] REVIEW_SKIPPED: {pdf.name} (front-matter article type)',flush=True)
            diagnostics.append(dict(file=pdf.name,record_type='REVIEW_SKIPPED',reason='Explicit review article; no original device measurements attributed'))
            continue
        calibrated_error=''
        matched=associations.get(str(pdf),[])
        print(f'[{index}/{len(pdfs)}] PROCESSING {pdf.name} | title={publication.title[:75]}',flush=True)
        for supplementary,evidence,score in matched:
            print(f'  MATCHED_SUPPLEMENT: {supplementary.path.name} | score={score} | evidence={evidence}',flush=True)
            diagnostics.append(dict(file=supplementary.path.name,record_type='MATCHED_SUPPLEMENT',
                                    matching_main=pdf.name,match_evidence=json.dumps(evidence),reason=f'match score={score}'))
        with TemporaryDirectory(prefix='psc_onepaper_') as tmp:
            # Discovery and matching are based on supplied document contents.
            # NO paper-specific calibrated recipes, DOI catalogs, coordinates or
            # handwritten measurements may enter the extraction path.
            canonical_doi=publication.doi or next((d.doi for d,ev,sc in matched if d.doi),'')
            if canonical_doi and not __import__('re').fullmatch(r'10\.\d{4,9}/[A-Za-z0-9._;()/:-]+',canonical_doi):
                canonical_doi=''
            target=Path(tmp)/'generic_extraction'
            paper_rows=[]
            from generic_stability_pipeline import auto_extract_publication
            supplement_paths=[d.path for d,ev,sc in matched]
            try:
                paper_rows,auto_diag=auto_extract_publication(
                    pdf,supplement_paths,canonical_doi or '',publication.title or '',
                    overlay_root=Path(tmp)/'auto_validation_overlays')
                diagnostics.extend(auto_diag)
                if paper_rows:
                    print(f'  GENERAL_AUTO_DIGITIZED: {len(paper_rows)} categories; '
                          f'{sum(r.get("stability_experiment_count",0) for r in paper_rows)} traces',flush=True)
            except Exception as exc:
                diagnostics.append(dict(file=pdf.name,record_type='GENERAL_READER_ERROR',reason=str(exc)))
                paper_rows=[]
            # The later native-axis and marker readers independently inspect
            # this document, even when the first reader produces no categories.
            if not paper_rows:
                from figure_discovery import stability_figures
                try:
                    found=stability_figures(pdf)
                    for fig in found:
                        diagnostics.append(dict(file=pdf.name,record_type='POTENTIAL_STABILITY_FIGURE',
                                                reason='Figure candidate; not yet numerically digitized',**fig))
                    if found:
                        print(f'  STABILITY_FIGURE_CANDIDATES: {len(found)}',flush=True)
                except Exception as exc:
                    diagnostics.append(dict(file=pdf.name,record_type='DISCOVERY_ERROR',reason=str(exc)))
            # v18.1: recover published PV tables/prose for STABILITY-TESTED categories only.
            # This does not alter curve digitization, fitting or selected T80 models.
            if paper_rows:
                from pv_table_recovery import recover as recover_category_pv
                recover_category_pv(paper_rows,pdf,[d.path for d,ev,sc in matched])
                # v18.6.3: category-aware J-V recovery from arbitrary publication
                # prose/table content. Additive, never edits stability or fits.
                from general_pv_measurement_recovery import recover as recover_general_pv
                n_pv_fields=recover_general_pv(paper_rows,pdf,[d.path for d,ev,sc in matched])
                if n_pv_fields:
                    print(f'  GENERAL_PV_RECOVERED: {n_pv_fields} missing PCE/Voc/Jsc/FF fields',flush=True)
                # v18.6.6: supplementary-table PV recovery for new layouts,
                # including Greek eta efficiency symbols and same-row J-V
                # calculation. Runs after existing extraction and fills only
                # missing fields; requires unambiguous contact/scan identity.
                from eta_supplement_pv_recovery import recover as recover_eta_supplement
                n_eta_fields=recover_eta_supplement(paper_rows,pdf,[d.path for d,ev,sc in matched])
                if n_eta_fields:
                    print(f'  ETA_SUPPLEMENT_JV_RECOVERED: {n_eta_fields} source-matched PV fields',flush=True)
                # v18.2: content-driven EQE current extraction.  This searches the
                # entire publication (main article + every paired supplement) for
                # explicit integrated-current values and printed figure endpoints.
                # It has NO DOI, filename, category-value or page-number allowlist.
                # Existing calibrated/previously extracted data are never replaced.
                from eqe_publication_recovery import recover as recover_eqe_currents
                found_eqe=recover_eqe_currents(paper_rows,pdf,[d.path for d,ev,sc in matched])
                if found_eqe:
                    print(f'  GENERAL_EQE_RECOVERED: {found_eqe} previously missing category currents',flush=True)
                # v18.3: FINISH this matched publication BEFORE advancing to next.
                # Promote existing JSON stability/fitting fields to the unified
                # Excel columns; recover absorber and optical gap from article/SI.
                from publication_completion import finalize_publication
                finalize_publication(paper_rows,pdf,[d.path for d,ev,sc in matched],
                                     {'doi':canonical_doi})
                # v18.4: numerical 50%-edge bandgap recovery from this
                # publication's calibrated EQE spectra, including shared edges.
                # Does not alter stability, fitted decay, or PV measurements.
                from eqe_bandgap_completion import recover_eqe_bandgaps
                completed_gaps=recover_eqe_bandgaps(
                    paper_rows,pdf,[d.path for d,ev,sc in matched],spectral_root=target)
                if completed_gaps:
                    print(f'  EQE_50PCT_BANDGAPS_RECOVERED: {completed_gaps}',flush=True)
                # v18.4 content-driven fallback: find any suitable EQE/IPCE
                # spectrum across THIS matched publication's PDFs, calibrate
                # printed wavelength/IPCE axes and detect overlapping red edges.
                from eqe_auto_cutoff import recover_shared_cutoff
                automated_gaps,auto_eqe_records=recover_shared_cutoff(
                    paper_rows,pdf,[d.path for d,ev,sc in matched])
                if automated_gaps:
                    print(f'  AUTO_SHARED_EQE_BANDGAPS: {automated_gaps}',flush=True)
                # v18.6.4: match written aging-initial PCE independently
                # from JV champion PCE; distinguish paired EQE-derived gaps.
                # Applied after generic PV, with source-type precedence.
                from contextual_pv_bandgap_recovery import recover as recover_source_context
                source_changes=recover_source_context(paper_rows,pdf,[d.path for d,ev,sc in matched])
                if source_changes:
                    print(f'  SOURCE_CONTEXTUAL_PV_BANDGAP_RECOVERED: {source_changes}',flush=True)
                # v18.6.4: optional raster-only J-V fallback. Runs AFTER the
                # source-matched J-V PCE becomes available. Accepts a trace
                # only if its calibrated axes, named colored legend, and
                # reconstructed maximum power agree with that J-V PCE.
                # Does not overwrite any existing PV metric or modify fits.
                from raster_jv_measurement_recovery import recover as recover_raster_jv
                raster_jv_fields=recover_raster_jv(paper_rows,pdf,[d.path for d,ev,sc in matched])
                if raster_jv_fields:
                    print(f'  RASTER_JV_RECOVERED: {raster_jv_fields} source-checked plotted Voc/Jsc/derived FF fields',flush=True)
                # v18.6.5: conservative *image-only* EQE/IPCE fallback.
                # Calibrate the printed wavelength and integrated-current axes,
                # resolve one explicitly labeled device, and recover its
                # actual cumulative-current plateau plus 50%-red-edge gap.
                # Never substitute J-V Jsc or fill unplotted formulations.
                from raster_eqe_measurement_recovery import recover as recover_raster_eqe
                n_raster_eqe=recover_raster_eqe(paper_rows,pdf,[d.path for d,ev,sc in matched])
                if n_raster_eqe:
                    print(f'  RASTER_EQE_RECOVERED: {n_raster_eqe} source-calibrated EQE current/band-gap fields',flush=True)
                # Additional category-specific published endpoint QC. Does not
                # modify or refit the curve or force fitted T80 to agree.
                from contextual_endpoint_check import verify as verify_named_endpoints
                verified=verify_named_endpoints(paper_rows,pdf,[d.path for d,ev,sc in matched])
                if verified:
                    print(f'  SOURCE_NAMED_ENDPOINT_CHECKS: {verified}',flush=True)

            # v18.6.7-FIX: scan EVERY source page, including end-appendix images,
            # and never gate discovery on an existing digitized stability row.
            try:
                from publication_evidence import scan as scan_all_stability
                source_figs,source_endpoints=scan_all_stability(
                    pdf,[d.path for d,ev,sc in matched])
                all_figures.extend(source_figs)
                all_reported_endpoints.extend(source_endpoints)
                if source_figs:
                    print(f'  FULL_DOCUMENT_FIGURE_SCAN: {len(source_figs)} source figure pages; '
                          f'{len(source_endpoints)} textual endpoint candidates',flush=True)
            except Exception as ex:
                diagnostics.append(dict(file=pdf.name,record_type='FULL_DOCUMENT_SCAN_ERROR',reason=str(ex)))
            # Marker-only extraction is independent of the continuous trace reader.
            # Keep its raw values out of fitted training data until identity and
            # error-bar/legend/artifact checks pass.
            try:
                from marker_stability_recovery import detect_pdf as read_markers
                for path in [pdf,*[d.path for d,ev,sc in matched]]:
                    if Path(path).suffix.lower()=='.pdf':
                        all_marker_candidates.extend(read_markers(path,source_name=Path(path).name))
            except Exception as ex:
                diagnostics.append(dict(file=pdf.name,record_type='MARKER_READER_ERROR',reason=str(ex)))
            # Newly calibrated continuous MPP measurements (not line-connected
            # marker/error-bar plots): recover with source-printed numeric axes.
            try:
                from vector_stability_promotion import obtain as obtain_continuous, to_rows as promote_continuous
                cont,verified=obtain_continuous(pdf,[d.path for d,ev,sc in matched])
                all_continuous_candidates.extend(cont)
                if verified and not paper_rows:
                    recovered=promote_continuous(verified,canonical_doi or publication.doi or '',
                                                publication.title,pdf.name)
                    if recovered:
                        paper_rows=recovered
                        print(f'  CONTINUOUS_MPP_RECOVERED: {len(recovered)} source-caption-identified categories, '
                              f'{sum(r["stability_experiment_count"] for r in recovered)} traces',flush=True)
                elif cont:
                    print(f'  CONTINUOUS_MPP_CANDIDATES: {len(cont)} source-calibrated panels',flush=True)
            except Exception as ex:
                diagnostics.append(dict(file=pdf.name,record_type='CONTINUOUS_MPP_ERROR',reason=str(ex)))
            # Marker glyphs with error bars: recover only the plotted symbol
            # centers. Never digitize connector lines or error-bar endpoints.
            # The preliminary rows explicitly disable automatic ML training.
            if not paper_rows:
                try:
                    from marker_stability_promotion import obtain as obtain_markers, to_rows as promote_markers
                    own_markers=[r for r in all_marker_candidates if r.get('source_file')==pdf.name]
                    if own_markers:
                        marker_audit,marker_verified=obtain_markers(pdf,raw_records=own_markers)
                        if marker_verified:
                            marker_rows=promote_markers(marker_verified,
                                canonical_doi or publication.doi or '',publication.title,pdf.name)
                            if marker_rows:
                                paper_rows=marker_rows
                                print(f'  MARKER_ONLY_PRELIMINARY: {len(marker_rows)} devices, '
                                      f'{sum(r["stability_experiment_count"] for r in marker_rows)} traces '
                                      '(NOT YET APPROVED FOR ML TRAINING)',flush=True)
                except Exception as ex:
                    diagnostics.append(dict(file=pdf.name,record_type='MARKER_PROMOTION_ERROR',reason=str(ex)))
            # Important: newly recovered fallback stability categories also
            # need their own PCE/Voc/Jsc/FF extraction pass. Previously the
            # primary PV reader ran only BEFORE these categories existed.
            if paper_rows and all(r.get('extraction_method') in (
                    'NATIVE_AXIS_CONTINUOUS_MPP_RASTER',
                    'RASTER_MARKER_CENTERS_PRELIMINARY') for r in paper_rows):
                supplemental_paths=[d.path for d,ev,sc in matched]
                try:
                    from pv_table_recovery import recover as source_pv
                    source_pv(paper_rows,pdf,supplemental_paths)
                    from general_pv_measurement_recovery import recover as source_general_pv
                    n_post=source_general_pv(paper_rows,pdf,supplemental_paths)
                    from eta_supplement_pv_recovery import recover as source_eta
                    n_post+=source_eta(paper_rows,pdf,supplemental_paths)
                    if n_post:
                        print(f'  FALLBACK_DEVICE_PV_RECOVERED: {n_post} matched numerical fields',flush=True)
                except Exception as ex:
                    diagnostics.append(dict(file=pdf.name,record_type='FALLBACK_DEVICE_PV_ERROR',reason=str(ex)))
            # v18.6.7: keep auditable evidence when text/raster compound
            # figures defeat the calibrated reader. This is NOT fake ML data.
            if not paper_rows:
                try:
                    from figure_rescue_audit import scan_paper
                    found=scan_paper(pdf,[d.path for d,ev,sc in matched],
                                     preview_dir=dest.parent/'stability_figure_review')
                    rescue_records.extend(found)
                    if found:
                        print(f'  STABILITY_EVIDENCE_RESCUED_FOR_REVIEW: {len(found)} source pages; '+
                              str(sum(z.get('calibrated_panels',0) for z in found))+
                              ' panel axes calibrated (NOT model-ready curves)',flush=True)
                        for record in found:
                            diagnostics.append(dict(file=record['source_file'],record_type='UNRESOLVED_STABILITY_GRAPH',
                                page=record['source_page'],reason=record['status'],
                                figure_labels=record['caption_evidence'][:200]))
                except Exception as ex:
                    diagnostics.append(dict(file=pdf.name,record_type='FIGURE_RESCUE_ERROR',reason=str(ex)))
                skipped.append((pdf.name,locals().get('calibrated_error','No automatically digitizable stability curve found')))
            for r in paper_rows:
                r['pdf_file']=str(pdf.name)
                if canonical_doi:r['doi']=canonical_doi
                # v15 and generic reader originally used different names for
                # the same signed/absolute Jsc discrepancy ML features.
                for legacy,canonical in (
                    ('jsc_discrepancy_signed_pct','jsc_difference_signed_pct'),
                    ('jsc_discrepancy_absolute_pct','jsc_difference_absolute_pct')):
                    if r.get(canonical) is None and r.get(legacy) is not None:
                        r[canonical]=r[legacy]
                    if r.get(legacy) is None and r.get(canonical) is not None:
                        r[legacy]=r[canonical]
                # v15 uses 'retention_last_pct' / 'stability_last_time_h';
                # v17 uses 'final_retention_pct' / 'final_aging_time_h'.
                # Normalize the export schema without changing either digitizer.
                if r.get('final_retention_pct') is None:
                    r['final_retention_pct']=r.get('retention_last_pct')
                if r.get('final_aging_time_h') is None:
                    r['final_aging_time_h']=r.get('stability_last_time_h')
                if r.get('stability_endpoint') is None and r.get('final_retention_pct') is not None and r.get('final_aging_time_h') is not None:
                    r['stability_endpoint']=str((r['final_retention_pct'],r['final_aging_time_h']))
                if not r.get('retention_reference') and r.get('primary_stability_metric'):
                    metric=str(r['primary_stability_metric']).lower()
                    if 'normalized' in metric:r['retention_reference']='PUBLISHED_NORMALIZED_RETENTION'
            
            # v18.6.7: J-V tables embedded within figure panels are not
            # necessarily labeled 'Table'. Read only physically aligned rows.
            # Also audit their candidates when stability digitization failed.
            try:
                from figure_pv_table_recovery import recover as recover_figure_pv
                n_fig_pv,candidates=recover_figure_pv(paper_rows,pdf,[d.path for d,ev,sc in matched])
                figure_pv_candidates.extend(candidates)
                if n_fig_pv:
                    print(f'  FIGURE_TABLE_PV_RECOVERED: {n_fig_pv} source-checked fields',flush=True)
                if not paper_rows and candidates:
                    print(f'  FIGURE_TABLE_PV_FOUND_AWAITING_STABILITY: {len(candidates)} rows',flush=True)
            except Exception as ex:
                diagnostics.append(dict(file=pdf.name,record_type='FIGURE_PV_ERROR',reason=str(ex)))

            # v18.6.2: apply the identical source-backed T80 comparison to
            # the calibrated reader (generic rows already carry this evidence).
            # Additive metadata only; never re-fit or edit measured curves.
            if paper_rows:
                from graph_endpoint_verification import annotate_calibrated_rows
                checked=annotate_calibrated_rows(paper_rows,pdf,[d.path for d,ev,sc in matched])
                if checked:print(f'  SOURCE_T80_COMPARISONS: {checked}',flush=True)

            # Hysteresis ML features: only source-matched F/R J-V PCE pairs.
            # Never confuse MPPT eta or a different active area with scan PCE.
            try:
                from hysteresis_recovery import recover as recover_hi
                own_figure_rows=[c for c in figure_pv_candidates if c.get('file')==pdf.name]
                hi_sources=recover_hi(paper_rows,pdf,[d.path for d,ev,sc in matched],
                                      extra_candidates=own_figure_rows)
                if hi_sources:
                    print(f'  HI_SCAN_PCE_CANDIDATES: {len(hi_sources)} source-aligned scan rows; '
                          f'{sum(r.get("hi_status")=="CALCULATED_FROM_SAME_SOURCE_JV_PAIR" for r in paper_rows)} paired device categories',flush=True)
                all_hi_candidates.extend(hi_sources)
            except Exception as ex:
                diagnostics.append(dict(file=pdf.name,record_type='HI_RECOVERY_ERROR',reason=str(ex)))
            # Shared canonical device name + scan-aware J-V selection, AFTER all
            # recovery routes, so forward values cannot overwrite reverse data.
            try:
                from device_identity_pv import apply as apply_scan_pv
                own_figure_rows=[c for c in figure_pv_candidates if c.get('file')==pdf.name]
                apply_scan_pv(paper_rows,pdf,[d.path for d,ev,sc in matched],own_figure_rows)
            except Exception as ex:
                diagnostics.append(dict(file=pdf.name,record_type='CANONICAL_JV_RECOVERY_ERROR',reason=str(ex)))
            # General complete J-V prose blocks can identify their device AFTER the
            # reported PCE/Voc/Jsc/FF measurements, even without a separate table.
            # Apply after reverse-first selection; never override a scan-qualified set.
            try:
                from narrative_jv_recovery import recover as recover_narrative_jv
                prose_rows=recover_narrative_jv(paper_rows,pdf,[d.path for d,ev,sc in matched])
                if prose_rows:
                    print(f'  COMPLETE_NARRATIVE_JV_CANDIDATES: {len(prose_rows)}; '+
                          str(sum(bool(r.get('narrative_jv_status','').startswith('SOURCE_LINKED')) for r in paper_rows))+
                          ' source-linked device categories',flush=True)
            except Exception as ex:
                diagnostics.append(dict(file=pdf.name,record_type='NARRATIVE_JV_RECOVERY_ERROR',reason=str(ex)))
            paper_rows,validation_stats=validate_publication(paper_rows,publication)
            if any(validation_stats.values()):
                print(f'  VALIDATION: {validation_stats}',flush=True)
            all_rows.extend(paper_rows)
            print(f'  COMPLETE: {len(paper_rows)} categories',flush=True)
    # Always export evidence and missing-graph diagnostics, including when
    # ALL source figures resisted unambiguous digitization.
    dest.parent.mkdir(parents=True,exist_ok=True)
    # Remove stale source-dependent audit files before writing this specific run.
    # Otherwise a no-hit run can display candidates from a previous publication.
    for audit_filename in ('figure_jv_table_candidates.csv','marker_only_stability_candidates.csv',
       'continuous_mpp_candidates.csv','stability_all_figure_pages.csv',
       'stability_reported_endpoints.csv','stability_figures_needing_review.csv',
       'hysteresis_scan_candidates.csv'):
        old=dest.parent/audit_filename
        if old.exists():old.unlink()
    from figure_pv_table_recovery import save_candidates
    save_candidates(figure_pv_candidates,dest.parent/'figure_jv_table_candidates.csv')
    from publication_evidence import save as save_evidence
    from marker_stability_recovery import save as save_markers
    save_markers(all_marker_candidates,dest.parent/'marker_only_stability_candidates.csv')
    save_evidence([{**{k:v for k,v in r.items() if k!='points'},
                    'points_json':json.dumps(r.get('points',[]))} for r in all_continuous_candidates],
        dest.parent/'continuous_mpp_candidates.csv',
        ['source_file','source_page','axis','color','number_of_points','t_end_h','retention_final_pct',
         'status','device_label','local_panel_text','plot_bounds','points_json'])
    save_evidence(all_figures,dest.parent/'stability_all_figure_pages.csv',
                  ['source_file','page','figure_captions','embedded_images','vector_drawings','evidence_state'])
    save_evidence(all_reported_endpoints,dest.parent/'stability_reported_endpoints.csv',
                  ['source_file','source_page','record_kind','measurement','reported_percent','reported_hours','source_context','status'])
    from figure_rescue_audit import save as save_figure_rescue
    save_figure_rescue(rescue_records,dest.parent/'stability_figures_needing_review.csv')
    import csv as _hi_csv
    with (dest.parent/'hysteresis_scan_candidates.csv').open('w',encoding='utf-8-sig',newline='') as _f:
        _fields=['device_label','scan_direction','pce_pct','source_file','source_page','methodology','evidence']
        _w=_hi_csv.DictWriter(_f,fieldnames=_fields,extrasaction='ignore');_w.writeheader()
        _w.writerows(all_hi_candidates)
    from extraction_detail_export import produce as write_details
    coverage=write_details(dest,all_rows,all_figures,all_reported_endpoints,
                           all_marker_candidates,all_continuous_candidates,
                           [p.path.name for p in pdfs])
    for c in coverage:
        print('  COVERAGE:',c['source_pdf'],c['coverage_status'],
              f"{c['exported_stability_traces']} traces; {c['exported_points']} points; "+
              c['unresolved_reason'],flush=True)
    review_path=save_extraction_audit(dest,diagnostics,skipped)
    if not all_rows:
        print('NO_VERIFIED_ML_CURVES: Source figures/PV evidence recorded in audit CSVs; '+
              'the pipeline has not invented degradation points or fitted T80.',flush=True)
        from openpyxl import Workbook
        wb=Workbook();wb.active.title='ML_dataset'
        wb.active.append(['status','note']);wb.active.append(['NO_VERIFIED_CURVE','See evidence sheets; no invented ML training data'])
        for tab,fname in [('PV_table_candidates','figure_jv_table_candidates.csv'),
            ('Stability_marker_review','marker_only_stability_candidates.csv'),
        ('Continuous_MPP','continuous_mpp_candidates.csv'),
            ('Stability_figures','stability_all_figure_pages.csv'),
            ('Reported_endpoints','stability_reported_endpoints.csv'),
            ('Stability_figure_review','stability_figures_needing_review.csv'),
            ('Stability_points','extracted_stability_points.csv'),
            ('PV_measured_values','extracted_pv_metrics.csv'),
            ('Extraction_coverage','extraction_coverage.csv'),
            ('HI_scan_sources','hysteresis_scan_candidates.csv')]:
            _append_audit_sheet(wb,tab,dest.parent/fname)
        wb.save(dest)
        print('EXTRACTION_REVIEW_QUEUE:',review_path,flush=True)
        print('SAVED EVIDENCE-ONLY WORKBOOK:',dest,flush=True)
        return
    n=write_excel(all_rows,dest,diagnostics=diagnostics)
    print(f'EXTRACTION_REVIEW_QUEUE: {review_path}',flush=True)
    if skipped:print('UNEXTRACTED MAIN PDFs:',skipped)
    if unmatched_supplements:print('UNPAIRED/AMBIGUOUS SUPPLEMENTS:',unmatched_supplements)
    print(f'SAVED ONE OUTPUT FILE: {dest}\nSTABILITY-TESTED CATEGORIES ONLY: {n} | STABILITY TRACES: {sum(r["stability_experiment_count"] for r in all_rows)}')



def save_extraction_audit(destination, diagnostics, skipped, documents=None):
    """Write machine-readable QC evidence alongside the single-sheet ML file.

    Never modifies extracted numeric data, reported PV values or fitted models.
    """
    import csv
    destination=Path(destination)
    review=destination.parent/'review_queue.csv'
    fields=('file','record_type','page','figure_labels','reason','matching_main','match_evidence')
    with review.open('w',newline='',encoding='utf-8-sig') as f:
        writer=csv.DictWriter(f,fieldnames=fields,extrasaction='ignore')
        writer.writeheader()
        for item in diagnostics:
            # Include graph calibration failures and rejected duplicate curves.
            writer.writerow({k:item.get(k,'') for k in fields})
        for name,reason in skipped:
            writer.writerow({'file':name,'record_type':'MAIN_ARTICLE_NOT_EXTRACTED','reason':reason})
    return review

if __name__=='__main__':main()
