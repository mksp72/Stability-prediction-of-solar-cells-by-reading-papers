"""Publication-local canonical device matching and reverse-scan-first J-V selection.

No DOI allowlists. No mixing independent scan measurements or JV and MPPT values.
"""
from __future__ import annotations
from pathlib import Path
from collections import defaultdict
import math
from hysteresis_recovery import canonical_scan_label, canonical_key

FIELDS=('pce0','voc0','jsc0','ff0')

def canonical_device(label):
    return canonical_scan_label(label)[0].strip()

def _number(value):
    try:
        value=float(value)
        return value if math.isfinite(value) else None
    except (TypeError,ValueError):return None

def _valid(vals):
    result={k:_number(vals.get(k)) for k in FIELDS}
    if result['ff0'] is not None and 0 < result['ff0'] <= 1.0:
        result['ff0']*=100
    bounds={'pce0':(0.5,45),'voc0':(0.2,3.2),'jsc0':(0.5,50),'ff0':(5,100)}
    for k,v in result.items():
        if v is not None and not(bounds[k][0]<v<=bounds[k][1]):return None
    if all(v is not None for v in result.values()):
        calculated=result['voc0']*result['jsc0']*result['ff0']/100
        if abs(calculated-result['pce0'])>max(1.2,0.09*result['pce0']):return None
    return result

def extract_candidates(main_pdf,supplements=(),figure_candidates=()):
    from eta_supplement_pv_recovery import _from_pymupdf_tables,_from_spatial,_from_lines
    import fitz
    out=[]
    for path in (Path(p) for p in [main_pdf,*supplements]):
        if not path.is_file() or path.suffix.lower()!='.pdf':continue
        try:
            with fitz.open(path) as doc:
                for pi in range(len(doc)):
                    page=doc[pi]
                    if not page.get_text('words'):continue
                    for reader in (_from_pymupdf_tables,_from_spatial,_from_lines):
                        try: entries=reader(doc,pi)
                        except Exception:entries=[]
                        for label,vals,number,evidence in entries:
                            base,direction=canonical_scan_label(label)
                            clean=_valid(vals)
                            if not canonical_key(base) or clean is None:continue
                            # Need at least two aligned metrics; isolated prose never used for replacement.
                            if sum(v is not None for v in clean.values())<2:continue
                            out.append(dict(device_label=label,device_canonical=base,
                                scan_direction=direction or 'unspecified',source_file=path.name,
                                source_page=number,method=reader.__name__,values=clean,
                                evidence=str(evidence)[:200]))
        except Exception:pass
    # Vector-outlined or rasterized SI tables: the caption is selectable but
    # numerical cells may have ZERO native PDF text. OCR is constrained to
    # source table pages and validated by same-row PCE closure.
    try:
        from raster_supplement_jv import extract as extract_raster_supplement
        for c in extract_raster_supplement(main_pdf,supplements):
            vals=_valid(c.get('values',{}))
            if vals is None or not all(vals[k] is not None for k in FIELDS):continue
            out.append(dict(device_label=c['device_label'],device_canonical=canonical_device(c['device_label']),
                scan_direction=c['scan_direction'],source_file=c['source_file'],
                source_page=c['source_page'],method=c['method'],values=vals,
                area_cm2=c.get('area_cm2'),evidence=c['evidence']))
    except Exception as err:
        # The ordinary text reader must still work on systems without OCR.
        print('  RASTER_SUPPLEMENT_JV_UNAVAILABLE:',str(err)[:130],flush=True)
    for c in figure_candidates:
        base,direction=canonical_scan_label(c.get('device_label',''))
        vals=_valid(c)
        if not canonical_key(base) or vals is None or sum(v is not None for v in vals.values())<2:continue
        out.append(dict(device_label=c.get('device_label',''),device_canonical=base,
            scan_direction=direction or 'unspecified',source_file=c.get('file',''),
            source_page=c.get('page',''),method='figure_table',values=vals,evidence=c.get('status','')))
    dedup={}
    for c in out:
        ident=(c['source_file'],str(c['source_page']),canonical_key(c['device_canonical']),
               c['scan_direction'],c.get('area_cm2'),tuple(round(c['values'][k],4) if c['values'][k] is not None else None for k in FIELDS))
        dedup.setdefault(ident,c)
    return list(dedup.values())

def apply(rows,main_pdf,supplements=(),figure_candidates=(),candidates=None):
    """Apply a complete source-aligned J-V set, reverse if both scans exist.

    Leaves unverified/no-pair existing metrics untouched with explicit status.
    Qualifier-only different labels are canonically named in every resulting row.
    """
    candidates=extract_candidates(main_pdf,supplements,figure_candidates) if candidates is None else candidates
    for row in rows:
        label=str(row.get('device_category',''))
        base,old_direction=canonical_scan_label(label)
        row['device_category_original']=label
        row['device_category']=base
        row['device_category_canonical']=base
        row['jv_selected_scan_direction']=None
        row['jv_selected_measurement_source']=None
        key=canonical_key(base)
        matched=[c for c in candidates if canonical_key(c['device_canonical'])==key]
        if not matched:
            row['jv_scan_selection_status']='NO_CANONICAL_DEVICE_JV_MATCH'
            continue
        # Never combine 0.25 and 1 cm2 scans. A matched explicit device
        # area is authoritative; otherwise use figure evidence where unique.
        def area_num(v):
            try:return float(str(v).replace('cm²','').replace('cm2','').strip())
            except (ValueError,TypeError):return None
        area=next((area_num(row.get(k)) for k in
                    ('device_area_cm2','active_area_cm2','area_cm2','cell_area_cm2')
                    if area_num(row.get(k)) is not None),None)
        if area is None:
            # Extracted stability categories from a specific figure may have
            # an area in the figure's textual device description.
            context=' '.join(str(row.get(k) or '') for k in
                ('primary_stability_conditions','primary_stability_figure','device_category_original',
                 'device_area','active_area','stability_condition'))
            import re
            area_hit=re.search(r'(\d+(?:\.\d+)?)\s*cm\s*(?:2|²)',context,re.I)
            if area_hit:area=float(area_hit.group(1))
        available=sorted({float(c['area_cm2']) for c in matched if c.get('area_cm2') is not None})
        if area is None and len(available)>1:
            # The SI explicitly says its representative J-V table is the one
            # shown in a numbered main-article figure panel. Resolve the area
            # from that particular panel caption, not from the smallest area
            # or an arbitrary first row. This is representative J-V area, not
            # a claim that every aging trace uses the same active area.
            try:
                import fitz
                import re
                refs=[]
                for sp in supplements:
                    with fitz.open(sp) as sd:
                        for page in list(sd)[:10]:
                            t=' '.join(page.get_text('text').split())
                            mm=re.search(r'(?i)Supplementary\s+Table\s+\d+.{0,210}?shown\s+in\s+Fig\.?\s*(\d+)\s*([a-z])',t)
                            if mm:refs.append((mm.group(1),mm.group(2).lower()))
                for figure_no,panel in refs:
                    with fitz.open(main_pdf) as md:
                        for page in md:
                            t=' '.join(page.get_text('text').split())
                            f=re.search(r'(?i)Fig\.?\s*'+re.escape(figure_no)+r'\s*[|:.]',t)
                            if not f:continue
                            cap=t[f.end():f.end()+1500]
                            mp=re.search(r'\b'+re.escape(panel)+r'\s*[,.)]\s*(.*?)(?=\s+[a-z]\s*[,.)]\s|$)',cap,re.I)
                            if mp:
                                am=re.search(r'(\d+(?:\.\d+)?)\s*[- ]?cm\s*(?:2|²)',mp.group(1),re.I)
                                if am:
                                    candidate_area=float(am.group(1))
                                    if any(abs(candidate_area-x)<.03 for x in available):
                                        area=candidate_area
                                        row['jv_area_match_status']='SI_TABLE_CROSSREFERENCED_REPRESENTATIVE_FIGURE_AREA'
                                        break
                        if area is not None:break
                    if area is not None:break
            except Exception as ex:
                row['jv_area_match_warning']=str(ex)[:160]
        if area is None and len(available)>1:
            row['jv_scan_selection_status']='MULTIPLE_ACTIVE_AREAS_REQUIRE_STABILITY_AREA_MATCH'
            row['jv_available_areas_cm2']='; '.join(map(str,available))
            continue
        if area is not None:
            matched=[c for c in matched if c.get('area_cm2') is None or abs(float(c['area_cm2'])-area)<0.03]
            if not matched:
                row['jv_scan_selection_status']='NO_JV_SET_FOR_MATCHED_ACTIVE_AREA'
                continue
        groups=defaultdict(list)
        for c in matched:groups[c['scan_direction']].append(c)
        complete={d:[c for c in arr if all(c['values'][k] is not None for k in FIELDS)]
                  for d,arr in groups.items()}
        reverse=complete.get('reverse',[])
        forward=complete.get('forward',[])
        unspecified=complete.get('unspecified',[])
        # Reversed/backward data win whenever paired scans are available;
        # also select reverse-only data where unambiguous.
        if reverse:direction='reverse';pool=reverse
        elif groups.get('reverse') and forward:
            row['jv_scan_selection_status']='REVERSE_PRESENT_BUT_INCOMPLETE_DO_NOT_SUBSTITUTE_FORWARD'
            continue
        elif forward:direction='forward';pool=forward
        elif unspecified:direction='unspecified';pool=unspecified
        else:
            row['jv_scan_selection_status']='NO_COMPLETE_SAME_ROW_JV_SET'
            continue
        sets={tuple(round(c['values'][k],5) for k in FIELDS) for c in pool}
        if len(sets)!=1:
            row['jv_scan_selection_status']='AMBIGUOUS_MULTIPLE_'+direction.upper()+'_JV_SETS'
            continue
        c=pool[0]
        source=f"{c['source_file']} p.{c['source_page']}; {c['device_label']}; {c['method']}"
        for field in FIELDS:
            previous=_number(row.get(field))
            new=round(c['values'][field],5)
            if previous is not None and abs(previous-new)>1e-5:
                row[field+'_previous_value']=previous
                row[field+'_previous_source']=row.get(field+'_source')
            row[field]=new
            row[field+'_source']=f'J-V {direction} scan, complete same-row metrics; {source}'
        row['pce_jv_pct']=row['pce0']
        row['pce_jv_source']=row['pce0_source']
        row['jv_selected_scan_direction']=direction
        row['jv_selected_area_cm2']=c.get('area_cm2') if c.get('area_cm2') is not None else area
        row['jv_selected_measurement_source']=source
        # Compute hysteresis from same-formulation, same-area, same-source
        # scan pairs. This also covers OCR-only supplementary tables.
        hi_groups=defaultdict(lambda:defaultdict(set))
        for obs in matched:
            if obs.get('scan_direction') not in ('reverse','forward'):continue
            v=obs['values'].get('pce0')
            if v is None:continue
            hi_groups[(obs['source_file'],str(obs['source_page']),obs.get('area_cm2'))][obs['scan_direction']].add(round(v,5))
        hi_pairs=[]
        for (file,page,test_area),two in hi_groups.items():
            if len(two['reverse'])==1 and len(two['forward'])==1:
                hi_pairs.append((file,page,test_area,next(iter(two['reverse'])),next(iter(two['forward']))))
        if len(hi_pairs)==1:
            from hysteresis_recovery import calculate_hi
            file,page,test_area,rev,fw=hi_pairs[0]
            signed,absolute=calculate_hi(rev,fw)
            row.update(hi_reverse_pce_pct=rev,hi_forward_pce_pct=fw,
                       hi_signed_pct=signed,hi_absolute_pct=absolute,
                       hi_status='CALCULATED_FROM_SAME_SOURCE_JV_PAIR',
                       hi_source=f'{file} p.{page}; area {test_area} cm2; reverse/forward J-V pair')
        elif len(hi_pairs)>1:
            row['hi_status']='MULTIPLE_PAIRED_JV_SETS_REQUIRE_AREA_SELECTION'
        row['jv_scan_selection_status']=('REVERSE_SELECTED_BOTH_DIRECTIONS_PRESENT' if direction=='reverse' and forward
            else 'REVERSE_SELECTED' if direction=='reverse'
            else 'FORWARD_ONLY_NO_REVERSE_AVAILABLE' if direction=='forward'
            else 'UNSPECIFIED_ONLY_NO_SCAN_QUALIFIED_SET')
    return candidates
