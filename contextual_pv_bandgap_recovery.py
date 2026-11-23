"""v18.6.4 source-grounded late-stage additions for publication rows.

Do not guess absent J-V values; do not borrow champion J-V measurements as
initial aged-device PCE. No DOI names, article-specific page maps, or manual
value lookup. Source precedence is enforced at measurement-type level.
"""
from __future__ import annotations
from pathlib import Path
import re
import fitz


def _clean(text):
    text=str(text).replace('−','-').replace('–','-').replace('‐','-')
    text=re.sub(r'(?<=[A-Za-z])-\s+(?=[a-z])','-',text)
    return ' '.join(text.split())


def _text_sources(main_pdf,supplements):
    for src in [main_pdf,*supplements]:
        if Path(src).suffix.lower()!='.pdf':continue
        try:
            with fitz.open(str(src)) as d:
                for i,p in enumerate(d):
                    yield (Path(src).name,i+1,_clean(p.get_text()))
        except Exception:continue


def _category_match(written,category):
    norm=lambda x:re.sub(r'[^a-z0-9]+','',x.lower())
    return norm(written)==norm(category)


def recover(rows,main_pdf,supplements=()):
    labels=[str(r.get('device_category') or '').strip() for r in rows]
    labels=[x for x in labels if len(x)>2]
    if not labels:return 0
    label_rx='|'.join(re.escape(x) for x in sorted(labels,key=len,reverse=True))
    counts=0
    for file,page,text in _text_sources(main_pdf,supplements):
        source=f'{file} p.{page}'
        # Initial performance of the actual MPP tracking population, not the
        # champion JV device: 'initial PCE values ... were 24.11% for A device'.
        for initial in re.finditer(r'\binitial\s+PCE\s+values\b',text,re.I):
            section=text[initial.end():initial.end()+450]
            paired=list(re.finditer(r'(\d{1,2}(?:\.\d+)?)\s*%\s+for\s+(?:the\s+)?('+label_rx+r')\s+(?:devices?|PSCs?|cells?)',section,re.I))
            if len(paired)<2:continue
            for match in paired:
                value=float(match.group(1));named=match.group(2)
                for row in rows:
                    if not _category_match(named,str(row.get('device_category') or '')):continue
                    if 1<value<40 and row.get('pce0') in (None,''):
                        row['pce0']=value
                        row['pce0_source']=source+' explicitly stated initial aging-device PCE'
                        counts+=1
        # Literature may state device-specific champion JV PCEs with the
        # name before or after the number. Keep distinct from aging pce0.
        for named in labels:
            esc=re.escape(named)
            matches=[]
            # X device delivers/achieves/reaches a (champion) PCE of NN%
            before=r'(?<![A-Za-z0-9])'+esc+r'\s+(?:devices?|PSCs?|cells?)\s+(?:achieves?|delivers?|reaches?|exhibits?|shows?)\s+(?:an?\s+)?(?:champion\s+)?PCE\s+of\s+(\d{1,2}(?:\.\d+)?)\s*%'
            matches.extend(re.finditer(before,text,re.I))
            # champion PCE of NN% for X device (one category).
            after=r'\b(?:champion\s+)?PCE\s+of\s+(\d{1,2}(?:\.\d+)?)\s*%\s+for\s+(?:the\s+)?'+esc+r'\s+(?:devices?|PSCs?|cells?)'
            matches.extend(re.finditer(after,text,re.I))
            # Wording like 'PCE increased to 26.02% for FACs@PIC device'.
            # Require the local same-sentence PCE/J-V context; the number
            # alone near a category must never be assumed to be its PCE.
            backwards=r'(\d{1,2}(?:\.\d+)?)\s*%\s+for\s+(?:the\s+)?'+esc+r'\s+(?:devices?|PSCs?|cells?)'
            for m in re.finditer(backwards,text,re.I):
                lead=text[max(0,m.start()-210):m.start()]
                sentence=re.split(r'[.!?]\s+(?=[A-Z])',lead)[-1]
                if re.search(r'\b(?:PCE|power conversion efficiency|efficiency)\b',sentence,re.I):
                    matches.append(m)
            for m in matches:
                value=float(m.group(1))
                if not 1<value<40:continue
                for row in rows:
                    if not _category_match(named,str(row.get('device_category') or '')):continue
                    if row.get('pce_jv_pct') in (None,''):
                        row['pce_jv_pct']=value
                        row['pce_jv_source']=source+' named champion/representative JV statement'
                        counts+=1
                    # Never overwrite experimental starting PCE with champion.
        # One figure's EQE-derived band gaps often show two DEVICE TYPES:
        # 'EQE-derived PV bandgaps ... for evaporated and solution-processed
        # PSCs are 1.67 and 1.685 eV, respectively'. A range such as
        # '1.65--1.72 eV' is NOT a device-specific measured gap.
        paired_gap=re.search(r'(?:EQE\)?[- ]derived|external quantum efficiency.{0,25}derived).{0,70}?band\s*gaps?.{0,105}?for\s+(?:the\s+)?([A-Za-z-]+)\s+and\s+([A-Za-z-]+)\s+(?:PSCs?|cells?|devices?)\s+(?:are|were)\s+(\d\.\d{1,3})\s+and\s+(\d\.\d{1,3})\s*eV\s*,?\s*respectively',text,re.I)
        if paired_gap:
            pairs=[(paired_gap.group(1),float(paired_gap.group(3))),
                   (paired_gap.group(2),float(paired_gap.group(4)))]
            for row in rows:
                cat=str(row.get('device_category') or '').lower()
                for device,gap in pairs:
                    root=device.lower().split('-')[0]
                    if not (root in ('evaporated','evaporation','solution') or len(root)>3):continue
                    if (root=='evaporated' and 'evaporation' not in cat and 'evaporat' not in cat):continue
                    if root!='evaporated' and root not in cat:continue
                    if not .7<=gap<=3:continue
                    current=row.get('bandgap_ev')
                    # Reported EQE band gap outranks the older first-number
                    # heuristic only if source explicitly distinguishes the
                    # plotted architectures. Retain previous in QC metadata.
                    if current not in (None,'') and abs(float(current)-gap)>0.002:
                        row['previous_bandgap_before_eqe_source_check_ev']=float(current)
                    if current in (None,'') or abs(float(current)-gap)>0.002:
                        row['bandgap_ev']=gap
                        row['bandgap_source']=source+' explicit paired EQE-derived band gaps, respectively'
                        counts+=1
    return counts
