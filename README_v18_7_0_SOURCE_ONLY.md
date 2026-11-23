# Perovskite Stability ML — v18.7.0 Source-Only

## What is different

This is a **publication-agnostic rebuild** of the user's v18.6.9.1 package. The former curated/legacy backend was removed from the production route. It contained embedded DOI-to-measurement mappings, predefined plot bounding boxes, manually entered stability points, fixed device composition values, and paper-specific rules. The removed code is **not shipped** in this ZIP.

The entry point `perovskite_onefile_final.py` processes the supplied file folder afresh. It no longer offers `--existing-results`, which could import previously curated or manually amended measurements.

### Evidence flow

1. Discover local PDFs/DOCX and extract publication metadata (article title, authors, DOI, year, journal). DOI is an **identifier only**, not a numeric lookup key.
2. Pair each article with its supplementary documents using printed DOI, title similarity, author overlap, publication year and other bibliographic evidence. Ambiguous matches remain unpaired.
3. For each publication bundle, scan plots, time axes, stability text and numerical endpoints throughout its documents. The generic plot reader may recover calibrated trajectories; the marker-only reader records preliminary symbol-center evidence; the continuous MPP reader processes continuous signals.
4. Recover J–V source rows from printed text, figure tables and matching supplementary files; use a generic OCR fallback on graphics-only PV tables. Canonicalize *standalone* F/Fw/FWD and R/Re/Rw/B/Bw/Backward scan labels without modifying chemical names.
5. Put source-aligned **reverse/backward** measurements into `pce0/voc0/jsc0/ff0` when available. Forward rows are retained for signed and absolute hysteresis, not used to override reverse values. Do not mix different active areas or unrelated champion and aging cells.
6. Derive normalized aging retention from the actual source axis and document the baseline. Keep unverified marker points in review status, outside ML-ready fits.
7. Write the main Excel workbook and detailed audit sheets with source-file/page evidence, ambiguities and missing fields.

### Run

```powershell
python -m pip install -r requirements.txt
python perovskite_onefile_final.py --pdf-dir "C:\path\to\papers" --output "C:\path\to\Perovskite_ML_dataset.xlsx"
```

Without `--pdf-dir`, the original folder picker appears. Install the Tesseract OCR **executable** and make it available on PATH for image-backed supplementary tables. DOCX graph digitization requires a usable Word-to-PDF converter: on Windows, Microsoft Word + `pywin32`; where LibreOffice is installed, its headless converter is supported. If conversion is unavailable the program reports this explicitly and cannot claim to digitize Word figures.

### Source-only rule

Never reintroduce `if doi == ...` to choose data, `IDXX` category catalogs, manually transcribed publication values, constant graph boxes tied to a paper, or historical workbook results as a recovery shortcut. Extend generic document parsers instead, using *new papers only as regression tests*.

### Honest limitations

This is **not a claim of complete numerical extraction for every visual figure**. The broad generic graph readers remain imperfect. Unknown legends, shared axes, overlapping error bars, truncated/rotated labels or raster quality can leave points unassigned or provisional. Missing values stay missing unless source evidence supports them; do not interpret a failed extractor as an experiment with zero measurements.

The previous version's higher category count may have relied on curated paper mappings; numerical output may decrease after removing them. That is intentional until a trustworthy general extraction algorithm recovers the source data independently.

## Local checks performed

- Core imports, Python syntax, CLI and source-only invariant tests passed.
- The Science Bulletin manuscript processed end-to-end using the new **generic-only** entry point: 2 categories, 4 **preliminary** marker traces (74 review points), 1 complete source-linked prose J–V set. The extracted modified device J–V columns were 22.0%, 1.16 V, 23.31 mA/cm², 81.2% FF; control remains blank rather than borrowing another device's parameters.
- Pairing test renamed both the Nature Materials main paper and its supplementary file to generic names: correct document-title/DOI pairing, **no** filename whitelist.
- The generic supplementary J–V reader recovered 16 complete source-derived measurement rows from the renamed paired documents, including separate scan directions and active areas.
- Synthetic regression verified reverse-scan priority, common abbreviations, preserving chemical names, and HI numerical formula.

**Not tested:** The complete 39-paper batch; full Nature Materials stability workbook end-to-end in this release; all PDF/Word graph layouts. Do not claim those are verified.
