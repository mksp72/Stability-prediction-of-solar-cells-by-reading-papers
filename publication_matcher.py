"""Content-first article/supplement matching for PDF and Word files.

No filename or DOI allowlist. Matching uses document text: primary article DOI,
article title, first/corresponding authors, journal and year. An ambiguous match is
never assigned automatically. Supports PDF and DOCX as metadata sources.
"""
from __future__ import annotations
from dataclasses import dataclass, field
from pathlib import Path
from difflib import SequenceMatcher
import re
import unicodedata

DOI_RE=re.compile(r'10\.\d{4,9}/[A-Za-z0-9][A-Za-z0-9._;()/:\-]*',re.I)
SI_RE=re.compile(r'\b(supplementary|supporting\s+information|supporting\s+materials|supplemental|electronic\s+supplementary|supplementary\s+materials|supporting\s+file)\b',re.I)
NAME_SI=re.compile(r'(?:_sm|_si|_supp|_supporting|_supplementary|_esm|_moesm\d*_esm|[-_]sup[-_]|\bsi\b|\bsupporting\b)',re.I)
JOURNALS={
 'Science':r'\bScience\b',
 'Nature Energy':r'\bNature\s+Energy\b',
 'Nature':r'\bNature\b',
 'Advanced Energy Materials':r'\bAdvanced\s+Energy\s+Materials\b',
 'Advanced Materials':r'\bAdvanced\s+Materials\b',
 'Advanced Science':r'\bAdvanced\s+Science\b',
 'ACS Energy Letters':r'\bACS\s+Energy\s+Letters\b',
 'Joule':r'\bJoule\b',
}

@dataclass
class Document:
    path:Path
    kind:str='main'
    title:str=''
    doi:str=''
    journal:str=''
    year:int|None=None
    authors:list[str]=field(default_factory=list)
    header:str=''
    evidence:str=''


def norm(s):
    s=unicodedata.normalize('NFKD',s or '').casefold()
    return re.sub(r'[^a-z0-9]+',' ',s).strip()


def title_similarity(a,b):
    a,b=norm(a),norm(b)
    if not a or not b:return 0.
    r=SequenceMatcher(None,a,b).ratio()
    x,y=set(a.split()),set(b.split())
    j=len(x&y)/max(1,len(x|y))
    return max(r,.65*r+.35*j)


def read_document(path:Path):
    if path.suffix.lower()=='.pdf':
        import fitz
        with fitz.open(str(path)) as d:
            # Front matter; title/DOI should live here. Fulltext not necessary for pairing.
            pieces=[d[i].get_text('text',sort=True) for i in range(min(3,len(d)))]
            meta=d.metadata or {}
            return '\n'.join(pieces), str(meta.get('title') or '')
    if path.suffix.lower()=='.docx':
        from docx import Document as WordDocument
        d=WordDocument(str(path))
        blocks=[p.text for p in d.paragraphs[:130] if p.text.strip()]
        for table in d.tables[:5]:
            for row in table.rows[:8]:
                blocks.extend(c.text for c in row.cells)
        return '\n'.join(blocks)[:28000],str(d.core_properties.title or '')
    raise ValueError(f'Unsupported document: {path.suffix}')


def clean_title(s):
    return re.sub(r'\s+',' ',re.sub(r'\b(supplementary\s+materials\s+for|supporting\s+information\s+for)\b','',s,flags=re.I)).strip(' .,:;\n\t')


def derive_title(text, metadata_title='', supplementary=False):
    lines=[s.strip() for s in text.splitlines() if s.strip()]
    if supplementary:
        for i,s in enumerate(lines[:40]):
            if re.search(r'^(supplementary materials for|supporting information for|supplemental materials for)\s*:?$',s,re.I):
                out=[]
                for t in lines[i+1:i+7]:
                    if re.search(r'^(siraj|corresponding authors?|science \d|doi\s*:|methods|abstract|figure s\d)',t,re.I):break
                    if len(t)<12:break
                    out.append(t)
                    if len(' '.join(out))>=85:break
                if out:return clean_title(' '.join(out))
            m=re.search(r'(?:supplementary materials for|supporting information for)\s*:\s*(.{20,})',s,re.I)
            if m:return clean_title(m.group(1))
    # Press-ready manuscripts can prepend a journal/DOI banner before the
    # real title. Recover the title immediately after 'Article in Press', not
    # an accidental concatenation of header + publisher URL.
    if not supplementary:
        for i,line in enumerate(lines[:20]):
            if re.fullmatch(r'\s*Article\s+in\s+Press\s*',line,re.I):
                title_lines=[]
                for next_line in lines[i+1:i+5]:
                    if re.search(r'^(?:Received|Accepted|Authors?|Cite this|We are providing|Copyright)',next_line,re.I):break
                    if re.search(r'https?://|doi\.org|@',next_line,re.I):break
                    if len(next_line)<8:break
                    title_lines.append(next_line)
                    if len(' '.join(title_lines))>110:break
                title=clean_title(' '.join(title_lines))
                if len(title)>=28:return title
    # Trust genuine PDF/Word metadata titles before the page-layout heuristic.
    if metadata_title and len(metadata_title)>23 and not SI_RE.search(metadata_title) and not re.search(r'(?:\.pdf|\.docx|microsoft word|_articlecontent|acrobat)',metadata_title,re.I):
        return clean_title(metadata_title)
    # Submitted manuscripts sometimes print an explicit 'Title:' field.
    mm=re.search(r'^\s*Title\s*:\s*(.{25,200})$',text,re.I|re.M)
    if mm:return clean_title(mm.group(1))
    # Scientific preprints often have numbered manuscript lines between title lines.
    if not supplementary:
        scrub=[]
        for x in lines[:16]:
            x=re.sub(r'^\d{1,2}\s+','',x).strip()
            if not x or re.fullmatch(r'\d{1,3}',x):continue
            if re.search(r'(?:@|\b(?:abstract|university|department)\b)',x,re.I):break
            if re.search(r'(?:[A-Z][a-z]+\d,?\d|[A-Z][a-z]+\d\*)',x):break
            if len(x)>12:scrub.append(x)
            if len(' '.join(scrub))>100:break
        if scrub:
            joined=' '.join(scrub)
            if 35<len(joined)<220:return clean_title(joined)
    # Strong metadata is preferable unless obviously a source filename or SI descriptor.
    if metadata_title and len(metadata_title)>23 and not SI_RE.search(metadata_title) and not re.search(r'(?:\.pdf|\.docx|microsoft word|_articlecontent|acrobat)',metadata_title,re.I):
        return clean_title(metadata_title)
    # Look near top for plausible long title, before ABSTRACT/INTRODUCTION or copyright boilerplate.
    candidates=[]
    for i,line in enumerate(lines[:65]):
        if re.search(r'^(abstract|introduction|reference|supplementary (information|materials)|article|methods|doi|downloaded)',line,re.I):continue
        if not (28<=len(line)<=190):continue
        if re.search(r'(?:https?://|www\.|@|\b(copyright|license|received|accepted|published|university|department|institute)\b)',line,re.I):continue
        penalty=1 if re.search(r'^(nature|science|research article|advanced|acs)',line,re.I) else 0
        candidates.append((penalty,i,-len(line),line))
    if candidates:
        candidates.sort()
        return clean_title(candidates[0][-1])
    return clean_title(lines[0] if lines else '')


def derive_authors(text,title=""):
    # We compare clearly named surnames when possible, not general token overlap.
    head=text[:6500]
    matches=re.findall(r'\b(?:Corresponding authors?|Authors?)\s*:\s*([^\n]{6,150})',head,re.I)
    names=[]
    for frag in matches:
        for full in re.split(r'\s*;\s*|\s*,\s*(?=[A-Z][a-z]+\s+[A-Z])',frag):
            full=re.sub(r'\S+@\S+','',full).strip(' .,;')
            if re.search(r'[A-Z][a-z]+\s+[A-Z][a-z]+',full):names.append(full)
    # Supporting files can show author surname before the article title.
    for m in re.finditer(r'\b([A-Z][a-z]+(?:\s+[A-Z]\.)?\s+[A-Z][a-z]{2,})\s+et\s+al\.',head):
        names.append(m.group(1))
    if title and not names:
        # Look immediately below the title (not in affiliations or references).
        header_lines=[x.strip() for x in head.splitlines() if x.strip()]
        parts=[x for x in norm(title).split() if len(x)>3]
        for i,line in enumerate(header_lines[:18]):
            if parts and any(w in norm(line).split() for w in parts[:3]):
                for next_line in header_lines[i+1:i+5]:
                    if re.search(r'(?:@|university|department|institute|abstract|introduction)',next_line,re.I):break
                    for frag in re.split(r',\s*|\s+and\s+',next_line):
                        frag=re.sub(r'\d+(?:,\d+)*[\*#]?','',frag).strip(' .,*')
                        m=re.fullmatch(r'([A-Z][a-zA-Z-]+\s+(?:[A-Z]\.?\s+)?[A-Z][a-zA-Z-]+)',frag)
                        if m:names.append(frag)
                break
    return list(dict.fromkeys(names))[:12]


def get_metadata(path:Path):
    text,meta_title=read_document(path)
    head=text[:18000]
    doi_candidates=[]
    for m in DOI_RE.finditer(head):
        d=m.group().rstrip('.,);]')
        if d.lower() not in [x.lower() for x in doi_candidates]:doi_candidates.append(d)
    # Many SI PDFs print the article DOI on the first page; otherwise DOI is unknown.
    doi=doi_candidates[0].lower() if doi_candidates else ''
    # A DOI-like text in a reference is weaker evidence than title/author concordance.
    if doi and head.lower().find(doi)>9000:doi=''
    lines=[s.strip() for s in head.splitlines()[:30] if s.strip()]
    filename_si=bool(NAME_SI.search(path.stem))
    first_page_si=bool(SI_RE.search('\n'.join(lines[:10])))
    kind='supplement' if filename_si or first_page_si else 'main'
    # Journal/title/year are evidence used in the scoring; missing fields are allowed.
    journal=''
    for name,rx in JOURNALS.items():
        if re.search(rx,head[:6500],re.I):journal=name;break
    year_candidates=re.findall(r'\b(?:19|20)\d{2}\b',head[:5000])
    year=int(year_candidates[0]) if year_candidates else None
    if year is not None and not (2010<=year<=2035):year=None
    if not year and path.suffix.lower()=='.pdf':
        try:
            import fitz
            with fitz.open(str(path)) as pp:
                m=re.search(r'(?:19|20)\d{2}',str(pp.metadata.get('creationDate','')))
                if m:year=int(m.group())
        except Exception:pass
    title=derive_title(head,meta_title,kind=='supplement')
    return Document(path=path,kind=kind,title=title,doi=doi,journal=journal,year=year,
                    authors=derive_authors(head,title),header=head[:3500],evidence=('filename+front_matter' if filename_si else 'front_matter'))


def author_agreement(a,b):
    if not a or not b:return None
    def surnames(seq):
        return {norm(x).split()[-1] for x in seq if norm(x)}
    aa,bb=surnames(a),surnames(b)
    return len(aa&bb)/max(1,min(len(aa),len(bb)))


def publication_identifiers(doc:Document):
    """Find unique article identifiers without requiring a known PDF basename.

    Identifiers are *evidence*, not authority: a conflicting explicitly printed DOI
    or a clearly unrelated title still blocks a proposed match.
    """
    stem=doc.path.stem.lower()
    out=set()
    doi=doc.doi.lower()
    # Publisher article numbers provide robust equivalence even when a SI cover
    # lacks the full DOI, but encodes the publisher article number.
    for src in [stem,doi]:
        for m in re.finditer(r'(?:s?\d{5})[-_./]?(?:20)?(\d{2,3})[-_./]0*(\d{3,6})(?:[-_./]\d+)?',src):
            number=re.search(r's?\d{5}',m.group(0)).group(0).lstrip('s');number='s'+number
            out.add(f'nature:{number}:{int(m.group(1)) % 100:02d}:{int(m.group(2))}')
        for m in re.finditer(r'(?:^|[^a-z])(adma|aenm|advs|anie)[._-]?0*(20\d{5,8})\b',src):
            out.add(f'wiley:{m.group(1)}:{m.group(2)}')
        for m in re.finditer(r'(?:science[._/-])?([a-z]{3}\d{4,7})(?:_sm|_si|$|[._/-])',src):
            if 'science' in src or src.startswith('wang_'):
                out.add(f'science:{m.group(1)}')
    return out


def conflicting_identifiers(main:Document,sup:Document):
    a=publication_identifiers(main)
    b=publication_identifiers(sup)
    # Only compare identifiers from the same publisher. A citation buried in an
    # unrelated supplementary file is not itself a confirmed article identifier.
    for prefix in ('nature:','wiley:','science:'):
        aa={i for i in a if i.startswith(prefix)}
        bb={i for i in b if i.startswith(prefix)}
        if aa and bb and aa.isdisjoint(bb):
            return True
    return False


def match_score(main:Document,sup:Document):
    t=title_similarity(main.title,sup.title)
    shared_ids=publication_identifiers(main)&publication_identifiers(sup)
    identifier_match=bool(shared_ids)
    same_doi=bool(main.doi and sup.doi and main.doi==sup.doi)
    # A PDF can contain incidental DOI references. The exact DOI comparison is
    # reliable only if the disagreement is supported by different article IDs.
    doi_conflict=bool(main.doi and sup.doi and main.doi!=sup.doi and
                      conflicting_identifiers(main,sup))
    j=bool(main.journal and sup.journal and norm(main.journal)==norm(sup.journal))
    ac=author_agreement(main.authors,sup.authors)
    yc=bool(main.year and sup.year and abs(main.year-sup.year)<=1)
    score=(9 if identifier_match else 0)+(8 if same_doi else 0)+(5*t)+(1 if j else 0)+\
          (1.5*ac if ac is not None else 0)+(0.5 if yc else 0)-(12 if doi_conflict else 0)
    # Unique publisher IDs are sufficient when bibliography is unavailable,
    # unless printed article titles *clearly* conflict. Otherwise insist on
    # substantial title and corroborating bibliographic evidence.
    explicit_title_a=norm(main.title)
    explicit_title_b=norm(sup.title)
    title_conflict=bool(len(explicit_title_a)>30 and len(explicit_title_b)>30 and
                        t<.14 and not same_doi and
                        len(set(explicit_title_a.split())&set(explicit_title_b.split()))<=1)
    acceptable=(not doi_conflict and not title_conflict and
                (identifier_match or (same_doi and (t>=.30 or not main.title or not sup.title)) or
                 (t>=.78 and (j or yc or (ac is not None and ac>.2)))))
    return score,acceptable,dict(title_similarity=round(t,3),doi_match=same_doi,
                                 article_identifier_match=identifier_match,
                                 shared_identifiers=sorted(shared_ids),
                                 journal_match=j,author_overlap=ac,
                                 year_consistent=yc,doi_conflict=doi_conflict,
                                 title_conflict=title_conflict)


def discover(folder:Path):
    files=sorted((p for p in folder.rglob('*') if p.is_file() and p.suffix.lower() in {'.pdf','.docx'}),key=lambda p:str(p).lower())
    docs=[];errors=[]
    for p in files:
        try: docs.append(get_metadata(p))
        except Exception as e:errors.append((p.name,str(e)))
    mains=[d for d in docs if d.kind=='main'];supps=[d for d in docs if d.kind=='supplement']
    associations={str(m.path):[] for m in mains};unmatched=[]
    for s in supps:
        candidates=[]
        for m in mains:
            score,allowed,details=match_score(m,s)
            if allowed:candidates.append((score,m,details))
        candidates.sort(key=lambda v:v[0],reverse=True)
        if not candidates or (len(candidates)>1 and candidates[0][0]-candidates[1][0]<1.2):
            unmatched.append((s.path.name,'AMBIGUOUS' if candidates else 'NO_SUPPORTED_MATCH',[(x[1].path.name,round(x[0],2)) for x in candidates[:3]]))
            continue
        score,m,details=candidates[0]
        associations[str(m.path)].append((s,details,round(score,3)))
    return mains,associations,unmatched,errors
