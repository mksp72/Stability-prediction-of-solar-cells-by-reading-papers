"""Optional DOCX-to-PDF figure bridge, preserving the original Word source.

On Windows requires Microsoft Word + pywin32 (optional), otherwise uses
LibreOffice on platforms where it is already installed. No conversion is
claimed if the necessary application is missing. Original files are unchanged.
"""
import os, shutil, subprocess
from pathlib import Path

def convert(path,outdir):
 path=Path(path).resolve();outdir=Path(outdir).resolve();outdir.mkdir(parents=True,exist_ok=True)
 dest=outdir/(path.stem+'.pdf')
 if dest.exists():
  raise FileExistsError(f'Conversion name collision: {dest.name}')
 if os.name=='nt':
  try:
   import win32com.client
  except ImportError as exc:
   raise RuntimeError('For DOCX figure digitization on Windows, install pywin32 and Microsoft Word') from exc
  app=None;doc=None
  try:
   app=win32com.client.DispatchEx('Word.Application')
   app.Visible=False;app.DisplayAlerts=0
   doc=app.Documents.Open(str(path),ReadOnly=True,AddToRecentFiles=False,Visible=False)
   doc.ExportAsFixedFormat(str(dest),17)
  finally:
   if doc is not None:doc.Close(False)
   if app is not None:app.Quit()
 else:
  executable=shutil.which('libreoffice') or shutil.which('soffice')
  if not executable:
   raise RuntimeError('DOCX image digitization requires a Word-to-PDF converter; no LibreOffice executable available')
  completed=subprocess.run([executable,'-env:UserInstallation=file:///tmp/psc_pdf_bridge_profile',
    '--headless','--convert-to','pdf','--outdir',str(outdir),str(path)],
    check=False,capture_output=True,text=True,timeout=180)
  if completed.returncode:
   raise RuntimeError('DOCX converter failed: '+completed.stderr[-500:])
 if not dest.is_file() or dest.stat().st_size<1000:
  raise RuntimeError('DOCX converter did not create a usable PDF')
 return dest


def replace_paths(documents,associations,outdir):
 """Return modified in-memory document paths with immutable source originals."""
 import hashlib
 warnings=[];mapping={}
 all_docs={id(d):d for d in documents}
 for items in associations.values():
  for supplement,ev,score in items:all_docs[id(supplement)]=supplement
 for d in all_docs.values():
  if d.path.suffix.lower()!='.docx':continue
  original=d.path
  try:
   unique=hashlib.sha1(str(original.resolve()).encode()).hexdigest()[:9]
   sub=Path(outdir)/unique
   target=convert(original,sub)
   d.path=target
   mapping[str(original)]=str(target)
  except Exception as exc:warnings.append((original.name,str(exc)))
 new_associations={}
 for orig,group in associations.items():
  new_associations[mapping.get(orig,orig)]=group
 return new_associations,warnings
