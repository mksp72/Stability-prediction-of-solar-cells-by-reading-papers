"""Reusable, publisher-independent stability axis semantics.

Never infer an axis solely from its numerical magnitude. Printed wording,
neighboring shared panel labels, and plot geometry jointly determine units.
Raw values are preserved; normalization is an additional derived column.
"""
from __future__ import annotations
import re, datetime as dt
import numpy as np

MONTHS={'jan':1,'january':1,'feb':2,'february':2,'mar':3,'march':3,
'apr':4,'april':4,'may':5,'jun':6,'june':6,'jul':7,'july':7,
'aug':8,'august':8,'sep':9,'sept':9,'september':9,'oct':10,'october':10,
'nov':11,'november':11,'dec':12,'december':12}

def date_value(s):
    s=str(s).strip().replace('–','-').replace('.','/').strip(',')
    for fmt in ('%d/%m/%Y','%d/%m/%y','%Y/%m/%d','%m/%Y','%b %Y','%B %Y', '%d %b %Y','%d %B %Y'):
        try:
            d=dt.datetime.strptime(s,fmt)
            # Month/year timepoints represent month starts; keep this explicit.
            return d.date().isoformat()
        except ValueError:pass
    m=re.fullmatch(r'(\d{1,2})/(\d{1,2})',s)
    if m and 1<=int(m[1])<=31 and 1<=int(m[2])<=12:
        return None # no year available: cannot convert across year boundary
    return None

def axis_kind(y_numbers, nearby='', caption=''):
    """Return kind/confidence, keeping ambiguous scales explicitly ambiguous."""
    a=np.asarray([float(v) for v in y_numbers],float)
    if not len(a):return 'UNKNOWN','UNRESOLVED'
    lo,hi=float(np.min(a)),float(np.max(a))
    txt=(nearby+' '+caption).lower()
    if re.search(r'\b(?:current density|jsc|ma\s*/?\s*cm|voltage|wavelength|eqe|ipce)\b',nearby,re.I):
        return 'NON_PCE_AXIS','HIGH'
    if re.search(r'\b(?:normalized|normalised|norm\.|relative|retention)\b',nearby,re.I):
        return ('NORMALIZED_FRACTION' if hi<=1.6 else 'RETENTION_PERCENT'),'HIGH'
    if re.search(r'\b(?:pce|efficiency|power conversion|s?po|output power)\b',nearby,re.I):
        if hi<=1.6:return 'NORMALIZED_FRACTION','MEDIUM'
        if 1.6<hi<=40:return 'ABSOLUTE_PCE_PERCENT','HIGH'
        if 40<hi<=110:return 'RETENTION_PERCENT','MEDIUM'
    if 0<=lo and hi<=1.6:return 'NORMALIZED_FRACTION','MEDIUM'
    if 1.6<hi<=35:return 'ABSOLUTE_PCE_PERCENT','MEDIUM'
    if 35<hi<=110:return 'RETENTION_PERCENT','MEDIUM'
    return 'UNKNOWN','UNRESOLVED'

def time_kind(tokens,shared_caption=''):
    text=' '.join(map(str,tokens))+' '+shared_caption
    if re.search(r'(?i)\b(?:ns|nanoseconds|microseconds|µs|ms|milliseconds|seconds|minutes)\b',text):
        return 'NON_AGING_TIME','HIGH'
    if re.search(r'(?i)\b(?:hours?|hrs?|aging\s*time\s*\(h\)|time\s*\(h\))\b',text):
        return 'HOURS','HIGH'
    if re.search(r'(?i)\b(?:days?|time\s*\(d\))\b',text):return 'DAYS','HIGH'
    if re.search(r'(?i)\b(?:months?)\b',text):return 'MONTHS','HIGH'
    if any(date_value(t) for t in tokens):return 'CALENDAR_DATE','HIGH'
    return 'UNRESOLVED','LOW'

def to_hours(raw_values,kind):
    """For calendar ticks use elapsed time from first observation; dates are retained separately."""
    if kind=='HOURS':return [float(x) for x in raw_values]
    if kind=='DAYS':return [float(x)*24 for x in raw_values]
    if kind=='CALENDAR_DATE':
        dates=[dt.date.fromisoformat(date_value(x)) for x in raw_values]
        return [(x-dates[0]).days*24. for x in dates]
    return None

def retention(raw,kind):
    a=np.asarray(raw,float)
    if kind=='NORMALIZED_FRACTION':return list(100*a),'PLOTTED_FRACTION_0_TO_1'
    if kind=='RETENTION_PERCENT':return list(a),'PLOTTED_PERCENT_RETENTION'
    if kind=='ABSOLUTE_PCE_PERCENT' and len(a) and a[0]>0:
        return list(100*a/a[0]),'FIRST_MEASURED_PCE_REVIEW_REFERENCE'
    return None,'UNRESOLVED'

def calendar_ticks(words, nearby_text=''):
    """Date labels with geometric coordinates; dates without year use explicit
    printed year from same caption/page only, not the current system year."""
    items=[]
    text_years=re.findall(r'\b(?:19|20)\d{2}\b',nearby_text)
    common_year=int(text_years[0]) if len(set(text_years))==1 else None
    sw=sorted(words,key=lambda w:(round(w[1]/5),w[0]))
    for i,w in enumerate(sw):
        token=str(w[4]);date=date_value(token)
        if date is None and common_year:
            m=re.fullmatch(r'(\d{1,2})/(\d{1,2})',token)
            if m:
                try:date=dt.date(common_year,int(m[2]),int(m[1])).isoformat()
                except ValueError:pass
        if date is None and i+1<len(sw):
            nxt=sw[i+1]
            if abs(nxt[1]-w[1])<5 and 0<nxt[0]-w[2]<20:
                date=date_value(token+' '+str(nxt[4]))
        if date:items.append((float((w[0]+w[2])/2),float((w[1]+w[3])/2),date,token))
    return items

def calendar_axis_candidates(words,caption=''):
    """Find adjacent-date X tick rows and convert actual dates to elapsed hours.
    Returns an audit candidate, not an automatic scientific curve certification."""
    ticks=calendar_ticks(words,caption)
    bins=[]
    for tick in ticks:
        g=next((b for b in bins if abs(np.median([x[1] for x in b])-tick[1])<6),None)
        if g is None:bins.append([tick])
        else:g.append(tick)
    results=[]
    for group in bins:
        group=sorted(group)
        if len(group)<3 or group[-1][0]-group[0][0]<45:continue
        dates=[dt.date.fromisoformat(t[2]) for t in group]
        hrs=np.array([(d-min(dates)).days*24 for d in dates],float)
        xs=np.array([t[0] for t in group],float)
        if len(set(hrs))<3 or np.any(np.diff(hrs)<=0):continue
        m,b=np.polyfit(xs,hrs,1)
        if m<=0:continue
        errors=np.abs(xs*m+b-hrs)
        if np.max(errors)>max(24,0.035*np.ptp(hrs)):continue
        results.append(dict(x0=float(min(xs)),x1=float(max(xs)),bottom=float(np.median([t[1] for t in group])-5),
            x_m=float(m),x_b=float(b),x_max=float(max(hrs)),tick_count=len(group),
            x_unit='CALENDAR_DATE',reference_date=min(dates).isoformat(),
            axis_evidence='PRINTED_CALENDAR_TICKS'))
    return results
