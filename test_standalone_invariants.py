"""Self-contained, synthetic checks: no dependency on any research article."""
from pathlib import Path
import ast, re
from hysteresis_recovery import canonical_scan_label, calculate_hi
from device_identity_pv import apply

BASE=Path(__file__).resolve().parent

def run():
    for p in BASE.glob('*.py'):
        if p.name.startswith('test_'):continue
        src=p.read_text(encoding='utf-8')
        ast.parse(src,filename=str(p))
        # No literal real research DOI or curated ID tags anywhere in the program.
        assert not re.search(r"10\.(?:1002|1016|1021|1038|1126)/[\w.\-]+",src),p.name
        assert not re.search(r"\bID(?:6|18|19|21|25|29)\b",src),p.name
    for suffix,direction in [('F','forward'),('Fw','forward'),('FWD','forward'),
                              ('R','reverse'),('Re','reverse'),('Rw','reverse'),
                              ('Bw','reverse'),('Backward','reverse')]:
        label,found=canonical_scan_label('Fabrication SnO2 '+suffix)
        assert label=='Fabrication SnO2' and found==direction,(label,found)
    assert canonical_scan_label('FAPbI3')[0]=='FAPbI3'
    signed,absolute=calculate_hi(18.,18.3)
    assert abs(signed+1.6666667)<0.0001 and abs(absolute-1.6666667)<0.0001
    # No publication fingerprint, filename, figure number, known study values.
    row={'device_category':'Fabrication SnO2','pce0':None,'voc0':None,'jsc0':None,'ff0':None}
    candidates=[]
    for direction,pce,voc,jsc,ff in [('forward',18.3,1.17,19.8,78.8),('reverse',18.,1.19,19.8,76.5)]:
        candidates.append({'device_label':'Fabrication SnO2 '+ ('F' if direction=='forward' else 'Bw'),
          'device_canonical':'Fabrication SnO2','scan_direction':direction,
          'source_file':'synthetic_supplement.pdf','source_page':3,'area_cm2':.25,
          'method':'same_row_table','values':dict(pce0=pce,voc0=voc,jsc0=jsc,ff0=ff)})
    apply([row],Path('synthetic_main.pdf'),(),candidates=candidates)
    assert row['pce0']==18. and row['voc0']==1.19 and row['ff0']==76.5,row
    assert row['jv_selected_scan_direction']=='reverse'
    assert row['hi_absolute_pct']==round(100*.3/18,6)
    print('PASS no article-specific constants, scan aliases, reverse priority, HI')
if __name__=='__main__':run()
