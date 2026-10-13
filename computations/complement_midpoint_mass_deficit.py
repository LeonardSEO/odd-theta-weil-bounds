"""Actual lower masses and deficit of the stored upward-binary midpoint masses."""
from pathlib import Path
from fractions import Fraction
from math import nextafter, inf
import json,hashlib,time
from flint import arb,acb,ctx
ctx.prec=192
root=Path(__file__).resolve().parents[1]
cellpath=root/'certificates/complement_signed_symbol_norm_cells.json'
rankpath=root/'certificates/complement_midpoint_rank.json'
cells=json.loads(cellpath.read_text())['kernel_cells']
rank=json.loads(rankpath.read_text())
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
assert sha(cellpath)==rank['cell_certificate_sha256']
rows=[(arb(n).log(),2*arb(p).log()/arb(n).sqrt())
      for n,p in [(2,2),(3,3),(4,2),(5,5),(7,7),(8,2)]]
pi=arb.pi();logpi=pi.log();scale=2**80
start=time.monotonic()
def delta_interval(k,den):
    lo=arb(k)/den;hi=arb(k+1)/den
    c=(lo+hi)/2;h=(hi-lo)/2;X=c+arb(0,h)
    L=arb(10) if k==0 else min(arb(10),2/(lo*lo)+pi/(2*lo))
    m=acb(arb(1)/4,c/2).digamma().real-logpi+arb(0,L*h)
    return sum((coef*(a*X).cos() for a,coef in rows),arb(0))-m

def arbf(q):
    return arb(q.numerator)/arb(q.denominator)

recordcells=[];deficit=Fraction(0);masssum=Fraction(0);lower_sum=Fraction(0)
for cell in cells:
    k=int(cell['k']);lower=arb(0);upper=arb(0)
    for j in range(4):
        delta=delta_interval(4*k+j,256)
        low=delta.lower();high=delta.upper()
        lower += (low if low>0 else arb(0))/(256*pi)
        upper += (high if high>0 else arb(0))/(256*pi)
    ln=int((lower.lower()*scale).floor().unique_fmpz())
    un=int((upper.upper()*scale).ceil().unique_fmpz())
    lf=Fraction(ln,scale)
    mf=Fraction(*nextafter(int(cell['measure_upper_dyadic_numerator'])/scale,inf).as_integer_ratio())
    assert mf>=lf
    masssum+=mf;lower_sum+=lf;deficit+=mf-lf
    recordcells.append({'k':k,'lower_mass_dyadic_numerator':str(ln),
                        'new_upper_mass_dyadic_numerator':str(un),
                        'midpoint_mass_exact':str(mf),'deficit_upper_exact':str(mf-lf)})
# A comfortable rational theorem endpoint, selected and compared exactly.
dden=10**8;dnum=(deficit.numerator*dden)//deficit.denominator+1
D=Fraction(dnum,dden)
assert deficit<D and masssum<14
print('PASS all2905cells true mass lower',float(lower_sum),
      'stored upward-binary mass',float(masssum),flush=True)
print('PASS whole actual mass deficit <',str(D),'diagnostic',float(deficit),flush=True)
record={'status':'PASS actual full-cell lower measures and upward-binary midpoint deficit',
 'source':'routes/complement-midpoint-mass-deficit.md','precision_bits':ctx.prec,
 'kernel_cells':len(cells),'new_fine_intervals':4*len(cells),
 'fine_width':'1/256','measure':'nu(dt)=r8(t)dt/pi on positive half-line',
 'lower_mass_exact':str(lower_sum),'stored_binary_mass_exact':str(masssum),
 'deficit_upper_exact':str(deficit),'deficit_theorem_upper_exact':str(D),
 'per_cell_dyadic_scale':str(scale),'cells':recordcells,
 'inherited_zero_support':'all omitted kernel cells and t>10000 have r8=0 by saved complete certificate/proof',
 'cell_certificate_sha256':sha(cellpath),'rank_certificate_sha256':sha(rankpath),
 'script_sha256':sha(Path(__file__)),'elapsed_seconds':time.monotonic()-start}
(root/'certificates/complement_midpoint_mass_deficit.json').write_text(json.dumps(record,indent=2)+'\n')
