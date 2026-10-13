"""New whole operator norm question for the actual midpoint mass deficit."""
from pathlib import Path
from fractions import Fraction
import json,hashlib,time
import numpy as np
from flint import arb,ctx
ctx.prec=192
root=Path(__file__).resolve().parents[1]
mpath=root/'certificates/complement_midpoint_mass_deficit.json'
fpath=root/'certificates/complement_fourier_weight_grid.json'
cpath=root/'certificates/complement_signed_symbol_norm_cells.json'
rpath=root/'certificates/complement_midpoint_rank.json'
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
r=json.loads(rpath.read_text());m=json.loads(mpath.read_text());f=json.loads(fpath.read_text())
assert sha(fpath)==r['Fourier_certificate_sha256']
assert sha(cpath)==r['cell_certificate_sha256']==m['cell_certificate_sha256']
scale=int(f['endpoint_scale']);unit=Fraction(1,2**53)
def gamma(n):return n*unit/(1-n*unit)
def upperfloat(x):return np.nextafter(float(x),np.inf)
def pairs(name):
    return (np.array([np.nextafter(int(a)/scale,-np.inf) for a,b in f[name]]),
            np.array([np.nextafter(int(b)/scale,np.inf) for a,b in f[name]]))
Wlo,Whi=pairs('W_outward_endpoint_pairs')
blo,bhi=pairs('score_overlap_outward_endpoint_pairs')
k=np.array([a['k'] for a in m['cells']],dtype=np.int64)
d_exact=[Fraction(a['deficit_upper_exact']) for a in m['cells']]
d=np.array([upperfloat(v) for v in d_exact])
N=len(k);assert N==2905 and np.all(d>0)
for di,exact in zip(d,d_exact):assert Fraction(*di.as_integer_ratio())>=exact
start=time.monotonic()
far=np.zeros(int(2*k.max()+2))
for j in range(6401,len(far)):
    far[j]=upperfloat((5*(-arb(j)/256).exp()).upper())
bslo=np.empty(N);bshi=np.empty(N)
for i,ki in enumerate(k):
    if ki<6400:bslo[i],bshi[i]=blo[ki],bhi[ki]
    else:
        t=arb(2*int(ki)+1)/128
        eb=upperfloat((5*t*(-t/4).exp()).upper())
        bslo[i],bshi[i]=-eb,eb
D=np.empty((N,N))
for i in range(N):
    diff=np.abs(k[i]-k);summ=k[i]+k+1
    def wb(idx):
        lo=np.zeros(N);hi=np.zeros(N);near=idx<=6400
        lo[near]=Wlo[idx[near]];hi[near]=Whi[idx[near]]
        lo[~near]=-far[idx[~near]];hi[~near]=far[idx[~near]]
        return lo,hi
    dl,dh=wb(diff);sl,sh=wb(summ)
    prod=np.stack([bslo[i]*bslo,bslo[i]*bshi,bshi[i]*bslo,bshi[i]*bshi])
    gl=(dl-sh)/2-prod.max(axis=0)/12
    gh=(dh-sl)/2-prod.min(axis=0)/12
    D[i,:]=np.maximum(np.abs(gl),np.abs(gh))+1e-12
assert np.isfinite(D).all() and np.max(D)<2
# Heuristic only selects positive exact dyadic Schur weights.
w=np.ones(N)
for step in range(120):
    y=D@(d*w);w=y/y.max()
w=np.maximum(w,1/1024)
nums=np.ceil(w*2**30).astype(np.int64)
w=nums.astype(np.float64)/2**30
assert np.all(nums>0)
rows=D@(d*w)/w
# gamma includes product, dense positive dot, and division. Relaxed
# rational endpoint absorbs final scalar bound evaluation as well.
g=float(gamma(2*N+20))
guarded=upperfloat(float(rows.max())/(1-g))
exactguard=Fraction(*guarded.as_integer_ratio())+Fraction(1,10**15)
den=10**8;num=(exactguard.numerator*den)//exactguard.denominator+1
endpoint=Fraction(num,den)
assert endpoint>exactguard
print('PASS full2905-frequency deficit norm <',str(endpoint),
      'diagnosticguarded',guarded,'maxrow',int(rows.argmax()),flush=True)
record={'status':'PASS complete midpoint mass-deficit operator Schur certificate',
 'source':'routes/complement-midpoint-mass-deficit.md MD4-MD6','precision_bits':ctx.prec,
 'dimension':N,'pole_deficit':0,'score_Gram_upper':'G0-bb*/12 >= actual G0-bb*/Q',
 'norm_upper_exact':str(endpoint),'guarded_float_upper':guarded,
 'exact_postscalar_allowance':'1/1000000000000000','IEEE_dot_gamma':g,
 'IEEE_operations_bound':2*N+20,'per_entry_IEEE_allowance':'1e-12',
 'weight_denominator':2**30,'weight_numerators':nums.tolist(),
 'mass_deficit_certificate_sha256':sha(mpath),'Fourier_certificate_sha256':sha(fpath),
 'cell_certificate_sha256':sha(cpath),'rank_certificate_sha256':sha(rpath),
 'all_coupling_retained':'all2905 squared frequency entries; pole difference exactlyzero',
 'all_frequency_tail':'5exp(-abs(nu)/4), score5abs(nu)exp(-abs(nu)/4), Arb outward',
 'script_sha256':sha(Path(__file__)),'elapsed_seconds':time.monotonic()-start}
(root/'certificates/complement_midpoint_deficit_norm.json').write_text(json.dumps(record,indent=2)+'\n')
