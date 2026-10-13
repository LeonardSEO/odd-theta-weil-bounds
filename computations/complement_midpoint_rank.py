"""New midpoint positive-price rank question; all source remainder retained."""
from pathlib import Path
import json
import hashlib
import time
from fractions import Fraction
import numpy as np
from flint import arb,ctx
from scipy.linalg import eigh

ctx.prec=160
root=Path(__file__).resolve().parents[1]
cellpath=root/'certificates/complement_signed_symbol_norm_cells.json'
fourierpath=root/'certificates/complement_fourier_weight_grid.json'
z=json.loads(cellpath.read_text());f=json.loads(fourierpath.read_text())
scale=int(f['endpoint_scale'])
def upperfloat(x):return np.nextafter(float(x),np.inf)
def pairs(name):
    p=f[name]
    lo=np.array([np.nextafter(int(a)/scale,-np.inf) for a,b in p])
    hi=np.array([np.nextafter(int(b)/scale,np.inf) for a,b in p])
    return lo,hi
Wlo,Whi=pairs('W_outward_endpoint_pairs')
vlo,vhi=pairs('pole_outward_endpoint_pairs')
blo,bhi=pairs('score_overlap_outward_endpoint_pairs')
cells=z['kernel_cells'];k=np.array([a['k'] for a in cells],dtype=np.int64)
mu=np.array([upperfloat(int(a['measure_upper_dyadic_numerator'])/(2**80)) for a in cells])
mu=np.append(mu,2.0);sm=np.sqrt(mu);n=len(k);size=n+1
start=time.monotonic()
far=np.zeros(int(2*k.max()+2))
for j in range(6401,len(far)):
    far[j]=upperfloat((5*(-arb(j)/256).exp()).upper())
bslo=np.zeros(n);bshi=np.zeros(n);vslo=np.zeros(n);vshi=np.zeros(n)
for i,ki in enumerate(k):
    if ki<6400:
        bslo[i],bshi[i]=blo[ki],bhi[ki]
        vslo[i],vshi[i]=vlo[ki],vhi[ki]
    else:
        t=arb(2*int(ki)+1)/128
        eb=upperfloat((5*t*(-t/4).exp()).upper())
        ev=upperfloat((5*(-t/4).exp()).upper())
        bslo[i],bshi[i]=-eb,eb
        vslo[i],vshi[i]=-ev,ev

A=np.empty((size,size));errors=np.empty((size,size))
for i in range(n):
    diff=np.abs(k[i]-k);summ=k[i]+k+1
    def wb(idx):
        lo=np.zeros(n);hi=np.zeros(n);near=idx<=6400
        lo[near]=Wlo[idx[near]];hi[near]=Whi[idx[near]]
        lo[~near]=-far[idx[~near]];hi[~near]=far[idx[~near]]
        return lo,hi
    dl,dh=wb(diff);sl,sh=wb(summ)
    prod=np.stack([bslo[i]*bslo,bslo[i]*bshi,bshi[i]*bslo,bshi[i]*bshi])
    low=((dl-sh)/2-prod.max(axis=0)/12)*sm[i]*sm[:n]
    high=((dh-sl)/2-prod.min(axis=0)/12)*sm[i]*sm[:n]
    A[i,:n]=(low+high)/2
    errors[i,:n]=(high-low)/2+1e-12
lo=(vslo+bslo/48)*sm[:n]*sm[n]
hi=(vshi+bshi/48)*sm[:n]*sm[n]
A[:n,n]=A[n,:n]=(lo+hi)/2
errors[:n,n]=errors[n,:n]=(hi-lo)/2+1e-12
A[n,n]=float(Fraction(3,500)-Fraction(1,192))*2
errors[n,n]=1e-12
# Symmetrization preserves the outward error norm price.
A=(A+A.T)/2;errors=(errors+errors.T)/2
delta=upperfloat(errors.sum(axis=1).max()/(1-(2*size+10)*2**-53))
print('Complete midpoint Gram entry-error upper',delta,'dimension',size,flush=True)
values,U=eigh(A,check_finite=False)
values=values[::-1];U=U[:,::-1]
base=11*Fraction(42,125*128**2)+Fraction(32237,10**12)
threshold=(Fraction(1,1000)-base)/Fraction(11,10)
count=int(np.count_nonzero(values>float(threshold)-delta))
print('UNVERIFIED eigenvalue locator: threshold',float(threshold),'count',count,
      'largest',values[0],'next',values[count],flush=True)
# A low-rank approximation can be used without assuming exact eigenvectors:
# its entire residual is bounded by the Frobenius norm, including all mixed
# retained/discarded coordinates. Rank stays at most chosen rank.
for rank in [count,count+10,count+25,count+50,count+100]:
    if rank>=size:continue
    approx=(U[:,:rank]*values[:rank])@U[:,:rank].T
    residual=A-approx
    fro=float(np.linalg.norm(residual,'fro'))
    # Conservative error for each matrix-product entry. All |A|<1 and
    # orthonormal locator entries are <=1; even a worst absolute sum
    # rank*max|lambda| is priced, independently of solver orthogonality.
    maxu=float(np.max(np.abs(U[:,:rank])))
    maxv=float(np.max(np.abs(values[:rank])))
    unit=2**-53;gamma=(2*rank+10)*unit/(1-(2*rank+10)*unit)
    roundprice=size*gamma*rank*maxu*maxu*maxv+1e-9
    remainder=upperfloat(fro+delta+roundprice)
    print('rank',rank,'full Frobenius residual',fro,'guarded remainder',remainder,flush=True)
    if remainder<float(threshold):
        selected=rank;certrem=remainder;selectedU=U[:,:rank];selectedv=values[:rank]
        break
else:
    raise RuntimeError('No tested rank certifies the complete target residual')

assert np.max(np.abs(A))<1 and np.max(np.abs(selectedU))<=1
assert np.all(selectedv>0)
# Preserve the actual binary64 factors as exact dyadic numbers in an npz.
factorpath=root/'certificates/complement_midpoint_rank_factors.npz'
np.savez_compressed(factorpath,U=selectedU,values=selectedv)
loss=upperfloat(float(base)+1.1*certrem)
assert loss<.001
record={'status':'PASS new actual midpoint rank residual calculation',
 'source':'routes/complement-midpoint-rank.md','dimension':size,
 'rank':selected,'gram_interval_error_upper':delta,
 'full_frobenius_residual_upper':certrem,'complete_loss_upper':loss,
 'base_price_exact':str(base),'midpoint_multiplier':'11/10',
 'threshold_exact':str(threshold),'locator_count':count,
 'all_coupling_retained':'full dense A-minus-rankN residual',
 'float_factor_interpretation':'stored binary64 factors are exact dyadic construction, not certified eigenvectors',
 'script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
 'cell_certificate_sha256':hashlib.sha256(cellpath.read_bytes()).hexdigest(),
 'Fourier_certificate_sha256':hashlib.sha256(fourierpath.read_bytes()).hexdigest(),
 'factor_file_sha256':hashlib.sha256(factorpath.read_bytes()).hexdigest(),
 'elapsed_seconds':time.monotonic()-start}
(root/'certificates/complement_midpoint_rank.json').write_text(json.dumps(record,indent=2)+'\n')
print('PASS candidate rank',selected,'complete loss',loss,flush=True)
