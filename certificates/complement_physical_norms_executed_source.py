"""Certify the actual stored 56-column physical frame, retaining Q and tau."""
from pathlib import Path
from fractions import Fraction
import hashlib
import json
import time
import numpy as np
from scipy.linalg import eigh
from flint import arb, acb, ctx

ctx.prec = 192
root = Path(__file__).resolve().parents[1]
paths = {
    'rank_script': root/'computations/complement_midpoint_rank.py',
    'rank_certificate': root/'certificates/complement_midpoint_rank.json',
    'cells': root/'certificates/complement_signed_symbol_norm_cells.json',
    'fourier': root/'certificates/complement_fourier_weight_grid.json',
    'factors': root/'certificates/complement_midpoint_rank_factors.npz',
}
sha = lambda p: hashlib.sha256(p.read_bytes()).hexdigest()
r = json.loads(paths['rank_certificate'].read_text())
expected = {'rank_script': r['script_sha256'], 'cells': r['cell_certificate_sha256'],
            'fourier': r['Fourier_certificate_sha256'], 'factors': r['factor_file_sha256']}
for key, digest in expected.items():
    assert sha(paths[key]) == digest, (key, 'input provenance mismatch')
f = json.loads(paths['fourier'].read_text())
z = json.loads(paths['cells'].read_text())
assert sha(root/'computations/complement_fourier_weight_grid.py') == f['script_sha256']
Z = np.load(paths['factors'])['U']
assert Z.shape == (2906,56) and np.isfinite(Z).all()
N, rank = Z.shape
start = time.monotonic()
unit = Fraction(1,2**53)
def gamma(operations):
    return operations*unit/(1-operations*unit)
def arbf(q):
    q = Fraction(q)
    return arb(q.numerator)/arb(q.denominator)
def upperfloat(x):
    return np.nextafter(float(x),np.inf)
def exactfloat(x):
    a,b = float(x).as_integer_ratio()
    return arb(a)/arb(b)
def upper_fraction(x):
    return Fraction(*float(x).as_integer_ratio())
def endpointpair(value, scale=2**100):
    return [int((value.lower()*scale).floor().unique_fmpz()),
            int((value.upper()*scale).ceil().unique_fmpz())]

def jet(x):
    phi = acb(0)
    dphi = acb(0)
    for n in range(1,9):
        y = arb.pi()*n*n*(2*x).exp()
        atom = (x/2-y).exp()
        phi += atom*(4*y*y-6*y)
        dphi += atom*(-8*y*y*y+30*y*y-15*y)
    return phi, dphi

def qintegrand(x, analytic):
    phi,dphi = jet(x)
    # A meromorphic callback: division gives a non-finite ball whenever
    # the complex evaluation domain contains a possible denominator zero.
    return dphi*dphi*(x/2).cosh()/phi

def tintegrand(x, analytic):
    phi,_ = jet(x)
    return phi*(x/2).sinh()**2/(x/2).cosh()

# The proof bounds the full theta correction on [0,2] by 2^-250
# and each two-sided physical tail beyond 2 by 2^-180.
assert arb(4).exp() > 54 and arb(4).exp() < 55
assert (-arb(240)).exp() < arbf(Fraction(1,2**340))
assert 2*(-arb(57)).exp() < arb(1)/2
poly = sum(arb(24)*arb(54)**(4-j)/arb([24,6,2,1,1][j])/arb(3)**(j+1)
           for j in range(5))
# Integral_54^infinity t^4 exp(-3t)dt, including both Q tails.
assert 512000*(-arb(162)).exp()*poly < arbf(Fraction(1,2**180))
full_norm_error = arbf(Fraction(1,2**250)+Fraction(1,2**180))
Qcalc = 2*acb.integral(qintegrand,0,2,abs_tol=arb(2)**-100,
                      rel_tol=arb(2)**-100,eval_limit=100000)
Tcalc = 2*acb.integral(tintegrand,0,2,abs_tol=arb(2)**-100,
                      rel_tol=arb(2)**-100,eval_limit=100000)
assert Qcalc.is_finite() and Tcalc.is_finite()
assert Qcalc.imag.contains(0) and Tcalc.imag.contains(0)
Q = Qcalc.real + arb(0,full_norm_error)
tau = Tcalc.real + arb(0,full_norm_error)
assert Q > 11 and Q < 12 and tau > arb(1)/200 and tau < arb(3)/500
print('PASS actual Q',Q,'actual tau',tau,flush=True)

# Certify the exact-binary coefficient Gram first, with no solver premise.
Wcenter = Z.T@Z
maxz = float(np.max(np.abs(Z)))
wentry = upperfloat(float(gamma(2*N+10))*N*maxz*maxz+1e-15)
Wdefect = upperfloat(float(np.max(np.sum(np.abs(Wcenter-np.eye(rank)),axis=1)))
                     /(1-float(gamma(2*rank+10)))+rank*wentry)
assert upper_fraction(Wdefect) < Fraction(1,10**8)
Wlower = Fraction(1,1)-Fraction(1,10**8)
Wupper = Fraction(1,1)+Fraction(1,10**8)
print('PASS coefficient Gram: (1-1e-8)I <= W <= (1+1e-8)I;',
      'actual guarded defect',Wdefect,flush=True)

scale = int(f['endpoint_scale'])
def pairs(name):
    p = f[name]
    return (np.array([np.nextafter(int(a)/scale,-np.inf) for a,b in p]),
            np.array([np.nextafter(int(b)/scale,np.inf) for a,b in p]))
Wlo,Whi = pairs('W_outward_endpoint_pairs')
vlo,vhi = pairs('pole_outward_endpoint_pairs')
blo,bhi = pairs('score_overlap_outward_endpoint_pairs')
cells = z['kernel_cells']
k = np.array([a['k'] for a in cells],dtype=np.int64)
mu = np.array([upperfloat(int(a['measure_upper_dyadic_numerator'])/(2**80))
               for a in cells]+[2.0])
sm = np.sqrt(mu)
n = len(k)
assert n == N-1 and len(k)==2905
# Mass interpretation exactly reproduces the original rank construction.
for m,a in zip(mu[:n],cells):
    assert upper_fraction(m) >= Fraction(int(a['measure_upper_dyadic_numerator']),2**80)
assert sum(map(upper_fraction,mu[:n])) < 14
qlo = np.nextafter(float(Q.lower()),-np.inf)
qhi = np.nextafter(float(Q.upper()),np.inf)
tlo = np.nextafter(float(tau.lower()),-np.inf)
thi = np.nextafter(float(tau.upper()),np.inf)
assert qlo > 11 and qhi < 12
far = np.zeros(int(2*k.max()+2))
for j in range(6401,len(far)):
    far[j] = upperfloat((5*(-arb(j)/256).exp()).upper())
bslo = np.zeros(n); bshi = np.zeros(n)
vslo = np.zeros(n); vshi = np.zeros(n)
for i,ki in enumerate(k):
    if ki < 6400:
        bslo[i],bshi[i] = blo[ki],bhi[ki]
        vslo[i],vshi[i] = vlo[ki],vhi[ki]
    else:
        t = arb(2*int(ki)+1)/128
        eb = upperfloat((5*t*(-t/4).exp()).upper())
        ev = upperfloat((5*(-t/4).exp()).upper())
        bslo[i],bshi[i] = -eb,eb
        vslo[i],vshi[i] = -ev,ev

A = np.empty((N,N)); errors = np.empty((N,N))
for i in range(n):
    diff = np.abs(k[i]-k); summ = k[i]+k+1
    def wb(idx):
        lo=np.zeros(n); hi=np.zeros(n); near=idx<=6400
        lo[near]=Wlo[idx[near]];hi[near]=Whi[idx[near]]
        lo[~near]=-far[idx[~near]];hi[~near]=far[idx[~near]]
        return lo,hi
    dl,dh=wb(diff);sl,sh=wb(summ)
    prod=np.stack([bslo[i]*bslo,bslo[i]*bshi,bshi[i]*bslo,bshi[i]*bshi])
    pmin=prod.min(axis=0);pmax=prod.max(axis=0)
    caplo=np.minimum(pmin/qlo,pmin/qhi)
    caphi=np.maximum(pmax/qlo,pmax/qhi)
    low=((dl-sh)/2-caphi)*sm[i]*sm[:n]
    high=((dh-sl)/2-caplo)*sm[i]*sm[:n]
    A[i,:n]=(low+high)/2
    errors[i,:n]=(high-low)/2+1e-12
lo=(vslo+np.minimum(bslo/(4*qlo),bslo/(4*qhi)))*sm[:n]*sm[n]
hi=(vshi+np.maximum(bshi/(4*qlo),bshi/(4*qhi)))*sm[:n]*sm[n]
A[:n,n]=A[n,:n]=(lo+hi)/2
errors[:n,n]=errors[n,:n]=(hi-lo)/2+1e-12
plow=2*(tlo-1/(16*qlo));phigh=2*(thi-1/(16*qhi))
A[n,n]=(plow+phigh)/2
errors[n,n]=(phigh-plow)/2+1e-12
A=(A+A.T)/2;errors=(errors+errors.T)/2
assert np.max(np.abs(A))<1 and np.max(sm)<2 and maxz<1
Gerror=upperfloat(errors.sum(axis=1).max()/(1-float(gamma(2*N+10))))

# Entire dense compression, including both frequency/pole mixed blocks.
AZ=A@Z
Hcenter=Z.T@AZ
Hcenter=(Hcenter+Hcenter.T)/2
# gamma_(2N+10) accounts for each dense dot; arbitrary product order
# obeys the conservative sum bound. The extra 1e-12 absorbs scalar
# guard evaluation, subtraction and symmetrization rounding.
first_entry=float(gamma(2*N+10))*N*float(np.max(np.abs(A)))*maxz
second_entry=float(gamma(2*N+10))*N*maxz*float(np.max(np.abs(AZ)))
compression_entry=N*maxz*first_entry+second_entry+1e-12
Herror=upperfloat(float(Wupper)*Gerror+rank*compression_entry)
print('PASS actual full Gram/compression error',Herror,flush=True)

# Eigensolver output only locates a preconditioner. Store that exact
# binary matrix, then certify its orthogonality and the complete residual.
values,V=eigh(Hcenter,check_finite=False)
assert np.isfinite(values).all() and np.max(np.abs(V))<1
Vm=float(np.max(np.abs(V)))
Vgram=V.T@V
ventry=float(gamma(2*rank+10))*rank*Vm*Vm+1e-15
odefect=upperfloat(np.max(np.sum(abs(Vgram-np.eye(rank)),axis=1))
                 /(1-float(gamma(2*rank+10)))+rank*ventry)
residual=Hcenter@V-V*values
rentry=float(gamma(2*rank+10))*rank*np.max(abs(Hcenter))*Vm+1e-14
rbound=upperfloat(float(np.linalg.norm(residual,'fro'))
                 /(1-float(gamma(2*rank*rank+10)))+rank*rentry)
assert odefect < 1
L=float(np.max(abs(values)))
# Polar-unitary residual bound, evaluated with Arb from exact binary data.
oarb=exactfloat(odefect);rarb=exactfloat(rbound);Larb=exactfloat(L)
spectral_error=exactfloat(Herror)+rarb/(1-oarb).sqrt()+2*Larb*(1+oarb).sqrt()*((1-oarb).rsqrt()-1)
Hlower=exactfloat(values[0])-spectral_error
Hupper=exactfloat(values[-1])+spectral_error
print('UNVERIFIED eigenvalue locators',values[0],values[-1],flush=True)
print('PASS actual certified H lower',Hlower,'H upper',Hupper,flush=True)
assert Hlower > 0
# A simple rational lower endpoint is selected downward, then independently
# compared using Arb; its relation to W is explicit rather than W=I.
scaleh=10**9
hnum=int((Hlower*scaleh/arbf(Wupper)).floor().unique_fmpz())
h=Fraction(hnum,scaleh)
assert h>0 and arbf(h)*arbf(Wupper) < Hlower
kappa=(1/arbf(h)).sqrt()
kappanum=int((kappa*10**6).upper().ceil().unique_fmpz())
kappabound=Fraction(kappanum,10**6)
assert kappa < arbf(kappabound)
assert exactfloat(values[0])+spectral_error > 0

outputfactor=root/'certificates/complement_physical_frame_locator.npz'
np.savez_compressed(outputfactor,V=V,values=values,Hcenter=Hcenter,Wcenter=Wcenter)
record={
 'status':'PASS actual physical rank56 and coefficient-metric conditioning',
 'source':'routes/complement-physical-frame-certificate.md',
 'dimension':N,'columns':rank,'kernel':'K=span(q)','physical_rank':rank,
 'coefficient_factor_interpretation':'all stored Z entries exact binary64 dyadic; no assumed orthogonality',
 'physical_Gram':'H=Z*D(G0-bb*/Q)DZ; tau is actual pole diagonal',
 'precision_bits':ctx.prec,'norm_integral_interval':'[0,2] doubled',
 'norm_integral_runtime':'existing Arb acb.integral, meromorphic callback',
 'norm_theta_error_exact':'1/2^250','norm_physical_tail_error_exact':'1/2^180',
 'norm_endpoint_scale':str(2**100),'Q_outward_endpoint_pair':endpointpair(Q),
 'tau_outward_endpoint_pair':endpointpair(tau),
 'Q_ball':str(Q),'tau_ball':str(tau),
 'coefficient_Gram_guarded_defect_upper':Wdefect,
 'coefficient_Gram_lower_exact':str(Wlower),'coefficient_Gram_upper_exact':str(Wupper),
 'physical_full_Gram_interval_error_upper':Gerror,
 'physical_compression_error_upper':Herror,
 'locator_orthogonality_error_upper':odefect,
 'locator_entire_residual_upper':rbound,
 'physical_spectral_error_ball':str(spectral_error),
 'physical_Gram_lower_ball':str(Hlower),'physical_Gram_upper_ball':str(Hupper),
 'generalized_lower_h_exact':str(h),'generalized_inequality':'H >= h W',
 'kappa_B_upper_exact':str(kappabound),'kappa_B_upper_ball':str(kappa),
 'all_coupling_retained':'entire dense frequency-frequency, frequency-pole and score-capture terms',
 'Fourier_far_tail':'5 exp(-|nu|/4); score 5|nu| exp(-|nu|/4), full theta',
 'pole_measure':2,'entry_IEEE_allowance':'1e-12',
 'source_hashes':{key:sha(path) for key,path in paths.items()},
 'Fourier_script_sha256':sha(root/'computations/complement_fourier_weight_grid.py'),
 'script_sha256':sha(Path(__file__)),
 'physical_locator_file_sha256':sha(outputfactor),
 'elapsed_seconds':time.monotonic()-start,
}
output=root/'certificates/complement_physical_frame_certificate.json'
output.write_text(json.dumps(record,indent=2)+'\n')
print('PASS rank',rank,'H >=',str(h),'W; kappa <=',str(kappabound),flush=True)
