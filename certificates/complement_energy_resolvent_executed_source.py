"""Shifted E8 resolvent certificate, using the entire physical residual.

Reuses the completed norm, Fourier, low-zero and exact-frame certificates.
Only the new zero-to-feature Gram and shifted trial question are evaluated.
"""
from pathlib import Path
from fractions import Fraction
import hashlib
import json
import math
import time
import numpy as np
from flint import arb, ctx

ctx.prec = 192
root = Path(__file__).resolve().parents[1]
started = time.monotonic()
unit = 2.0**-53
def up(x):
    return float(np.nextafter(float(x), np.inf))
def gamma(n):
    return up(n*unit/(1-n*unit))
def fro(x):
    return up((np.linalg.norm(x, 'fro')/(1-gamma(2*x.size+10))+1e-15)*(1+1e-12))
def product_error(x, y):
    g = gamma(2*x.shape[1]+10)
    return fro(g*(np.abs(x)@np.abs(y))/(1-g))
def sha(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()
def arbf(x):
    f = Fraction(x)
    return arb(f.numerator)/arb(f.denominator)
def boundnorm(a):
    """Certified spectral norm of a symmetric binary matrix via full locator."""
    symround = gamma(3)*fro(a)
    a = (a+a.T)/2
    values, v = np.linalg.eigh(a)
    o = fro(v.T@v-np.eye(len(a)))+product_error(v.T, v)
    r = fro(a@v-v*values)+product_error(a, v)
    r += gamma(3)*(fro(a@v)+2*fro(v*values))+1e-12
    assert o < 1e-6
    ob,rb,lb = arb(up(o+1e-12)),arb(up(r+1e-12)),arb(float(max(abs(values))))
    error = (rb/(1-ob).sqrt()+2*lb*(1+ob).sqrt()*((1-ob).rsqrt()-1)
             +arb(up(symround)))
    return up((lb+error).upper()),up((arb(float(values[-1]))+error).upper())

paths = {name:root/'certificates'/file for name,file in {
    'cells':'complement_signed_symbol_norm_cells.json',
    'fourier':'complement_fourier_weight_grid.json',
    'factor':'complement_midpoint_rank_factors.npz',
    'physical':'complement_physical_frame_certificate.json',
    'physical_locator':'complement_physical_frame_locator.npz',
    'norms':'complement_physical_norms.json',
    'sampling':'complement_low_zero_frame_certificate.json',
    'sampling_matrices':'complement_low_zero_frame_matrices.npz',
    'taylor':'complement_low_zero_taylor_grid.npz',
    'taylor_record':'complement_low_zero_taylor_grid_certificate.json',
    'mass_deficit':'complement_midpoint_mass_deficit.json',
    'deficit_norm':'complement_midpoint_deficit_norm.json',
}.items()}
cachepath = root/'certificates/complement_energy_resolvent_matrices.npz'
cacherecordpath = root/'certificates/complement_energy_resolvent.json'
cache_record = json.loads(cacherecordpath.read_text()) if cacherecordpath.exists() else None
cached = (cache_record is not None and cachepath.exists()
          and sha(cachepath)==cache_record['matrix_sha256']
          and all(sha(paths[name])==cache_record['input_sha256'][name]
                  for name in ['cells','fourier','factor','norms','taylor']))
cache = np.load(cachepath) if cached else None
physical = json.loads(paths['physical'].read_text())
sampling = json.loads(paths['sampling'].read_text())
assert sha(paths['factor']) == physical['source_hashes']['factors']
assert sha(paths['cells']) == physical['source_hashes']['cells']
assert sha(paths['fourier']) == physical['source_hashes']['fourier']
assert sha(paths['physical_locator']) == physical['physical_locator_file_sha256']
assert sha(paths['sampling_matrices']) == sampling['matrix_file_sha256']
assert sha(paths['taylor']) == sampling['Taylor_grid_sha256']
assert sha(paths['factor']) == sampling['factor_file_sha256']
assert sha(paths['cells']) == sampling['cell_certificate_sha256']
Z = np.load(paths['factor'])['U']
loc = np.load(paths['physical_locator'])
low = np.load(paths['sampling_matrices'])
H, N = loc['Hcenter'], loc['Wcenter']
sample, Gzero = low['A_center'], low['Gzero_center']
zeros, zero_error = low['zero_center'], low['zero_error']
assert Z.shape == (2906,56) and sample.shape == (2469,56)
cells = json.loads(paths['cells'].read_text())['kernel_cells']
k = np.array([c['k'] for c in cells], dtype=np.int64)
mu = np.array([np.nextafter(int(c['measure_upper_dyadic_numerator'])/2**80,np.inf)
               for c in cells]+[2.0])
sm = np.sqrt(mu)
smerror = np.array([up(abs(arb(float(m)).sqrt()-arb(float(s))).upper())
                    for m,s in zip(mu,sm)])
times = (2*k+1)/128
f = json.loads(paths['fourier'].read_text())
norms = json.loads(paths['norms'].read_text())
Q, tau = arb(norms['Q_ball']), arb(norms['tau_ball'])
qlo,qhi = np.nextafter(float(Q.lower()),-np.inf),up(Q.upper())
tlo,thi = np.nextafter(float(tau.lower()),-np.inf),up(tau.upper())
scale = int(f['endpoint_scale'])
def pairs(key):
    return (np.array([np.nextafter(int(a)/scale,-np.inf) for a,b in f[key]]),
            np.array([up(int(b)/scale) for a,b in f[key]]))
wlo,whi = pairs('W_outward_endpoint_pairs')
vlo,vhi = pairs('pole_outward_endpoint_pairs')
blo,bhi = pairs('score_overlap_outward_endpoint_pairs')
far = np.zeros(2*int(k.max())+2)
for j in range(6401,6401 if cached else len(far)):
    far[j] = up((5*(-arb(j)/256).exp()).upper())
bl,bh,vl,vh = [np.zeros(len(k)) for _ in range(4)]
for j,kj in enumerate([] if cached else k):
    if kj < 6400:
        bl[j],bh[j],vl[j],vh[j] = blo[kj],bhi[kj],vlo[kj],vhi[kj]
    else:
        t = arb(2*int(kj)+1)/128
        eb,ev = up((5*t*(-t/4).exp()).upper()),up((5*(-t/4).exp()).upper())
        bl[j],bh[j],vl[j],vh[j] = -eb,eb,-ev,ev
GB, GBerror = np.empty((2906,2906)),np.empty((2906,2906))
for j,kj in enumerate([] if cached else k):
    def wb(index):
        lo,hi = -far[index],far[index]
        near = index<=6400
        lo[near],hi[near] = wlo[index[near]],whi[index[near]]
        return lo,hi
    dl,dh = wb(abs(kj-k)); sl,sh = wb(kj+k+1)
    prod = np.stack([bl[j]*bl,bl[j]*bh,bh[j]*bl,bh[j]*bh])
    pmin,pmax = prod.min(axis=0),prod.max(axis=0)
    lo = ((dl-sh)/2-np.maximum(pmax/qlo,pmax/qhi))*sm[j]*sm[:-1]
    hi = ((dh-sl)/2-np.minimum(pmin/qlo,pmin/qhi))*sm[j]*sm[:-1]
    GB[j,:-1],GBerror[j,:-1] = (lo+hi)/2,(hi-lo)/2+1e-12
lo = (vl+np.minimum(bl/(4*qlo),bl/(4*qhi)))*sm[:-1]*sm[-1]
hi = (vh+np.maximum(bh/(4*qlo),bh/(4*qhi)))*sm[:-1]*sm[-1]
GB[:-1,-1] = GB[-1,:-1] = (lo+hi)/2
GBerror[:-1,-1] = GBerror[-1,:-1] = (hi-lo)/2+1e-12
lo,hi = 2*(tlo-1/(16*qlo)),2*(thi-1/(16*qhi))
GB[-1,-1],GBerror[-1,-1] = (lo+hi)/2,(hi-lo)/2+1e-12
if cached:
    GB = cache['GB']
    gerror = cache_record['B_Gram_error_operator_upper']
    gnorm = cache_record['B_Gram_operator_upper']
else:
    GB = (GB+GB.T)/2
    gerror = up(GBerror.sum(axis=1).max()/(1-gamma(2*len(GB)+10)))
    gnorm = up(np.abs(GB).sum(axis=1).max()/(1-gamma(2*len(GB)+10)))
assert gerror < 3.3e-9

grid = np.load(paths['taylor'])
wc,we,vc,ve = [grid[key] for key in ['W_center','W_error','pole_center','pole_error']]
far2 = np.zeros(math.ceil((3000+times.max()+1)*8)+2)
for j in range(1022,1022 if cached else len(far2)):
    far2[j] = up((5*(-arb(j)/32).exp()).upper())
def fourier(nu, ae, pole=False):
    absolute = abs(nu)
    index = np.minimum(np.rint(absolute*8).astype(np.int64),1024)
    h = absolute-index/8
    he = gamma(3)*(absolute+index/8)
    cg,rg = (vc,ve) if pole else (wc,we)
    val,mag,rad = [np.zeros_like(absolute) for _ in range(3)]
    for d in range(19,-1,-1):
        val = val*h+cg[index,d]
        mag = mag*abs(h)+abs(cg[index,d])
        rad = rad*abs(h)+rg[index,d]
    err = (rad*(1+gamma(50))+120*(2*(abs(h)+he))**20
           +gamma(50)*mag+(.012 if pole else .11)*(ae+he))*(1+1e-12)+1e-300
    fi = np.maximum(np.floor(absolute*8).astype(np.int64)-1,1022)
    val,err = np.where(absolute<=128,val,0),np.where(absolute<=128,err,far2[fi])
    return val*np.sign(nu) if pole else val,np.nextafter(err,np.inf)

# Qb contains every zero-to-feature correlation before compression by Z.
Qb,Qberror = np.empty((2469,2906)),np.empty((2469,2906))
for j,(z,ze) in enumerate([] if cached else zip(zeros,zero_error)):
    ae = ze+gamma(3)*(abs(z)+times)
    dm,de = fourier(z-times,ae); sp,se = fourier(z+times,ae)
    vp,vep = fourier(np.array([z]),np.array([ze]),pole=True)
    un = np.append((dm-sp)/2,vp[0]); er = np.append((de+se)/2+gamma(3)*(abs(dm)+abs(sp)),vep[0])
    Qb[j] = sm*un
    Qberror[j] = (sm*er+smerror*(abs(un)+er)+gamma(3)*abs(Qb[j]))*(1+1e-12)
if cached:
    Qb = cache['Qb']
    qerror = cache_record['zero_B_error_Frobenius_upper']
    print('Reuse completed full B/zero Gram; only changed finite guards evaluated',flush=True)
else:
    qerror = fro(Qberror)
gzerror = float(sampling['Gzero_error_operator_upper'])
gznorm = float(sampling['Gzero_center_operator_upper'])
sample_error = float(sampling['A_evaluation_rounding_Frobenius_error'])
herror = float(physical['physical_compression_error_upper'])
GBZ = cache['GBZ'] if cached else GB@Z
GBZerror = gerror*math.sqrt(1+1e-8)+product_error(GB,Z)
mid = cache['mid'] if cached else GBZ.T@GBZ
miderror = 2*fro(GBZ)*GBZerror+GBZerror**2+product_error(GBZ.T,GBZ)
L = low['L_center']; lerror = float(sampling['L_error_operator_upper'])
print('PASS full B/zero physical correlations; errors',gerror,qerror,flush=True)

# A genuine operator error between the actual continuous Rq and Rmid.
mass = json.loads(paths['mass_deficit'].read_text())
assert mass['cell_certificate_sha256'] == sha(paths['cells'])
assert [c['k'] for c in mass['cells']] == list(k)
deficit = Fraction(mass['deficit_theorem_upper_exact'])
delta0 = Fraction(21,1024000)
# The full weighted Gram gives the entire missing-mass operator price.
# Its row bound retains every frequency correlation; the pole deficit is zero.
deficitnorm = json.loads(paths['deficit_norm'].read_text())
dnorm_exact = Fraction(deficitnorm['norm_upper_exact'])
dnorm = up(arbf(dnorm_exact))
# Accepted dense rank residual: GB<=Gamma, ||Gamma-Z Lambda Z*||<=3/12500.
largest = float(np.load(paths['factor'])['values'].max())
bgram_upper = (arbf(Fraction(100000001,100000000))*arb(largest)
               +arbf(Fraction(3,12500)))
# ||Rq-Rtrue_mid|| <= 2 sqrt(||Rtrue_mid|| delta0)+delta0.
# Rtrue_mid<=Rmid, whose norm is <=||GB||. Mass deficit costs ||w||^2<=1/2.
eta_ball = 2*(bgram_upper*arbf(delta0)).sqrt()+arbf(delta0)+arbf(dnorm_exact)
eta_ball += arbf(Fraction(32237,10**12)+Fraction(14,10**19))
eta = up(eta_ball.upper())
records = []
trial_matrices = {}
for mu_exact in [Fraction(1,20),Fraction(1,10),Fraction(3,20),Fraction(1,5),Fraction(1,4)]:
    delta_exact = mu_exact-Fraction(979236349,2000000000000)
    delta = float(delta_exact)
    # Binary X is a trial definition; numerical solve accuracy is no premise.
    T = L+mid+delta*H
    X = np.linalg.solve(T,H)
    assert np.isfinite(X).all()
    Xgram = X.T@X
    Xnorm = up(math.sqrt(boundnorm(Xgram)[1]+product_error(X.T,X)))
    Ygram = X.T@H@X
    Yerror = herror*Xnorm**2+product_error(X.T,H@X)+Xnorm*product_error(H,X)
    Ynorm = math.sqrt(max(0,boundnorm(Ygram)[1]+Yerror))
    # Entire physical residual, including every B/zero mixed pairing.
    DX = np.eye(56)-delta*X
    D1 = Z@DX-GBZ@X
    D2 = -2*sample@X
    e1 = (product_error(Z,DX)+product_error(GBZ,X)+GBZerror*Xnorm
          +gamma(3)*(fro(Z@DX)+fro(GBZ@X))+gamma(3)*delta*fro(X)*math.sqrt(1+1e-8))
    e2 = 2*(sample_error*Xnorm+product_error(sample,X))+gamma(3)*fro(D2)
    residue = D1.T@(GB@D1)+2*D1.T@(Qb.T@D2)+D2.T@(Gzero@D2)
    d1,d2 = fro(D1),fro(D2)
    re = gerror*d1*d1+2*qerror*d1*d2+gzerror*d2*d2
    re += (product_error(D1.T,GB@D1)+d1*product_error(GB,D1)
           +2*product_error(D1.T,Qb.T@D2)+2*d1*product_error(Qb.T,D2)
           +product_error(D2.T,Gzero@D2)+d2*product_error(Gzero,D2))
    re += gamma(5)*(fro(residue)+d1*d1*gnorm+2*d1*d2*fro(Qb)+d2*d2*gznorm)
    rnorm = (math.sqrt(max(0,boundnorm(residue)[1]+re))
             +math.sqrt(gnorm+gerror)*e1+math.sqrt(gznorm+gzerror)*e2)
    # The exact actual K differs from this finite quadratic by eta*Y*Y.
    K = H@X+X.T@H-X.T@T@X
    ke = (2*herror*Xnorm+(lerror+miderror+delta*herror)*Xnorm**2
          +2*product_error(H,X)+product_error(X.T,T@X)+Xnorm*product_error(T,X)
          +gamma(5)*(2*fro(H@X)+fro(X.T@T@X)))
    kupper = boundnorm(K+eta*Ygram)[1]+ke+eta*Yerror
    residual_upper = rnorm+eta*Ynorm
    priceball = arbf(Fraction(11,10))*(arb(up(kupper+1e-6))
                  +arb(up(residual_upper+1e-6))**2/arbf(delta_exact))
    price = up(priceball.upper())
    # Includes coefficient target N >= (1-1e-8)I and scalar rounding reserve.
    record = {'mu_exact':str(mu_exact),'delta_exact':str(delta_exact),'eta_upper':eta,
              'trial_X_operator_upper':Xnorm,'trial_Y_physical_operator_upper':Ynorm,
              'midpoint_physical_residual_upper':rnorm,'actual_physical_residual_upper':residual_upper,
              'captured_K_upper':kupper,'resolvent_price_upper':price,
              'passes':price<1-1e-8,'finite_K_error_upper':ke,'finite_residual_Gram_error_upper':re}
    records.append(record)
    key = str(mu_exact).replace('/','_')
    trial_matrices['X_'+key] = X
    trial_matrices['Ygram_'+key] = Ygram
    trial_matrices['K_'+key] = K
    trial_matrices['residual_Gram_'+key] = residue
    print('Shift',mu_exact,'full resolvent price',price,'PASS' if record['passes'] else 'FAIL',flush=True)
assert any(r['passes'] and Fraction(r['mu_exact'])<Fraction(258019032237,10**12) for r in records)
chosen = next(r for r in records if r['mu_exact']=='1/10')
assert arb(chosen['captured_K_upper'])+arb('0.000001') < arbf(Fraction(101,200))
assert arb(chosen['actual_physical_residual_upper'])+arb('0.000001') < arbf(Fraction(77,500))
theta = Fraction(41,50)
final_price_exact = Fraction(11,10)*(Fraction(101,200)
                  +Fraction(77,500)**2/(Fraction(1,10)-Fraction(979236349,2000000000000)))
final_price = arbf(final_price_exact)
assert final_price < arbf(theta*Fraction(99999999,100000000))
full_loss = Fraction(979236349,2000000000000)+theta*(Fraction(1,10)-Fraction(979236349,2000000000000))
assert full_loss < Fraction(8209,100000)
factorpath = root/'certificates/complement_energy_resolvent_matrices.npz'
np.savez_compressed(factorpath,GB=GB,Qb=Qb,GBZ=GBZ,mid=mid,**trial_matrices)
record = {'status':'PASS at least one improved full odd quotient floor',
          'source':'routes/complement-energy-resolvent-advance.md','shifts':records,
          'continuous_midpoint_error_ball':str(eta_ball),
          'mass_deficit_upper_exact':str(deficit),'B_Gram_operator_upper':gnorm,
          'B_Gram_spectral_upper_ball':str(bgram_upper),
          'weighted_deficit_operator_upper':dnorm,
          'weighted_deficit_operator_upper_exact':str(dnorm_exact),
          'chosen_mu_exact':'1/10','chosen_K_upper_exact':'101/200',
          'chosen_residual_upper_exact':'77/500','chosen_resolvent_upper_exact':str(theta),
          'chosen_price_exact':str(final_price_exact),
          'retained_energy_coefficient_exact':'9/50','full_floor_exact':str(full_loss),
          'chosen_price_Arb_ball':str(final_price),
          'scalar_rounding_reserve_exact':'1/1000000',
          'B_Gram_error_operator_upper':gerror,'zero_B_error_Frobenius_upper':qerror,
          'coefficient_target_lower_exact':'99999999/100000000',
          'high_zero_operator_error_exact':'7/5000000000000000000',
          'prime_tail_exact':'32237/1000000000000',
          'all_couplings':'all2906Bfeatures/all2469zeros, physical mixed blocks, actualQ/tau/score/pole',
          'input_sha256':{name:sha(path) for name,path in paths.items()},
          'full_Gram_cache_reused':cached,
          'matrix_sha256':sha(factorpath),'script_sha256':sha(Path(__file__)),
          'elapsed_seconds':time.monotonic()-started}
(root/'certificates/complement_energy_resolvent.json').write_text(json.dumps(record,indent=2)+'\n')
