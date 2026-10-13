"""Outward whole-action certificate for an exact binary enriched trial.

RE2--RE5 own the analytical bridge. Numerical solves select trial coefficients
only. Existing completed theta/zero/frame certificates are reused by hash.
The convolution uses native acb.dft on exact binary inputs, not scipy FFT.
"""
from pathlib import Path
from fractions import Fraction
import ast
import hashlib
import json
import math
import time
import numpy as np
from flint import arb, acb, arb_mat, ctx

ctx.prec = 128
start = time.monotonic()
root = Path(__file__).resolve().parents[1]
p = root / 'certificates'
sha = lambda path: hashlib.sha256(path.read_bytes()).hexdigest()
owner = root / 'computations/complement_energy_resolvent.py'
main = json.loads((p/'complement_energy_resolvent.json').read_text())
assert sha(owner) == main['script_sha256']
names = {'up', 'gamma', 'fro', 'product_error', 'arbf', 'boundnorm'}
definitions = [n for n in ast.parse(owner.read_text()).body
               if isinstance(n, ast.FunctionDef) and n.name in names]
assert {n.name for n in definitions} == names
unit = 2.0**-53
exec(compile(ast.Module(body=definitions, type_ignores=[]), str(owner), 'exec'))
paths = {n: p/n for n in [
    'complement_residual_enrichment_pilot.json',
    'complement_residual_enrichment_pilot_matrices.npz',
    'complement_energy_resolvent.json', 'complement_energy_resolvent_matrices.npz',
    'complement_energy_resolvent_tighter_matrices.npz',
    'complement_midpoint_rank_factors.npz',
    'complement_physical_frame_certificate.json', 'complement_physical_frame_locator.npz',
    'complement_physical_norms.json', 'complement_signed_symbol_norm_cells.json',
    'complement_low_zero_frame_certificate.json', 'complement_low_zero_frame_matrices.npz',
    'complement_low_zero_taylor_grid_certificate.json', 'complement_low_zero_taylor_grid.npz',
    'complement_fourier_weight_grid.json',
]}
pilot = json.loads(paths['complement_residual_enrichment_pilot.json'].read_text())
for n, hsh in pilot['input_sha256'].items():
    assert sha(paths[n]) == hsh
for n, hsh in pilot['additional_input_sha256_postprocessing'].items():
    assert sha(paths[n]) == hsh
assert sha(paths['complement_residual_enrichment_pilot_matrices.npz']) == pilot['matrix_sha256']
physical = json.loads(paths['complement_physical_frame_certificate.json'].read_text())
sampling = json.loads(paths['complement_low_zero_frame_certificate.json'].read_text())
assert sha(paths['complement_low_zero_frame_matrices.npz']) == sampling['matrix_file_sha256']
assert sha(paths['complement_low_zero_taylor_grid.npz']) == sampling['Taylor_grid_sha256']
assert sha(paths['complement_physical_frame_locator.npz']) == physical['physical_locator_file_sha256']
assert sha(paths['complement_fourier_weight_grid.json']) == physical['source_hashes']['fourier']
assert sha(paths['complement_midpoint_rank_factors.npz']) == physical['source_hashes']['factors']
saved = np.load(paths['complement_residual_enrichment_pilot_matrices.npz'])
Z = np.load(paths['complement_midpoint_rank_factors.npz'])['U']
# These rounded products themselves define the exact dyadic trial Y.
# No assertion about their equality to a real matrix product is needed.
b = saved['Bcoef'] @ saved['X']
v = saved['Vcoef'] @ saved['X']
assert b.shape == (2906, 56) and v.shape == (2469, 56)
cache = np.load(paths['complement_energy_resolvent_matrices.npz'])
GB, Qb = cache['GB'], cache['Qb']
low = np.load(paths['complement_low_zero_frame_matrices.npz'])
Gz, zeros, ze = low['Gzero_center'], low['zero_center'], low['zero_error']
assert np.all(ze < 1e-9) and zeros[-1]+ze[-1] < 3000
loc = np.load(paths['complement_physical_frame_locator.npz'])
H, N = loc['Hcenter'], loc['Wcenter']
norms = json.loads(paths['complement_physical_norms.json'].read_text())
Q = arb(norms['Q_ball'])
gerr = main['B_Gram_error_operator_upper']
qerr = main['zero_B_error_Frobenius_upper']
gzerr = sampling['Gzero_error_operator_upper']
bop = Fraction('0.255985513')
herror = physical['physical_compression_error_upper']

def reserve(x):
    # Accepted float guard has a relative scalar reserve, as in its owner.
    return arb(float(x))*(1+arbf(Fraction(1,10**10))) + arbf(Fraction(1,10**25))

def ballfloat(x):
    c = float(x.mid())
    return c, up(abs(x-arb(c)).upper())

def exact(a):
    return arb_mat([[arb(float(x)) for x in row] for row in a])

def afrob(a):
    return sum((a[i,j]**2 for i in range(a.nrows())
                for j in range(a.ncols())), arb(0)).sqrt()

def spectral_upper(a, key):
    # Eigensolver is a locator; complete native residual/polar correction pays it.
    values, V = np.linalg.eigh((a+a.T)/2)
    av = exact(V)
    ident = arb_mat(np.eye(56).astype(int).tolist())
    diag = arb_mat([[arb(float(values[i])) if i == j else arb(0)
                     for j in range(56)] for i in range(56)])
    aa = (exact(a)+exact(a).transpose())/2
    o = afrob(av.transpose()*av-ident)
    r = afrob(aa*av-av*diag)
    assert o < arbf(Fraction(1,10**6))
    l = arb(float(np.max(abs(values))))
    correction = r/(1-o).sqrt()+2*l*(1+o).sqrt()*((1-o).rsqrt()-1)
    record['locators'][key] = {'orthogonality':str(o), 'whole_residual':str(r),
                              'polar_correction':str(correction)}
    return arb(float(values[-1]))+correction

record = {'status':'INCOMPLETE', 'precision_bits':ctx.prec, 'locators':{},
          'trial':'exact binary b and v computed from saved pilot Bcoef,Vcoef,X; no rank or solve premise',
          'input_sha256':{n:sha(path) for n,path in paths.items()},
          'source_sha256':sha(Path(__file__)),
          'route_sha256':sha(root/'routes/complement-residual-enrichment.md')}
cells = json.loads(paths['complement_signed_symbol_norm_cells.json'].read_text())['kernel_cells']
times = np.array([(2*c['k']+1)/128 for c in cells])
mu = np.array([np.nextafter(int(c['measure_upper_dyadic_numerator'])/2**80,np.inf)
               for c in cells]+[2.0])
sm = np.sqrt(mu)
sme = np.array([up(abs(arb(float(m)).sqrt()-arb(float(s))).upper())
                for m,s in zip(mu,sm)])
fourier = json.loads(paths['complement_fourier_weight_grid.json'].read_text())
scale = int(fourier['endpoint_scale'])
bb, bbe = np.zeros(2906), np.zeros(2906)
for j, t in enumerate(times):
    k = cells[j]['k']
    if k < 6400:
        lo,hi = fourier['score_overlap_outward_endpoint_pairs'][k]
        bb[j],bbe[j] = ballfloat((arb(lo)+arb(hi))/(2*scale))
        bbe[j] += up((arb(hi)-arb(lo))/(2*scale))
    else:
        bbe[j] = up((5*arb(float(t))*(-arb(float(t))/4).exp()).upper())
bb[-1] = -.25

# Exact finite Rmid minus the retained pole: B D, with its input error paid.
gbb, qbv = GB@b, Qb.T@v
d = gbb+qbv
d[-1] = 0
d_error = (gerr*fro(b)+qerr*fro(v)+product_error(GB,b)
           +product_error(Qb.T,v)+gamma(3)*(fro(gbb)+fro(qbv)))
finite_map_error = reserve(d_error)*arbf(bop).sqrt()
blocks = [Z,b,d]
vb = [np.zeros((2469,56)),v,np.zeros((2469,56))]
PG, PE = {}, {}
for i in range(3):
    for j in range(i,3):
        a1,a2,a3 = GB@blocks[j],Qb.T@vb[j],Gz@vb[j]
        p1,p2,p3,p4 = blocks[i].T@a1,blocks[i].T@a2,vb[i].T@(Qb@blocks[j]),vb[i].T@a3
        center = p1+p2+p3+p4
        err = gerr*fro(blocks[i])*fro(blocks[j])
        err += qerr*(fro(blocks[i])*fro(vb[j])+fro(vb[i])*fro(blocks[j]))
        err += gzerr*fro(vb[i])*fro(vb[j])
        err += product_error(blocks[i].T,a1)+fro(blocks[i])*product_error(GB,blocks[j])
        err += product_error(blocks[i].T,a2)+fro(blocks[i])*product_error(Qb.T,vb[j])
        err += product_error(vb[i].T,Qb@blocks[j])+fro(vb[i])*product_error(Qb,blocks[j])
        err += product_error(vb[i].T,a3)+fro(vb[i])*product_error(Gz,vb[j])
        err += gamma(5)*(fro(p1)+fro(p2)+fro(p3)+fro(p4))
        PG[i,j], PE[i,j] = center, reserve(err)
PG[0,0],PE[0,0] = H,reserve(herror)
print('Closed exact trial/finite physical Grams',time.monotonic()-start,flush=True)

Bcoeff = np.concatenate(blocks,axis=1)
Vcoeff = np.concatenate(vb,axis=1)
wb = sm[:,None]*Bcoeff
wbe = sme[:,None]*abs(Bcoeff)+gamma(3)*abs(wb)
overlap = sm*bb
overlaperr = sm*bbe+sme*(abs(bb)+bbe)+gamma(3)*abs(overlap)
score_center = -(overlap@Bcoeff)/float(Q.mid())
score_err = ((overlaperr@abs(Bcoeff))/(1-gamma(6000))
             +gamma(6000)*(abs(overlap)@abs(Bcoeff))/(1-gamma(6000)))/float(Q.lower())
score_err += abs(score_center)*up(abs(Q-arb(float(Q.mid()))).upper())/float(Q.lower())
score_err += gamma(5)*abs(score_center)
score_err = np.nextafter(score_err*(1+1e-10)+1e-30,np.inf)

den,T,R = 20,3300,3
K = den*T
t = np.arange(K+1,dtype=float)/den
te = gamma(3)*t
h = arb(1)/den
pi = arb.pi()
grid = np.load(paths['complement_low_zero_taylor_grid.npz'])
wc,we,pc,pe = [grid[n] for n in ['W_center','W_error','pole_center','pole_error']]
farsize = math.ceil((T+3000+1)*8)+4
far = np.array([up((5*(-arb(j)/32).exp()).upper()) for j in range(farsize)])

def transform(nu, ae, pole=False):
    absolute = abs(nu)
    ix = np.minimum(np.rint(absolute*8).astype(np.int64),1024)
    z = absolute-ix/8
    he = gamma(3)*(absolute+ix/8)
    cg,eg = (pc,pe) if pole else (wc,we)
    value,mag,err = np.zeros_like(z),np.zeros_like(z),np.zeros_like(z)
    for k in range(19,-1,-1):
        value = value*z+cg[ix,k]
        mag = mag*abs(z)+abs(cg[ix,k])
        err = err*abs(z)+eg[ix,k]
    err = (err*(1+gamma(50))+120*(2*(abs(z)+he))**20
           +gamma(50)*mag+(.012 if pole else .11)*(ae+he))*(1+1e-10)+1e-300
    return value*np.sign(nu) if pole else value,np.nextafter(err,np.inf)

# New score samples only, using the existing proven whole-theta GL bound.
score = np.zeros(K+1)
scoree = np.zeros(K+1)
rho = Fraction(23,20)
egl = 24*240*Fraction(3)**32*rho**(-640)/(1-1/rho)
baseerror = arbf(egl+32*Fraction(1,2)**128+Fraction(1,2)**1100)
nodes=[]
for j in range(320):
    x,qweight = arb.legendre_p_root(320,j,weight=True)
    x = (x+1)*3/2
    qweight *= 3/2
    phi=arb(0)
    for n in range(1,9):
        y=pi*n*n*(2*x).exp()
        phi += (x/2-y).exp()*(4*y*y-6*y)
    nodes.append((x,2*qweight*phi))
for k in range(128*den+1):
    arg=arb(k)/den
    value=sum((qweight*(arg*x).cos() for x,qweight in nodes),arb(0))
    score[k],scoree[k]=ballfloat(-arg*(value+arb(0,baseerror)))
for k in range(128*den+1,K+1):
    scoree[k]=up((5*(arb(k)/den)*(-arb(k)/(4*den)).exp()).upper())
print('Closed full Xi/score frequency samples',time.monotonic()-start,flush=True)

fh = np.zeros((K+1,168))
rowerror = np.zeros((K+1,3))
rowcoeff = [np.linalg.norm(wb[:-1,i*56:(i+1)*56],axis=1)*(1+1e-10)
            +np.linalg.norm(wbe[:-1,i*56:(i+1)*56],axis=1)*(1+1e-10) for i in range(3)]
vrowcoeff = [np.linalg.norm(Vcoeff[:,i*56:(i+1)*56],axis=1)*(1+1e-10) for i in range(3)]
for first in range(0,K+1,256):
    last=min(first+256,K+1)
    tt=t[first:last,None]
    tterr=te[first:last,None]
    # Every omitted far sample receives its complete exponential enclosure.
    bc=np.where((times>=tt[0,0]-128)&(times<=tt[-1,0]+128))[0]
    vc=np.where((zeros>=tt[0,0]-128)&(zeros<=tt[-1,0]+128))[0]
    bm,be1=transform(tt-times[bc],tterr+gamma(3)*(tt+times[bc]))
    bp,be2=transform(tt+times[bc],tterr+gamma(3)*(tt+times[bc]))
    bn=(bm-bp)/2
    bne=(be1+be2)/2+gamma(3)*(abs(bm)+abs(bp))
    # The positive argument may exceed 128; use its complete envelope there.
    afar=tt+times[bc]>128
    bp[afar]=0
    be2[afar]=far[np.maximum(np.floor((tt+times[bc])[afar]*8).astype(int)-1,1022)]
    bn=(bm-bp)/2
    bne=(be1+be2)/2+gamma(3)*(abs(bm)+abs(bp))
    am=abs(tt-times[bc]); farminus=am>128
    bm[farminus]=0
    be1[farminus]=far[np.maximum(np.floor(am[farminus]*8).astype(int)-1,1022)]
    bn=(bm-bp)/2
    bne=(be1+be2)/2+gamma(3)*(abs(bm)+abs(bp))
    vm,ve1=transform(tt-zeros[vc],tterr+ze[vc]+gamma(3)*(tt+zeros[vc]))
    vp,ve2=transform(tt+zeros[vc],tterr+ze[vc]+gamma(3)*(tt+zeros[vc]))
    for vals,ers,args in [(vm,ve1,abs(tt-zeros[vc])),(vp,ve2,tt+zeros[vc])]:
        beyond=args>128
        vals[beyond]=0
        ers[beyond]=far[np.maximum(np.floor(args[beyond]*8).astype(int)-1,1022)]
    vn=(vm-vp)/2
    vne=(ve1+ve2)/2+gamma(3)*(abs(vm)+abs(vp))
    pv,pve=transform(tt[:,0],te[first:last],pole=True)
    farpole=tt[:,0]>128
    pv[farpole]=0
    pve[farpole]=far[np.maximum(np.floor(tt[farpole,0]*8).astype(int)-1,1022)]
    pb = bn@wb[bc]
    pz = vn@Vcoeff[vc]
    pp = pv[:,None]*wb[-1]
    ps = score[first:last,None]*score_center
    fh[first:last]=pb+pz+pp+ps
    eb=(bne@abs(wb[bc])+abs(bn)@wbe[bc])/(1-gamma(2*len(bc)+10))
    eb+=gamma(2*len(bc)+10)*(abs(bn)@abs(wb[bc]))/(1-gamma(2*len(bc)+10))
    ev=vne@abs(Vcoeff[vc])/(1-gamma(2*len(vc)+10))
    ev+=gamma(2*len(vc)+10)*(abs(vn)@abs(Vcoeff[vc]))/(1-gamma(2*len(vc)+10))
    ep=pve[:,None]*(abs(wb[-1])+wbe[-1])+abs(pv[:,None])*wbe[-1]
    es=scoree[first:last,None]*(abs(score_center)+score_err)+abs(score[first:last,None])*score_err
    ef=(eb+ev+ep+es+gamma(5)*(abs(pb)+abs(pz)+abs(pp)+abs(ps)))*(1+1e-10)
    omitb=np.ones(len(times),dtype=bool);omitb[bc]=False
    omitv=np.ones(len(zeros),dtype=bool);omitv[vc]=False
    obarg=abs(tt-times[omitb])
    ovarg=abs(tt-zeros[omitv])
    ober=far[np.maximum(np.floor(obarg*8).astype(int)-1,1022)]
    ober+=far[np.maximum(np.floor((tt+times[omitb])*8).astype(int)-1,1022)]
    over=far[np.maximum(np.floor(ovarg*8).astype(int)-1,1022)]
    over+=far[np.maximum(np.floor((tt+zeros[omitv])*8).astype(int)-1,1022)]
    for i in range(3):
        rowerror[first:last,i]=(np.linalg.norm(ef[:,i*56:(i+1)*56],axis=1)
            +ober@rowcoeff[i][omitb]/2+over@vrowcoeff[i][omitv]/2)*(1+1e-10)+1e-30
    if first % 16384 == 0:
        print('Closed sample rows',first,'/',K,time.monotonic()-start,flush=True)
assert np.isfinite(fh).all() and np.isfinite(rowerror).all()

ell,elle=np.empty(K+1),np.empty(K+1)
logs=[arb(n).log() for n in [2,3,4,5,7,8]]
weights=[2*arb(n).log()/arb(m).sqrt() for n,m in zip([2,3,2,5,7,2],[2,3,4,5,7,8])]
for k in range(K+1):
    arg=arb(k)/den
    symbol=acb(arb(1)/4,arg/2).digamma().real-pi.log()
    symbol-=sum((w*(arg*l).cos() for w,l in zip(weights,logs)),arb(0))
    ell[k],elle[k]=ballfloat(symbol)
action=ell[:,None]*fh[:,56:112]
actionerror=(abs(ell)*rowerror[:,1]+elle*(np.linalg.norm(fh[:,56:112],axis=1)*(1+1e-10)+rowerror[:,1])
             +gamma(3)*np.linalg.norm(action,axis=1)*(1+1e-10))*(1+1e-10)
node=reserve(float(np.sum(actionerror))/den/math.pi)

# Complete analytic alias and omitted frequency/physical tails, native outward.
cn=np.array([up(np.linalg.norm(wb[j,56:112])*(1+1e-10)+np.linalg.norm(wbe[j,56:112])*(1+1e-10))
             for j in range(2906)])
vn=np.array([up(np.linalg.norm(v[j])*(1+1e-10)) for j in range(2469)])
an=reserve(np.linalg.norm(score_center[56:112])*(1+1e-10)+np.linalg.norm(score_err[56:112])*(1+1e-10))
A=20*(sum(((4480+160*arb(float(a)))*arb(float(c)) for a,c in zip(times,cn[:-1])),arb(0))
      +sum(((4480+160*(arb(float(a))+arb(float(e))))*arb(float(c)) for a,e,c in zip(zeros,ze,vn)),arb(0))
      +4480*arb(float(cn[-1]))+47200*an)
P=2*pi/h
alias=A/pi*(arb(R)/4).exp()*(-P/4).exp()/(1-(-P/4).exp())
rr=(-h/4).exp()
first=arb(T)+h
sum0=1/(1-rr)
sum1=first*sum0+h*rr/(1-rr)**2
sum2=first**2*sum0+2*first*h*rr/(1-rr)**2+h**2*rr*(1+rr)/(1-rr)**3
weighted=arb(0)
for a,c in list(zip(times,cn[:-1]))+list(zip(zeros+ze,vn)):
    weighted += arb(float(c))*((arb(float(a))/4).exp()+(-arb(float(a))/4).exp())/2
weighted += arb(float(cn[-1]))
tail=h/pi*(-first/4).exp()*5*(weighted*(200*sum0+10*sum1)+an*(200*sum1+10*sum2))
L1nodes=reserve(float(np.linalg.norm(action,axis=1).sum())/den/math.pi)
B_h=L1nodes+node+tail
spatial=(A/(2*pi)+B_h)*arbf(Fraction(1,2)**550)
egamma=(alias+tail+node)/arb(2).sqrt()+spatial+finite_map_error
record['action_prices']={name:str(val) for name,val in {
    'strip_L1_A':A,'alias_before_w':alias,'frequency_tail_before_w':tail,
    'node_before_w':node,'physical_tail':spatial,'finite_Rmid_map_error':finite_map_error,
    'whole_map_error':egamma}.items()}
record.update(frequency_step_exact='1/20',frequency_limit=T,physical_radius=R)
print('Closed whole action error',egamma,time.monotonic()-start,flush=True)

# Whole physical cross Grams of the binary action, with all sample uncertainties.
SA,SE={},{}
afrow=np.linalg.norm(action,axis=1)*(1+1e-10)
for i in range(3):
    f=fh[:,i*56:(i+1)*56]
    SA[i]=f.T@action/den/math.pi
    SE[i]=reserve(float(np.dot(rowerror[:,i],afrow))/den/math.pi
                  +product_error(f.T,action)/den/math.pi+gamma(5)*fro(SA[i]))
scorecross=score@action/den/math.pi
scorecross_error=reserve(float(np.dot(scoree,afrow))/den/math.pi
    +product_error(score.reshape(1,-1),action)/den/math.pi+gamma(5)*fro(scorecross.reshape(1,-1)))

# Native outward convolution on exact binary Lhat and What samples.
limit=128*den
nu=np.arange(limit+1,dtype=float)/den
kernel,kerre=transform(nu,gamma(3)*nu)
kernelfull=np.concatenate((kernel[:0:-1],kernel))
kernelerr=reserve(float(kerre[0]+2*kerre[1:].sum()))
kernelerr += 10*(-(arb(128)+h)/4).exp()/(1-(-h/4).exp())
full=np.concatenate((-action[:0:-1],action),axis=0)
length=len(full)
nfft=1<<int((length+len(kernelfull)-2).bit_length())
kvec=kernelfull.tolist()+[0.0]*(nfft-len(kernelfull))
kfft=acb.dft(kvec)
convolution=np.empty_like(full)
dft_radius=np.zeros(56)
for j in range(56):
    vec=full[:,j].tolist()+[0.0]*(nfft-length)
    aft=acb.dft(vec)
    result=acb.dft([a*k for a,k in zip(aft,kfft)],inverse=True)
    radii=np.empty(length)
    for i,z in enumerate(result[limit:limit+length]):
        convolution[i,j],radii[i]=ballfloat(z.real)
    dft_radius[j]=fro(radii.reshape(-1,1))
    if j%8==0:
        print('Closed native DFT convolution columns',j,'/56',time.monotonic()-start,flush=True)
del kfft,aft,result
c=1/(2*den*math.pi)
Jraw=c*c*(full.T@convolution)
Jraw_error=(h/(2*pi))**2*(kernelerr*reserve(fro(full))**2
    +reserve(fro(full))*reserve(fro(dft_radius.reshape(-1,1))))
Jraw_error+=reserve(c*c*product_error(full.T,convolution)+gamma(5)*fro(Jraw))
J=Jraw-np.outer(scorecross,scorecross)/float(Q.mid())
Jerror=Jraw_error+(2*reserve(fro(scorecross.reshape(1,-1)))*scorecross_error+scorecross_error**2)/Q
Jerror+=reserve(fro(scorecross.reshape(1,-1)))**2*abs(1/Q-1/arb(float(Q.mid())))
Jerror+=reserve(gamma(5)*(fro(Jraw)+fro(np.outer(scorecross,scorecross))/float(Q.lower())))
record['DFT']={'native':'acb.dft and inverse; all inputs exact binary; all output radii paid',
    'length':nfft,'kernel_halfwidth':limit,'kernel_l1_error':str(kernelerr),
    'extracted_convolution_radius_Frobenius':float(fro(dft_radius.reshape(-1,1))),
    'whole_smooth_action_Gram_error':str(Jerror)}

# Etilde action is smooth action plus B d; all output components are retained.
GY,GE=PG[1,1],PE[1,1]
CY,CYe=PG[0,1],PE[0,1]
YA=SA[1]+PG[1,2]
YAe=SE[1]+PE[1,2]+reserve(gamma(3)*(fro(SA[1])+fro(PG[1,2])))
CA=SA[0]+PG[0,2]
CAe=SE[0]+PE[0,2]+reserve(gamma(3)*(fro(SA[0])+fro(PG[0,2])))
AG=J+PG[2,2]+SA[2]+SA[2].T
AGe=Jerror+PE[2,2]+2*SE[2]+reserve(gamma(5)*(fro(J)+fro(PG[2,2])+2*fro(SA[2])))
s=Fraction(1,20)-Fraction(979236349,2000000000000)
sf=float(s)
shifterr=abs(arbf(s)-arb(sf))
Kh=CY+CY.T-sf*GY-(YA+YA.T)/2
Ke=2*CYe+arbf(s)*GE+YAe+shifterr*reserve(fro(GY))
Ke+=reserve(gamma(6)*(2*fro(CY)+abs(sf)*fro(GY)+fro(YA)))
Rh=H-CA-CA.T-sf*(CY+CY.T)+AG+sf*(YA+YA.T)+sf*sf*GY
Re=reserve(herror)+2*CAe+2*arbf(s)*CYe+AGe+2*arbf(s)*YAe+arbf(s)**2*GE
Re+=shifterr*reserve(2*fro(CY)+2*fro(YA)+3*abs(sf)*fro(GY))
Re+=reserve(gamma(12)*(fro(H)+2*fro(CA)+2*abs(sf)*fro(CY)+fro(AG)
                     +2*abs(sf)*fro(YA)+sf*sf*fro(GY)))
for mat in [GY,Kh,Rh]:
    assert np.isfinite(mat).all()
epsilon=2*(arbf(Fraction('0.255985513'))*arbf(Fraction(21,1024000))).sqrt()
epsilon+=arbf(Fraction(21,1024000))+arbf(Fraction('0.00855266'))
ynorm=(spectral_upper(GY,'Y_Gram')+GE).sqrt()
choices=[]
for beta in [Fraction(2),Fraction(5,2),Fraction(3)]:
    # tau=epsilon*||Y||/egamma minimizes the scalar combined error price.
    tau=epsilon*ynorm/egamma
    bfac=1+arbf(beta)
    dfac=1+1/arbf(beta)
    ycoef=epsilon+dfac*(1+tau)*epsilon**2/arbf(s)
    eyeprice=egamma*ynorm+dfac*(1+1/tau)*egamma**2/arbf(s)
    cycoef=float(ycoef.mid())
    crcoef=float((bfac/arbf(s)).mid())
    upper=1.1*(Kh+cycoef*GY+crcoef*Rh)
    assembly=reserve(gamma(6)*fro(upper))
    matrixerror=arbf(Fraction(11,10))*(Ke+ycoef*GE+bfac/arbf(s)*Re+eyeprice)
    matrixerror+=arbf(Fraction(11,10))*(abs(ycoef-arb(cycoef))*reserve(fro(GY))
        +abs(bfac/arbf(s)-arb(crcoef))*reserve(fro(Rh)))+assembly
    top=spectral_upper(upper,'coupled_beta_'+str(beta)) + matrixerror
    choices.append({'beta_exact':str(beta),'tau_ball':str(tau),
                    'central_largest_eigenvalue':float(np.linalg.eigvalsh((upper+upper.T)/2)[-1]),
                    'finite_matrix_error_ball':str(matrixerror),
                    'whole_price_upper_ball':str(top), 'pass':bool(top<arbf(Fraction(99999999,100000000))),
                    '_upper':top})
best=min(choices,key=lambda r:float(r['_upper'].upper()))
record['finite_Gram_errors']={name:str(value) for name,value in {
    'Y_Gram':GE,'C_Y_cross':CYe,'Y_action_cross':YAe,'C_action_cross':CAe,
    'whole_action_Gram':AGe,'K':Ke,'whole_residual_Gram':Re,'Y_norm_upper':ynorm}.items()}
record['choices']=[{k:v for k,v in r.items() if k!='_upper'} for r in choices]
record['status']='PASS full RE5 enriched-trial certificate' if best['pass'] else 'NOT PASS: explicit paid upper bound exceeds target'
record['chosen_beta_exact']=best['beta_exact']
record['chosen_price_upper_ball']=str(best['_upper'])
record['elapsed_seconds']=time.monotonic()-start
out=p/'complement_residual_enrichment_action_matrices.npz'
np.savez_compressed(out,b=b,v=v,d=d,GY=GY,CY=CY,YA=YA,CA=CA,AG=AG,K=Kh,residual_Gram=Rh)
record['matrix_sha256']=sha(out)
(p/'complement_residual_enrichment_action.json').write_text(json.dumps(record,indent=2)+'\n')
print(json.dumps({k:record[k] for k in ['status','chosen_price_upper_ball','finite_Gram_errors','elapsed_seconds']},indent=2),flush=True)
