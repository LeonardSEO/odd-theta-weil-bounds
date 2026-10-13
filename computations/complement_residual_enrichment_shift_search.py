"""Adaptive fixed-enriched-trial shifts; reuse all completed physical/action data."""
from pathlib import Path
from fractions import Fraction
import ast
import hashlib
import json
import time
import numpy as np
from flint import arb, arb_mat, ctx

ctx.prec=192
start=time.monotonic()
root=Path(__file__).resolve().parents[1]
p=root/'certificates'
sha=lambda path:hashlib.sha256(path.read_bytes()).hexdigest()
stagepath=p/'complement_residual_enrichment_action.json'
matrixpath=p/'complement_residual_enrichment_action_matrices.npz'
stage=json.loads(stagepath.read_text())
snapshot=p/'complement_residual_enrichment_action_executed_source.py'
assert sha(snapshot)==stage['source_sha256']
assert sha(matrixpath)==stage['matrix_sha256']
assert all(sha(p/n)==h for n,h in stage['input_sha256'].items())
assert stage['DFT']['native'].startswith('acb.dft and inverse')
matrix=np.load(matrixpath)
GY,Kh,Rh=[matrix[n] for n in ['GY','K','residual_Gram']]
record={'status':'INCOMPLETE','precision_bits':ctx.prec,'locators':{},
    'action_record_sha256':sha(stagepath),'action_matrix_sha256':sha(matrixpath),
    'action_executed_source_sha256':sha(snapshot),'source_sha256':sha(Path(__file__)),
    'route_sha256':sha(root/'routes/complement-residual-enrichment.md'),
    'calculation':'RE5 assembly only; completed point, whole-action, cross/Gram and native DFT stages reused',
    'input_sha256':stage['input_sha256']}
owner=root/'computations/complement_energy_resolvent.py'
main=json.loads((p/'complement_energy_resolvent.json').read_text())
assert sha(owner)==main['script_sha256']
names={'up','gamma','fro','arbf'}
unit=2.0**-53
defs=[n for n in ast.parse(owner.read_text()).body if isinstance(n,ast.FunctionDef) and n.name in names]
exec(compile(ast.Module(body=defs,type_ignores=[]),str(owner),'exec'))
names={'reserve','exact','afrob','spectral_upper'}
defs=[n for n in ast.parse(snapshot.read_text()).body if isinstance(n,ast.FunctionDef) and n.name in names]
exec(compile(ast.Module(body=defs,type_ignores=[]),str(snapshot),'exec'))
GE,Ke,Re,ynorm=[arb(stage['finite_Gram_errors'][n]) for n in ['Y_Gram','K','whole_residual_Gram','Y_norm_upper']]
egamma=arb(stage['action_prices']['whole_map_error'])
# The stage tail used binary (zero_center+zero_error) as its envelope shift.
# Recompute that small scalar price with the exact sum of the two dyadics;
# this pays its addition rounding without repeating any successful node/DFT work.
saved=np.load(p/'complement_residual_enrichment_pilot_matrices.npz')
low=np.load(p/'complement_low_zero_frame_matrices.npz')
cells=json.loads((p/'complement_signed_symbol_norm_cells.json').read_text())['kernel_cells']
times=np.array([(2*c['k']+1)/128 for c in cells])
mu=np.array([np.nextafter(int(c['measure_upper_dyadic_numerator'])/2**80,np.inf)
             for c in cells]+[2.0])
sm=np.sqrt(mu)
sme=np.array([up(abs(arb(float(m)).sqrt()-arb(float(t))).upper()) for m,t in zip(mu,sm)])
trialb,trialv=matrix['b'],matrix['v']
wb=sm[:,None]*trialb
wbe=sme[:,None]*abs(trialb)+gamma(3)*abs(wb)
cn=np.array([up(np.linalg.norm(wb[j])*(1+1e-10)+np.linalg.norm(wbe[j])*(1+1e-10)) for j in range(2906)])
vn=np.array([up(np.linalg.norm(trialv[j])*(1+1e-10)) for j in range(2469)])
weighted=sum((arb(float(c))*((arb(float(a))/4).exp()+(-arb(float(a))/4).exp())/2
              for a,c in zip(times,cn[:-1])),arb(0))+arb(float(cn[-1]))
for a,e,c in zip(low['zero_center'],low['zero_error'],vn):
    shift=arb(float(a))+arb(float(e))
    weighted+=arb(float(c))*((shift/4).exp()+(-shift/4).exp())/2
# A/20 >= 47200 ||alpha||; this avoids rebuilding the completed overlap grid.
alpha_upper=arb(stage['action_prices']['strip_L1_A'])/(20*47200)
h=arb(1)/20
r=(-h/4).exp()
first=arb(3300)+h
ss0=1/(1-r)
ss1=first*ss0+h*r/(1-r)**2
ss2=first**2*ss0+2*first*h*r/(1-r)**2+h**2*r*(1+r)/(1-r)**3
tail=h/arb.pi()*(-first/4).exp()*5*(weighted*(200*ss0+10*ss1)+alpha_upper*(200*ss1+10*ss2))
egamma+=tail/arb(2).sqrt()
record['additional_exact_shift_frequency_tail_ball']=str(tail)

rho=Fraction(979236349,2000000000000)
nlo=Fraction(99999999,100000000)
epsilon=2*(arbf(Fraction('0.255985513'))*arbf(Fraction(21,1024000))).sqrt()
epsilon+=arbf(Fraction(21,1024000))+arbf(Fraction('0.00855266'))
CY,YA,CA,AG=[matrix[n] for n in ['CY','YA','CA','AG']]
H=np.load(p/'complement_physical_frame_locator.npz')['Hcenter']
physical=json.loads((p/'complement_physical_frame_certificate.json').read_text())
herror=physical['physical_compression_error_upper']
CYe,YAe,CAe,AGe=[arb(stage['finite_Gram_errors'][n]) for n in
    ['C_Y_cross','Y_action_cross','C_action_cross','whole_action_Gram']]
tau=egamma/(epsilon*ynorm)
rows={}
scale=10**12

def ceiling(value):
    binary=Fraction(up(value.upper()))
    out=Fraction(-(-binary.numerator*scale//binary.denominator),scale)
    assert value < arbf(out)
    return out

def evaluate(mu,beta):
    if (mu,beta) in rows: return rows[mu,beta]
    s=mu-rho
    assert s>0 and beta>0
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
    bfac=1+arbf(beta)
    dfac=1+1/arbf(beta)
    ycoef=epsilon+dfac*(1+tau)*epsilon**2/arbf(s)
    eyeprice=egamma*ynorm+dfac*(1+1/tau)*egamma**2/arbf(s)
    cycoef=float(ycoef.mid())
    crcoef=float((bfac/arbf(s)).mid())
    upper=1.1*(Kh+cycoef*GY+crcoef*Rh)
    assembly=afrob(exact(upper)-arbf(Fraction(11,10))*(exact(Kh)
        +arb(cycoef)*exact(GY)+arb(crcoef)*exact(Rh)))
    err=arbf(Fraction(11,10))*(Ke+ycoef*GE+bfac/arbf(s)*Re+eyeprice)
    err+=arbf(Fraction(11,10))*(abs(ycoef-arb(cycoef))*reserve(fro(GY))
        +abs(bfac/arbf(s)-arb(crcoef))*reserve(fro(Rh)))+assembly
    top=spectral_upper(upper,'mu_'+str(mu)+'_beta_'+str(beta))+err
    theta=ceiling(top/arbf(nlo))
    alpha=rho+theta*s
    row=dict(mu_exact=str(mu),beta_exact=str(beta),theta_exact=str(theta),
        retained_energy_exact=str(1-theta),alpha_exact=str(alpha),
        theta_display=float(theta),retained_energy_display=float(1-theta),
        alpha_display=float(alpha),passes=theta<1,price_ball=str(top),
        matrix_error_ball=str(err),K_error_ball=str(Ke),residual_Gram_error_ball=str(Re))
    rows[mu,beta]=row
    return row

# No monotone scalar-certificate premise; all candidate points are independently certified.
shifts=[Fraction(n,1000) for n in [20,25,30,35,40,45,50,60,75,100,150]]
betas=[Fraction(1,8),Fraction(1,4),Fraction(3,8)]+[Fraction(n,4) for n in [2,3,4,5,6,7,8,9,10,12]]
for mu in shifts:
    for beta in betas: evaluate(mu,beta)
# Independently refine scalar, energy, and mu=1/20 tradeoffs with rational neighborhoods.
# This search uses no bisection or exclusion assertion: every new point pays full RE5.
for iteration in range(9):
    valid=[key for key,row in rows.items() if row['passes']]
    chosen=[min(valid,key=lambda key:Fraction(rows[key]['alpha_exact'])),
            min(valid,key=lambda key:Fraction(rows[key]['theta_exact'])),
            min([key for key in rows if key[0]==Fraction(1,20)],
                key=lambda key:Fraction(rows[key]['theta_exact']))]
    dm=Fraction(1,400*2**iteration)
    dbeta=Fraction(1,8*2**iteration)
    for index,(mu,beta) in enumerate(chosen):
        for candidate_beta in [beta-dbeta,beta+dbeta]:
            if candidate_beta>0: evaluate(mu,candidate_beta)
        if index<2:
            for candidate_mu in [mu-dm,mu+dm]:
                if candidate_mu>rho: evaluate(candidate_mu,beta)
valid=[row for row in rows.values() if row['passes']]
assert valid
bestloss=min(valid,key=lambda row:Fraction(row['alpha_exact']))
bestenergy=min(valid,key=lambda row:Fraction(row['theta_exact']))
pareto=[row for row in valid if not any(
    Fraction(o['alpha_exact'])<=Fraction(row['alpha_exact']) and
    Fraction(o['theta_exact'])<=Fraction(row['theta_exact']) and
    (Fraction(o['alpha_exact'])<Fraction(row['alpha_exact']) or
     Fraction(o['theta_exact'])<Fraction(row['theta_exact'])) for o in valid)]
record.update(status='PASS cached fixed enriched-trial adaptive shift certificates',
    shifts=list(rows.values()),best_tested_scalar=bestloss,best_tested_energy=bestenergy,
    pareto=sorted(pareto,key=lambda row:Fraction(row['alpha_exact'])),
    search_scope='Finite rational candidates and local refinement only; no excluded intervals or proved global optimum',
    rho_exact=str(rho),epsilon_ball=str(epsilon),tau_ball=str(tau),
    lattice_exact=str(Fraction(1,scale)),whole_action_error_ball=str(egamma),
    trial='Same exact dyadic Y as action certificate; no new domain, rank, or inverse assertion',
    elapsed_seconds=time.monotonic()-start)
(p/'complement_residual_enrichment_shift_search.json').write_text(json.dumps(record,indent=2)+'\n')
print(json.dumps(dict(count=len(rows),bestloss=bestloss,bestenergy=bestenergy,pareto=record['pareto'],seconds=record['elapsed_seconds']),indent=2))
