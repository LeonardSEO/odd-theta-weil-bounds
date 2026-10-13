"""RE5 assembly only; reuse the completed outward enriched action calculation."""
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
s=Fraction(1,20)-Fraction(979236349,2000000000000)
epsilon=2*(arbf(Fraction('0.255985513'))*arbf(Fraction(21,1024000))).sqrt()
epsilon+=arbf(Fraction(21,1024000))+arbf(Fraction('0.00855266'))
choices=[]
for beta in [Fraction(2),Fraction(5,2),Fraction(3)]:
    tau=egamma/(epsilon*ynorm)
    bfac=1+arbf(beta)
    dfac=1+1/arbf(beta)
    ycoef=epsilon+dfac*(1+tau)*epsilon**2/arbf(s)
    eyeprice=egamma*ynorm+dfac*(1+1/tau)*egamma**2/arbf(s)
    cycoef=float(ycoef.mid())
    crcoef=float((bfac/arbf(s)).mid())
    upper=1.1*(Kh+cycoef*GY+crcoef*Rh)
    assembly=afrob(exact(upper)-arbf(Fraction(11,10))*(exact(Kh)
        +arb(cycoef)*exact(GY)+arb(crcoef)*exact(Rh)))
    matrixerror=arbf(Fraction(11,10))*(Ke+ycoef*GE+bfac/arbf(s)*Re+eyeprice)
    matrixerror+=arbf(Fraction(11,10))*(abs(ycoef-arb(cycoef))*reserve(fro(GY))
        +abs(bfac/arbf(s)-arb(crcoef))*reserve(fro(Rh)))+assembly
    top=spectral_upper(upper,'coupled_beta_'+str(beta))+matrixerror
    choices.append({'beta_exact':str(beta),'tau_ball':str(tau),
        'central_largest_eigenvalue':float(np.linalg.eigvalsh((upper+upper.T)/2)[-1]),
        'finite_matrix_error_ball':str(matrixerror),'whole_price_upper_ball':str(top),
        '_upper':top})
best=min(choices,key=lambda x:float(x['_upper'].upper()))
target=Fraction(79,100)*Fraction(99999999,100000000)
record['status']='PASS full RE5 enriched-trial certificate' if best['_upper']<arbf(target) else 'NOT PASS chosen rational target'
record['choices']=[{k:v for k,v in c.items() if k!='_upper'} for c in choices]
record.update(chosen_beta_exact=best['beta_exact'],chosen_price_upper_ball=str(best['_upper']),
    rational_resolvent_coefficient_exact='79/100',coefficient_metric_lower_exact='99999999/100000000',
    exact_target=str(target),mu_exact='1/20',rho_exact='979236349/2000000000000',
    shift_exact=str(s),whole_action_error_ball=str(egamma),elapsed_seconds=time.monotonic()-start)
record['retained_energy_coefficient_exact']='21/100'
record['full_floor_exact']=str(Fraction(979236349,2000000000000)+Fraction(79,100)*s)
record['analytical_scope']='Entire original odd form domain under RE1-RE5 and existing original-operator/theta/score/frame ER premises; no even or Li exhaustion claim'
(p/'complement_residual_enrichment_action_certificate.json').write_text(json.dumps(record,indent=2)+'\n')
print(json.dumps(record,indent=2),flush=True)
