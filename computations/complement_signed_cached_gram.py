"""Signed whole-Gram assembly from existing exact binary action data only."""
from pathlib import Path
from fractions import Fraction
import hashlib,json,time
import numpy as np
from flint import arb,arb_mat,ctx
ctx.prec=192
start=time.monotonic();root=Path(__file__).resolve().parents[1];p=root/'certificates'
sha=lambda f:hashlib.sha256(f.read_bytes()).hexdigest()
action=json.loads((p/'complement_residual_enrichment_action.json').read_text())
final=json.loads((p/'complement_residual_enrichment_action_certificate.json').read_text())
physical=json.loads((p/'complement_physical_frame_certificate.json').read_text())
assert sha(p/'complement_residual_enrichment_action.json')==final['action_record_sha256']
assert sha(p/'complement_residual_enrichment_action_matrices.npz')==final['action_matrix_sha256']
assert sha(p/'complement_residual_enrichment_action_executed_source.py')==final['action_executed_source_sha256']
assert all(sha(p/n)==h for n,h in action['input_sha256'].items())
a=np.load(p/'complement_residual_enrichment_action_matrices.npz')
l=np.load(p/'complement_physical_frame_locator.npz')
H,N=l['Hcenter'],l['Wcenter']; CY,YA,CA,AG,GY=[a[n] for n in ['CY','YA','CA','AG','GY']]
# Float solves/eigenvectors select exact dyadic maps, never proof endpoints.
X=np.linalg.solve(N,CY)
_,V=np.linalg.eigh((GY+GY.T)/2);W=V[:,-5:]
rec={'status':'INCOMPLETE','precision_bits':ctx.prec,'locators':{},
 'scope':'Complete signed action and whole Grams; certified 5-dimensional enriched positive subspace, no global sign',
 'input_sha256':{n:sha(p/n) for n in ['complement_residual_enrichment_action.json','complement_residual_enrichment_action_matrices.npz','complement_residual_enrichment_action_certificate.json','complement_physical_frame_certificate.json','complement_physical_frame_locator.npz']},
 'source_sha256':sha(Path(__file__)),'reuse':'No integration, zero validation, point grid or DFT replay; exact-input finite matrix assembly only'}
def rat(x):
 f=Fraction(x);return arb(f.numerator)/f.denominator
def exact(A):return arb_mat([[arb(float(x)) for x in row] for row in A])
def norm(A):return sum((abs(A[i,j]).upper()**2 for i in range(A.nrows()) for j in range(A.ncols())),arb(0)).sqrt()
def sym(A):return (A+A.transpose())/2
def extrema(A,key):
 ac=np.array([[float(A[i,j].mid()) for j in range(A.ncols())] for i in range(A.nrows())]);vals,v=np.linalg.eigh((ac+ac.T)/2)
 av=exact(v);I=arb_mat(np.eye(len(vals),dtype=int).tolist());di=arb_mat([[arb(float(vals[i])) if i==j else arb(0) for j in range(len(vals))] for i in range(len(vals))])
 o=norm(av.transpose()*av-I);r=norm(A*av-av*di);assert o<rat('0.000001')
 m=arb(float(max(abs(vals))));cor=r/(1-o).sqrt()+2*m*(1+o).sqrt()*((1-o).rsqrt()-1)
 rec['locators'][key]={'orthogonality':str(o),'whole_residual':str(r),'polar_correction':str(cor)}
 return arb(float(vals[0]))-cor,arb(float(vals[-1]))+cor
E={k:arb(v) for k,v in action['finite_Gram_errors'].items()}
herr=arb(float(physical['physical_compression_error_upper']))
n0=rat('0.99999999');nerr=rat('0.00000002');eps=rat(Fraction(3,12500));t8=rat(Fraction(32237,10**12));yn=rat('2.159290');cn=rat('0.507');eg=rat('0.00001421821')
ne,xe,cye,he,yae,cae,age,gye,we=map(exact,[N,X,CY,H,YA,CA,AG,GY,W]);I=arb_mat(np.eye(56,dtype=int).tolist())
assert norm(ne-I)<rat('0.00000001')
xnorm=norm(xe);solve_res=norm(cye-ne*xe)
selector_error=(E['C_Y_cross']+solve_res+nerr*xnorm)/n0
error=eg+(eps+t8)*yn+cn*selector_error
print('action error',error,'selector',selector_error,'residual',solve_res,'X',xnorm,flush=True)
assert error<rat('0.00056'); e=rat('0.00056')
# Exact Arb assembly of center matrices: no floating product is an endpoint.
dcraw=yae-cye.transpose()*xe;dc=sym(dcraw);zc=cae-he*xe
ac=sym(age-cae.transpose()*xe-xe.transpose()*cae+xe.transpose()*he*xe)
derr=E['Y_action_cross']+E['C_Y_cross']*xnorm
zerr=E['C_action_cross']+herr*xnorm
aerr=E['whole_action_Gram']+2*E['C_action_cross']*xnorm+herr*xnorm**2
lo,hi=extrema(dc,'signed_central_56'); alo,ahi=extrema(ac,'signed_action_56')
# Restrict to exact dyadic W; pay its true norm by certified orthogonality.
wdef=norm(we.transpose()*we-arb_mat(np.eye(5,dtype=int).tolist()));assert wdef<rat('0.000001')
wn2=1+wdef;wn=wn2.sqrt();g5=sym(we.transpose()*gye*we);d5=sym(we.transpose()*dc*we);a5=sym(we.transpose()*ac*we)
raw5=we.transpose()*dcraw*we
glo,ghi=extrema(g5,'restricted_physical_Gram');dlo,dhi=extrema(d5,'restricted_signed_central');a5lo,a5hi=extrema(a5,'restricted_signed_action')
physical_lower=glo-E['Y_Gram']*wn2
signed_lower=dlo-derr*wn2-e*yn*wn2
assert physical_lower>rat('0.5')
assert signed_lower>rat('0.001')
# Full residual against a dyadic projection selector; no inverse is a premise.
mid=lambda A:np.array([[float(A[i,j].mid()) for j in range(A.ncols())] for i in range(A.nrows())])
L=np.linalg.solve(mid(g5),mid(raw5));le=exact(L);ln=norm(le)
r5=sym(a5-raw5.transpose()*le-le.transpose()*raw5+le.transpose()*g5*le)
rerr=wn2*(aerr+2*derr*ln+E['Y_Gram']*ln**2)
rlo,rhi=extrema(r5,'whole_signed_trial_residual')
# The actual P_F Zhat coefficient is Gtrue^-1 rawtrue; price its selector error.
coeferror=(norm(raw5-g5*le)+wn2*(derr+E['Y_Gram']*ln))/physical_lower
# Certified lower witness of actual outside-trial coupling, coefficient metric.
rmaxlower=rhi-2*arb(rec['locators']['whole_signed_trial_residual']['polar_correction'])
coupling_lower=(rmaxlower-rerr).lower().sqrt()-yn*wn*coeferror-e*wn
assert coupling_lower>rat('0.22')
beta=rat('0.01')
coupling_upper=(1+beta)*(rhi+rerr)+(1+1/beta)*e**2*wn2
# ||YW x||^2 <= ||Y||^2 ||W||^2 ||x||^2; retain exact reserve.
assert yn**2*wn2<rat('5')
# Hence a(YW x)> (1/1000)||x||^2 >= (1/5000)||YW x||^2.
rec.update(status='PASS cached whole signed Grams and enriched 5-dimensional positivity',
 selector_frobenius_ball=str(xnorm),selector_native_residual_ball=str(solve_res),coefficient_center_defect_ball=str(norm(ne-I)),
 coefficient_center_uncertainty_exact='1/50000000',capture_error_exact='3/12500',
 whole_signed_action_error_ball=str(error),whole_signed_action_error_endpoint_exact='7/12500',
 signed_central_Gram_error_ball=str(derr),C_signed_action_cross_error_ball=str(zerr),whole_signed_action_Gram_error_ball=str(aerr),
 restricted_selector_orthogonality_ball=str(wdef),restricted_physical_Gram_lower_ball=str(physical_lower),
 restricted_signed_form_coefficient_lower_ball=str(signed_lower),
 restricted_rank=5,restricted_signed_norm_lower_exact='1/5000',
 whole_signed_trial_residual_error_ball=str(rerr),whole_signed_trial_residual_upper_ball=str(rhi+rerr),
 outside_trial_coupling_norm_lower_ball=str(coupling_lower),outside_trial_coupling_norm_lower_exact='11/50',
 outside_trial_coupling_squared_upper_ball=str(coupling_upper),
 coupling_metric='input coefficient Euclidean norm of exact dyadic YW; no normalization or tail gap assumed',
 central56_lower_before_action_error_ball=str(lo-derr),
 full56_rank_unestablished='Y Gram uncertainty exceeds minimum center eigenvalue; no enlarged inverse used',
 conditioning='restricted exact physical Gram >=1/2 I; signed action error after orthonormalization <=sqrt(2)*(7/12500)*sqrt(1+wdef)',
 elapsed_seconds=time.monotonic()-start)
# Persist exact binary selectors and finite centers, paying each downcast.
cent=lambda A:np.array([[float(A[i,j].mid()) for j in range(A.ncols())] for i in range(A.nrows())])
items={'Dcenter':(dc,derr),'Cactioncenter':(zc,zerr),'actionGramcenter':(ac,aerr),
 'restrictedGcenter':(g5,E['Y_Gram']*wn2),'restrictedDcenter':(d5,derr*wn2),
 'restrictedAGcenter':(a5,aerr*wn2),'restrictedRawCross':(raw5,derr*wn2),
 'restrictedResidualcenter':(r5,rerr)}
centers={k:cent(A) for k,(A,err) in items.items()}
rec['saved_center_operator_error_balls']={k:str(err+norm(A-exact(centers[k]))) for k,(A,err) in items.items()}
rec['saved_center_conversion_error_balls']={k:str(norm(A-exact(centers[k]))) for k,(A,err) in items.items()}
out=p/'complement_signed_cached_gram_matrices.npz'
np.savez_compressed(out,X=X,W=W,L=L,**centers)
rec['matrix_sha256']=sha(out)
(p/'complement_signed_cached_gram_certificate.json').write_text(json.dumps(rec,indent=2)+'\n')
print(json.dumps({k:rec[k] for k in ['status','whole_signed_action_error_ball','restricted_physical_Gram_lower_ball','restricted_signed_form_coefficient_lower_ball','outside_trial_coupling_norm_lower_ball','outside_trial_coupling_squared_upper_ball','elapsed_seconds']},indent=2))
