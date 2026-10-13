"""One new shift, mu=2/25, using the completed full physical Gram data."""
from pathlib import Path
from fractions import Fraction
import ast
import hashlib
import json
import math
import time
import numpy as np
from flint import arb, ctx

ctx.prec = 192
started = time.monotonic()
root = Path(__file__).resolve().parents[1]
owner = root/'computations/complement_energy_resolvent.py'
mainpath = root/'certificates/complement_energy_resolvent.json'
main = json.loads(mainpath.read_text())
digest = lambda p: hashlib.sha256(p.read_bytes()).hexdigest()
assert digest(owner) == main['script_sha256']
# Directly reuse only the owning numerical guards; do not execute its run.
names = {'up','gamma','fro','product_error','arbf','boundnorm'}
parsed = ast.parse(owner.read_text())
definitions = [node for node in parsed.body
               if isinstance(node,ast.FunctionDef) and node.name in names]
assert {node.name for node in definitions} == names
unit = 2.0**-53
exec(compile(ast.Module(body=definitions,type_ignores=[]),str(owner),'exec'))
paths = {key:root/'certificates'/filename for key,filename in {
    'frame':'complement_midpoint_rank_factors.npz',
    'physical':'complement_physical_frame_certificate.json',
    'physical_locator':'complement_physical_frame_locator.npz',
    'sampling':'complement_low_zero_frame_certificate.json',
    'sampling_matrices':'complement_low_zero_frame_matrices.npz',
    'completed_full_Gram':'complement_energy_resolvent_matrices.npz',
}.items()}
assert digest(paths['completed_full_Gram']) == main['matrix_sha256']
assert digest(paths['frame']) == main['input_sha256']['factor']
assert digest(paths['physical_locator']) == main['input_sha256']['physical_locator']
assert digest(paths['sampling_matrices']) == main['input_sha256']['sampling_matrices']
physical,sampling = [json.loads(paths[key].read_text()) for key in ['physical','sampling']]
assert digest(paths['frame']) == physical['source_hashes']['factors']
assert digest(paths['physical_locator']) == physical['physical_locator_file_sha256']
assert digest(paths['sampling_matrices']) == sampling['matrix_file_sha256']
Z = np.load(paths['frame'])['U']
physical_data = np.load(paths['physical_locator'])
low = np.load(paths['sampling_matrices'])
completed = np.load(paths['completed_full_Gram'])
H,N = physical_data['Hcenter'],physical_data['Wcenter']
sample,Gzero,L = low['A_center'],low['Gzero_center'],low['L_center']
GB,Qb,GBZ,mid = [completed[key] for key in ['GB','Qb','GBZ','mid']]
gerror,gnorm,qerror = [main[key] for key in ['B_Gram_error_operator_upper','B_Gram_operator_upper','zero_B_error_Frobenius_upper']]
gzerror,gznorm,sample_error,lerror = [sampling[key] for key in ['Gzero_error_operator_upper','Gzero_center_operator_upper','A_evaluation_rounding_Frobenius_error','L_error_operator_upper']]
herror = physical['physical_compression_error_upper']
GBZerror = gerror*math.sqrt(1+1e-8)+product_error(GB,Z)
miderror = 2*fro(GBZ)*GBZerror+GBZerror**2+product_error(GBZ.T,GBZ)
eta = up(arb(main['continuous_midpoint_error_ball']).upper())
rho = Fraction(979236349,2000000000000)
mu_exact = Fraction(2,25)
delta_exact = mu_exact-rho
delta = float(delta_exact)
T = L+mid+delta*H
X = np.linalg.solve(T,H)
assert np.isfinite(X).all()
Xnorm = up(math.sqrt(boundnorm(X.T@X)[1]+product_error(X.T,X)))
Ygram = X.T@H@X
Yerror = herror*Xnorm**2+product_error(X.T,H@X)+Xnorm*product_error(H,X)
Ynorm = math.sqrt(max(0,boundnorm(Ygram)[1]+Yerror))
DX = np.eye(56)-delta*X
D1,D2 = Z@DX-GBZ@X,-2*sample@X
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
K = H@X+X.T@H-X.T@T@X
ke = (2*herror*Xnorm+(lerror+miderror+delta*herror)*Xnorm**2
      +2*product_error(H,X)+product_error(X.T,T@X)+Xnorm*product_error(T,X)
      +gamma(5)*(2*fro(H@X)+fro(X.T@T@X)))
kupper = boundnorm(K+eta*Ygram)[1]+ke+eta*Yerror
residual_upper = rnorm+eta*Ynorm
priceball = arbf(Fraction(11,10))*(arb(up(kupper+1e-6))
              +arb(up(residual_upper+1e-6))**2/arbf(delta_exact))
price = up(priceball.upper())
print('PASS new mu',mu_exact,'K',kupper,'residual',residual_upper,'full price',price,flush=True)
theta = Fraction(97,100)
captured_exact,residual_exact = Fraction(563,1000),Fraction(79,500)
assert arb(kupper)+arbf(Fraction(1,10**6)) < arbf(captured_exact)
assert arb(residual_upper)+arbf(Fraction(1,10**6)) < arbf(residual_exact)
price_exact = Fraction(11,10)*(captured_exact+residual_exact**2/delta_exact)
assert price_exact < theta*Fraction(99999999,100000000)
loss = rho+theta*delta_exact
assert loss < Fraction(2,25)
record = {'status':'PASS shift 2/25' if price<1-1e-8 else 'FAIL sufficient criterion at shift 2/25',
          'source':'routes/complement-energy-resolvent-tighter.md',
          'mu_exact':str(mu_exact),'delta_exact':str(delta_exact),'eta_upper':eta,
          'trial_X_operator_upper':Xnorm,'trial_Y_physical_operator_upper':Ynorm,
          'midpoint_physical_residual_upper':rnorm,'actual_physical_residual_upper':residual_upper,
          'captured_K_upper':kupper,'resolvent_price_upper':price,
          'captured_price_part':1.1*kupper,'residual_price_part':1.1*residual_upper**2/delta,
          'finite_K_error_upper':ke,'finite_residual_Gram_error_upper':re,
          'scalar_reserve_exact':'1/1000000','coefficient_target_lower_exact':'99999999/100000000',
          'all_couplings':'unchanged full physical B/zero Gram and inherited complete continuous/prime/zero errors',
          'owning_main_certificate_sha256':digest(mainpath),'owning_utilities_sha256':digest(owner),
          'input_sha256':{key:digest(path) for key,path in paths.items()},
          'main_result_unchanged':True,'no_zero_grid_Gram_loop':True}
record.update(chosen_K_upper_exact=str(captured_exact),chosen_residual_upper_exact=str(residual_exact),
              chosen_resolvent_upper_exact=str(theta),chosen_price_exact=str(price_exact),
              retained_energy_coefficient_exact=str(1-theta),full_floor_exact=str(loss))
output = root/'certificates/complement_energy_resolvent_tighter_matrices.npz'
np.savez_compressed(output,X=X,Ygram=Ygram,K=K,residual_Gram=residue,D1=D1,D2=D2)
record['matrix_sha256'] = digest(output)
record['script_sha256'] = digest(Path(__file__))
record['elapsed_seconds'] = time.monotonic()-started
(root/'certificates/complement_energy_resolvent_tighter.json').write_text(json.dumps(record,indent=2)+'\n')
