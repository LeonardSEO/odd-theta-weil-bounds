"""Diagnostic direct-d8,+ residual enrichment; no proof from quadrature samples."""
from pathlib import Path
from fractions import Fraction
import hashlib
import json
import time
import numpy as np
from scipy.special import digamma
from scipy.signal import fftconvolve
from numpy.polynomial.legendre import leggauss

root = Path(__file__).resolve().parents[1]
p = root / 'certificates'
start = time.monotonic()
main = json.loads((p / 'complement_energy_resolvent.json').read_text())
et = json.loads((p / 'complement_energy_resolvent_tighter.json').read_text())
digest = lambda x: hashlib.sha256(x.read_bytes()).hexdigest()
assert digest(p / 'complement_energy_resolvent_matrices.npz') == main['matrix_sha256']
assert digest(p / 'complement_energy_resolvent_tighter_matrices.npz') == et['matrix_sha256']
assert digest(p / 'complement_midpoint_rank_factors.npz') == main['input_sha256']['factor']
assert digest(p / 'complement_low_zero_frame_matrices.npz') == main['input_sha256']['sampling_matrices']
Z = np.load(p / 'complement_midpoint_rank_factors.npz')['U']
cache = np.load(p / 'complement_energy_resolvent_matrices.npz')
trial = np.load(p / 'complement_energy_resolvent_tighter_matrices.npz')
low = np.load(p / 'complement_low_zero_frame_matrices.npz')
phys = np.load(p / 'complement_physical_frame_locator.npz')
grid = np.load(p / 'complement_low_zero_taylor_grid.npz')
norms = json.loads((p / 'complement_physical_norms.json').read_text())
from flint import arb
Q = float(arb(norms['Q_ball']).mid())
cells = json.loads((p / 'complement_signed_symbol_norm_cells.json').read_text())['kernel_cells']
times = np.array([(2 * c['k'] + 1) / 128 for c in cells])
mu = np.array([np.nextafter(int(c['measure_upper_dyadic_numerator']) / 2**80, np.inf)
               for c in cells] + [2.0])
sm = np.sqrt(mu)
zeros = low['zero_center']
GB, Qb, H = cache['GB'], cache['Qb'], phys['Hcenter']
D1, D2 = trial['D1'], trial['D2']
# Exact projected residual is the analytical target; this center projection
# is exploratory only. A certificate must enclose its H^-1 and mixed Grams.
cross = Z.T @ (GB @ D1 + Qb.T @ D2)
projection = np.linalg.solve(H, cross)
Bcoef = np.concatenate((Z, D1 - Z @ projection), axis=1)
Vcoef = np.concatenate((np.zeros((len(zeros), 56)), D2), axis=1)
G = np.empty((112, 112))
G[:56, :56] = H
G[:56, 56:] = cross - H @ projection
G[56:, :56] = G[:56, 56:].T
G[56:, 56:] = (trial['residual_Gram'] - cross.T @ projection
               - projection.T @ cross + projection.T @ H @ projection)
G = (G + G.T) / 2

def transform(arg, pole=False):
    absolute = np.abs(arg)
    index = np.minimum(np.rint(absolute * 8).astype(np.int64), 1024)
    h = absolute - index / 8
    coeff = grid['pole_center' if pole else 'W_center']
    result = np.zeros_like(arg)
    for d in range(19, -1, -1):
        result = result * h + coeff[index, d]
    result = np.where(absolute <= 128, result, 0)
    return result * np.sign(arg) if pole else result

nodes, weights = leggauss(256)
x, weights = 1.5 * (nodes + 1), 3 * weights
y = np.pi * np.arange(1, 9)[None, :]**2 * np.exp(2 * x[:, None])
phi = np.sum(np.exp(x[:, None]/2 - y) * (4*y*y - 6*y), axis=1)
def score_transform(arg):
    result = np.zeros_like(arg)
    near = np.abs(arg) <= 128
    result[near] = -arg[near] * (np.cos(arg[near, None] * x) @ (weights * phi))
    return result

delta_t = 0.25
limit = 3100.0
t = np.arange(0, limit + delta_t/2, delta_t)
fh = np.empty((len(t), 112))
bscore = score_transform(times)
scorecoef = -(np.append(bscore, -0.25) * sm) @ Bcoef / Q
weighted_B = sm[:-1, None] * Bcoef[:-1]
for first in range(0, len(t), 256):
    tt = t[first:first+256, None]
    bf = (transform(tt-times) - transform(tt+times)) / 2
    vf = (transform(tt-zeros) - transform(tt+zeros)) / 2
    fh[first:first+len(tt)] = (bf @ weighted_B + vf @ Vcoef
        + transform(tt[:, 0], pole=True)[:, None] * sm[-1] * Bcoef[-1]
        + score_transform(tt[:, 0])[:, None] * scorecoef)
print('Fourier trial built', time.monotonic()-start, flush=True)
primes = np.array([2, 3, 4, 5, 7, 8], dtype=float)
vonmangoldt = np.log(np.array([2, 3, 2, 5, 7, 2], dtype=float))
k8 = np.cos(t[:, None] * np.log(primes)) @ (2 * vonmangoldt / np.sqrt(primes))
mGamma = np.real(digamma(0.25 + 0.5j*t)) - np.log(np.pi)
d = np.maximum(mGamma - k8, 0)
full = np.concatenate((-fh[:0:-1], fh), axis=0)
dfull = np.concatenate((d[:0:-1], d))
action = dfull[:, None] * full
A = delta_t / (2*np.pi) * (full.T @ action)
kernel_index = np.arange(-512, 513)
kernel = transform(kernel_index * delta_t)
conv = fftconvolve(action, kernel[:, None], mode='same', axes=0)
J = delta_t**2 / (4*np.pi**2) * (action.T @ conv)
bfull = np.concatenate((-score_transform(t)[:0:-1], score_transform(t)))
v = delta_t / (2*np.pi) * (bfull @ action)
J -= np.outer(v, v)/Q
A, J = (A+A.T)/2, (J+J.T)/2
# A second operator has the global ER2--ER3 error bridge on this new span:
# S L8 S + Rmid - 2(Sp_s) tensor (Sp_s). Its smooth multiplier is ell8.
ell = mGamma-k8
smooth_action = np.concatenate((ell[:0:-1], ell))[:, None]*full
smooth_A = delta_t/(2*np.pi)*(full.T@smooth_action)
smooth_conv = fftconvolve(smooth_action, kernel[:, None], mode='same', axes=0)
smooth_J = delta_t**2/(4*np.pi**2)*(smooth_action.T@smooth_conv)
smooth_v = delta_t/(2*np.pi)*(bfull@smooth_action)
smooth_J -= np.outer(smooth_v, smooth_v)/Q
smooth_B = np.zeros((2906, 112))
Bscorecoef = -np.append(bscore, -0.25)*sm/Q
for first in range(0, len(t), 256):
    tt = t[first:first+256, None]
    bh = np.empty((len(tt), 2906))
    bh[:, :-1] = (transform(tt-times)-transform(tt+times))/2*sm[:-1]
    bh[:, -1] = transform(tt[:, 0], pole=True)*sm[-1]
    bh += score_transform(tt[:, 0])[:, None]*Bscorecoef
    smooth_B += delta_t/np.pi*(bh.T@(ell[first:first+len(tt), None]*fh[first:first+len(tt)]))
BU = GB@Bcoef+Qb.T@Vcoef
finite_action_coef = BU.copy()
finite_action_coef[-1] = 0  # Rmid minus its exactly retained pole column.
bridge_A = smooth_A+BU.T@finite_action_coef
bridge_J = (smooth_J+finite_action_coef.T@GB@finite_action_coef
    +smooth_B.T@finite_action_coef+finite_action_coef.T@smooth_B)
bridge_A, bridge_J = (bridge_A+bridge_A.T)/2, (bridge_J+bridge_J.T)/2
# Whitening only chooses the trial. The actual selected binary coefficients
# need no numerical-inverse premise in the whole residual identity.
ge, gv = np.linalg.eigh(G)
select = ge > 1e-11
W = gv[:, select] / np.sqrt(ge[select])
s = float(Fraction(1, 20) - Fraction(979236349, 2000000000000))
M = A + s*G
X = W @ np.linalg.solve(W.T @ M @ W, W.T @ G[:, :56])
K = G[:56] @ X + X.T @ G[:, :56] - X.T @ M @ X
cross_action = A[:56] + s*G[:56]
residue = H - cross_action @ X - X.T @ cross_action.T
residue += X.T @ (J + 2*s*A + s*s*G) @ X
residue = (residue+residue.T)/2
price_matrix = 1.1 * (K + residue / s)
bridge_M = bridge_A+s*G
bridge_X = W@np.linalg.solve(W.T@bridge_M@W, W.T@G[:, :56])
bridge_K = G[:56]@bridge_X+bridge_X.T@G[:, :56]-bridge_X.T@bridge_M@bridge_X
bridge_cross = bridge_A[:56]+s*G[:56]
bridge_R = H-bridge_cross@bridge_X-bridge_X.T@bridge_cross.T
bridge_R += bridge_X.T@(bridge_J+2*s*bridge_A+s*s*G)@bridge_X
bridge_R = (bridge_R+bridge_R.T)/2
bridge_Y = bridge_X.T@G@bridge_X
epsilon_R = 2*np.sqrt(0.255985513*21/1024000)+21/1024000+0.00855266
from scipy.optimize import minimize_scalar
def bridged_price(beta):
    upper = bridge_K+epsilon_R*bridge_Y
    upper += ((1+beta)*bridge_R+(1+1/beta)*epsilon_R**2*bridge_Y)/s
    return float(np.linalg.eigvalsh((upper+upper.T)*0.55)[-1])
optimal = minimize_scalar(bridged_price, bounds=(0.01, 100), method='bounded')
out = {
    'status': 'DIAGNOSTIC ONLY: uncertified symbol/Fourier quadrature',
    'mu_exact': '1/20', 'frequency_step': delta_t, 'frequency_limit': limit,
    'trial_columns': 112, 'selected_center_Gram_rank': int(select.sum()),
    'center_Gram_minimum': float(ge[0]),
    'projected_midpoint_residual_cross_norm': float(np.linalg.norm(cross-H@projection)),
    'whole_action_Gram_minimum': float(np.linalg.eigvalsh(J)[0]),
    'whole_residual_Gram_minimum': float(np.linalg.eigvalsh(residue)[0]),
    'K_maximum': float(np.linalg.eigvalsh((K+K.T)/2)[-1]),
    'whole_residual_squared_maximum': float(np.linalg.eigvalsh(residue)[-1]),
    'matrix_resolvent_price_maximum': float(np.linalg.eigvalsh((price_matrix+price_matrix.T)/2)[-1]),
    'separate_norm_price': 1.1*(float(np.linalg.eigvalsh((K+K.T)/2)[-1])
        + float(np.linalg.eigvalsh(residue)[-1])/s),
    'global_midpoint_difference_upper': epsilon_R,
    'smooth_plus_midpoint_robust_beta': float(optimal.x),
    'smooth_plus_midpoint_robust_price': float(optimal.fun),
    'smooth_plus_midpoint_whole_residual_squared_maximum': float(np.linalg.eigvalsh(bridge_R)[-1]),
    'smooth_plus_midpoint_Gram_minimum': float(np.linalg.eigvalsh(bridge_J)[0]),
    'elapsed_seconds': time.monotonic()-start,
    'input_sha256': {name: digest(p/name) for name in [
        'complement_energy_resolvent_matrices.npz',
        'complement_energy_resolvent_tighter_matrices.npz',
        'complement_low_zero_frame_matrices.npz',
        'complement_midpoint_rank_factors.npz',
        'complement_low_zero_taylor_grid.npz']},
    'source_sha256': digest(Path(__file__)),
}
assert np.isfinite(X).all() and np.isfinite(price_matrix).all()
(p/'complement_residual_enrichment_pilot.json').write_text(json.dumps(out, indent=2)+'\n')
np.savez_compressed(p/'complement_residual_enrichment_pilot_matrices.npz',
    Bcoef=Bcoef, Vcoef=Vcoef, G=G, A=A, J=J, X=X, K=K, residual_Gram=residue)
print(json.dumps(out, indent=2), flush=True)
