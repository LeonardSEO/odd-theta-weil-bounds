"""Actual complete low-zero sampling of the stored exact 56-column frame.

Uses the already installed Python/FLINT runtime; no dependency installation.
Every numpy value used to define a matrix is interpreted as exact binary64.
All finite evaluation, zero-location, and matrix-arithmetic prices are retained.
"""
from pathlib import Path
from fractions import Fraction
import hashlib
import json
import math
import time
import numpy as np
from flint import arb, acb, ctx

ctx.prec = 192
root = Path(__file__).resolve().parents[1]
started = time.monotonic()
cutoff = 3000
order = 20
frequency_limit = 128
grid_denominator = 8
quadrature_order = 320
unit = 2.0**-53


def up(x):
    return float(np.nextafter(float(x), np.inf))


def gamma(n):
    return up(n*unit/(1-n*unit))


def ball_float(x):
    center = float(x.mid())
    radius = up(abs(x-arb(center)).upper())
    return center, radius


def norm_upper(x):
    # The 1e-12 relative guard includes square/sum/sqrt and this scalar price.
    return up(float(np.linalg.norm(x, 'fro'))*(1+1e-12)+1e-300)


def product_error(x, y):
    absolute = np.abs(x) @ np.abs(y)
    return gamma(2*x.shape[1]+10)*absolute/(1-gamma(2*x.shape[1]+10))


def hashfile(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


cellpath = root/'certificates/complement_signed_symbol_norm_cells.json'
factorpath = root/'certificates/complement_midpoint_rank_factors.npz'
cells = json.loads(cellpath.read_text())['kernel_cells']
k = np.array([c['k'] for c in cells], dtype=np.int64)
# This is exactly the upward binary mass choice in complement_midpoint_rank.py.
mu = np.array([np.nextafter(int(c['measure_upper_dyadic_numerator'])/(2**80),
                            np.inf) for c in cells])
mu = np.append(mu, 2.0)
times = (2*k+1)/128
Z = np.load(factorpath)['U']
assert Z.shape == (len(mu), 56)
assert np.max(np.abs(Z)) <= 1 and np.min(mu) > 0
sqrtmu = np.sqrt(mu)
sqrtmu_error = np.array([ball_float(arb(float(m)).sqrt()-arb(float(s)))[1]
                         + abs(ball_float(arb(float(m)).sqrt()-arb(float(s)))[0])
                         for m, s in zip(mu, sqrtmu)])
sqrtmu_error = np.nextafter(sqrtmu_error*(1+1e-12), np.inf)

# Native rigorous zero isolation plus the complete multiplicity count.
zero_count = int(arb(cutoff).zeta_nzeros().unique_fmpz())
roots = acb.zeta_zeros(1, zero_count)
assert all(r.real == arb(1)/2 for r in roots)
assert roots[0].imag.lower() > 0 and roots[-1].imag.upper() < cutoff
assert all(a.imag.upper() < b.imag.lower() for a, b in zip(roots, roots[1:]))
# There are zero_count distinct actual roots and zero_count total multiplicities.
# Hence this complete list has multiplicity one at every returned ordinate.
root_scale = 2**128
root_pairs = [[int((r.imag.lower()*root_scale).floor().unique_fmpz()),
               int((r.imag.upper()*root_scale).ceil().unique_fmpz())]
              for r in roots]
zerofloat = np.array([ball_float(r.imag)[0] for r in roots])
zeroerror = np.array([ball_float(r.imag)[1] for r in roots])
print('PASS actual zero intervals/count', zero_count, 'through', cutoff, flush=True)

# New arbitrary-frequency Taylor evaluator, not a replay of the old grid loop.
# Derivative j receives 4^j times the analytic ellipse/theta price and the
# common real-tail price. The complete Taylor remainder is 120(2|h|)^order.
rho = Fraction(23, 20)
egl = (24*240*Fraction(3)**(frequency_limit//4)*rho**(-2*quadrature_order)
       /(1-1/rho))
theta_price = 32*Fraction(1, 2)**128
real_tail = Fraction(1, 2)**1100
errors = [arb(str((4**j*egl+4**j*theta_price+real_tail).numerator))
          /arb(str((4**j*egl+4**j*theta_price+real_tail).denominator))
          for j in range(order)]
nodes = []
pi = arb.pi()
for j in range(quadrature_order):
    x, q = arb.legendre_p_root(quadrature_order, j, weight=True)
    x = (x+1)*3/2
    q *= 3/2
    phi = arb(0)
    for n in range(1, 9):
        y = pi*n*n*(2*x).exp()
        phi += (x/2-y).exp()*(4*y*y-6*y)
    W = phi/(x/2).cosh()
    weight = 2*q*W
    wp = weight*(x/2).sinh()
    powers = [arb(1)]
    for d in range(1, order):
        powers.append(powers[-1]*x/d)
    nodes.append((x, [weight*p for p in powers], [wp*p for p in powers]))

grid_size = frequency_limit*grid_denominator+1
wc = np.empty((grid_size, order)); we = np.empty_like(wc)
vc = np.empty_like(wc); ve = np.empty_like(wc)
for j in range(grid_size):
    frequency = arb(j)/grid_denominator
    wvals = [arb(0) for _ in range(order)]
    vvals = [arb(0) for _ in range(order)]
    for x, wpower, vpower in nodes:
        sine = (frequency*x).sin(); cosine = (frequency*x).cos()
        ctrig = [cosine, -sine, -cosine, sine]
        strig = [sine, cosine, -sine, -cosine]
        for d in range(order):
            wvals[d] += wpower[d]*ctrig[d % 4]
            vvals[d] += vpower[d]*strig[d % 4]
    for d in range(order):
        allowance = errors[d]/math.factorial(d)
        wc[j, d], we[j, d] = ball_float(wvals[d]+arb(0, allowance))
        vc[j, d], ve[j, d] = ball_float(vvals[d]+arb(0, allowance))
    if j % 256 == 0:
        print('Certified Taylor grid', j, '/', grid_size-1, flush=True)

assert up(float(egl)) < 1e-15
gridpath = root/'certificates/complement_low_zero_taylor_grid.npz'
np.savez_compressed(gridpath, W_center=wc, W_error=we,
                    pole_center=vc, pole_error=ve)

# A dyadic lower argument of the exponential is used for every far entry;
# no floating exponential is a certificate premise.
far_size = math.ceil((cutoff+float(times.max())+1)*grid_denominator)+2
far = np.zeros(far_size)
for j in range(frequency_limit*grid_denominator-2, far_size):
    far[j] = up((5*(-arb(j)/(4*grid_denominator)).exp()).upper())


def fourier(nu, argument_error, pole=False):
    absolute = np.abs(nu)
    near = absolute <= frequency_limit
    index = np.minimum(np.rint(absolute*grid_denominator).astype(np.int64),
                       grid_size-1)
    h = absolute-index/grid_denominator
    herror = gamma(3)*(absolute+index/grid_denominator)
    center_grid, radius_grid = (vc, ve) if pole else (wc, we)
    value = np.zeros_like(absolute)
    magnitude = np.zeros_like(absolute)
    interval = np.zeros_like(absolute)
    for d in range(order-1, -1, -1):
        value = value*h+center_grid[index, d]
        magnitude = magnitude*np.abs(h)+np.abs(center_grid[index, d])
        interval = interval*np.abs(h)+radius_grid[index, d]
    # |W'|<.11; |v'|<.012, both proved in the original Fourier route.
    derivative_bound = .012 if pole else .11
    remainder = 120*(2*(np.abs(h)+herror))**order
    error = (interval*(1+gamma(2*order+10))+remainder
             +gamma(2*order+10)*magnitude
             +derivative_bound*(argument_error+herror))*(1+1e-12)+1e-300
    far_index = np.maximum(np.floor(absolute*grid_denominator).astype(np.int64)-1,
                           frequency_limit*grid_denominator-2)
    value = np.where(near, value, 0.0)
    error = np.where(near, error, far[far_index])
    if pole:
        value *= np.sign(nu)
    return value, np.nextafter(error, np.inf)


# The exact score subtraction vanishes at these actual zeros, and only here.
samples = np.empty((zero_count, 56)); sample_error = np.empty_like(samples)
absZ = np.abs(Z)
for j, (g, ge) in enumerate(zip(zerofloat, zeroerror)):
    am = g-times; ap = g+times
    ae = ge+gamma(3)*(abs(g)+times)
    wm, em = fourier(am, ae); wp, ep = fourier(ap, ae)
    vp, ev = fourier(np.array([g]), np.array([ge]), pole=True)
    unscaled = np.append((wm-wp)/2, vp[0])
    unscaled_error = np.append((em+ep)/2+gamma(3)*(np.abs(wm)+np.abs(wp)), ev[0])
    row = sqrtmu*unscaled
    erow = (sqrtmu*unscaled_error
            +sqrtmu_error*(np.abs(unscaled)+unscaled_error)
            +gamma(3)*np.abs(row))*(1+1e-12)
    samples[j] = row @ Z
    sample_error[j] = ((erow @ absZ)/(1-gamma(2*len(mu)+10))
                       +gamma(2*len(mu)+10)*(np.abs(row) @ absZ)
                       /(1-gamma(2*len(mu)+10)))*(1+1e-12)+1e-300
sampling_error = norm_upper(sample_error)
sampling_norm = norm_upper(samples)
print('Complete A evaluation/rounding Frobenius error', sampling_error, flush=True)

# Zero-representer Gram on the same odd Hilbert space (measure 2 dx).
Gzero = np.empty((zero_count, zero_count)); Gerror = np.empty_like(Gzero)
for j, (g, ge) in enumerate(zip(zerofloat, zeroerror)):
    argument_error = ge+zeroerror+gamma(3)*(abs(g)+np.abs(zerofloat))
    dm, de = fourier(g-zerofloat, argument_error)
    sm, se = fourier(g+zerofloat, argument_error)
    Gzero[j] = (dm-sm)/2
    Gerror[j] = (de+se)/2+gamma(3)*(np.abs(dm)+np.abs(sm))
Gzero = (Gzero+Gzero.T)/2
Gerror = (Gerror+Gerror.T)/2
gerror = up(float(Gerror.sum(axis=1).max())/(1-gamma(2*zero_count+10)))
gnorm = up(float(np.abs(Gzero).sum(axis=1).max())/(1-gamma(2*zero_count+10)))

# Full finite matrices, including every low-zero/low-zero coupling.
L = 2*samples.T @ samples
lround = 2*norm_upper(product_error(samples.T, samples))
lerror = up(2*(2*sampling_norm*sampling_error+sampling_error**2)+lround)
GA = Gzero @ samples
GAround = norm_upper(product_error(Gzero, samples))
O = 4*samples.T @ GA
oround = (4*sampling_norm*GAround
          +4*norm_upper(product_error(samples.T, GA)))
oerror = up(4*((2*sampling_norm*sampling_error+sampling_error**2)*(gnorm+gerror)
              +sampling_norm**2*gerror)+oround)
L = (L+L.T)/2; O = (O+O.T)/2

# Certify the rectangular singular factor directly: the tiny eigenvalue is
# never obtained by subtracting a coarse normal-Gram error from an eigensolve.
U, singular, Vt = np.linalg.svd(samples, full_matrices=False)
identity = np.eye(56)
udefect = U.T @ U-identity
vdefect = Vt @ Vt.T-identity
uerror = up(norm_upper(udefect)+norm_upper(product_error(U.T, U))+1e-14)
verror = up(norm_upper(vdefect)+norm_upper(product_error(Vt, Vt.T))+1e-14)
US = U*singular
approximation = US @ Vt
residual = samples-approximation
resprice = (norm_upper(residual)+norm_upper(product_error(US, Vt))
            +gamma(3)*(norm_upper(samples)+norm_upper(approximation)))
sigma_lower_ball = (arb(float(singular[-1]))*(1-arb(uerror)).sqrt()
                    *(1-arb(verror)).sqrt()-arb(up(resprice))-arb(sampling_error))
assert sigma_lower_ball > 0
sigma_lower = float(sigma_lower_ball.lower())
lambda_lower = 2*sigma_lower_ball**2
assert lambda_lower > arb('1e-15')
print('PASS complete actual sampling L >= 1e-15 I; sigma lower',
      sigma_lower_ball, flush=True)

arraypath = root/'certificates/complement_low_zero_frame_matrices.npz'
np.savez_compressed(arraypath, A_center=samples, A_error=sample_error,
                    zero_center=zerofloat, zero_error=zeroerror,
                    Gzero_center=Gzero, Gzero_error=Gerror,
                    L_center=L, O_center=O,
                    svd_U=U, svd_singular=singular, svd_Vt=Vt)
record = {
    'status': 'PASS complete actual low-zero sampling, rank 56',
    'cutoff': cutoff, 'positive_zero_count': zero_count,
    'multiplicity': 'all one: distinct certified roots exhaust total count',
    'zero_endpoint_scale': str(root_scale), 'zero_outward_endpoint_pairs': root_pairs,
    'root_precision_bits': ctx.prec,
    'zero_method': 'native acb.zeta_zeros(1,N), arb(3000).zeta_nzeros()',
    'exact_frame': 'stored binary64 U=Z; upward binary masses from owning source; pole mass 2',
    'fourier_convention': 'nonunitary + sign; FT(wC alpha)(gamma)=i A_gamma alpha',
    'both_ordinate_signs': 'L=2 A* A and O=4 A* Gzero A',
    'quadrature_nodes': quadrature_order, 'theta_atoms': 8,
    'Taylor_order': order, 'Taylor_grid_denominator': grid_denominator,
    'Taylor_frequency_limit': frequency_limit,
    'ellipse_GL_error_exact': str(egl),
    'theta_derivative_price': '4^j 32/2^128',
    'real_derivative_tail_price': '2^-1100 for j=0,...,19',
    'Taylor_remainder': '120 (2 abs(h))^20, including h rounding',
    'far_Fourier_envelope': '5 exp(-abs(nu)/4), dyadic lower argument, native Arb',
    'A_evaluation_rounding_Frobenius_error': sampling_error,
    'A_center_Frobenius_upper': sampling_norm,
    'Gzero_error_operator_upper': gerror,
    'Gzero_center_operator_upper': gnorm,
    'L_error_operator_upper': lerror, 'O_error_operator_upper': oerror,
    'SVD_U_orthogonality_defect_upper': uerror,
    'SVD_V_orthogonality_defect_upper': verror,
    'SVD_whole_residual_upper': up(resprice),
    'sampling_sigma_lower': sigma_lower,
    'coefficient_sampling_L_lower_exact': '1/1000000000000000',
    'all_coupling_retained': 'all 2469 positive zeros; full zero Gram; both ordinate signs',
    'scope': 'S_3000 only; all omitted zeros require the separate actual-band tail',
    'script_sha256': hashfile(Path(__file__)),
    'cell_certificate_sha256': hashfile(cellpath),
    'factor_file_sha256': hashfile(factorpath),
    'Taylor_grid_sha256': hashfile(gridpath),
    'matrix_file_sha256': hashfile(arraypath),
    'elapsed_seconds': time.monotonic()-started,
}
(root/'certificates/complement_low_zero_frame_certificate.json').write_text(
    json.dumps(record, indent=2)+'\n')
print('PASS saved complete low-zero certificate; L/O errors', lerror, oerror,
      'elapsed', record['elapsed_seconds'], flush=True)
