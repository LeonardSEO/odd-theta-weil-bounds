"""Certified whole-theta Fourier W/pole grids for the new R8 norm question."""
from flint import arb, ctx
from fractions import Fraction
from pathlib import Path
import hashlib
import json
import time

ctx.prec = 192
root = Path(__file__).resolve().parents[1]
rho = Fraction(23, 20)
error = (24 * 240 * 3**25 * rho**(-512) / (1 - 1/rho)
         + 32 * Fraction(1, 2)**128 + Fraction(1, 2)**1160)
assert error < Fraction(1, 10**12)
radius = arb(str(error.numerator)) / arb(str(error.denominator))
pi = arb.pi()
nodes = []
for k in range(256):
    z, weight = arb.legendre_p_root(256, k, weight=True)
    x = (z + 1) * 3 / 2
    weight = weight * 3 / 2
    phi = arb(0)
    for n in range(1, 9):
        y = pi * n*n * (2*x).exp()
        phi += (x/2-y).exp() * (4*y*y-6*y)
    W = phi / (x/2).cosh()
    nodes.append((x, 2*weight*W, 2*weight*W*(x/2).sinh(), 2*weight*phi))

scale = 2**80
def endpoints(value):
    value += arb(0, radius)
    lower = int((value.lower()*scale).floor().unique_fmpz())
    upper = int((value.upper()*scale).ceil().unique_fmpz())
    return [lower, upper]

started = time.monotonic()
W_grid = []
W_derivative_grid = []
for k in range(6401):
    nu = arb(k)/64
    value = arb(0)
    derivative = arb(0)
    for x,a,_,_ in nodes:
        value += a*(nu*x).cos()
        derivative -= x*a*(nu*x).sin()
    W_grid.append(endpoints(value))
    W_derivative_grid.append(endpoints(derivative))
print('PASS W Fourier grid 6401 outward pairs', flush=True)
pole_grid = []
score_grid = []
for k in range(6400):
    nu = arb(2*k+1)/128
    value = arb(0)
    phi_value = arb(0)
    for x,_,a,p in nodes:
        value += a*(nu*x).sin()
        phi_value += p*(nu*x).cos()
    pole_grid.append(endpoints(value))
    # Add the complete Fourier error before multiplying by the frequency.
    phi_value += arb(0, radius)
    score_grid.append(endpoints(-nu*phi_value))
print('PASS pole Fourier grid 6400 outward pairs', flush=True)

record = {
    'status': 'PASS (actual certified execution)',
    'precision_bits': 192,
    'runtime': 'existing python-flint 0.9.0 Arb native certified Legendre nodes',
    'source': 'routes/complement-fourier-weight-grid-proof.md',
    'quadrature': '256-node native outward Gauss-Legendre on [0,3]',
    'theta_atoms': 8,
    'full_error_fraction': str(error),
    'full_error_upper_decimal': str(radius),
    'error_includes': ['analytic GL remainder', 'entire theta lattice after atom 8',
                      'both physical tails outside [-3,3]'],
    'endpoint_scale': str(scale),
    'W_frequency_grid': 'nu=k/64, k=0,...,6400',
    'W_outward_endpoint_pairs': W_grid,
    'W_derivative_outward_endpoint_pairs': W_derivative_grid,
    'pole_frequency_grid': 'nu=(2k+1)/128, k=0,...,6399',
    'pole_outward_endpoint_pairs': pole_grid,
    'score_overlap_outward_endpoint_pairs': score_grid,
    'W_large_frequency_envelope': '5 exp(-abs(nu)/4)',
    'pole_large_frequency_envelope': '5 exp(-abs(nu)/4)',
    'score_large_frequency_envelope': '5 abs(nu) exp(-abs(nu)/4)',
    'script_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
    'elapsed_seconds': time.monotonic()-started,
}
path = root/'certificates/complement_fourier_weight_grid.json'
path.write_text(json.dumps(record, separators=(',', ':'))+'\n')
print('PASS complete Fourier certificate', path, 'elapsed', record['elapsed_seconds'], flush=True)
