"""Coupled whole-residual resolvent certificate; immutable completed Gram inputs."""
from pathlib import Path
from fractions import Fraction
import ast
import hashlib
import json
import math
import time
import numpy as np
from flint import arb, arb_mat, ctx

ctx.prec = 192
started = time.monotonic()
root = Path(__file__).resolve().parents[1]
owner = root / 'computations/complement_energy_resolvent.py'
mainpath = root / 'certificates/complement_energy_resolvent.json'
sha = lambda p: hashlib.sha256(p.read_bytes()).hexdigest()
main = json.loads(mainpath.read_text())
assert sha(owner) == main['script_sha256']
# The accepted owner supplies its complete absolute-product/Frobenius guards only.
names = {'up', 'gamma', 'fro', 'product_error', 'arbf'}
definitions = [n for n in ast.parse(owner.read_text()).body
               if isinstance(n, ast.FunctionDef) and n.name in names]
assert {n.name for n in definitions} == names
unit = 2.0 ** -53
exec(compile(ast.Module(body=definitions, type_ignores=[]), str(owner), 'exec'))
paths = {key: root / 'certificates' / filename for key, filename in {
    'factor': 'complement_midpoint_rank_factors.npz',
    'physical': 'complement_physical_frame_certificate.json',
    'physical_locator': 'complement_physical_frame_locator.npz',
    'sampling': 'complement_low_zero_frame_certificate.json',
    'sampling_matrices': 'complement_low_zero_frame_matrices.npz',
}.items()}
cachepath = root / 'certificates/complement_energy_resolvent_matrices.npz'
assert sha(cachepath) == main['matrix_sha256']
assert all(sha(p) == main['input_sha256'][k] for k, p in paths.items())
physical, sampling = [json.loads(paths[k].read_text()) for k in ['physical', 'sampling']]
Z = np.load(paths['factor'])['U']
loc = np.load(paths['physical_locator'])
H, N = loc['Hcenter'], loc['Wcenter']
low = np.load(paths['sampling_matrices'])
sample, Gzero, L = [low[k] for k in ['A_center', 'Gzero_center', 'L_center']]
cache = np.load(cachepath)
GB, Qb, GBZ, mid = [cache[k] for k in ['GB', 'Qb', 'GBZ', 'mid']]
assert Z.shape == (2906, 56) and sample.shape == (2469, 56)
assert GB.shape == (2906, 2906) and Qb.shape == (2469, 2906)
assert all(np.isfinite(a).all() for a in [Z, H, N, sample, Gzero, L, GB, Qb, GBZ, mid])

rho = Fraction(979236349, 2000000000000)
eta = Fraction('0.013155646076930')
nlo = Fraction(physical['coefficient_Gram_lower_exact'])
nhi = Fraction(physical['coefficient_Gram_upper_exact'])
assert arb(main['continuous_midpoint_error_ball']) < arbf(eta)
gerror, gnorm, qerror = [main[k] for k in [
    'B_Gram_error_operator_upper', 'B_Gram_operator_upper', 'zero_B_error_Frobenius_upper']]
gzerror, gznorm, sample_error, lerror = [sampling[k] for k in [
    'Gzero_error_operator_upper', 'Gzero_center_operator_upper',
    'A_evaluation_rounding_Frobenius_error', 'L_error_operator_upper']]
herror = physical['physical_compression_error_upper']

def exact(a):
    return arb_mat([[arb(float(x)) for x in row] for row in a])

def sym(a):
    return (a + a.transpose()) / 2

def afrob(a):
    return sum((a[i, j] ** 2 for i in range(a.nrows())
                for j in range(a.ncols())), arb(0)).sqrt()

def rational_upper(a, denominator=10 ** 12):
    f = Fraction.from_float(up(a.upper()))
    numerator = (f.numerator * denominator + f.denominator - 1) // f.denominator
    q = Fraction(numerator, denominator)
    assert a < arbf(q)
    return q

def guard(x):
    # This reserve encloses scalar evaluation of each accepted binary guard.
    return arb(float(x)) + arbf(Fraction(1, 10 ** 12))

def pe(x, y):
    return guard(product_error(x, y))

def ff(x):
    return guard(fro(x))

def locate(a, center, key):
    """Full exact-input Arb residual and polar correction; solver is a locator."""
    values, v = np.linalg.eigh((center + center.T) / 2)
    av = exact(v)
    identity = arb_mat(np.eye(len(values)).astype(int).tolist())
    diagonal = arb_mat([[arb(float(values[i])) if i == j else arb(0)
                         for j in range(len(values))] for i in range(len(values))])
    o = afrob(av.transpose() * av - identity)
    r = afrob(a * av - av * diagonal)
    assert o < arbf(Fraction(1, 1000000))
    lb = arb(float(max(abs(values))))
    correction = (r / (1 - o).sqrt()
                  + 2 * lb * (1 + o).sqrt() * ((1 - o).rsqrt() - 1))
    upper = arb(float(values[-1])) + correction
    matrices['locator_values_' + key] = values
    matrices['locator_vectors_' + key] = v
    return upper, {'orthogonality_Frobenius_ball': str(o),
                   'whole_residual_Frobenius_ball': str(r),
                   'polar_correction_ball': str(correction),
                   'largest_eigenvalue_upper_ball': str(upper)}

gbzerror = guard(gerror) * arbf(nhi).sqrt() + pe(GB, Z)
miderror = 2 * ff(GBZ) * gbzerror + gbzerror ** 2 + pe(GBZ.T, GBZ)
epsilons = [Fraction(1, 10), Fraction(1, 5), Fraction(1, 4), Fraction(1, 3),
            Fraction(1, 2), Fraction(2, 3), Fraction(1), Fraction(3, 2),
            Fraction(2), Fraction(3)]
records, matrices = [], {}
AH, AL, AM = [exact(a) for a in [H, L, mid]]
for mu in [Fraction(1, 20), Fraction(3, 100), Fraction(1, 40), Fraction(1, 25)]:
    s = mu - rho
    sf = float(s)
    T = L + mid + sf * H
    X = np.linalg.solve(T, H)
    assert np.isfinite(X).all()
    key = str(mu).replace('/', '_')
    AX = exact(X)
    xn2, xn_locator = locate(AX.transpose() * AX, X.T @ X, 'X_' + key)
    xn = xn2.sqrt()
    Y = X.T @ H @ X
    K = H @ X + X.T @ H - X.T @ T @ X
    AY, AK = sym(exact(Y)), sym(exact(K))
    exact_Ycenter = AX.transpose() * AH * AX
    exact_Kcenter = AH * AX + AX.transpose() * AH - AX.transpose() * (AL + AM + arbf(s) * AH) * AX
    ye = guard(herror) * xn2 + afrob(AY - exact_Ycenter)
    ke = (2 * guard(herror) * xn
          + (guard(lerror) + miderror + arbf(s) * guard(herror)) * xn2
          + afrob(AK - exact_Kcenter))
    DX = np.eye(56) - sf * X
    first, second = Z @ DX, GBZ @ X
    D1, D2 = first - second, -2 * sample @ X
    shift_error = abs(arbf(Fraction.from_float(sf) - s))
    e1 = (pe(Z, DX) + pe(GBZ, X) + gbzerror * xn
          + guard(gamma(3)) * (ff(first) + ff(second))
          + (guard(gamma(3)) * arb(sf) + shift_error) * ff(X) * arbf(nhi).sqrt())
    e2 = 2 * (guard(sample_error) * xn + pe(sample, X)) + guard(gamma(3)) * ff(D2)
    ed = (guard(gnorm + gerror).sqrt() * e1
          + guard(gznorm + gzerror).sqrt() * e2)
    S = D1.T @ (GB @ D1) + 2 * D1.T @ (Qb.T @ D2) + D2.T @ (Gzero @ D2)
    AS = sym(exact(S))
    d1, d2 = ff(D1), ff(D2)
    re = guard(gerror) * d1 ** 2 + 2 * guard(qerror) * d1 * d2 + guard(gzerror) * d2 ** 2
    re += (pe(D1.T, GB @ D1) + d1 * pe(GB, D1)
           + 2 * pe(D1.T, Qb.T @ D2) + 2 * d1 * pe(Qb.T, D2)
           + pe(D2.T, Gzero @ D2) + d2 * pe(Gzero, D2))
    re += guard(gamma(5)) * (ff(S) + d1 ** 2 * guard(gnorm)
                            + 2 * d1 * d2 * ff(Qb) + d2 ** 2 * guard(gznorm))
    sn, sn_locator = locate(AS, S, 'residual_' + key)
    rcoeff = re + 2 * (sn + re).sqrt() * ed + ed ** 2
    finite = {name: rational_upper(value) for name, value in {
        'X_operator_upper': xn, 'Y_Gram_error_upper': ye, 'K_error_upper': ke,
        'D1_coefficient_error_upper': e1, 'D2_coefficient_error_upper': e2,
        'physical_coefficient_error_upper': ed, 'residual_Gram_error_upper': re,
        'full_midpoint_residual_Gram_error_upper': rcoeff,
        'GBZ_error_upper': gbzerror, 'mid_error_upper': miderror}.items()}
    choices = []
    for epsilon in epsilons:
        ay = eta + (1 + 1 / epsilon) * eta ** 2 / s
        ass = (1 + epsilon) / s
        center = 1.1 * (K + float(ay) * Y + float(ass) * S)
        beta = Fraction(11, 10) * (finite['K_error_upper']
               + ay * finite['Y_Gram_error_upper']
               + ass * finite['full_midpoint_residual_Gram_error_upper'])
        choices.append((np.linalg.eigvalsh((center + center.T) / 2)[-1] + float(beta),
                        epsilon, ay, ass, beta, center))
    _, epsilon, ay, ass, beta, center = min(choices, key=lambda item: item[0])
    coupled = arbf(Fraction(11, 10)) * (AK + arbf(ay) * AY + arbf(ass) * AS)
    top, coupled_locator = locate(coupled, center, 'coupled_' + key)
    # The final scalar reserve is exact and additional to all listed errors.
    reserve = Fraction(1, 1000000)
    price = top + arbf(beta + reserve)
    p = rational_upper(price)
    theta = Fraction((p / nlo * 1000).__ceil__() + 1, 1000)
    margin = theta * nlo - p
    assert margin > 0 and price < arbf(p)
    loss = rho + theta * s
    passes = theta < 1
    record = {'mu_exact': str(mu), 'shift_exact': str(s), 'epsilon_exact': str(epsilon),
              'passes': passes, 'coupled_Y_coefficient_exact': str(ay),
              'coupled_residual_coefficient_exact': str(ass),
              'finite_error_guards_exact': {k: str(v) for k, v in finite.items()},
              'matrix_error_price_exact': str(beta), 'final_scalar_reserve_exact': str(reserve),
              'complete_price_upper_exact': str(p), 'complete_price_Arb_ball': str(price),
              'theta_exact': str(theta), 'coefficient_target_margin_exact': str(margin),
              'retained_energy_coefficient_exact': str(1 - theta) if passes else None,
              'full_floor_exact': str(loss) if passes else None,
              'locators': {'X_Gram': xn_locator, 'whole_residual_Gram': sn_locator,
                           'coupled_matrix': coupled_locator}}
    records.append(record)
    for name, value in {'X': X, 'Ygram': Y, 'K': K, 'residual_Gram': S,
                        'D1': D1, 'D2': D2, 'coupled_locator_center': center}.items():
        matrices[name + '_' + key] = value
    print('PASS' if passes else 'FAIL sufficient criterion', 'mu', mu,
          'epsilon', epsilon, 'price', float(p), 'theta', theta,
          'floor', float(loss) if passes else None, flush=True)

assert records[0]['passes']
output = root / 'certificates/complement_matrix_residual_certificate_matrices.npz'
np.savez_compressed(output, **matrices)
record = {'status': 'PASS coupled whole-residual certificate on exact original quotient',
          'source': 'routes/complement-matrix-residual-certificate.md',
          'method': 'Young matrix residual plus complete coefficient-map and Gram errors',
          'rho_exact': str(rho), 'eta_upper_exact': str(eta),
          'coefficient_metric': 'exact N=Z*Z; N_lower exact from accepted physical certificate',
          'coefficient_N_lower_exact': str(nlo), 'coefficient_N_upper_exact': str(nhi),
          'shifts': records, 'tested_shift_set_exact': [r['mu_exact'] for r in records],
          'epsilon_locator_set_exact': [str(e) for e in epsilons],
          'optimality_claim': False, 'main_result_and_inputs_unchanged': True,
          'all_couplings': 'all 2906 B features, all 2469 V0 features and complete rectangular mixed Gram',
          'no_zero_Gram_grid_or_integral_loop': True, 'Arb_precision_bits': ctx.prec,
          'owning_main_certificate_sha256': sha(mainpath), 'owning_guards_sha256': sha(owner),
          'completed_full_Gram_sha256': sha(cachepath),
          'input_sha256': {k: sha(p) for k, p in paths.items()},
          'matrix_sha256': sha(output), 'script_sha256': sha(Path(__file__)),
          'elapsed_seconds': time.monotonic() - started}
(root / 'certificates/complement_matrix_residual_certificate.json').write_text(json.dumps(record, indent=2) + '\n')
