"""Exact endpoint postprocessing only; no matrix calculation is repeated."""
from pathlib import Path
from fractions import Fraction
import hashlib
import json
from flint import arb, ctx

ctx.prec = 192
root = Path(__file__).resolve().parents[1]
path = root / 'certificates/complement_matrix_residual_certificate.json'
record = json.loads(path.read_text())
rho = Fraction(record['rho_exact'])
nlo = Fraction(record['coefficient_N_lower_exact'])
def aq(q):
    return arb(q.numerator) / q.denominator
results = []
for previous in record['shifts']:
    if not previous['passes']:
        continue
    price = Fraction(previous['complete_price_upper_exact'])
    assert arb(previous['complete_price_Arb_ball']) < aq(price)
    theta = Fraction((price / nlo * 10 ** 6).__ceil__(), 10 ** 6)
    margin = theta * nlo - price
    assert margin > 0 and 0 < theta < 1
    loss = rho + theta * Fraction(previous['shift_exact'])
    results.append({'mu_exact': previous['mu_exact'],
                    'epsilon_exact': previous['epsilon_exact'],
                    'unchanged_complete_price_upper_exact': str(price),
                    'theta_exact': str(theta),
                    'coefficient_target_margin_exact': str(margin),
                    'retained_energy_coefficient_exact': str(1 - theta),
                    'full_floor_exact': str(loss)})
    print(previous['mu_exact'], 'theta', theta, 'retained', 1 - theta,
          'floor', loss, float(loss), 'margin', margin)
record['tight_exact_postprocessing'] = {
    'denominator': 10 ** 6, 'matrix_or_trial_run_repeated': False,
    'method': 'ceil(unchanged exact guarded price/Nlower * 10^6)/10^6',
    'source': 'computations/complement_matrix_residual_tight_endpoints.py',
    'source_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
    'results': results}
path.write_text(json.dumps(record, indent=2) + '\n')
