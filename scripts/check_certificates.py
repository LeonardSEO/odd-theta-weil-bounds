"""Read-only portable identity, finite spectral guard and exact endpoint check.

This checks the final 56-dimensional matrix comparison using saved full-error
prices. It does not regenerate upstream integrals, zeros or large Gram inputs.
"""
from pathlib import Path
from fractions import Fraction
import ast
import runpy
import sys
import hashlib
import json
import numpy as np
from flint import arb, arb_mat, ctx

if not __debug__:
    raise RuntimeError('Run without -O: certificate assertions must be enabled')
ctx.prec = 192
root = Path(__file__).resolve().parents[1]
manifest = json.loads((root / 'CERTIFICATE_MANIFEST.json').read_text())
for name, item in manifest['files'].items():
    data = (root / name).read_bytes()
    assert len(data) == item['bytes'], name
    assert hashlib.sha256(data).hexdigest() == item['sha256'], name

# Historical hashes are aliases only for explicitly verified provenance adaptations.
provenance_by_path = {}
for adaptation in manifest.get('provenance_adaptations', []):
    path = root / adaptation['release_file']
    assert hashlib.sha256(path.read_bytes()).hexdigest() == adaptation['release_sha256']
    payload = json.loads(path.read_text())
    assert adaptation['omitted_fields'] == ['execution_provenance.command', 'execution_provenance.cwd']
    for key in ['command', 'cwd']:
        assert payload['execution_provenance'].pop(key) == adaptation['release_values'][key]
    canonical = json.dumps(payload, sort_keys=True, separators=(',', ':'), ensure_ascii=False).encode()
    assert hashlib.sha256(canonical).hexdigest() == adaptation['canonical_payload_sha256']
    provenance_by_path[path.resolve()] = adaptation

def historical_digest(path):
    path = Path(path).resolve()
    digest = hashlib.sha256(path.read_bytes()).hexdigest()
    adaptation = provenance_by_path.get(path)
    if adaptation is not None:
        assert digest == adaptation['release_sha256']
        return adaptation['original_sha256']
    return digest

# Resolve embedded digests to included bytes or declared provenance-only adaptations.
def hashes(value, field=''):
    if isinstance(value, dict):
        for key, child in value.items():
            yield from hashes(child, field + '.' + key)
    elif isinstance(value, list):
        for index, child in enumerate(value):
            yield from hashes(child, field + '[' + str(index) + ']')
    elif isinstance(value, str) and len(value) == 64 and all(c in '0123456789abcdef' for c in value):
        yield field, value

links = manifest['certificate_hash_links']
for relative, entries in links.items():
    actual = dict(hashes(json.loads((root / relative).read_text())))
    assert set(actual) == set(entries), relative
    for field, target in entries.items():
        assert actual[field] == historical_digest(root / target), (relative, field)
assert set(links) == {str(p.relative_to(root)) for folder in ['certificates', 'routes'] for p in (root / folder).glob('*.json')}
for mapping in manifest['source_adaptations']:
    assert hashlib.sha256((root / mapping['original_file']).read_bytes()).hexdigest() == mapping['original_sha256']
    assert hashlib.sha256((root / mapping['release_file']).read_bytes()).hexdigest() == mapping['release_sha256']
print('PASS complete certificate input/source/owner hash closure and provenance payload identity')
if sys.argv[1:] == ['--hash-only']:
    print('PASS hash-only check; numerical checks were not rerun')
    raise SystemExit(0)
if sys.argv[1:]:
    raise SystemExit('Usage: check_certificates.py [--hash-only]')


def aq(value):
    value = Fraction(value)
    return arb(value.numerator) / value.denominator

def exact(array):
    return arb_mat([[arb(float(x)) for x in row] for row in array])

def symmetric(array):
    value = exact(array)
    return (value + value.transpose()) / 2

def frobenius(value):
    return sum((value[i, j] ** 2 for i in range(value.nrows())
                for j in range(value.ncols())), arb(0)).sqrt()

record = json.loads((root / 'certificates/complement_matrix_residual_certificate.json').read_text())
matrices = np.load(root / 'certificates/complement_matrix_residual_certificate_matrices.npz')
nlo = Fraction(record['coefficient_N_lower_exact'])
rho = Fraction(record['rho_exact'])
eta = Fraction(record['eta_upper_exact'])
for item in record['tight_exact_postprocessing']['results']:
    old = next(x for x in record['shifts'] if x['mu_exact'] == item['mu_exact'])
    key = item['mu_exact'].replace('/', '_')
    s = Fraction(old['shift_exact'])
    epsilon = Fraction(old['epsilon_exact'])
    ay = eta + (1 + 1 / epsilon) * eta ** 2 / s
    ass = (1 + epsilon) / s
    assert ay == Fraction(old['coupled_Y_coefficient_exact'])
    assert ass == Fraction(old['coupled_residual_coefficient_exact'])
    coupled = aq(Fraction(11, 10)) * (symmetric(matrices['K_' + key])
               + aq(ay) * symmetric(matrices['Ygram_' + key])
               + aq(ass) * symmetric(matrices['residual_Gram_' + key]))
    vectors = exact(matrices['locator_vectors_coupled_' + key])
    values = matrices['locator_values_coupled_' + key]
    diagonal = arb_mat([[arb(float(values[i])) if i == j else arb(0)
                        for j in range(len(values))] for i in range(len(values))])
    identity = arb_mat(np.eye(len(values)).astype(int).tolist())
    defect = frobenius(vectors.transpose() * vectors - identity)
    residual = frobenius(coupled * vectors - vectors * diagonal)
    assert defect < aq(Fraction(1, 1000000))
    correction = (residual / (1 - defect).sqrt()
                  + 2 * arb(float(max(abs(values)))) * (1 + defect).sqrt()
                    * ((1 - defect).rsqrt() - 1))
    upper = arb(float(values[-1])) + correction
    guards = {k: Fraction(v) for k, v in old['finite_error_guards_exact'].items()}
    beta = Fraction(11, 10) * (guards['K_error_upper']
           + ay * guards['Y_Gram_error_upper']
           + ass * guards['full_midpoint_residual_Gram_error_upper'])
    assert beta == Fraction(old['matrix_error_price_exact'])
    price = Fraction(item['unchanged_complete_price_upper_exact'])
    assert upper + aq(beta + Fraction(old['final_scalar_reserve_exact'])) < aq(price)
    theta = Fraction(item['theta_exact'])
    assert 0 < theta < 1
    assert theta * nlo - price == Fraction(item['coefficient_target_margin_exact']) > 0
    assert 1 - theta == Fraction(item['retained_energy_coefficient_exact'])
    assert rho + theta * s == Fraction(item['full_floor_exact'])
    print('PASS coupled finite guard', item['mu_exact'], 'retained', 1 - theta,
          'loss', item['full_floor_exact'])

barrier = json.loads((root / 'routes/complement-positive-rest-barrier-certificate.json').read_text())
assert arb(barrier['a_ball']) > 6 and arb(barrier['D_ball']) < 55000
assert arb(barrier['partial_midpoint_ball_with_a_lower_6']) > 1875
r = Fraction(3, 88)
delta = Fraction(21, 1024000)
margin = aq(r / 10 - delta - Fraction(11, 10) * Fraction(3, 12500)
            - Fraction(32237, 10 ** 12)) - 2 * aq(delta * r).sqrt()
assert margin > aq(Fraction(7, 5000))
assert rho + Fraction(7, 5000) == Fraction('0.0018896181745')
print('PASS fixed-price scalar barrier and', len(manifest['files']), 'file fingerprints')

search = json.loads((root / 'certificates/complement_residual_enrichment_shift_search.json').read_text())
for filename, expected in search['input_sha256'].items():
    assert historical_digest(root / 'certificates' / filename) == expected
for filename, field in [
    ('certificates/complement_residual_enrichment_action.json', 'action_record_sha256'),
    ('certificates/complement_residual_enrichment_action_matrices.npz', 'action_matrix_sha256'),
    ('certificates/complement_residual_enrichment_action_executed_source.py', 'action_executed_source_sha256'),
    ('computations/complement_residual_enrichment_shift_search.py', 'source_sha256'),
    ('routes/complement-residual-enrichment.md', 'route_sha256'),
]:
    assert hashlib.sha256((root / filename).read_bytes()).hexdigest() == search[field]
assert len(search['shifts']) == 231
for row in search['shifts']:
    mu, beta, theta = [Fraction(row[k]) for k in ['mu_exact', 'beta_exact', 'theta_exact']]
    assert mu > rho and beta > 0
    assert Fraction(row['retained_energy_exact']) == 1 - theta
    assert Fraction(row['alpha_exact']) == rho + theta * (mu - rho)
    assert row['passes'] == (theta < 1)
    assert arb(row['price_ball']) < aq(theta * nlo)
passed = [x for x in search['shifts'] if x['passes']]
assert min(passed, key=lambda x: Fraction(x['alpha_exact'])) == search['best_tested_scalar']
assert max(passed, key=lambda x: Fraction(x['retained_energy_exact'])) == search['best_tested_energy']
assert Fraction(search['best_tested_scalar']['alpha_exact']) == Fraction(31324476155738682727829, 10**24)
for mu, beta, retention, ceiling in [('1/25', '3/2', '0.132584532412', '0.034762'),
                                     ('1/20', '1277/512', '0.212458845389', '0.039482')]:
    row = next(x for x in passed if x['mu_exact'] == mu and x['beta_exact'] == beta)
    assert Fraction(row['retained_energy_exact']) == Fraction(retention)
    assert Fraction(row['alpha_exact']) < Fraction(ceiling)
signed_price = (Fraction(71, 100) * (Fraction(1532843, 10**11)
                + Fraction(1797275, 10**14) + Fraction(1, 10**88))
                + Fraction(1, 10**157) + Fraction(32237, 10**12) * Fraction(2159290, 10**6))
assert signed_price < Fraction(11, 10**6)
print('PASS all 231 enriched endpoint/hash checks, theorem transfer and signed-action rational endpoint')

# New portable finite signed-center checks; upstream action remains accepted.
signed = json.loads((root / 'certificates/complement_signed_cached_gram_certificate.json').read_text())
for filename, expected in signed['input_sha256'].items():
    assert historical_digest(root / 'certificates' / filename) == expected
assert hashlib.sha256((root / 'computations/complement_signed_cached_gram.py').read_bytes()).hexdigest() == signed['source_sha256']
assert hashlib.sha256((root / 'certificates/complement_signed_cached_gram_matrices.npz').read_bytes()).hexdigest() == signed['matrix_sha256']
sm = np.load(root / 'certificates/complement_signed_cached_gram_matrices.npz')
se = {k: arb(v) for k, v in signed['saved_center_operator_error_balls'].items()}

def spectrum(array):
    target = symmetric(array)
    values, vectors = np.linalg.eigh((array + array.T) / 2)
    av = exact(vectors)
    identity = arb_mat(np.eye(len(values)).astype(int).tolist())
    diagonal = arb_mat([[arb(float(values[i])) if i == j else arb(0)
                        for j in range(len(values))] for i in range(len(values))])
    defect = frobenius(av.transpose() * av - identity)
    residual = frobenius(target * av - av * diagonal)
    assert defect < aq(Fraction(1, 1000000))
    correction = residual / (1 - defect).sqrt() + 2 * arb(float(max(abs(values)))) * (1 + defect).sqrt() * ((1 - defect).rsqrt() - 1)
    return arb(float(values[0])) - correction, arb(float(values[-1])) + correction

w = exact(sm['W'])
wd = frobenius(w.transpose() * w - arb_mat(np.eye(5).astype(int).tolist()))
wn2 = 1 + wd
yn = aq(Fraction('2.159290'))
e = aq(Fraction(7, 12500))
assert arb(signed['whole_signed_action_error_ball']) < e
physical_lower = spectrum(sm['restrictedGcenter'])[0] - se['restrictedGcenter']
signed_lower = spectrum(sm['restrictedDcenter'])[0] - se['restrictedDcenter'] - e * yn * wn2
assert physical_lower > aq(Fraction(1, 2))
assert signed_lower > aq(Fraction(1, 1000))
assert yn**2 * wn2 < 5
reslo, reshi = spectrum(sm['restrictedResidualcenter'])
raw, gram, selector = [exact(sm[k]) for k in ['restrictedRawCross', 'restrictedGcenter', 'L']]
selector_norm = frobenius(selector)
selector_error = (frobenius(raw - gram * selector) + se['restrictedRawCross']
                  + se['restrictedGcenter'] * selector_norm) / physical_lower
# The first spectral endpoint is the least eigenvalue; certify a separate
# binary Rayleigh witness for the largest residual eigenvalue below.
rv = np.linalg.eigh((sm['restrictedResidualcenter'] + sm['restrictedResidualcenter'].T) / 2)[1][:, -1:]
ray = exact(rv)
rayleigh = (ray.transpose() * symmetric(sm['restrictedResidualcenter']) * ray)[0, 0] / (ray.transpose() * ray)[0, 0]
coupling_lower = (rayleigh - se['restrictedResidualcenter']).sqrt() - yn * wn2.sqrt() * selector_error - e * wn2.sqrt()
assert coupling_lower > aq(Fraction(11, 50))
coupling_squared_upper = aq(Fraction(101, 100)) * (reshi + se['restrictedResidualcenter']) + 101 * e**2 * wn2
assert coupling_squared_upper < aq(Fraction('0.062842088203'))

st = json.loads((root / 'certificates/complement_enriched_signed_tail_capture_certificate.json').read_text())
for filename, expected in st['input_sha256'].items():
    assert historical_digest(root / 'certificates' / filename) == expected
assert hashlib.sha256((root / 'certificates/complement_enriched_signed_tail_capture_executed_source.py').read_bytes()).hexdigest() == st['source_sha256']
assert arb(st['verified_lower_witness_ball']) > aq(Fraction(7, 50))
print('PASS portable 5D physical/signed spectral guards, whole outside coupling, ST2 recorded barrier and lineage')

# Reuse the immutable finite assembly only; stop before the parameter-search loop.
source = root / 'computations/complement_residual_enrichment_shift_search.py'
tree = ast.parse(source.read_text())
setup = []
for node in tree.body:
    # Preserve all math; only route historical input hashes through the declared alias.
    if isinstance(node, ast.Assign) and len(node.targets) == 1 and isinstance(node.targets[0], ast.Name) and node.targets[0].id == 'sha':
        assert isinstance(node.value, ast.Lambda)
        continue
    setup.append(node)
    if isinstance(node, ast.FunctionDef) and node.name == 'evaluate':
        break
assert isinstance(setup[-1], ast.FunctionDef) and setup[-1].name == 'evaluate'
namespace = {'__file__': str(source), '__name__': '__certificate_finite_replay__', 'sha': historical_digest}
exec(compile(ast.Module(body=setup, type_ignores=[]), str(source), 'exec'), namespace)
published = [search['best_tested_scalar'], search['best_tested_energy']]
published += [next(row for row in passed if row['mu_exact'] == mu and row['beta_exact'] == beta)
              for mu, beta in [('1/25', '3/2'), ('1/20', '1277/512')]]
for row in published:
    replay = namespace['evaluate'](Fraction(row['mu_exact']), Fraction(row['beta_exact']))
    # Reconstructed price, not the stored price_ball, must meet the published theta.
    assert arb(replay['price_ball']) < aq(Fraction(row['theta_exact']) * nlo)
    print('PASS reconstructed published enriched row', row['mu_exact'], row['beta_exact'])
runpy.run_path(str(root / 'computations/complement_enriched_signed_tail_capture_check.py'), run_name='__certificate_st2__')
print('PASS bounded certificate verification; upstream cached integral/zero/DFT/Gram error claims remain premises')
