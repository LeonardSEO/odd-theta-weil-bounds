"""Small scalar certificate for the actual positive midpoint comparison barrier.

Only two new jet integrals and a nonnegative partial midpoint sum are evaluated.
The complete theta-tail allowances are proved in the companion route note.
"""
from fractions import Fraction
from pathlib import Path
import hashlib
import json
import math
from flint import arb, acb, ctx

ctx.prec = 192
root = Path(__file__).resolve().parents[1]
cells_path = root / 'certificates/complement_signed_symbol_norm_cells.json'
norms_path = root / 'certificates/complement_physical_norms.json'
cells = json.loads(cells_path.read_text())
norms = json.loads(norms_path.read_text())
assert hashlib.sha256(cells_path.read_bytes()).hexdigest() == 'cf4e6b03f62d5f787f41fb104f7912a7fcadd6333da3032ac086b36f5020df30'
Q = arb(norms['Q_ball'])
assert Q > 11 and Q < 12

def rational(value):
    value = Fraction(value)
    return arb(value.numerator) / value.denominator

def jet(x):
    p, d, h = acb(0), acb(0), acb(0)
    for n in range(1, 9):
        y = arb.pi() * n * n * (2*x).exp()
        atom = (x/2-y).exp()
        p += atom*(4*y*y-6*y)
        d += atom*(-8*y**3+30*y*y-15*y)
        h += atom*(-32*y**5+360*y**4-1056*y**3+810*y*y-90*y)
    return p, d, h

def cross_integrand(x, analytic):
    p, d, h = jet(x)
    return d*h*(x/2).cosh()/p

def square_integrand(x, analytic):
    p, d, h = jet(x)
    return h*h*(x/2).cosh()/p

# Full lattice and physical-tail estimates, with exact scalar comparisons.
assert rational(Fraction(90,81)+Fraction(810,27)+Fraction(1056,9)+120+32) < 301
assert rational(Fraction(1024,8000)/(1-Fraction(58,20**5))) < 1
assert (arb(4)).exp() > 54 and arb(4).exp() < 55
assert (-arb(240)).exp() < rational(Fraction(1,2**340))
assert 3*(-arb(57)).exp() < arb(1)/2
assert arb(301)*220**3*9**10 < arb(2)**65
assert 20000*55**3 < 2**32
assert 12*arb(2)**85*arb(2)**-270 < arb(2)**-170
tail_poly = sum(arb(math.factorial(8)//math.factorial(j))*54**j/arb(3)**(9-j)
                for j in range(9))
assert arb(10)**11*(-arb(162)).exp()*tail_poly < arb(2)**-140
assert arb(2)**-140+arb(2)**-170 < arb(2)**-130
error = arb(2)**-130
cross = 2*acb.integral(cross_integrand, 0, 2,
                      abs_tol=arb(2)**-80, rel_tol=arb(2)**-80,
                      eval_limit=100000)
square = 2*acb.integral(square_integrand, 0, 2,
                       abs_tol=arb(2)**-80, rel_tol=arb(2)**-80,
                       eval_limit=100000)
assert cross.is_finite() and square.is_finite()
assert cross.imag.contains(0) and square.imag.contains(0)
Hq = cross.real + arb(0, error)
Hh = square.real + arb(0, error)
assert Hq > 72 and Hq < 77 and Hh > 0 and Hh < 55000
a = Hq/Q
D = Hh-Hq*Hq/Q
assert a > 6 and a < 7 and D > 0 and D < 55000

partial = arb(0)
count = 0
for cell in cells['kernel_cells']:
    t = arb(2*cell['k']+1)/128
    if t > 20:
        continue
    # Same exact upward binary64 mass as B in MP/PF/ER.
    mu = math.nextafter(int(cell['measure_upper_dyadic_numerator'])/2**80,
                        math.inf)
    num, den = mu.as_integer_ratio()
    s = acb(arb(1)/2, t)
    xi = s*(s-1)/2*arb.pi()**(-s/2)*(s/2).gamma()*s.zeta()
    assert xi.is_finite() and xi.imag.contains(0)
    partial += arb(num)/den*t*t*(t*t+arb(25)/4)**2*xi.real*xi.real
    count += 1
assert count == 635 and partial > 1875
assert Fraction(1875,55000) > Fraction(3,100)

delta0 = Fraction(21,1024000)
epsilon = Fraction(3,12500)
t8 = Fraction(32237,10**12)
margin = rational(Fraction(3,1000)-delta0-Fraction(11,10)*epsilon-t8) - 2*rational(delta0*Fraction(3,100)).sqrt()
assert margin > rational(Fraction(1,1000))
# The stronger endpoint follows from r > 3/88, without using a float locator.
strong_margin = rational(Fraction(3,880)-delta0-Fraction(11,10)*epsilon-t8) - 2*rational(delta0*Fraction(3,88)).sqrt()
assert strong_margin > rational(Fraction(7,5000))

record = {
    'status': 'PASS actual score-projected first equality-jet barrier',
    'precision_bits': ctx.prec,
    'Hq_ball': str(Hq), 'Hh_ball': str(Hh),
    'a_ball': str(a), 'D_ball': str(D),
    'complete_each_integral_allowance': '1/2^130',
    'new_norm_method': '2 acb.integral_0^2; meromorphic callbacks; abs_tol=rel_tol=2^-80',
    'partial_midpoint_cells': count,
    'partial_midpoint_ball_with_a_lower_6': str(partial),
    'r_strict_lower': '3/88',
    'general_margin_at_r_3_100': str(margin),
    'strong_margin_at_r_3_88': str(strong_margin),
    'strict_barrier': '(R56-E_Y)(k)>7/5000 ||k||^2',
    'cells_sha256': hashlib.sha256(cells_path.read_bytes()).hexdigest(),
    'existing_norm_cache_sha256': hashlib.sha256(norms_path.read_bytes()).hexdigest(),
    'script_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
    'owning_note': 'routes/complement-positive-rest-barrier.md',
}
out = Path(__file__).with_name('complement-positive-rest-barrier-certificate.json')
out.write_text(json.dumps(record, indent=2)+'\n')
print(json.dumps(record, indent=2))
