#!/usr/bin/env python3
"""Require complete #print axioms coverage and only standard logical axioms."""
from pathlib import Path
import re
import sys

ALLOWED = {'propext', 'Classical.choice', 'Quot.sound'}
ROOT = Path(__file__).resolve().parents[1]


def check(output, expected):
    seen = {}
    pattern = r"'([^']+)' (?:depends on axioms:\s*\[([^\]]*)\]|does not depend on any axioms)"
    for name, body in re.findall(pattern, output):
        if name in seen:
            raise ValueError(f'Duplicate axiom report: {name}')
        axioms = {a.strip() for a in body.split(',') if a.strip()}
        if axioms - ALLOWED:
            raise ValueError(f'Unsupported axioms for {name}: {sorted(axioms - ALLOWED)}')
        seen[name] = axioms
    if not expected or set(seen) != set(expected):
        raise ValueError(f'Axiom coverage mismatch: missing {set(expected) - set(seen)}, extra {set(seen) - set(expected)}')
    return len(seen)


def self_test():
    assert check("'A' depends on axioms: [propext, Classical.choice, Quot.sound]", ['A']) == 1
    assert check("'A' does not depend on any axioms", ['A']) == 1
    for bad in ["'A' depends on axioms: [sorryAx]", "'A' depends on axioms: [Lean.ofReduceBool]", "'A' depends on axioms: [customClaim]", '', "'A' depends on axioms: []\n'A' depends on axioms: []"]:
        try:
            check(bad, ['A'])
        except ValueError:
            continue
        raise AssertionError(f'Accepted invalid log: {bad}')
    print('PASS: axiom parser rejection and coverage checks')


if __name__ == '__main__':
    if sys.argv[1:] == ['--self-test']:
        self_test()
    else:
        if len(sys.argv) != 2:
            raise SystemExit('Usage: python3 scripts/check_axioms.py LOG | --self-test')
        expected = re.findall(r'^#print axioms ([\w.]+)$', (ROOT / 'lean/Axioms.lean').read_text(), re.M)
        count = check(Path(sys.argv[1]).read_text(), expected)
        print(f'PASS: {count} complete transitive axiom reports; standard logical axioms only')
