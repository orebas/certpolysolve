#!/usr/bin/env python3
"""Independent exact product/distinct-root checks for prepared isolation inputs."""
from fractions import Fraction as Q
from functools import reduce
from math import gcd, lcm
from pathlib import Path
import hashlib
import json


def product(roots):
    coefficients = [(Q(1), Q(0))]
    for real, imag in roots:
        next_coefficients = [(Q(0), Q(0)) for _ in range(len(coefficients) + 1)]
        for j, (x, y) in enumerate(coefficients):
            a, b = next_coefficients[j]
            next_coefficients[j] = (a - x * real + y * imag, b - x * imag - y * real)
            a, b = next_coefficients[j + 1]
            next_coefficients[j + 1] = (a + x, b + y)
        coefficients = next_coefficients
    denominator = lcm(*(value.denominator for pair in coefficients for value in pair))
    integers = [(int(x * denominator), int(y * denominator)) for x, y in coefficients]
    content = reduce(gcd, (abs(value) for pair in integers for value in pair))
    return [(x // content, y // content) for x, y in integers]


def main():
    directory = Path(__file__).resolve().parent / 'isolation-fixtures'
    manifest = json.loads((directory / 'manifest.json').read_text())
    for case in manifest['cases']:
        path = directory / case['input']
        body = path.read_bytes()
        assert hashlib.sha256(body).hexdigest() == case['sha256'], case['id']
        lines = body.decode().splitlines()
        assert lines[0] == 'PSCASE 1'
        assert int(lines[1]) == case['degree'] + 1 == len(lines) - 2
        roots = [tuple(map(Q, pair)) for pair in case['roots']]
        assert len(roots) == case['degree'] == len(set(roots)), case['id']
        rows = [tuple(map(int, line.split())) for line in lines[2:]]
        assert all(len(row) == 4 and row[2:] == (0, 0) for row in rows)
        assert [(row[0], row[1]) for row in rows] == product(roots), case['id']
        print(case['id'], 'exact product and distinct roots verified')


if __name__ == '__main__':
    main()
