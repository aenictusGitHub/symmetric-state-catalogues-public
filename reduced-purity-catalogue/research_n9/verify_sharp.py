"""Verify the sharp four-body purity certificate with SymPy and JSON data."""
import json
from collections import defaultdict
from math import comb
from pathlib import Path

import sympy as s
from sympy.polys.matrices import DomainMatrix

root = Path(__file__).resolve().parent
r = s.Symbol("r")
p = 6174*r**3 + 252*r**2 - 173*r - 4
alpha = s.CRootOf(p, 2)
field = s.QQ.algebraic_field(alpha)
a = field.unit
lo, hi = s.Rational(0), s.Rational(1)
assert p.subs(r, lo) < 0 < p.subs(r, hi)
# Descartes' rule isolates the unique positive root.
for _ in range(160):
    mid = (lo+hi)/2
    if p.subs(r, mid) > 0:
        hi = mid
    else:
        lo = mid

def decode(coefficients):
    value = field.zero
    for coefficient in coefficients:
        value = value*a + field.convert(s.Rational(coefficient))
    return value

def positive(value):
    lower = upper = s.Rational(0)
    for coefficient in value.to_list():
        corners = (lower*lo, lower*hi, upper*lo, upper*hi)
        c = s.QQ.to_sympy(coefficient)
        lower, upper = min(corners)+c, max(corners)+c
    return lower > 0

def multiply(left, right):
    result = defaultdict(lambda: field.zero)
    for (u,v), c in left.items():
        for (w,x), d in right.items():
            result[tuple(sorted(u+w)), tuple(sorted(v+x))] += c*d
    return {key:value for key,value in result.items() if value}

norm = {((k,),(k,)):field.convert(comb(9,k)) for k in range(10)}
target = defaultdict(lambda: field.zero)
for i in range(5):
    for h in range(5):
        for j in range(6):
            for l in range(6):
                key = tuple(sorted((i+j,h+l))), tuple(sorted((h+j,i+l)))
                target[key] += field.convert(comb(4,i)*comb(4,h)*comb(5,j)*comb(5,l))
num = 777924*r**4+12348*r**2-1680*r+1159
den = 3528*(21*r**2+1)**2
beta = field.from_sympy((num/den).subs(r,alpha))
for key,value in multiply(norm,norm).items():
    target[key] -= beta*value
lhs = multiply(target,norm)

basis = {(block["charge"],block["twice_spin"]):block["orbits"]
         for block in json.loads((root/"cubic_covariants.json").read_text())}
rhs = defaultdict(lambda:field.zero)
certificate = json.loads((root/"sharp_certificate_t4.json").read_text())
for block in certificate:
    J = block["twice_spin"]
    orbits = basis[block["charge"],J]
    kernel = [[decode(v) for v in row] for row in block["kernel"]]
    gram = [[decode(v) for v in row] for row in block["gram"]]
    size = len(kernel)
    assert len(gram) == size and all(len(row)==size for row in gram)
    assert all(len(row)==len(orbits) for row in kernel)
    assert all(gram[i][j]==gram[j][i] for i in range(size) for j in range(size))
    for k in range(1,size+1):
        assert positive(DomainMatrix([row[:k] for row in gram[:k]],(k,k),field).det())
    for level in range(J+1):
        polynomials = []
        for row in kernel:
            f = defaultdict(lambda:field.zero)
            for weight,orbit in zip(row,orbits):
                if weight:
                    for (u,v),coefficient in orbit[level]:
                        f[tuple(u),tuple(v)] += weight*field.convert(s.Rational(coefficient))
            polynomials.append({key:value for key,value in f.items() if value})
        for i,f in enumerate(polynomials):
            for j,g in enumerate(polynomials):
                weight = gram[i][j]/field.convert(comb(J,level))
                for key,value in multiply(f,{(v,u):c for (u,v),c in g.items()}).items():
                    rhs[key] += weight*value
assert lhs == {key:value for key,value in rhs.items() if value}

# Independently evaluate the candidate marginal to check attainment.
c = [s.sqrt(21)*r,0,0,1,0,0,-1,0,0,-s.sqrt(21)*r]
matrix = s.Matrix(5,6,lambda i,j:c[i+j]*s.sqrt(s.Rational(comb(4,i)*comb(5,j),comb(9,i+j))))
assert s.cancel(s.trace((matrix*matrix.T)**2)/(2+42*r**2)**2-num/den)==0
print(f"Sharp identity and positivity of {len(certificate)} blocks verified exactly. Attainment verified.")
