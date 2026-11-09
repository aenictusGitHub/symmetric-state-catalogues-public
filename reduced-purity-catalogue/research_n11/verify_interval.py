"""Degree-eight lower-bound certificate via invariant interpolation and interval bounds.

Run with Sage's Python. The verifier uses a rational Gram seed and proves an exact positive correction.
"""
import sys
import gzip
import json
from collections import defaultdict
from math import comb
from pathlib import Path

from sage.all import GF, QQ, AA, RealIntervalField, NumberField, PolynomialRing, matrix, vector

root = Path(__file__).resolve().parent
t = int(sys.argv[1])
assert t in (4,5)
with gzip.open(root/f'quartic_input_t{t}.json.gz','rt',encoding='utf-8') as stream:
    source = json.load(stream)
certificate = json.loads((root/f"interval_certificate_t{t}.json").read_text())
beta = QQ(certificate['bound'])
points = [source['points'][i] for i in certificate['rows']]

def multiply(u, v):
    x,y = u
    z,w = v
    return x*z-y*w, x*w+y*z

def prepare(record):
    J = record['J']
    polys = [{(tuple(u),tuple(v)): QQ(c) for (u,v),c in poly} for poly in record['orbits']]
    return [polys[i:i+J+1] for i in range(0,len(polys),J+1)]

def derivative(poly, up):
    result = defaultdict(lambda: QQ.zero())
    for (u,v), c in poly.items():
        for side, indices in enumerate((u,v)):
            step = -1 if up != bool(side) else 1
            for k in set(indices):
                if 0 <= k+step <= 11:
                    shifted = list(indices)
                    shifted.remove(k)
                    shifted.append(k+step)
                    shifted.sort()
                    factor = (k if step == -1 else 11-k)*indices.count(k)*(-1 if side else 1)
                    key = (tuple(shifted),v) if side == 0 else (u,tuple(shifted))
                    result[key] += factor*c
    return {key:value for key,value in result.items() if value}

def check_orbits(orbits, charge, J):
    for orbit in orbits:
        assert all(len(u)==(4+charge)//2 and len(v)==(4-charge)//2 and
                   11*charge-2*(sum(u)-sum(v))==J for u,v in orbit[0])
        assert not derivative(orbit[0], True)
        for level in range(1,J+1):
            assert derivative(orbit[level-1],False) == {
                key:level*c for key,c in orbit[level].items()}
        assert not derivative(orbit[-1],False)

# Precompute integer Gaussian monomial values once per point.
cache = {}
def monomial(u,v):
    key = u,v
    if key not in cache:
        data = []
        for z in points:
            q = (1,0)
            for i in u:
                q = multiply(q,z[i])
            for i in v:
                q = multiply(q,(z[i][0],-z[i][1]))
            data.append(q)
        cache[key] = list(zip(*data))
    return cache[key]

def evaluate(orbits,J):
    for level in range(J+1):
        mons = sorted(set().union(*(o[level] for o in orbits)))
        C = matrix(QQ,[[o[level].get(m,0) for m in mons] for o in orbits])
        data = [monomial(*m) for m in mons]
        yield tuple(C*matrix(QQ,[d[k] for d in data]) for k in (0,1))

# The real, conjugation-even invariant subspace has symmetric matrices
# on the multiplicity spaces of the holomorphic quartic representation.
F = GF(1000000009)
weight_basis = defaultdict(list)
dimension = 0
for block in source['holomorphic']:
    J = block['J']
    orbits = prepare(block)
    check_orbits(orbits,4,J)
    dimension += len(orbits)*(len(orbits)+1)//2
    for orbit in orbits:
        for level,poly in enumerate(orbit):
            weight_basis[J-2*level].append(poly)
assert sum(map(len,weight_basis.values())) == comb(15,4)
for polys in weight_basis.values():
    mons = sorted(set().union(*polys))
    assert matrix(F,[[F(poly.get(mon,0)) for mon in mons] for poly in polys]).rank()==len(polys)
assert dimension==177
def target(z):
    norm = sum(comb(11,k)*(u*u+v*v) for k,(u,v) in enumerate(z))
    purity = 0
    for i in range(t+1):
        for h in range(t+1):
            x=y=0
            for j in range(12-t):
                u,v=multiply(z[i+j],(z[h+j][0],-z[h+j][1]))
                x+=comb(11-t,j)*u
                y+=comb(11-t,j)*v
            purity+=comb(t,i)*comb(t,h)*(x*x+y*y)
    return norm**2*(purity-beta*norm**2)

interval = RealIntervalField(160)
def check_positive(gram, perturbation):
    size = gram.nrows()
    A = [[interval(gram[i,j]) for j in range(size)] for i in range(size)]
    for i in range(size):
        A[i][i] -= size*perturbation
    for k in range(size):
        pivot=A[k][k]
        assert pivot.lower()>0
        for i in range(k+1,size):
            for j in range(k+1,size):
                A[i][j]-=A[i][k]*A[k][j]/pivot

selected={i:k for k,i in enumerate(certificate['cols'])}
assert len(selected)==len(points)==dimension
B=[[interval(0) for _ in points] for _ in points]
seed=[QQ(x)/certificate['scale'] for x in certificate['seed']]
total=[interval(0)]*len(points)
offset=0
grams=[]
for block in source['blocks']:
    J=block['J']
    orbits=prepare(block)
    check_orbits(orbits,block['charge'],J)
    transform=matrix(QQ,block['transform']).change_ring(interval)
    size=transform.nrows()
    gram=matrix(QQ,size)
    pivots=[]
    for i in range(size):
        for j in range(i,size):
            gram[i,j]=gram[j,i]=seed[offset]
            if offset in selected:
                pivots.append((selected[offset],i,j))
            offset+=1
    grams.append(gram)
    H=transform.transpose()*gram.change_ring(interval)*transform
    sums=[interval(0)]*len(points)
    for level,(real,imag) in enumerate(evaluate(orbits,J)):
        weight=interval(QQ.one()/comb(J,level))
        for Z in (real.change_ring(interval),imag.change_ring(interval)):
            product=H*Z
            sums=[old+weight*u.dot_product(v)
                  for old,u,v in zip(sums,Z.columns(),product.columns())]
            if pivots:
                X=transform*Z
                for k,i,j in pivots:
                    factor=weight*(1 if i==j else 2)
                    for m in range(len(points)):
                        B[m][k]+=factor*X[i,m]*X[j,m]
    total=[old+val for old,val in zip(total,sums)]
    print('Checked block',block['charge'],J,size,flush=True)
assert offset==len(seed)
scales=[QQ(sum(comb(11,k)*(u*u+v*v) for k,(u,v) in enumerate(z)))**4 for z in points]
for i in range(len(points)):
    B[i]=[v/interval(scales[i]) for v in B[i]]
residual=[interval((target(z)-val)/scale) for z,val,scale in zip(points,total,scales)]
C=[[interval(QQ(x)/certificate['scale']) for x in row] for row in certificate['inverse']]
assert all(len(row)==len(points) for row in C) and len(C)==len(points)
n=len(points)
error=max(sum(abs(interval(int(i==j))-sum(C[i][k]*B[k][j] for k in range(n)))
              for j in range(n)).upper() for i in range(n))
assert error<1
# Recenter the rational Gram seed using the approximate inverse. The
# remaining correction is then much smaller than the original residual.
updates = [QQ(sum(C[i][j]*residual[j] for j in range(n)).center()) for i in range(n)]
offset = 0
for gram in grams:
    for i in range(gram.nrows()):
        for j in range(i, gram.ncols()):
            if offset in selected:
                gram[i,j] += updates[selected[offset]]
                gram[j,i] = gram[i,j]
            offset += 1
residual = [residual[i]-sum(B[i][j]*interval(updates[j]) for j in range(n)) for i in range(n)]
correction=max(abs(sum(C[i][j]*residual[j] for j in range(n))).upper() for i in range(n))
delta=(interval(correction)/(1-interval(error))).upper()
for gram in grams:
    check_positive(gram,interval(delta))
print('Inverse residual bound:',error)
print('Gram-entry correction bound:',delta)
print('Degree-eight lower-bound certificate verified by exact arithmetic and rigorous intervals.')
