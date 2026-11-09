"""Sharp degree-eight certificate via invariant interpolation and interval bounds.

Run with Sage's Python. The verifier uses a rational Gram seed and proves an exact positive correction.
"""
import gzip
import json
from collections import defaultdict
from math import comb
from pathlib import Path

from sage.all import GF, QQ, AA, RealIntervalField, NumberField, PolynomialRing, matrix, vector

root = Path(__file__).resolve().parent
with gzip.open(root/'quartic_input.json.gz','rt',encoding='utf-8') as stream:
    source = json.load(stream)
certificate = json.loads((root/"interval_certificate_t4.json").read_text())
R = PolynomialRing(QQ, 'r')
r = R.gen()
p = 435*r**2 - 17
K = NumberField(p, 'a')
a = K.gen()
beta = QQ(113)/515
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
                if 0 <= k+step <= 8:
                    shifted = list(indices)
                    shifted.remove(k)
                    shifted.append(k+step)
                    shifted.sort()
                    factor = (k if step == -1 else 8-k)*indices.count(k)*(-1 if side else 1)
                    key = (tuple(shifted),v) if side == 0 else (u,tuple(shifted))
                    result[key] += factor*c
    return {key:value for key,value in result.items() if value}

def check_orbits(orbits, charge, J):
    for orbit in orbits:
        assert all(len(u)==(4+charge)//2 and len(v)==(4-charge)//2 and
                   8*charge-2*(sum(u)-sum(v))==J for u,v in orbit[0])
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
assert sum(map(len,weight_basis.values())) == comb(12,4)
for polys in weight_basis.values():
    mons = sorted(set().union(*polys))
    assert matrix(F,[[F(poly.get(mon,0)) for mon in mons] for poly in polys]).rank()==len(polys)
assert dimension==62
# On the real family (z0,z4,z8)=(1,a,-1), evaluate the numerator
# of the derivative of normalized purity when z4 varies.
witness = {0:K.one(),4:a,8:-K.one()}
def purity_slope(t):
    norm = 2+70*a*a
    norm_slope = 140*a
    purity = K.zero()
    slope = K.zero()
    for i in range(t+1):
        for h in range(t+1):
            v = dv = K.zero()
            for j in range(9-t):
                x,y = witness.get(i+j,K.zero()),witness.get(h+j,K.zero())
                v += comb(8-t,j)*x*y
                dv += comb(8-t,j)*(int(i+j==4)*y+x*int(h+j==4))
            weight = comb(t,i)*comb(t,h)
            purity += weight*v*v
            slope += 2*weight*v*dv
    return slope*norm-2*purity*norm_slope

# Value and derivative along the normalized family are independent,
# as witnessed by n^4 and n^2 F2.
assert purity_slope(2) != 0
assert purity_slope(4) == 0

def decode(entries):
    return sum((QQ(x)*a**i for i,x in enumerate(reversed(entries))), K.zero())

def target(z):
    norm = sum(comb(8,k)*(u*u+v*v) for k,(u,v) in enumerate(z))
    purity = 0
    for i in range(5):
        for h in range(5):
            x=y=0
            for j in range(5):
                u,v=multiply(z[i+j],(z[h+j][0],-z[h+j][1]))
                x+=comb(4,j)*u
                y+=comb(4,j)*v
            purity+=comb(4,i)*comb(4,h)*(x*x+y*y)
    return norm**2*(purity-beta*norm**2)

interval = RealIntervalField(160)
root_interval = interval(max(p.roots(AA,multiplicities=False)))
def enclose(value):
    bound = interval(0)
    for coefficient in reversed(K(value).list()):
        bound = bound*root_interval+interval(coefficient)
    return bound

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

def at(poly):
    return sum((c*prod(witness[k] for k in u+v) for (u,v),c in poly.items()
                if all(k in witness for k in u+v)),K.zero())

from math import prod
selected={i:k for k,i in enumerate(certificate['cols'])}
assert len(selected)==len(points)==dimension-2
B=[[interval(0) for _ in points] for _ in points]
seed=[QQ(x)/certificate['scale'] for x in certificate['seed']]
total=[K.zero()]*len(points)
offset=0
grams=[]
for block in source['blocks']:
    J=block['J']
    orbits=prepare(block)
    check_orbits(orbits,block['charge'],J)
    face=matrix(K,[[decode(v) for v in row] for row in block['kernel']])
    for level in range(J+1):
        assert not any(face*vector(K,[at(orbit[level]) for orbit in orbits]))
    transform=matrix(QQ,block['transform']).change_ring(K)*face
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
    H=transform.transpose()*gram.change_ring(K)*transform
    components=[matrix(QQ,[[(v.list()+[QQ.zero()]*3)[k] for v in row] for row in H.rows()])
                for k in range(3)]
    transform_parts=[matrix(QQ,[[(v.list()+[QQ.zero()]*3)[k] for v in row] for row in transform.rows()])
                     for k in range(3)] if pivots else []
    sums=[[QQ.zero()]*len(points) for _ in range(3)]
    for level,(real,imag) in enumerate(evaluate(orbits,J)):
        weight=QQ.one()/comb(J,level)
        for k,part in enumerate(components):
            for Z in (real,imag):
                product=part._multiply_pari(Z)
                sums[k]=[old+weight*u.dot_product(v)
                         for old,u,v in zip(sums[k],Z.columns(),product.columns())]
        if pivots:
            projected=[]
            for Z in (real,imag):
                values=[T._multiply_pari(Z) for T in transform_parts]
                projected.append([[sum(interval(values[k][i,j])*root_interval**k for k in range(3))
                                   for j in range(len(points))] for i in range(size)])
            X,Y=projected
            for k,i,j in pivots:
                factor=weight*(1 if i==j else 2)
                for m in range(len(points)):
                    B[m][k]+=factor*(X[i][m]*X[j][m]+Y[i][m]*Y[j][m])
    total=[old+sum((v*a**k for k,v in enumerate(values)),K.zero())
           for old,values in zip(total,zip(*sums))]
    print('Checked block',block['charge'],J,size,flush=True)
assert offset==len(seed)
scales=[QQ(sum(comb(8,k)*(u*u+v*v) for k,(u,v) in enumerate(z)))**4 for z in points]
for i in range(len(points)):
    B[i]=[v/interval(scales[i]) for v in B[i]]
residual=[enclose((target(z)-val)/scale) for z,val,scale in zip(points,total,scales)]
C=[[interval(QQ(x)/certificate['scale']) for x in row] for row in certificate['inverse']]
assert all(len(row)==len(points) for row in C) and len(C)==len(points)
n=len(points)
error=max(sum(abs(interval(int(i==j))-sum(C[i][k]*B[k][j] for k in range(n)))
              for j in range(n)).upper() for i in range(n))
assert error<1
correction=max(abs(sum(C[i][j]*residual[j] for j in range(n))).upper() for i in range(n))
delta=(interval(correction)/(1-interval(error))).upper()
for gram in grams:
    check_positive(gram,interval(delta))
assert target([(witness.get(i,K.zero()),0) for i in range(9)])==0
print('Inverse residual bound:',error)
print('Gram-entry correction bound:',delta)
print('Sharp degree-eight certificate verified by exact arithmetic and rigorous intervals.')
