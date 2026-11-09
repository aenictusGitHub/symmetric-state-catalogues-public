"""Degree-eight lower-bound certificate via invariant interpolation and interval bounds.

Run with Sage's Python. The verifier uses a rational Gram seed and proves an exact positive correction.
"""
import sys
import gzip
import json
from collections import defaultdict
from math import comb
from fractions import Fraction
from pathlib import Path

from sage.all import GF, QQ, AA, RealIntervalField, NumberField, PolynomialRing, matrix, vector

root = Path(__file__).resolve().parent
t = 5
with gzip.open(root/'sharp_input_t5.json.gz','rt',encoding='utf-8') as stream:
    source = json.load(stream)
certificate = json.loads((root/'sharp_certificate_t5.json').read_text())
points = [source['points'][i] for i in certificate['rows']]
interval = RealIntervalField(180)

def invert(A):
    """Enclose the inverse of every matrix represented by A."""
    n = A.nrows()
    augmented = [list(A.row(i))+[interval(int(i==j)) for j in range(n)] for i in range(n)]
    for i in range(n):
        pivot = max(range(i,n),key=lambda j:abs(augmented[j][i].center()))
        augmented[i],augmented[pivot]=augmented[pivot],augmented[i]
        divisor = augmented[i][i]
        assert 0 not in divisor
        augmented[i] = [v/divisor for v in augmented[i]]
        for j in range(n):
            if i!=j:
                factor=augmented[j][i]
                augmented[j]=[u-factor*v for u,v in zip(augmented[j],augmented[i])]
    return matrix(interval,[row[n:] for row in augmented])

def bound_norm(A):
    return max(sum(abs(x) for x in row).upper() for row in A.rows())

# Isolate the stationary state in the real even-Dicke chart z0=1.
R=PolynomialRing(QQ,5,'q')
q=R.gens()
z={0:R.one(), **{2*k+2:v for k,v in enumerate(q)}}
norm=sum(comb(11,k)*v*v for k,v in z.items())
purity=R.zero()
for i in range(6):
    for h in range(6):
        inner=sum(comb(6,j)*z.get(i+j,0)*z.get(h+j,0) for j in range(7))
        purity+=comb(5,i)*comb(5,h)*inner**2
gradient=[norm*purity.derivative(v)-2*purity*norm.derivative(v) for v in q]
jacobian=[[g.derivative(v) for v in q] for g in gradient]

def evaluate_poly(poly,box):
    powers=[[v**k for k in range(7)] for v in box]
    return sum(interval(c)*prod(powers[i][e] for i,e in enumerate(m)) for m,c in poly.dict().items())

from math import prod
center=[QQ(Fraction(row[0])) for row in json.loads((root/'witness_t5.json').read_text())['z'][2::2]]
radius=QQ(1)/10**36
box=[interval(c-radius,c+radius) for c in center]
center_box=list(map(interval,center))
J0=matrix(interval,[[evaluate_poly(poly,center_box) for poly in row] for row in jacobian])
approx_inverse=matrix(interval,[[interval(QQ(v.center())) for v in row] for row in invert(J0).rows()])
Jbox=matrix(interval,[[evaluate_poly(poly,box) for poly in row] for row in jacobian])
identity=matrix(interval,5,5,lambda i,j:int(i==j))
contraction=bound_norm(identity-approx_inverse*Jbox)
center_step=approx_inverse*vector(interval,[evaluate_poly(poly,center_box) for poly in gradient])
assert (interval(max(abs(v).upper() for v in center_step))+interval(contraction)*interval(radius)).upper()<radius
beta=evaluate_poly(purity,box)/evaluate_poly(norm,box)**2
witness={0:interval(1), **{2*k+2:v for k,v in enumerate(box)}}
print('Stationary state isolated by interval contraction.',flush=True)

def jet(poly):
    result=[interval(0) for _ in range(6)]
    for (u,v),c in poly.items():
        if any(k not in witness for k in u+v):continue
        indices=u+v
        result[0]+=c*prod(witness[k] for k in indices)
        for k in set(indices)-{0}:
            rest=list(indices)
            rest.remove(k)
            result[k//2]+=c*indices.count(k)*prod(witness[j] for j in rest)
    return result

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
jets=[]
for block in source['holomorphic']:
    J = block['J']
    orbits = prepare(block)
    check_orbits(orbits,4,J)
    dimension += len(orbits)*(len(orbits)+1)//2
    values=[[jet(poly) for poly in orbit] for orbit in orbits]
    for i in range(len(orbits)):
        for j in range(i,len(orbits)):
            column=[interval(0) for _ in range(6)]
            for level,(u,v) in enumerate(zip(values[i],values[j])):
                weight=QQ.one()/comb(J,level)
                column[0]+=weight*u[0]*v[0]
                for k in range(1,6):
                    column[k]+=weight*(u[0]*v[k]+u[k]*v[0])
            jets.append(column)
    for orbit in orbits:
        for level,poly in enumerate(orbit):
            weight_basis[J-2*level].append(poly)
assert sum(map(len,weight_basis.values())) == comb(15,4)
for polys in weight_basis.values():
    mons = sorted(set().union(*polys))
    assert matrix(F,[[F(poly.get(mon,0)) for mon in mons] for poly in polys]).rank()==len(polys)
assert dimension==177
# Verify independence of value plus the five derivatives at the root.
jets=matrix(interval,jets).transpose()
import numpy as np
import scipy.linalg as la
numeric=np.array([[float(v.center()) for v in row] for row in jets.rows()])
pivots=list(map(int,la.qr(numeric,pivoting=True,mode='economic')[2][:6]))
invert(jets.matrix_from_columns(pivots))
print('Six independent invariant conditions verified.',flush=True)
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
assert len(selected)==len(points)==dimension-6
B=[[interval(0) for _ in points] for _ in points]
seed=[QQ(x)/certificate['scale'] for x in certificate['seed']]
total=[QQ.zero()]*len(points)
offset=0
grams=[]
for block in source['blocks']:
    J=block['J']
    orbits=prepare(block)
    check_orbits(orbits,block['charge'],J)
    levels=[l for l in range(J+1) if (J-2*l-block['charge']*11)%4==0]
    piv=block['pivots']
    free=[i for i in range(len(orbits)) if i not in piv]
    assert len(piv)==len(levels)
    kernel=matrix(interval,len(free),len(orbits))
    for i,j in enumerate(free):kernel[i,j]=1
    if levels:
        E=matrix(interval,[[jet(o[l])[0] for o in orbits] for l in levels])
        solution=-invert(E.matrix_from_columns(piv))*E.matrix_from_columns(free)
        for i,j in enumerate(piv):
            for k in range(len(free)):kernel[k,j]=solution[i,k]
    transform=matrix(QQ,block['transform']).change_ring(interval)*kernel
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
correction=max(abs(sum(C[i][j]*residual[j] for j in range(n))).upper() for i in range(n))
delta=(interval(correction)/(1-interval(error))).upper()
for gram in grams:
    check_positive(gram,interval(delta))
print('Inverse residual bound:',error)
print('Gram-entry correction bound:',delta)
print('Certified minimum enclosure:',beta)
print('Sharp degree-eight certificate verified for N=11, t=5.')
