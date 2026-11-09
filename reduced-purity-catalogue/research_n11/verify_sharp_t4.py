"""Rigorous degree-ten equality certificate for N=11, t=4.

Use Python with python-flint available and pass Sage's Python as the second argument. The input is produced by
export_quintic_t4.py from a numerical fit and contains only untrusted seeds.
"""
import gzip
import json
import sys
import subprocess
from pathlib import Path
from collections import Counter, defaultdict
from fractions import Fraction
from functools import lru_cache
from itertools import combinations_with_replacement
from math import comb

from flint import acb, acb_mat, arb, arb_mat, ctx, fmpq

from quintic_gradient import gradient_highest


ctx.prec = 180
ctx.threads = 4
with gzip.open(sys.argv[1], 'rt') as stream:
    source = json.load(stream)

def rational(x):
    x = Fraction(x)
    return arb(fmpq(x.numerator,x.denominator))

def enclose(bounds):
    lo, hi = map(rational,bounds)
    return arb((lo+hi)/2, ((hi-lo)/2).abs_upper())

stationary_script = Path(__file__).with_name('verify_stationary_t4.py')
root_log = subprocess.check_output([sys.argv[2],str(stationary_script),'--json'],text=True)
*messages, root_data = root_log.splitlines()
print(*messages,sep='\n',flush=True)
root_data = json.loads(root_data)
even, box, beta = (root_data[key] for key in ('even','box','beta'))
witness = [acb(int(k==0)) for k in range(12)]
for (k,part), value in zip(even,box):
    witness[k] += acb(0,enclose(value)) if part else acb(enclose(value))
beta = enclose(beta)

def ladder(p, raising):
    out = defaultdict(lambda: 0)
    for (a,b), coeff in p.items():
        for ket, inds in ((True,a),(False,b)):
            step = (-1 if ket else 1) if raising else (1 if ket else -1)
            for k in set(inds):
                if 0 <= k+step <= 11:
                    factor = (k if step==-1 else 11-k)*(1 if ket else -1)*inds.count(k)
                    moved = list(inds)
                    moved.remove(k)
                    moved.append(k+step)
                    shifted = tuple(sorted(moved))
                    out[(shifted,b) if ket else (a,shifted)] += factor*coeff
    return {m:c for m,c in out.items() if c != 0}

def orbit(p,J,charge):
    assert p and all(len(a)+len(b)==5 and len(a)-len(b)==charge and
                     charge*11-2*(sum(a)-sum(b))==J for a,b in p)
    assert not ladder(p, True)
    yield p
    for level in range(1,J+1):
        p = {m:c*Fraction(1,level) for m,c in ladder(p,False).items()}
        yield p
    assert not ladder(p, False)

def multiply(seq):
    result = acb(1)
    for term in seq:
        result *= term
    return result

@lru_cache(maxsize=None)
def witness_monomial(a,b):
    return multiply(witness[k] for k in a)*multiply(witness[k].conjugate() for k in b)

def evaluate_at_witness(p):
    return sum(rational(c)*witness_monomial(a,b) for (a,b),c in p.items())

def jet(p):
    vals = [evaluate_at_witness(p)]
    for k,part in even:
        u = v = acb(0)
        for (a,b),c in p.items():
            if k in a:
                i = a.index(k)
                u += rational(c*a.count(k))*witness_monomial(a[:i]+a[i+1:],b)
            if k in b:
                i = b.index(k)
                v += rational(c*b.count(k))*witness_monomial(a,b[:i]+b[i+1:])
        vals.append(acb(0,1)*(u-v) if part else u+v)
    return vals

# Count real, conjugation-even invariants in bidegree (5,5) from
# the weights of Sym^5(V_{11/2}).
weight_counts = Counter(map(sum,combinations_with_replacement(range(12),5)))
dimension = sum(comb(weight_counts[i]-weight_counts[i-1]+1,2) for i in range(28))
assert dimension==992
points = source['points']
size = len(points)
assert size==dimension-10
selected = {col:i for i,col in enumerate(source['cols'])}
assert len(selected)==size
scale = source['scale']
B = arb_mat(size,size)
total = [arb(0) for _ in points]
offset = 0
grams = []
samples = []
target = []
for point in points:
    z = [acb(x,y) for x,y in point]
    n = sum(comb(11,k)*(x*x+y*y) for k,(x,y) in enumerate(point))
    F = sum(comb(4,i)*comb(4,h)*abs(sum(comb(7,j)*z[i+j]*z[h+j].conjugate()
             for j in range(8)))**2 for i in range(5) for h in range(5))
    target.append(F/(n*n)-beta)
    samples.append([v/arb(n).sqrt() for v in z])

jets_certified = False
jets = []
for record in source['blocks']:
    charge,J = record['charge'],record['J']
    levels = range((J-charge*11)%4//2,J+1,2)
    if 'gradient_scale' in record:
        assert charge==3 and J in (31,33)
        # Check the highest-weight and ladder relations coefficientwise in beta.
        A = list(orbit(gradient_highest(J,Fraction(0)),J,charge))
        A1 = list(orbit(gradient_highest(J,Fraction(1)),J,charge))
        L = rational(record['gradient_scale'])
        orbits = [[{m:L*(rational(p.get(m,0))+beta*rational(q.get(m,0)-p.get(m,0)))
                    for m in p.keys()|q.keys()} for p,q in zip(A,A1)]]
        transform = acb_mat([[1]])
    else:
        raw = [{(tuple(a),tuple(b)):Fraction(c) for a,b,c in p} for p in record['highest']]
        orbits = [list(orbit(p,J,charge)) for p in raw]
        values = [[evaluate_at_witness(o[l]) for o in orbits] for l in levels]
        E = [[getattr(x,part) for x in row] for part in ('real','imag') for row in values]
        pivots = record['pivots']
        free = [j for j in range(len(raw)) if j not in pivots]
        assert len(pivots)==len(set(pivots))==len(E)
        invertible = arb_mat([[row[j] for j in pivots] for row in E])
        kernel_pivots = -invertible.solve(arb_mat([[row[j] for j in free] for row in E]))
        K = arb_mat(len(free),len(raw))
        for i,j in enumerate(free):
            K[i,j]=1
        for i in range(len(free)):
            for k,j in enumerate(pivots):
                K[i,j]=kernel_pivots[k,i]
        transform = acb_mat(arb_mat([[rational(v) for v in row] for row in record['transform']])*K)

        if not jets_certified:
            raw_jets = [[jet(p) for p in o] for o in orbits]
            for i in range(len(raw)):
                for j in range(i,len(raw)):
                    column = [arb(0) for _ in range(10)]
                    for l,(u,v) in enumerate(zip(raw_jets[i],raw_jets[j])):
                        column[0] += (u[0].conjugate()*v[0]).real/comb(J,l)
                        for k in range(1,10):
                            column[k] += (u[k].conjugate()*v[0]+u[0].conjugate()*v[k]).real/comb(J,l)
                    jets.append(column)
            if len(jets)>=10:
                import numpy as np
                from scipy.linalg import qr
                order = qr(np.array([[float(v.mid()) for v in col] for col in jets]).T,
                           pivoting=True,mode='economic')[2][:10]
                minor = arb_mat([[jets[int(i)][k] for i in order] for k in range(10)])
                jets_certified = not minor.det().contains(0)
                if jets_certified:
                    print('Verified ten independent value and first derivative conditions.',flush=True)

    m = transform.nrows()
    G = arb_mat(m,m)
    positions = []
    index = 0
    for i in range(m):
        for j in range(i,m):
            G[i,j]=G[j,i]=rational(Fraction(record['gram'][index],scale))
            if offset in selected:
                positions.append((selected[offset],i,j))
            index += 1
            offset += 1
    assert index==len(record['gram'])
    grams.append(G)
    for l in range(J+1):
        polys = [o[l] for o in orbits]
        mons = sorted(set().union(*polys))
        coeff = acb_mat([[c if isinstance(c,arb) else rational(c)
                          for c in (p.get(mon,0) for mon in mons)] for p in polys])
        evaluated = acb_mat([[multiply(z[k] for k in a)*multiply(z[k].conjugate() for k in b)
                              for z in samples] for a,b in mons])
        X = transform*coeff*evaluated
        GX = acb_mat(G)*X
        w = rational(Fraction(1,comb(J,l)))
        for i in range(m):
            row = [X[i,k].conjugate() for k in range(size)]
            for k, value in enumerate(row):
                total[k] += w*(value*GX[i,k]).real
        for col,i,j in positions:
            factor = w if i==j else 2*w
            for k in range(size):
                B[k,col] += factor*(X[i,k].conjugate()*X[j,k]).real
    print('Checked Gram block:',charge,J,m,flush=True)

assert jets_certified
assert offset>max(selected)
C = arb_mat([[rational(Fraction(v,scale)) for v in row] for row in source['inverse']])
residual = arb_mat([[a-b] for a,b in zip(target,total)])
error = C*B
norm = max(sum((error[i,j]-int(i==j)).abs_upper() for j in range(size)).upper() for i in range(size))
assert norm < 1
updates = C*residual
correction = max(updates[i,0].abs_upper() for i in range(size))/(1-norm)
for G in grams:
    m = G.nrows()
    for i in range(m):
        G[i,i] -= m*correction
    for k in range(m):
        assert G[k,k]>0
        for i in range(k+1,m):
            for j in range(k+1,m):
                G[i,j] -= G[i,k]*G[k,j]/G[k,k]

print('Inverse error bound:',norm)
print('Gram entry correction bound:',correction)
print('Verified the sharp degree-ten certificate for N=11, t=4.')
