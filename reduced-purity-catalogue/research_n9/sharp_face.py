"""Test the cubic SOS face forced by vanishing at the four-body candidate."""
import sys,pickle
from pathlib import Path
sys.path.insert(0,'tmp/logconvex_mixed_N6_20261013/theory/deps')
import sympy as s
from collections import Counter
from math import comb,factorial
from sympy.polys.matrices import DomainMatrix
t=int(sys.argv[1])
assert t in (3,4)
d=pickle.loads(Path(__file__).with_name(f'basis{t}.pkl').read_bytes())
r=s.symbols('r')
p=(6174*r**3+252*r**2-173*r-4 if t==4 else 53802*r**3+1260*r**2-1439*r-20)
K=s.QQ.algebraic_field(s.CRootOf(p,2))
a=K.unit
z={0:42*a,3:K.one,6:-K.one,9:-42*a}
def at(poly):
    result=K.zero
    for (u,v),c in poly.items():
        if any(k not in z for k in u+v): continue
        term=K.convert(c)
        for k in u+v:term*=z[k]
        result+=term
    return result
def dm(rows,n):
    return DomainMatrix(rows,(len(rows),n),K)
def norm3(key):
    u,v=key
    if u!=v:return 0
    return factorial(3)//s.prod(factorial(m) for m in Counter(u).values())*s.prod(comb(9,k) for k in u)
purity=((777924*r**4+12348*r**2-1680*r+1159)/3528 if t==4 else (1555848*r**4+40572*r**2-1680*r+2405)/7056)/(21*r**2+1)**2
beta=K.from_sympy(purity.subs(r,K.ext.as_expr()))
A=d['A']
offset=0;blocks=[];columns=[]
for (charge,J),orbits in d['blocks']:
    n=len(orbits)
    ker=dm([[at(o[l]) for o in orbits] for l in range(J+1)],n).nullspace().to_list()
    size=len(ker)
    if size: blocks.append((charge,J,ker))
    print((charge,J),n,size,flush=True)
    # Gram coefficients in original block coordinates.
    for i in range(size):
        for j in range(i,size):
            col={}
            index=offset
            for u in range(n):
                for v in range(u,n):
                    val=ker[i][u]*ker[j][v]
                    if i!=j:val+=ker[j][u]*ker[i][v]
                    if val:col[index]=val
                    index+=1
            columns.append(col)
    offset+=n*(n+1)//2
T=dm([[col.get(i,K.zero) for col in columns] for i in range(A.cols)],len(columns))
rows=d["rows"]
M=DomainMatrix.from_Matrix(A[rows,:]).convert_to(K)*T
b=dm([[K.convert(v)+(K.convert(s.Rational(1,4))-beta)*K.convert(norm3(k))] for v,k in ((d["b"][i],d["keys"][i]) for i in rows)],1)
# The full coefficient system must be consistent before checking positivity.
pivots=M.hstack(b).rref()[1]
print('columns',len(columns),'augmented pivots',pivots,flush=True)
if len(columns) in pivots:
    print('Exact coefficient system inconsistent on equality face')
    sys.exit(0)
print('Exact coefficient system consistent')
import numpy as np
import cvxpy as cp
import scipy.linalg as la
def floating(v):
    return float(K.to_sympy(v).evalf(30))
An=np.array([[floating(v) for v in row] for row in M.to_list()])
bn=np.array([floating(row[0]) for row in b.to_list()])
sc=[]
for _,_,ker in blocks:
    weights=[2.**round(np.log2(np.linalg.norm([floating(v) for v in row]))) for row in ker]
    sc.extend(weights[i]*weights[j] for i in range(len(ker)) for j in range(i,len(ker)))
An/=np.array(sc)
row_scale=np.linalg.norm(An,axis=1)
An/=row_scale[:,None]
bn/=row_scale
ix=la.qr(An.T,pivoting=True)[2][:len(pivots)]
gs=[cp.Variable((len(ker),len(ker)),symmetric=True) for _,_,ker in blocks]
x=cp.hstack([g[i,j] for g in gs for i in range(g.shape[0]) for j in range(i,g.shape[0])])
eps=cp.Variable()
prob=cp.Problem(cp.Maximize(eps),[g>>eps*np.eye(g.shape[0]) for g in gs]+[An[ix]@x==bn[ix]])
prob.solve(solver='CLARABEL',tol_gap_abs=1e-10,tol_feas=1e-10,max_iter=300)
print(prob.status,eps.value,flush=True)
if x.value is not None:
    print('residual',max(abs(An@x.value-bn)),flush=True)
# Round free variables, then enforce all coefficient equalities in Q(r).
raw=[K.convert(s.Rational(round(v*10**12),10**12)/s.Rational(scale)) for v,scale in zip(x.value,sc)]
rr,piv=M.hstack(b).rref()
for row,j in zip(rr.to_list(),piv):
    raw[j]=row[-1]-sum((row[k]*raw[k] for k in range(len(raw)) if k not in piv),K.zero)
original=T*dm([[v] for v in raw],1)
full_rhs=dm([[K.convert(v)+(K.convert(s.Rational(1,4))-beta)*K.convert(norm3(k))] for v,k in zip(d['b'],d['keys'])],1)
assert (DomainMatrix.from_Matrix(A).convert_to(K)*original-full_rhs).is_zero_matrix
print('Full coefficient identity verified exactly',flush=True)

# Rational interval arithmetic in the positive embedding of Q(r).
lo,hi=s.Rational(0),s.Rational(1)
for _ in range(160):
    m=(lo+hi)/2
    if p.subs(r,m)>0:hi=m
    else:lo=m
def positive(v):
    l=h=s.Rational(0)
    for c in v.to_list():
        corners=(l*lo,l*hi,h*lo,h*hi)
        c=s.QQ.to_sympy(c)
        l,h=min(corners)+c,max(corners)+c
    return l>0
import json
def encode(v):
    return [str(s.QQ.to_sympy(c)) for c in v.to_list()]
cert=[];offset=0
for charge,J,ker in blocks:
    n=len(ker)
    g=[[K.zero for _ in range(n)] for _ in range(n)]
    for i in range(n):
        for j in range(i,n):
            g[i][j]=g[j][i]=raw[offset];offset+=1
    assert all(positive(dm([row[:k] for row in g[:k]],k).det()) for k in range(1,n+1))
    cert.append(dict(charge=charge,twice_spin=J,kernel=[[encode(v) for v in row] for row in ker],
                     gram=[[encode(v) for v in row] for row in g]))
Path(__file__).with_name(f'sharp_certificate_t{t}.json').write_text(json.dumps(cert,indent=2))
print('Exact positivity verified in all reduced blocks',flush=True)
