"""Solver-free rational verification of the saved degree-six compatibility certificate.
The local pickle supplies only the explicitly constructed polynomial covariants;
this verifier rebuilds both sides independently and checks positive Gram matrices.
"""
from pathlib import Path
from collections import defaultdict,Counter
from math import comb,factorial
import sympy as s,pickle,json
root=Path(__file__).resolve().parent
D=json.loads((root/'cubic_covariants.json').read_text());C=json.loads((root/'lower_certificate.json').read_text());N=8

def mul(p,q):
 out=defaultdict(lambda:s.S.Zero)
 for(a,b),v in p.items():
  for(c,d),w in q.items():out[(tuple(sorted(a+c)),tuple(sorted(b+d)))]+=v*w
 return {k:v for k,v in out.items() if v}
def conj(p):return {(b,a):v for(a,b),v in p.items()}
P=defaultdict(lambda:s.S.Zero)
for i in range(5):
 for h in range(5):
  for j in range(5):
   for l in range(5):P[(tuple(sorted((i+j,h+l))),tuple(sorted((h+j,i+l))))]+=comb(4,i)*comb(4,h)*comb(4,j)*comb(4,l)
norm={((k,),(k,)):s.Integer(comb(N,k))for k in range(N+1)}
for m,v in mul(norm,norm).items():P[m]-=s.Rational(C['bound'])*v
lhs=mul(P,norm);rhs=defaultdict(lambda:s.S.Zero)
for block,record in zip(D,C['blocks']):
 charge,J=block['charge'],block['twice_spin']
 orb=[[{(tuple(k[0]),tuple(k[1])):s.Rational(v) for k,v in poly}for poly in o]for o in block['orbits']]
 assert charge==record['charge'] and J==record['twice_spin']
 G=s.Matrix(record['gram']);assert G==G.T
 _,dd=G.LDLdecomposition(hermitian=False);assert all(dd[i,i]>0 for i in range(dd.rows))
 for i in range(G.rows):
  for j in range(G.cols):
   for level in range(J+1):
    for m,v in mul(orb[i][level],conj(orb[j][level])).items():rhs[m]+=G[i,j]*v/s.binomial(J,level)
assert all(lhs.get(k,0)==rhs.get(k,0)for k in set(lhs)|set(rhs))
# Attainment of the exact upper bound; this is not an optimality assertion.
p=s.Rational(87,412);c=[s.sqrt(p),0,0,0,s.sqrt(1-2*p),0,0,0,-s.sqrt(p)]
A=s.Matrix(5,5,lambda a,b:c[a+b]*s.sqrt(s.binomial(4,a)*s.binomial(4,b)/s.binomial(8,a+b)))
assert s.simplify(s.trace((A*A.H)**2))==s.Rational(113,515)
# Exact minimization of the remaining C4 residue classes.
u=s.symbols('u',nonnegative=True)
for i,j,bound in [(1,5,s.Rational(1,4)),(2,6,s.Rational(25,98)),(3,7,s.Rational(1,4))]:
 z=[0]*9;z[i]=s.sqrt(u);z[j]=s.sqrt(1-u)
 M=s.Matrix(5,5,lambda a,b:z[a+b]*s.sqrt(s.binomial(4,a)*s.binomial(4,b)/s.binomial(8,a+b)))
 F=s.factor(s.trace((M*M.H)**2)).subs(s.conjugate(s.sqrt(1-u)),s.sqrt(1-u))
 # Use the real polynomial continuation on 0<=u<=1.
 F=s.factor(s.trace((M*M.T)**2));uc=s.solve(s.diff(F,u),u)[0]
 assert 0<=uc<=1 and s.simplify(F.subs(u,uc))==bound
report={'exact_global_lower_bound':C['bound'],'exact_attained_upper_bound':'113/515','full_polynomial_identity':True,'exact_positive_Gram_blocks':len(C['blocks']),'global_minimum_proved':False}
(root/'verification.json').write_text(json.dumps(report,indent=2));print(report)
