import sympy as s
from sympy.physics.wigner import clebsch_gordan
from scipy.optimize import linprog
n=7;j=s.Rational(7,2)
def K(t):
 B=lambda a,b:s.sqrt(s.binomial(t,a)*s.binomial(n-t,b)/s.binomial(n,a+b))
 A=s.zeros(64)
 for a in range(t+1):
  for c in range(t+1):
   for b in range(n-t+1):
    for d in range(n-t+1):
     i,h,k,l=a+b,c+d,c+b,a+d
     A[8*i+h,8*k+l]+=B(a,b)*B(c,d)*B(c,b)*B(a,d)
 return A
V={}
for J in [1,3,5,7]:
 v=s.zeros(64,1)
 for i in range(8):
  for k in range(8):
   if n-i-k==J:v[8*i+k]=clebsch_gordan(j,j,J,j-i,j-k,J)
 V[J]=v
rows=[]
for t in range(1,4):
 A=K(t);r=[s.simplify((v.T*A*v)[0]) for v in V.values()];rows.append(r);print('P',t,r)
# Linear rank inequalities P1>=1/2,P2>=1/3, normalization and p>=0
res=linprog(list(map(float,rows[2])),A_ub=[[-float(v) for v in row] for row in rows[:2]],b_ub=[-.5,-1/3],A_eq=[[1]*4],b_eq=[1],bounds=(0,None),method='highs');print(res.fun,res.x,res.ineqlin,res.eqlin)
c=s.zeros(8,1);c[0]=s.sqrt(30)/10;c[5]=s.sqrt(70)/10
for t in [1,2,3]:
 A=s.Matrix(t+1,8-t,lambda a,b:c[a+b]*s.sqrt(s.binomial(t,a)*s.binomial(7-t,b)/s.binomial(7,a+b)));print('AC purity',t,s.simplify(s.trace((A*A.H)**2)))
P=[s.ones(1,4)]+[s.Matrix([r]) for r in rows]+[s.Matrix([r]) for r in rows[::-1]]+[s.ones(1,4)]
S=[s.ones(1,4)]
for t in range(1,8):S.append(s.simplify(s.binomial(7,t)*(2**t*P[t]-sum((s.binomial(t,k)/s.binomial(7,k)*S[k] for k in range(t)),s.zeros(1,4)))))
a=[]
for t in range(8):
 f=lambda l:s.factorial(7+l+1)*s.factorial(7-l)*s.factorial(t)**2/(s.factorial(t+l+1)*s.factorial(t-l)*s.factorial(7)**2)
 a.append(s.simplify((P[t]-sum((f(l)*a[l] for l in range(t)),s.zeros(1,4)))/f(t)))
cons=[*S[1:],*a[1:]];names=[f'S{k}' for k in range(1,8)]+[f'a{k}' for k in range(1,8)]
res=linprog(list(map(float,rows[2])),A_ub=[[-float(v) for v in r] for r in cons],b_ub=[0]*len(cons),A_eq=[[1]*4],b_eq=[1],bounds=(0,None),method='highs')
print('FULL LP',res.fun,res.x)
for name,r,m in zip(names,cons,res.ineqlin.marginals):
 if abs(m)>1e-8:print(name,r,'dual',m)
print('eq dual',res.eqlin.marginals)
for i,r in enumerate(S):print('S',i,list(r))
for i,r in enumerate(a):print('a',i,list(r))
import json
from pathlib import Path
D=json.loads(Path('output/github_publication_20261106/catalogues/symmetric-sector-catalogue/data/catalogue.js').read_text().split('=',1)[1].strip().rstrip(';'))
cs=[list(map(float,r)) for r in cons];bs=[0.]*len(cs);nn=names.copy()
for e in D['entries']:
 if e['n']!=7 or e['status']!='proven':continue
 row=S[e['k']];sgn=1 if e['mode']=='min' else -1;v=e.get('numericValue',e['decimal']);cs.append([sgn*float(z) for z in row]);bs.append(sgn*v);nn.append(f"known S{e['k']} {e['mode']} {v}")
res=linprog(list(map(float,rows[2])),A_ub=[[-z for z in r] for r in cs],b_ub=[-b for b in bs],A_eq=[[1]*4],b_eq=[1],bounds=(0,None));print('With known bounds',res.fun)
for name,m in zip(nn,res.ineqlin.marginals):
 if abs(m)>1e-8:print(name,m)
pp=s.Matrix.vstack(s.ones(1,4),S[1],S[3],a[2]).inv()*s.Matrix([1,0,0,0]);print('Exact relaxed point',list(pp));assert (s.Matrix([rows[2]])*pp)[0]==s.Rational(1,4)
assert all((r*pp)[0]>=0 for r in cons)
