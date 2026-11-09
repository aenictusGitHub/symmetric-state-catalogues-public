import sympy as s
from sympy.physics.wigner import clebsch_gordan
from scipy.optimize import linprog
n=8;j=s.Integer(4)
def K(t):
 B=lambda a,b:s.sqrt(s.binomial(t,a)*s.binomial(n-t,b)/s.binomial(n,a+b))
 A=s.zeros((n+1)**2)
 for a in range(t+1):
  for c in range(t+1):
   for b in range(n-t+1):
    for d in range(n-t+1):
     i,h,k,l=a+b,c+d,c+b,a+d
     A[(n+1)*i+h,(n+1)*k+l]+=B(a,b)*B(c,d)*B(c,b)*B(a,d)
 return A
V={}
for J in [0,2,4,6,8]:
 v=s.zeros((n+1)**2,1)
 for i in range(n+1):
  for k in range(n+1):
   if n-i-k==J:v[(n+1)*i+k]=clebsch_gordan(j,j,J,j-i,j-k,J)
 V[J]=v
rows=[]
for t in range(1,5):
 A=K(t);r=[s.simplify((v.T*A*v)[0]) for v in V.values()];rows.append(r);print('P',t,r)

P=[s.ones(1,5)]+[s.Matrix([r]) for r in rows]+[s.Matrix([r]) for r in rows[-2::-1]]+[s.ones(1,5)]
S=[s.ones(1,5)]
for t in range(1,n+1):
 S.append(s.simplify(s.binomial(n,t)*(2**t*P[t]-sum((s.binomial(t,k)/s.binomial(n,k)*S[k] for k in range(t)),s.zeros(1,5)))))
a=[]
for t in range(n+1):
 f=lambda l:s.factorial(n+l+1)*s.factorial(n-l)*s.factorial(t)**2/(s.factorial(t+l+1)*s.factorial(t-l)*s.factorial(n)**2)
 a.append(s.simplify((P[t]-sum((f(l)*a[l] for l in range(t)),s.zeros(1,5)))/f(t)))
cons=[*S[1:],*a[1:]]; names=[f'S{k}' for k in range(1,n+1)]+[f'a{k}' for k in range(1,n+1)]
res=linprog(list(map(float,rows[-1])),A_ub=[[-float(v) for v in r] for r in cons],b_ub=[0]*len(cons),A_eq=[[1]*5],b_eq=[1],bounds=(0,None),method='highs')
print('FULL LP',res.fun,res.x)
for name,r,m in zip(names,cons,res.ineqlin.marginals):
 if abs(m)>1e-8:print(name,list(r),'dual',m)
print('eq dual',res.eqlin.marginals)
for i,r in enumerate(S):print('S',i,list(r))
for i,r in enumerate(a):print('a',i,list(r))
import json
from pathlib import Path
Path(__file__).with_name('rows.json').write_text(json.dumps({'P':[[str(x) for x in r] for r in P],'S':[[str(x) for x in r] for r in S],'a':[[str(x) for x in r] for r in a]}))
