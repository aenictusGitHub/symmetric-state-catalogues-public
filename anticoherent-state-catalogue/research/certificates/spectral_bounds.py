"""Exact two-copy SU(2) spectral bounds, using fractions only."""
from fractions import Fraction as F
from math import comb, factorial
import json
from pathlib import Path

def C(n,k):
 return comb(n,k) if 0<=k<=n else 0

def purity_eigenvalue(n,d,t):
 # d=n-J. Use unnormalized Dicke sums U_a, ||U_a||²=C(n,a).
 # Highest-weight vector has coefficient v_c=(-1)^c (n-d+c)!(n-c)!/((n-d)!n!).
 return sum((F((-1)**c*factorial(n-d+c)*factorial(n-c),factorial(n-d)*factorial(n))*C(t,c)*C(n-t,d-c) for c in range(d+1)),F(0))/C(n,d)

def sector_eigenvalue(n,k,d):
 return C(n,k)*sum(((-1)**(k+t)*2**t*C(k,t)*purity_eigenvalue(n,d,t) for t in range(k+1)),F(0))

def record(n,k):
 ev={n-d:sector_eigenvalue(n,k,d) for d in range(0,n+1,2)}
 lo=min(ev.values());hi=max(ev.values())
 return {'n':n,'k':k,'eigenvalues':{str(j):str(x) for j,x in ev.items()},'lower':str(lo),'upper':str(hi),'lowerFloat':float(lo),'upperFloat':float(hi),'minJ':[j for j,v in ev.items() if v==lo],'maxJ':[j for j,v in ev.items() if v==hi]}

if __name__=='__main__':
 records=[record(n,k) for n in range(2,13) for k in range(1,n+1)]
 out=Path(__file__).with_name('spectral-bounds.json');out.write_text(json.dumps(records,indent=2))
 for n in range(2,13):
  print('n=',n)
  print('  '+'; '.join(f'S{k}: [{record(n,k)["lower"]}, {record(n,k)["upper"]}]' for k in range(1,n+1)))
