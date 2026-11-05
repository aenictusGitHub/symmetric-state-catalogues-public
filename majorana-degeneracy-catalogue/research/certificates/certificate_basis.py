"""Exact rational basis for sector certificates. Standard library only."""
from spectral_bounds import *

def multipole_eigenvalue(n,L,d):
 return F((2*L+1)*C(n,L),(L+1)*C(n+L+1,n))*sum(((-1)**(t+L)*C(L,t)*C(t+L,t)*purity_eigenvalue(n,d,t) for t in range(L+1)),F(0))
def constraints(n):
 ds=list(range(0,n+1,2));dim=len(ds);g=[]
 for idx,d in enumerate(ds):g.append((f'p_{n-d}',[F(int(i==idx)) for i in range(dim)]))
 for L in range(1,n+1):g.append((f'a_{L}',[multipole_eigenvalue(n,L,d) for d in ds]))
 for t in range(1,n//2+1):
  g.append((f'r_{t}-1/{t+1}',[purity_eigenvalue(n,d,t)-F(1,t+1) for d in ds]))
  g.append((f'1-r_{t}',[1-purity_eigenvalue(n,d,t) for d in ds]))
 # Established published sector bounds, separate from the bounds being optimized.
 if n>=3:
  g.append((f'C({n},2)-S_2',[F(C(n,2))-sector_eigenvalue(n,2,d) for d in ds]))
 if n>=9:
  g.append((f'C({n},4)-S_4',[F(C(n,4))-sector_eigenvalue(n,4,d) for d in ds]))
 if n==9:
  g.append(('S_8-9',[sector_eigenvalue(n,8,d)-F(9) for d in ds]))
 g.append((f'full_sector_max-S_{n}',[F(2**(n-1)+(1 if n%2==0 else 0))-sector_eigenvalue(n,n,d) for d in ds]))
 # Reduced-state spin-flip and sector inequalities: source Eq.(27).
 for t in range(1,n+1):
  alt=[sum(((-1)**m*C(n-m,t-m)*(F(1) if m==0 else sector_eigenvalue(n,m,d)) for m in range(t+1)),F(0)) for d in ds]
  g.append((f'shadow_reduced_{t}',alt))
  g.append((f'shadow_reduced_upper_{t}',[F(2**t*C(n,t))-x for x in alt]))
  even=[sum((C(n-2*m,t-2*m)*(F(1) if m==0 else sector_eigenvalue(n,2*m,d)) for m in range((t-1)//2+1)),F(0)) for d in ds]
  upper=F(2**t,4)*(3+(-1)**t)*C(n,t)
  g.append((f'even_reduced_upper_{t}',[upper-x for x in even]))
 return g

