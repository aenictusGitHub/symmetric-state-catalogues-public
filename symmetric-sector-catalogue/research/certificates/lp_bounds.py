from spectral_bounds import *
import numpy as np
from scipy.optimize import linprog
import sympy as sp

from certificate_basis import multipole_eigenvalue, constraints

def solve_bound(n,k,mode):
 gs=constraints(n);G=np.array([[float(x) for x in row] for name,row in gs]);s=[sector_eigenvalue(n,k,d) for d in range(0,n+1,2)]
 obj=np.array([float(x) for x in s])*(1 if mode=='min' else -1)
 scales=np.maximum(np.max(abs(G),axis=1),1e-12);Gn=G/scales[:,None]
 lp=linprog(obj,A_ub=-Gn,b_ub=np.zeros(len(gs)),A_eq=[np.ones(len(s))],b_eq=[1],bounds=[(None,None)]*len(s),method='highs')
 if not lp.success:raise ValueError(lp.message)
 # Exactify primal vertex from nonzero coordinates and active constraints.
 active=[i for i,row in enumerate(Gn) if abs(row@lp.x)<1e-8]
 M=sp.Matrix([[sp.Rational(x) for x in gs[i][1]] for i in active]+[[1]*len(s)]);rhs=sp.Matrix([0]*len(active)+[1])
 sol=sp.linsolve((M,rhs));vec=list(sol)[0]
 if any(x.free_symbols for x in vec):raise ValueError(('underdetermined',n,k,mode,vec))
 bound=sum(sp.Rational(x)*p for x,p in zip(s,vec));B=F(int(bound.p),int(bound.q))
 # Recover sparse nonnegative dual exactly.
 target=[(B-x if mode=='max' else x-B) for x in s]
 dual=linprog(np.zeros(len(gs)),A_eq=Gn.T,b_eq=np.array([float(x) for x in target]),bounds=[(0,None)]*len(gs),method='highs')
 if not dual.success:raise ValueError(('dual',n,k,mode,dual.message))
 support=[i for i,x in enumerate(dual.x) if x>1e-8]
 if support:
  D=sp.Matrix([[sp.Rational(gs[i][1][j]) for i in support] for j in range(len(s))]);b=sp.Matrix([sp.Rational(x) for x in target]);vals=list(sp.linsolve((D,b)))[0]
  if any(x.free_symbols or x<0 for x in vals):raise ValueError(('bad dual',n,k,mode,vals))
  cert={gs[i][0]:str(v) for i,v in zip(support,vals) if v!=0}
 else:cert={}
 return {'n':n,'k':k,'mode':mode,'value':str(B),'decimal':float(B),'certificate':cert,'pJ_vertex':{str(n-d):str(v) for d,v in zip(range(0,n+1,2),vec)}}

if __name__=='__main__':
 records=[]
 for n in range(2,13):
  print('n',n)
  for k in range(1,n+1):
   lo=solve_bound(n,k,'min');hi=solve_bound(n,k,'max');records +=[lo,hi]
   print(f' S{k}: [{lo["value"]}, {hi["value"]}]')
 Path(__file__).with_name('lp-bounds.json').write_text(json.dumps(records,indent=2))
