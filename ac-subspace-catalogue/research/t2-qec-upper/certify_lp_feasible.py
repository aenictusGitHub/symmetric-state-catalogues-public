"""Rational feasible enumerators show the stated relaxation cannot exclude k.
This does NOT prove a physical code exists.
"""
from pathlib import Path
from fractions import Fraction as F
import json,importlib.util
import numpy as np
HERE=Path(__file__).resolve().parent
spec=importlib.util.spec_from_file_location('intrinsic_lp',HERE/'intrinsic_lp.py');mod=importlib.util.module_from_spec(spec);spec.loader.exec_module(mod)

def system(n,k):
 d=n+1;M=[[mod.racah(n,a,b)for a in range(d)]for b in range(d)];I=[[F(i==j)for j in range(d)]for i in range(d)]
 eq=[([F(1)]*d,F(1)),(I[0],F(k,d)),(I[1],F(0)),(I[2],F(0)),(M[1],F(0)),(M[2],F(0))]
 # All inequalities are r dot q <= rhs.
 iq=[([-v for v in row],F(0))for row in I]
 def add(row,b=0):iq.append((list(map(F,row)),F(b)))
 for l in range(d):add([I[l][a]-k*M[l][a]for a in range(d)]);add([-v for v in M[l]])
 if n%2 and k%2:add([(-1)**l for l in range(d)],F(k-1,k))
 for l in range(d):
  row=[M[l][a]*(-1)**a for a in range(d)]
  add([-v for v in row]);add(row,F(min(2*l+1,F(k*(k+(-1)**(n-l)),2)),k))
 c=F(n*(n+2),12);r=d-4*k
 fac=[[F(1)]*d,[3-F(l*(l+1),2)/c for l in range(d)],[-4+F(l*(l+1),2)/c for l in range(d)]]
 fac[1][0]=F(3);fac[2][0]=F(r,k);dims=[k,3*k,r]
 for a in range(3):
  for b in range(a,3):
   f=[u*v for u,v in zip(fac[a],fac[b])]
   for l in range(d):
    row=[M[l][i]*f[i]for i in range(d)];sh=[v*(-1)**i for i,v in enumerate(row)]
    add([-v for v in row]);add([-v for v in sh]);add(sh,F(2*l+1,k))
    if a==b:add([I[l][i]*f[i]-dims[a]*row[i]for i in range(d)])
 return eq,iq

def rationalize(n,k):
 eq,iq=system(n,k);A=np.array([[float(v)for v in r]for r,b in iq]);b=np.array([float(b)for r,b in iq]);E=np.array([[float(v)for v in r]for r,b in eq]);e=np.array([float(b)for r,b in eq]);
 from scipy.optimize import linprog
 sol=linprog(np.zeros(n+1),A_ub=A,b_ub=b,A_eq=E,b_eq=e,bounds=(None,None),method='highs')
 assert sol.success
 active=[i for i,v in enumerate(b-A@sol.x)if abs(v)<1e-8]
 equations=list(eq)+[iq[i]for i in active]
 # Pick independent rows at numerical point; exact rank and solution afterward.
 from scipy.linalg import qr
 rows=np.array([[float(v)for v in row]for row,val in equations]);scale=np.maximum(np.max(abs(rows),axis=1),1e-20)
 _,_,perm=qr((rows/scale[:,None]).T,pivoting=True,mode='economic')
 selected=perm[:n+1];basis=[equations[i][0]for i in selected];rhs=[equations[i][1]for i in selected]
 rows=[list(row)+[val]for row,val in zip(basis,rhs)];d=n+1
 for col in range(d):
  pivot=min((r for r in range(col,d)if rows[r][col]),key=lambda r:rows[r][col].numerator.bit_length()+rows[r][col].denominator.bit_length())
  rows[col],rows[pivot]=rows[pivot],rows[col];v=rows[col][col];rows[col]=[z/v for z in rows[col]]
  for r in range(d):
   if r!=col and rows[r][col]:
    v=rows[r][col];rows[r]=[a-v*b for a,b in zip(rows[r],rows[col])]
 x=[row[-1]for row in rows]
 assert all(sum(a*v for a,v in zip(row,x))==val for row,val in eq)
 slacks=[val-sum(a*v for a,v in zip(row,x))for row,val in iq]
 assert min(slacks)>=0,(n,k,float(min(slacks)))
 out={'passed':True,'n':n,'k':k,'point':[str(v)for v in x],'equalitiesChecked':len(eq),'inequalitiesChecked':len(iq),'minimumExactSlack':str(min(slacks)),'method':'Exact rational reconstruction and direct verification of every enumerator, Cauchy, shadow, projector/error/remainder and applicable skew-time-reversal inequality.','scope':'Feasible point of necessary-condition relaxation only; not an AC2 code, and no claim of physical projector existence.'}
 (HERE/f'lp-feasible-N{n}-k{k}.json').write_text(json.dumps(out,indent=2)+'\n')
 print(n,k,'passed',len(iq),flush=True)
 return out
if __name__=='__main__':
 for n,k in [(35,8),(41,9),(42,9)]:rationalize(n,k)
