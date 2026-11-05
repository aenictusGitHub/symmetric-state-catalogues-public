"""Build an exact algebraic existence certificate from the exploratory C5 candidate."""
from pathlib import Path
from collections import defaultdict
from math import comb,prod
import json
import numpy as np
import sympy as s
HERE=Path(__file__).resolve().parent
OUT=HERE/'n24-k5';OUT.mkdir(exist_ok=True)
N=24;k=5
slots=[(a,r,p) for a in range(3) for r in range(a,25 if a<2 else 13,5) for p in (0,1)]
index={v:i for i,v in enumerate(slots)}
fixed=[1,11,21];unknown=[i for i in range(len(slots)) if i not in fixed]
coefficient_map={}
for a in range(k):
 for r in range(a,N+1,k):
  aa,rr=(a,r) if a<3 else (4-a,24-r)
  if aa==2 and rr>12:rr=24-rr
  coefficient_map[a,r]=(index[aa,rr,0],index[aa,rr,1])

def polynomial(i,j,q,p,part):
 terms=defaultdict(int)
 for r in range(q,N+1):
  if (i,r-q) not in coefficient_map or (j,r) not in coefficient_map:continue
  ar,ai=coefficient_map[i,r-q];br,bi=coefficient_map[j,r]
  m=N-2*r
  w=comb(N,r)*prod(range(r-q+1,r+1))*(3*m*m-N*(N+2) if p==2 else m**p)
  pairs=[(ar,br,1),(ai,bi,1)] if part==0 else [(ar,bi,1),(ai,br,-1)]
  for a,b,sign in pairs:
   if a not in fixed and b not in fixed:terms[tuple(sorted((a,b)))]+=sign*w
 return {'basisPair':[i,j],'raisingPower':q,'zPower':p,'part':['real','imaginary'][part],'constant':-int(q==p==part==0 and i==j),'terms':[[a,b,v] for (a,b),v in sorted(terms.items()) if v]}
full=[]
for p in range(3):
 for i in range(k):
  for j in range(i,k):
   for part in range(1 if i==j else 2):
    e=polynomial(i,j,0,p,part)
    if e['terms'] or e['constant']:full.append(e)
for q,p in [(1,0),(1,1),(2,0)]:
 for i in range(k):
  for j in range(k):
   for part in range(2):
    e=polynomial(i,j,q,p,part)
    if e['terms']:full.append(e)
mons=[(-1,-1)]+sorted({(a,b) for e in full for a,b,v in e['terms']});mi={v:i for i,v in enumerate(mons)}
rows=[]
for e in full:
 row=[0]*len(mons);row[0]=e['constant']
 for a,b,v in e['terms']:row[mi[a,b]]=v
 rows.append(row)
M=s.Matrix(rows)
_,selected=M.T.rref()
print('slots',len(slots),'unknown',len(unknown),'physical eqs',len(full),'poly rank',len(selected),flush=True)
base=M[list(selected),:]
# Recover every discarded equation as a rational linear combination.
_,pivot=base.rref();square=base[:,list(pivot)];reconstruction=M[:,list(pivot)]*square.inv()
assert reconstruction*base==M
system={'n':N,'dimension':k,'order':2,'slots':slots,'coefficientMap':{f'{a},{r}':list(v) for (a,r),v in coefficient_map.items()},'fixedSlotIndices':fixed,'unknownSlotIndices':unknown,'fullEquations':full,'independentEquationIndices':selected,'equations':[full[i] for i in selected],'fullEquationReconstruction':[[str(v) for v in row] for row in reconstruction.tolist()]}
(OUT/'system.json').write_text(json.dumps(system,indent=2)+'\n')
# Rotate and rephase the numerical candidate into the exact reflection ansatz.
d=json.loads((HERE/'cyclic-n24-k5.json').read_text());z=np.array([complex(*v) for v in d['components']])
theta=-np.angle(z[17]/z[7])/10;z*=np.exp(1j*np.arange(25)*theta)
for a in range(3):z[a::k]*=np.exp(-1j*np.angle(z[a]))
x=[]
for a,r,p in slots:x.append((z[r].real if p==0 else z[r].imag)/np.sqrt(comb(N,r)))
(OUT/'initial.json').write_text(json.dumps(x,indent=2)+'\n')
if len(selected)!=len(unknown):raise SystemExit('Not square; inspect before certifying')

import mpmath as mp
mp.mp.dps=120
ix={u:i for i,u in enumerate(unknown)}
eqs=system['equations']
def f(*x):
 return tuple(mp.mpf(e['constant'])+sum(v*x[ix[a]]*x[ix[b]] for a,b,v in e['terms']) for e in eqs)
def jac(*x):
 out=mp.matrix(len(x))
 for ei,e in enumerate(eqs):
  for a,b,v in e['terms']:out[ei,ix[a]]+=v*x[ix[b]];out[ei,ix[b]]+=v*x[ix[a]]
 return out
sol=mp.findroot(f,[mp.mpf(str(x[i])) for i in unknown],J=jac,tol=mp.mpf('1e-110'),maxsteps=30)
print('sol residual',mp.nstr(max(abs(v) for v in f(*sol)),8),flush=True)
D=10**80;AD=10**70
center=[0]*len(slots)
for u,v in zip(unknown,sol):center[u]=int(mp.nint(v*D))
inv=jac(*sol)**-1;A=[[int(mp.nint(inv[i,j]*AD)) for j in range(len(unknown))] for i in range(len(unknown))]
from fractions import Fraction as F
n=len(unknown);fnum=[];jnum=[]
for e in eqs:
 value=e['constant']*D*D;grad=[0]*len(slots)
 for a,b,v in e['terms']:value+=v*center[a]*center[b];grad[a]+=v*center[b];grad[b]+=v*center[a]
 fnum.append(value);jnum.append([grad[i] for i in unknown])
defect=F(max(sum(abs((AD*D if i==j else 0)-sum(A[i][p]*jnum[p][j] for p in range(n))) for j in range(n)) for i in range(n)),AD*D)
correction=F(max(abs(sum(A[i][p]*fnum[p] for p in range(n))) for i in range(n)),AD*D*D)
anorm=F(max(sum(abs(v) for v in row) for row in A),AD)
hessian=max(2*sum(abs(v) for a,b,v in e['terms']) for e in eqs)
radius=F(1,10**60);contraction=defect+anorm*hessian*radius;selfmap=correction+contraction*radius
assert contraction<1 and selfmap<radius
print('bounds', {k:float(v) for k,v in {'defect':defect,'correction':correction,'anorm':anorm,'hessian':hessian,'contraction':contraction,'selfmap':selfmap}.items()},flush=True)
system.update({'centerDenominator':str(D),'inverseDenominator':str(AD),'fullCenterNumerators':list(map(str,center)),'inverseNumerators':[[str(v) for v in row] for row in A],'radius':str(radius),'fixedDenominator':'1','fixedNumerators':{str(u):'0' for u in fixed},'bounds':{key:str(val) for key,val in [('inverseDefect',defect),('newtonCorrection',correction),('inverseInfinityNorm',anorm),('jacobianLipschitzBound',F(hessian)),('contractionFactor',contraction),('selfMapRadiusBound',selfmap)]}})
(OUT/'certificate.json').write_text(json.dumps(system,indent=2)+'\n')
# Full-precision Dicke coefficients at the certified root centre.
basis=[]
for a in range(k):
 coeff=[]
 for r in range(N+1):
  if (a,r) not in coefficient_map:coeff.append(['0','0'])
  else:
   ri,ii=coefficient_map[a,r];coeff.append([mp.nstr(mp.sqrt(comb(N,r))*mp.mpf(center[u])/D,85) for u in [ri,ii]])
 basis.append(coeff)
(OUT/'numerical-basis.json').write_text(json.dumps(basis,indent=2)+'\n')
names={u:f'x{i+1}' for i,u in enumerate(unknown)}
atom=lambda u:names.get(u,'0')
analytical={'id':'algebraic-n24-k5-t2','n':24,'t':2,'dimension':5,'label':'Exact algebraic five-dimensional second-order AC subspace','evidence':'exact','exactnessKind':'isolated-algebraic-root','basis':[],
 'algebraicDefinition':{'kind':'unique-real-root-in-rational-box','variables':[names[u] for u in unknown],'coordinateSlots':slots,'unknownSlotIndices':unknown,'fixedCoordinates':{str(u):'0' for u in fixed},'parameterIntervals':{names[u]:[str(F(center[u],D)-radius),str(F(center[u],D)+radius)] for u in unknown},'equations':['+'.join(([str(e['constant'])] if e['constant'] else [])+[f'({v})*({atom(a)})*({atom(b)})' for a,b,v in e['terms']]) or '0' for e in eqs],'certificate':'certificate.json','definition':'The unique real zero of the 23 integer quadratic equations within the stated rational box. Exact rational contraction proves existence and uniqueness; all Gram and rank-one/two equations follow by exact rational linear combinations.'},
 'derivationLatex':r'\begin{aligned}u_1&=\sum_{\ell=0}^4a_\ell|D_{24}^{(5\ell)}\rangle,&u_5&=Ru_1,\\u_2&=\sum_{\ell=0}^4b_\ell|D_{24}^{(1+5\ell)}\rangle,&u_4&=Ru_2,\\u_3&=c_2(|D_{24}^{(2)}\rangle+|D_{24}^{(22)}\rangle)+c_7(|D_{24}^{(7)}\rangle+|D_{24}^{(17)}\rangle)+c_{12}|D_{24}^{(12)}\rangle,\\a_\ell&=c_{5\ell},\qquad b_\ell=c_{1+5\ell},\\R|D_{24}^{(r)}\rangle&=|D_{24}^{(24-r)}\rangle.\end{aligned}',
 'parameterLatex':r'\begin{aligned}(\Re z_0,\Re z_5,\Im z_5,\ldots,\Re z_{20},\Im z_{20})&=(x_1,\ldots,x_9),\\(\Re z_1,\Re z_6,\Im z_6,\ldots,\Re z_{21},\Im z_{21})&=(x_{10},\ldots,x_{18}),\\(\Re z_2,\Re z_7,\Im z_7,\Re z_{12},\Im z_{12})&=(x_{19},\ldots,x_{23}),\\\Im z_0=\Im z_1=\Im z_2&=0,\qquad c_r=\sqrt{\binom{24}{r}}z_r.\end{aligned}',
 'provenance':'New cyclic-residue search followed by unitary-reflection reduction and an exact isolated algebraic root certificate. This construction attains the general upper bound k<=floor((N+5)/5)=5 at N=24. A short radical expression is not asserted.'}
for a in range(k):
 exact=[];numeric=[]
 for r in range(N+1):
  if (a,r) not in coefficient_map:exact.append('0');numeric.append([0,0])
  else:
   ri,ii=coefficient_map[a,r];exact.append(f'sqrt({comb(N,r)})*(({atom(ri)})+I*({atom(ii)}))');numeric.append([float(mp.sqrt(comb(N,r))*mp.mpf(center[u])/D) for u in [ri,ii]])
 analytical['basis'].append({'id':f'algebraic-n24-k5-t2-basis-{a+1}','n':24,'label':f'Basis vector {a+1}','exactCoefficients':exact,'coefficients':numeric,'support':list(range(a,25,5)),'stateLatex':rf'|u_{{{a+1}}}\rangle=\sum_{{\ell=0}}^4 c_{{{a}+5\ell}}|D_{{24}}^{{({a}+5\ell)}}\rangle','coefficientDefinitionLatex':r'c_r=\sqrt{\binom{24}{r}}z_r'})
(OUT/'basis.json').write_text(json.dumps(analytical,indent=2)+'\n')
