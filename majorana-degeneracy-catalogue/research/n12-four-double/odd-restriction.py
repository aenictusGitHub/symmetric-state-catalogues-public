from pathlib import Path
import sympy as s,json,time
h,z=s.symbols('h z',real=True)
Q=z**4+2*h*z**2+1
RB=[z**4+1,z**2,s.I*(z**4-1),z**3-z,s.I*(z**3+z)]
p=[[s.Poly(s.expand(Q**2*R),z).nth(r) for r in range(13)] for R in RB]
pairs=[(i,j)for i in range(5)for j in range(i+1,5)]
rows=[('diag',a,a)for a in range(4)]+[(kind,a,c)for a in range(4)for c in range(a+1,4)for kind in['re','im']]
M=[]
for kind,a,c in rows:
 line=[]
 for i,j in pairs:
  v=s.I*(-1)**(a+c)*sum(s.binomial(9,b)*(p[i][a+b]*s.conjugate(p[j][c+b])-p[j][a+b]*s.conjugate(p[i][c+b]))/(s.binomial(12,a+b)*s.binomial(12,c+b)) for b in range(10))
  v=s.expand(s.im(v) if kind=='im' else s.re(v));line.append(v)
 # Remove rational row denominators only (do not change any h-dependent ranks).
 scale=s.ilcm(*[s.denom(q)for v in line for q in s.Poly(v,h).all_coeffs()]);M.append([s.expand(v*scale)for v in line])
M=s.Matrix(M);_,indices=M.subs(h,s.Rational(1,3)).T.rref();print('generic rank',len(indices),'rows',indices,flush=True)
selected=M[list(indices),:];print('computing determinant',flush=True);det=s.factor(selected.det(method='domain-ge'));print('det',det,flush=True)
print('rank h0',M.subs(h,0).rank(),flush=True)
Path(__file__).with_suffix('.json').write_text(json.dumps({'basis':list(map(str,RB)),'pairs':pairs,'rowDescriptions':rows,'matrix':[[str(v)for v in M.row(i)]for i in range(M.rows)],'selectedRows':indices,'determinant':str(det),'rankAtZero':int(M.subs(h,0).rank())},indent=2)+'\n')
