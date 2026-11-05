import sympy as s,json
from pathlib import Path
R=Path(__file__).resolve().parent
A,B,T=s.symbols('A B T',real=True)
d=json.loads((R/'reduced-polynomials.json').read_text())
loc={'A':A,'B':B,'T':T};h1=s.sympify(d['h1'],locals=loc);h2=s.sympify(d['h2'],locals=loc)
subs=s.subresultants(h1,h2,A)
out=[]
for p in subs:
 p=s.factor(p)
 degree=s.degree(p,A)
 print('degreeA',degree,'length',len(str(p)),flush=True)
 if degree==1:
  p=s.Poly(p,A).primitive()[1].as_expr()
  rat=s.factor(-p.subs(A,0)/s.diff(p,A))
  num,den=s.fraction(rat)
  print('A=',rat,flush=True)
  out={'Aformula':str(rat),'Anumerator':str(num),'Adenominator':str(den)}
res=s.sympify(d['resultantB'],locals=loc)
pol=[f for f,m in s.factor_list(res)[1] if s.degree(f,B)==11][0]
out['Bpolynomial']=str(pol)
out['Binterval']=['2498461150300332175000977146562933123120/10^40','2498461150300332175000977146562933123121/10^40']
out['Tnumerator']=d['Tnumerator'];out['Tdenominator']=d['Tdenominator']
(R/'exact-algebraic.json').write_text(json.dumps(out,indent=2)+'\n')
