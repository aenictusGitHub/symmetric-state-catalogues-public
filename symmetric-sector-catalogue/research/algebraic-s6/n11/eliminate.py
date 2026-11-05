import json,sympy as s,time
from pathlib import Path
R=Path(__file__).resolve().parent
raw=json.loads((R/'polynomials.json').read_text())
a,b,c,d=s.symbols('a b c d',real=True)
F=s.sympify(raw['sectorPolynomial'],locals={'a':a,'b':b,'c':c,'d':d})
Fd=s.diff(F,d).subs(d,1)
G=[55*s.diff(F,a).subs(d,1)-a*Fd,s.diff(F,b).subs(d,1)-3*b*Fd,5*s.diff(F,c).subs(d,1)-42*c*Fd]
G=[s.Poly(g,a,b,c).primitive()[1].as_expr() for g in G]
print('quartic equations',G,flush=True)
(R/'quartic-equations.json').write_text(json.dumps([str(g) for g in G],indent=2)+'\n')
start=time.time()
gb=s.groebner(G,a,b,c,order='grevlex',method='f5b')
print('grevlex done',time.time()-start,'length',len(gb.polys),'zero?',gb.is_zero_dimensional,flush=True)
(R/'groebner-grevlex.json').write_text(json.dumps([str(p.as_expr()) for p in gb.polys]))
if gb.is_zero_dimensional:
 lex=gb.fglm('lex')
 print('lex done',time.time()-start,flush=True)
 (R/'groebner-lex.json').write_text(json.dumps([str(p.as_expr()) for p in lex.polys]))
 for p in lex.polys:print(s.factor(p.as_expr()),flush=True)
