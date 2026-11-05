from polynomials import *
from pathlib import Path
import time,json
u,z=s.symbols('u z')
G=[s.Poly(s.expand((F[0]/b).subs(d,z*b)),b).as_dict(),s.Poly(s.expand(F[1].subs(d,z*b)),b).as_dict(),s.Poly(s.expand((F[2]/b).subs(d,z*b)),b).as_dict()]
H=[]
for g in G:
 assert set(k[0] for k in g)<= {0,2}
 H.append(s.expand(g.get((2,),0)*u+g.get((0,),0)))
print('H=',H,flush=True)
A=s.diff(H[0],u);B=H[0].subs(u,0)
res12=s.factor(s.resultant(H[0],H[1],u));res13=s.factor(s.resultant(H[0],H[2],u))
print('r12=',res12,flush=True);print('r13=',res13,flush=True)
r12=s.Poly(res12,c,z).primitive()[1].as_expr();r13=s.Poly(res13,c,z).primitive()[1].as_expr()
start=time.time();res=s.factor(s.resultant(r12,r13,z));print('resultant_c',res,'seconds',time.time()-start,flush=True)
(Path(__file__).resolve().parent / 'elimination.txt').write_text('H='+str(H)+'\nA='+str(A)+'\nB='+str(B)+'\nr12='+str(r12)+'\nr13='+str(r13)+'\nresultant_c='+str(res)+'\n')
# Groebner lex in z,c should yield linear rational-expression z(c) on generic branch.
start=time.time();gb=s.groebner([r12,r13],z,c,order='lex');print('GB secs',time.time()-start,flush=True)
for p in gb.polys:print(p.as_expr(),flush=True)
(Path(__file__).resolve().parent / 'groebner.txt').write_text(str([p.as_expr() for p in gb.polys]))
