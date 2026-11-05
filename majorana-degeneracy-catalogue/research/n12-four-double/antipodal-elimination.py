from pathlib import Path
import json,time
import sympy as s
h,b,c,d,e,f,t,y=s.symbols('h b c d e f t y',real=True)
D=json.loads(Path(__file__).with_name('antipodal-reduction.json').read_text());E=[s.sympify(v,locals={str(k):k for k in[h,b,c,d,e,f]}) for v in D['conditions']]
sol=s.solve([E[1].subs({d:1,e:1,f:t}),E[2].subs({d:1,e:1,f:t})],[b,c]);print('b',s.factor(sol[b]));print('c',s.factor(sol[c]),flush=True)
v=s.Symbol('v',real=True)
# d=1, e²=v, f²=t²v, ef=tv.
P=[]
for j in[0,3,4]:
 expr=s.Poly(E[j],e,f).terms();out=0
 for (ep,fp),coeff in expr:
  assert(ep+fp)%2==0
  out+=coeff*t**fp*v**((ep+fp)//2)
 out=s.factor(out.subs({b:sol[b],c:sol[c],d:1}))
 P.append(out)
 print('P',j,out,flush=True)
vv=s.solve(P[-1],v)[0];print('v',s.factor(vv),flush=True)
R=[s.factor(s.together(p.subs(v,vv)).as_numer_denom()[0]) for p in P[:2]]
for r in R:print('R',r,flush=True)
RR=[]
for r in R:
 p=s.Poly(r,t);q=0
 for (i,),coef in p.terms():assert i%2==0;q+=coef*y**(i//2)
 RR.append(s.Poly(q,y).primitive()[1].as_expr())
res=s.factor(s.resultant(*RR,y));print('resultant',res,flush=True)
Path(__file__).with_suffix('.json').write_text(json.dumps({'genericAssumptions':'d,e,f nonzero and stated reconstruction denominators nonzero','b':str(s.factor(sol[b])),'c':str(s.factor(sol[c])),'v':str(s.factor(vv)),'residuals':list(map(str,RR)),'resultant':str(res)},indent=2)+'\n')
