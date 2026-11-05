from pathlib import Path
import sympy as s,json
h,b,c,d,e,f=s.symbols('h b c d e f',real=True)
E=[s.sympify(v,locals={str(k):k for k in[h,b,c,d,e,f]}) for v in json.loads(Path(__file__).with_name('antipodal-reduction.json').read_text())['conditions']]
# e=0,f !=0 => d=0. c from remaining offdiagonal.
C=s.solve(E[2].subs({e:0,d:0,f:1}),c)[0]
print('e0 c=',s.factor(C))
polys=[s.factor(s.together(E[j].subs({e:0,d:0,c:C})).as_numer_denom()[0]) for j in[0,3]]
M=s.Matrix([[s.expand(p).coeff(b,2),s.expand(p).coeff(f,2)] for p in polys]);print('e0 det',s.factor(M.det()),flush=True)
# e=f=0, d!=0 => c=-2bh/(h²+2)
C2=-2*b*h/(h*h+2)
polys2=[s.factor(s.together(E[j].subs({e:0,f:0,c:C2})).as_numer_denom()[0]) for j in[0,3]]
M2=s.Matrix([[s.expand(p).coeff(b,2),s.expand(p).coeff(d,2)] for p in polys2]);print('ef0 dnonzero det',s.factor(M2.det()),flush=True)
# e=f=d=0; b can be set1 unlesszero.
polys3=[s.factor(E[j].subs({e:0,f:0,d:0,b:1})) for j in[0,3]]
res3=s.factor(s.resultant(*polys3,c));print('even real resultant',res3,flush=True)
print('even real h0',list(map(lambda P:s.factor(P.subs(h,0)),polys3)))
for val in[-1,0,1]:
 print('even real',val,s.solve([P.subs(h,val) for P in polys3],c),flush=True)
Path(__file__).with_suffix('.json').write_text(json.dumps({'e_zero_branch':{'c':str(s.factor(C)),'determinant':str(s.factor(M.det())),'matrix':str(M)},'e_f_zero_d_nonzero_branch':{'c':str(C2),'determinant':str(s.factor(M2.det())),'matrix':str(M2)},'even_real_branch':{'polynomials':list(map(str,polys3)),'resultant':str(res3)}},indent=2)+'\n')
