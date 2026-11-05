from pathlib import Path
import json,sys,importlib.util
import sympy as s
here=Path(__file__).resolve().parent
root=here.parents[1]
here.mkdir(exist_ok=True)
spec=importlib.util.spec_from_file_location('orders',root/'research/exact-orders/build.py')
mod=importlib.util.module_from_spec(spec);spec.loader.exec_module(mod)
sources=[(21,4,{0:s.sqrt(915),21:s.sqrt(915),4:13*s.sqrt(14),17:-13*s.sqrt(14),8:17*s.sqrt(7),13:17*s.sqrt(7)},4*s.sqrt(663),r'\frac{\sqrt{915}\bigl(|D_{21}^{(0)}\rangle+|D_{21}^{(21)}\rangle\bigr)+13\sqrt{14}\bigl(|D_{21}^{(4)}\rangle-|D_{21}^{(17)}\rangle\bigr)+17\sqrt7\bigl(|D_{21}^{(8)}\rangle+|D_{21}^{(13)}\rangle\bigr)}{4\sqrt{663}}'),(22,5,{1:3*s.sqrt(11),6:s.sqrt(154),11:s.sqrt(119),16:-s.sqrt(154),21:3*s.sqrt(11)},s.Integer(25),r'\frac{3\sqrt{11}\bigl(|D_{22}^{(1)}\rangle+|D_{22}^{(21)}\rangle\bigr)+\sqrt{154}\bigl(|D_{22}^{(6)}\rangle-|D_{22}^{(16)}\rangle\bigr)+\sqrt{119}|D_{22}^{(11)}\rangle}{25}')]
states=[]
for n,t,raw,den,formula in sources:
 c=[s.simplify(raw.get(r,s.S.Zero)/den) for r in range(n+1)]
 assert s.simplify(sum(x*s.conjugate(x) for x in c))==1
 checks=[];next_witness=None
 for q in range(1,t+2):
  delta=(mod.reduction(c,n,q)-s.eye(q+1)/(q+1)).applyfunc(s.simplify)
  zero=delta==s.zeros(q+1);assert zero==(q<=t),(n,q,delta)
  dist=s.simplify(sum(x*s.conjugate(x) for x in delta))
  checks.append({'t':q,'identityVerified':zero,'frobeniusDistanceSquared':str(dist),'purity':str(s.simplify(dist+s.Rational(1,q+1)))})
  if q==t+1:
   a,b=(0,q) if delta[0,q]!=0 else next((a,b) for a in range(q+1) for b in range(q+1) if delta[a,b]!=0)
   next_witness={'order':q,'row':a,'column':b,'differenceFromIdentity':str(delta[a,b]),'differenceLatex':s.latex(delta[a,b])}
 states.append({'id':f'radical-n{n}-t{t}','n':n,'achievedOrder':t,'witnessesOrders':list(range(1,t+1)),'label':f'Exact signed radical state · order {t}','status':'proven-existence','stateLatex':formula,'support':[r for r,x in enumerate(c) if x!=0],'exactCoefficients':[str(x) for x in c],'coefficients':[[float(s.re(x)),float(s.im(x))] for x in c],'majoranaCoordinatePolicy':'numeric-roots-retain-explicit-poles','normalization':'sum_r |c_r|^2 = 1, verified exactly','coefficientConvention':'Normalized Dicke basis; r counts excitations.','references':[{'title':'Twenty-one and twenty-two qubits: exact states and rigorous bounds','url':'./research/n21-n22/index.html'},{'title':'Independent exact reduced-state verification','url':'./research/n21-n22/independent-verification.json'}],'provenance':'Independent analytical construction by rational moment matching and signed coherence cancellation.','notes':f'This state is exactly {t}-AC and is not {t+1}-AC. Optimality is recorded separately from exact attainment.','verification':{'method':'Exact Dicke-basis partial traces through the next order, independently cross-checked using a Schmidt matrix.','exactNormalization':True,'reductions':checks,'exactOrderVerified':True,'nextOrderWitness':next_witness,'globalMaximumOrderProven':False}})
 print(n,t,next_witness,flush=True)
(here/'states.json').write_text(json.dumps({'states':states},indent=2,ensure_ascii=False)+'\n')
