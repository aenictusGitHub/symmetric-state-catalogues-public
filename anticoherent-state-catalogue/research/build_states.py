"""Extract/source and exactly verify full anticoherence witnesses.

All matrices use the normalized Dicke basis. Full t-AC means rho_t=I/(t+1),
not vanishing even multipoles alone. Sources use r=j-m to map spin to Dicke.
"""
from pathlib import Path
import json,math
import sympy as s
R=Path(__file__).resolve().parent

def reduced(c,n,t):
 M=s.Matrix(t+1,n-t+1,lambda a,b:c.get(a+b,0)*s.sqrt(s.binomial(t,a)*s.binomial(n-t,b)/s.binomial(n,a+b)))
 return (M*M.conjugate().T).applyfunc(s.simplify)
def exact_order(c,n,claimed):
 norm=s.simplify(sum(v*s.conjugate(v) for v in c.values()))
 assert norm==1,(n,norm)
 checks=[]
 for t in range(1,min(n//2,claimed+1)+1):
  rho=reduced(c,n,t);delta=(rho-s.eye(t+1)/s.Integer(t+1)).applyfunc(s.simplify)
  residual=s.simplify(sum(x*s.conjugate(x) for x in delta))
  checks.append({'t':t,'identityVerified':residual==0,'frobeniusDistanceSquared':str(residual),'purity':str(s.simplify(s.trace(rho*rho)))})
  if t<=claimed:assert residual==0,(n,t,residual)
  else:assert residual!=0,(n,'has higher order',t)
 return checks
def ref(file,page,location):return {'file':file,'page':page,'location':location}
def state(n,t,label,c,references,provenance='published exact construction',notes='',id=None):
 checks=exact_order(c,n,t)
 coeff=[c.get(r,s.S(0)) for r in range(n+1)]
 ket=lambda r:rf'|D_{{{n}}}^{{({r})}}\rangle'
 def term(r,v):
  coefficient=s.latex(v)
  if v.is_Add:coefficient=r'\left('+coefficient+r'\right)'
  return coefficient+ket(r)
 formula=' + '.join(term(r,v) for r,v in sorted(c.items())).replace('+ -','- ')
 out={'id':id or f'ac-n{n}-t{t}','n':n,'achievedOrder':t,'witnessesOrders':list(range(1,t+1)),'label':label,'status':'proven-existence','stateLatex':formula,'support':sorted(c),'exactCoefficients':[str(v) for v in coeff],'coefficients':[[float(s.re(v)),float(s.im(v))] for v in coeff],'normalization':'sum_r |c_r|^2 = 1, verified exactly','coefficientConvention':'Normalized Dicke basis |D_n^(r)>; r counts excitations. Source spin states are mapped with r=j-m.','references':references,'provenance':provenance,'notes':notes,'verification':{'method':'Exact SymPy evaluation of each t-body reduced matrix in Sym^t(C^2); rho_t=I_(t+1)/(t+1).','exactNormalization':True,'reductions':checks,'globalMaximumOrderProven':False}}
 print(n,t,label,'exactly verified',flush=True)
 return out

def main():
 states=[]
 paper17='PRA96_032304(2017).pdf';paper15='PRA92_052333(2015).pdf';paper20='Quantum4_285(2020).pdf'
 for n in (2,3,5):states.append(state(n,1,'GHZ / spin-cat state',{0:1/s.sqrt(2),n:1/s.sqrt(2)},[ref(paper17,2,'Eq. (1)')]))
 states.append(state(4,2,'Tetrahedral state',{0:s.Rational(1,2),2:s.I/s.sqrt(2),4:s.Rational(1,2)},[ref(paper17,2,'Eq. (2)'),ref('PRA90_032314(2014).pdf',4,'Eq. (25) at mu=i sqrt(2); uniqueness discussion continues on PDF p5')]))
 states.append(state(6,3,'Octahedral state',{1:1/s.sqrt(2),5:1/s.sqrt(2)},[ref(paper17,10,'Eq. (61)'),ref(paper15,9,'C4-symmetric six-qubit construction')]))
 states.append(state(7,2,'Two-term cyclic construction',{0:s.sqrt(s.Rational(3,10)),5:s.sqrt(s.Rational(7,10))},[ref(paper17,9,'Eq. (60), j=7/2')]))
 for n in (8,10):
  states.append(state(n,3,'Even-qubit cyclic construction',{0:s.sqrt(s.Rational(n+2,6*n)),n//2:s.sqrt(s.Rational(2*(n-1),3*n)),n:s.sqrt(s.Rational(n+2,6*n))},[ref(paper15,9,'Eq. (48); paragraph immediately after gives order 3 for even N>=8')]))
 states.append(state(9,2,'Small-angle rotosensor, chi=0',{1:-s.sqrt(21)/8,5:s.sqrt(30)/8,9:s.sqrt(13)/8},[ref(paper20,11,'Eq. (65), chi=0')],notes='The same source explicitly reports no full t>=3 anticoherent states at j=9/2. The q=3 label for a different N=9 state in the 2022 Table IV is inconsistent with its nonzero rank-3 multipole.'))
 endpoint=(-35*s.sqrt(10)+s.I*s.sqrt(10022))/288
 states.append(state(11,3,'Four-term cyclic construction',{1:endpoint,4:5*s.sqrt(3)/18,7:5*s.sqrt(3)/18,10:endpoint},[ref(paper15,4,'Dicke-coefficient characterization, Eqs. (16)–(17)')],provenance='Independently derived exact construction using the published characterization',notes='This radical witness is derived here, not copied from the rounded numerical N=11 vector in the 2017 Appendix. The two outer probabilities are 29/108 and the two inner probabilities are 25/108. Their phases cancel the only possible off-diagonal third-order constraint exactly.'))
 states.append(state(12,4,'Four-anticoherent cyclic family, theta=0',{1:s.sqrt(7)/5,6:s.sqrt(11)/5,11:s.sqrt(7)/5},[ref(paper17,10,'Eq. (62), theta=0')],id='ac-n12-t4-family',notes='Eq. (62) supplies a continuous family. Its theta=pi/2 endpoint is 5-anticoherent.'))
 states.append(state(12,5,'Icosahedral state',{1:s.sqrt(7)/5,6:s.I*s.sqrt(11)/5,11:s.sqrt(7)/5},[ref(paper17,10,'Eq. (62), theta=pi/2; Sec. IV.B.5')]))
 # The local mixed-state AC manuscript explicitly supplies these pure inputs.
 local={
  4:(2,{0:1,2:s.I*s.sqrt(2),4:1},'C2'),
  5:(1,{1:1,4:1},'C3'),
  6:(3,{1:1,5:1},'C4'),
  7:(2,{0:2,3:-s.sqrt(7),6:-s.sqrt(7)},'C5'),
  8:(3,{0:s.sqrt(5),4:-s.sqrt(14),8:-s.sqrt(5)},'C6'),
  9:(2,{0:s.sqrt(13),4:s.sqrt(30),8:-s.sqrt(21)},'C7'),
  10:(3,{1:s.sqrt(5),5:s.sqrt(6),9:-s.sqrt(5)},'C8')}
 for n,(t,un,equation) in local.items():
  if n in (4,6):
   continue
  norm=s.sqrt(sum(v*s.conjugate(v) for v in un.values()));c={r:s.simplify(v/norm) for r,v in un.items()}
  states.append(state(n,t,'Alternative exact representative',c,[{'title':'Independent exact reduced-state verification','url':'./research/states.json','location':f'ac-n{n}-t{t}-manuscript'}],provenance='Normalization and anticoherence are verified independently from the saved coefficients.',id=f'ac-n{n}-t{t}-manuscript',notes='The saved coefficients are divided by their exact norm; every reduction through the achieved AC order is checked independently.'))
 for st in states:
  if st['id'].endswith('-manuscript'):st['internalProvenance']='Exact construction explicitly listed in the supplied local manuscript'
 # Two exact larger source constructions make the source's arbitrary-order
 # existence discussion concrete without claiming an exhaustive larger table.
 n=24
 un={0:s.sqrt((n+2)*(n+4)*(7*n-4)),n//4:4*s.sqrt(2*(n-2)*(n-1)*(n+2)),n//2:2*s.sqrt(3*(n-1)*(n*n+16)),3*n//4:4*s.sqrt(2*(n-2)*(n-1)*(n+2)),n:s.sqrt((n+2)*(n+4)*(7*n-4))}
 norm=s.sqrt(sum(v*v for v in un.values()));c={r:s.simplify(v/norm) for r,v in un.items()}
 states.append(state(n,5,'Five-anticoherent cyclic family at N=24',c,[ref(paper15,10,'Eqs. (51)–(52), m=5')]))
 c={0:s.sqrt(7062)/343,7:s.sqrt(29315)/343,14:3*s.sqrt(451)/343,21:s.sqrt(36777)/343,28:-3*s.sqrt(451)/343,35:s.sqrt(29315)/343,42:-s.sqrt(7062)/343}
 states.append(state(42,7,'Seven-anticoherent D7d construction',c,[ref(paper15,10,'Figure 5 caption')]))
 states.sort(key=lambda e:(e['n'],e['achievedOrder']))
 output={'definition':'A pure symmetric N-qubit state is t-anticoherent iff rho_t is the maximally mixed state I_(t+1)/(t+1) on the symmetric t-qubit subspace. Equivalently all spin moments through order t are independent of direction. This requires every multipole rank1,...,t to vanish, not only the even ranks.','definitionReferences':[ref(paper15,3,'Eq. (10) and anticoherence definition'),ref(paper17,3,'Sec. II.B.2, Eq. (13)')],'states':states,'sourceCaveats':[{'file':'PhysRevResearch4_013178(2022).pdf','page':17,'location':'Table IV, j=9/2 row','issue':'The printed q=3 label conflicts with both the listed nonzero rank3 multipole and the explicit state. The normalized printed state has achieved order2, not3: rho_3[0,3]=(10-sqrt(42))/126 is nonzero.'},{'file':'PhysRevResearch4_013178(2022).pdf','page':17,'location':'Table IV, j=5 row','issue':'Printed endpoint amplitudes1/sqrt(11) with middle amplitude sqrt(3/5) have norm squared43/55. Use the normalized exact construction from PRA92(2015) Eq48; do not silently quote the erroneous row.'}]}
 (R/'states.json').write_text(json.dumps(output,indent=2)+'\n')
 print('Wrote',len(states),'exactly verified witnesses.',flush=True)

if __name__=='__main__':main()
