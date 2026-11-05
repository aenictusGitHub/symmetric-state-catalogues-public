#!/usr/bin/env python3
from pathlib import Path
import json,math
from copy import deepcopy
import sympy as s

BASE=Path(__file__).resolve().parent

def verify_reduction(n,t,c):
 M=s.zeros(t+1,n-t+1)
 for a in range(t+1):
  for b in range(n-t+1):M[a,b]=c[a+b]*s.sqrt(s.Rational(math.comb(t,a)*math.comb(n-t,b),math.comb(n,a+b)))
 residual=(M*M.conjugate().T-s.eye(t+1)/s.Integer(t+1)).applyfunc(s.simplify)
 assert residual==s.zeros(t+1)

def record(n,t,g,prob,source,location,construction):
 c=[s.sqrt(prob.get(k,s.Integer(0))) for k in range(n+1)]
 assert s.simplify(sum(x*s.conjugate(x) for x in c))==1
 for tt in range(1,t+1):verify_reduction(n,tt,c)
 low=min(prob);high=max(prob);pole=max(low,n-high)
 assert pole==g
 # Pole multiplicity is exact; the matching universal upper bound excludes
 # a larger multiplicity at any other point, independently of factorization.
 formula=' + '.join(s.latex(c[k])+f'|D_{{{n}}}^{{({k})}}\\rangle' for k in prob)
 achieved=3 if (n,t)==(6,2) else t
 if achieved>t:verify_reduction(n,achieved,c)
 return {'id':f'gmax-n{n}-t{t}','n':n,'t':t,'achievedOrder':achieved,'maxMultiplicity':g,'maxMultiplicityVerified':'An exact pole multiplicity g and matching universal moment upper bound.','support':list(prob),'probabilities':{str(k):str(v) for k,v in prob.items()},'exactCoefficients':[str(x) for x in c],'coefficients':[[float(s.re(x)),float(s.im(x))] for x in c],'stateLatex':formula,'construction':construction,'references':[{'file':source,'page':10 if (n,t)==(6,2) else 9 if source.startswith('PRA96') else 4,'location':location}],'verification':{'normalization':'exact','rho_t':'I/(t+1), verified exactly for every order1,...,achievedOrder','poleMultiplicities':{'zero':low,'infinity':n-high}}}

def main():
 states=[]
 for n in range(2,13):
  g=n//2
  if n%2==0:prob={g:s.Integer(1)};loc='Eq. (58)';description='Balanced Dicke state'
  else:prob={0:s.Rational(1,n+1),(n+1)//2:s.Rational(n,n+1)};loc='Eq. (59), r=j-m';description='Two-term odd-N family'
  states.append(record(n,1,g,prob,'PRA96_032304(2017).pdf',loc,description))
 for n in [4,6,7,8,9,10,11,12]:
  g=(n-1)//3
  if n==4:prob={0:s.Rational(1,3),3:s.Rational(2,3)};loc='Eq. (60), g=1';source='PRA96_032304(2017).pdf';description='Two-term tetrahedral representative'
  elif n==6:prob={1:s.Rational(1,2),5:s.Rational(1,2)};loc='Eq. (61); also 3-anticoherent';source='PRA96_032304(2017).pdf';description='Octahedral representative'
  else:
   delta=n-3*g-1;L=n-g
   prob={0:s.Rational(9*(n-2*g)-n*delta,18*L),3:s.Rational(n*delta,18*(L-3)),L:s.Rational(n*(3*(L-3)-delta),6*L*(L-3))};prob={k:v for k,v in prob.items() if v}
   loc='Eqs. (16)–(17), construction independently derived here';source='PRA92_052333(2015).pdf';description='General extremal-degeneracy family for all N≥7'
   if delta==0:loc='Eq. (60), g='+str(g);source='PRA96_032304(2017).pdf';description='Published two-term family (also a special case of the general construction)'
  states.append(record(n,2,g,prob,source,loc,description))
 real=record(12,3,2,{2:s.Rational(7,16),6:s.Rational(1,8),10:s.Rational(7,16)},'PRA92_052333(2015).pdf','Eqs. (16)–(17), construction independently derived here','Double octahedron: six distinct points, each of multiplicity two')
 real['majoranaPolynomialLatex']='\\frac{\\sqrt{462}}{4}z^2(z^4+1)^2';real['multiplicityPartition']=[2]*6
 z=s.symbols('z');c=list(map(s.sympify,real['exactCoefficients']));poly=s.expand(sum((-1)**r*s.sqrt(s.binomial(12,r))*c[r]*z**r for r in range(13)))
 assert s.simplify(poly-s.sqrt(462)*z**2*(z**4+1)**2/4)==0
 assert s.gcd(s.Poly((z**4+1)**2,z),s.Poly(s.diff((z**4+1)**2,z),z)).monic().as_expr()==z**4+1
 real['verification']['exactPolynomialFactorization']=True;states.append(real)
 imaginary=deepcopy(real);imaginary['id']+='-phase-i';c[6]*=s.I
 for tt in range(1,4):verify_reduction(12,tt,c)
 imaginary['exactCoefficients']=[str(x) for x in c];imaginary['coefficients']=[[float(s.re(x)),float(s.im(x))] for x in c]
 imaginary['stateLatex']=' + '.join(s.latex(c[k])+f'|D_{{12}}^{{({k})}}\\rangle' for k in [2,6,10]);imaginary['construction']='Phase variant with two double points and eight simple points'
 imaginary['majoranaPolynomialLatex']='\\frac{\\sqrt{462}}{4}z^2(z^8+2iz^4+1)';imaginary['multiplicityPartition']=[2,2]+[1]*8
 q=s.Poly(z**8+2*s.I*z**4+1,z,extension=s.I);assert s.gcd(q,q.diff()).degree()==0
 imaginary['verification']['nonzeroRootPolynomialSquarefree']=True;states.append(imaginary)
 result={'scope':'Exact representatives attaining maximal Majorana-point multiplicity at orders1 and2 for all applicable n2..12, plus two N12 order3 examples with different multiplicity partitions. No 2-AC state exists at n2,3,5.','basisConvention':'Normalized Dicke coefficients; r counts excitations. Source spin coefficients use r=j-m. Majorana polynomial P(z)=sum_r (-1)^r sqrt(C(n,r)) c_r z^r, with n-degree(P) roots at infinity. Multiplicity is unaffected by the pole convention.','states':states}
 (BASE/'maximal-multiplicity-witnesses.json').write_text(json.dumps(result,indent=2)+'\n')
 print('Exact reduced-matrix checks passed for',len(states),'witnesses.')

if __name__=='__main__':main()
