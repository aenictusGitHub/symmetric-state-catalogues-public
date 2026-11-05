#!/usr/bin/env python3
"""Exact D2 exclusion for n=12 and full partition(2^4,1^4).

Uses the separately verified antipodal-doubled-point theorem for two size-two
orbits; this file covers every size-four doubled orbit and every residual D2
character, including the vanishing-leading-coefficient branch. SymPy only.
"""
from pathlib import Path
import json
import sympy as s
HERE=Path(__file__).resolve().parent
A,B,X,Y,U,V,z=s.symbols('A B X Y U V z',real=True)
a=A+s.I*B;b=X+s.I*Y
Q=z**4+a*z*z+1
L=2*a+b;M=a*a+2*a*b+3;N=a*a*b+4*a+2*b
F=s.expand(s.re(33*L+3*M*s.conjugate(L)+N*s.conjugate(M)))
H=s.expand(s.im(165*L+9*M*s.conjugate(L)+N*s.conjugate(M)))
G=s.expand(8712+12*L*s.conjugate(L)-8*M*s.conjugate(M)-3*N*s.conjugate(N))
# Independently check these as necessary spin-moment equations.
p=s.Poly(s.expand(Q**2*(z**4+b*z*z+1)),z)
weights=[s.expand(p.nth(r)*s.conjugate(p.nth(r)))/s.binomial(12,r)for r in range(13)]
zz=s.expand(sum(((r-6)**2-14)*weights[r]for r in range(13)))
pp=s.expand(sum(s.Rational((r+1)*(r+2),s.binomial(12,r))*p.nth(r+2)*s.conjugate(p.nth(r))for r in range(11)))
odd=s.expand(sum((r-5)*s.Rational((r+1)*(r+2),s.binomial(12,r))*p.nth(r+2)*s.conjugate(p.nth(r))for r in range(11)))
assert s.expand(198*zz-G)==0
assert s.expand(s.Rational(33,4)*pp-F)==0
assert s.expand(-s.Rational(33,4)/s.I*odd-H)==0
# Residual two-dimensional density: u=|beta|²,v=|gamma|²,
# w=Re(beta*conjugate(gamma)),t with the matching imaginary convention.
C=s.Matrix([[e.subs({X:0,Y:0}),e.coeff(X,2),e.coeff(X,1).subs(Y,0),e.coeff(Y,1).subs(X,0)]for e in[F,H,G]])
for row,e in enumerate([F,H,G]):
 assert s.expand(e-(C[row,0]+C[row,1]*(X*X+Y*Y)+C[row,2]*X+C[row,3]*Y))==0
# beta=0: u=w=t=0 andv>0. Its firsttwo equations forceA=B=0.
assert s.factor(C[0,1])==2*A*(A*A+B*B+5)
assert s.factor(C[1,1])==2*B*(A*A+B*B+7)
assert C[:,1].subs({A:0,B:0})==s.zeros(3,1)
# beta!=0: normalize beta=1. a=0 is immediately inconsistent.
assert G.subs({A:0,B:0})==8640
K=C[:,1:];D=s.factor(K.det())
DP=s.Poly(D,A,B)
assert all(i%2==j%2==0 and coefficient<0 for(i,j),coefficient in DP.terms())
assert DP.coeff_monomial(A*A)<0 and DP.coeff_monomial(B*B)<0
# Hence D<0 unless A=B=0, already handled.
nums=[]
for index in range(3):
 T=K.copy();T[:,index]=-C[:,0];nums.append(s.factor(T.det()))
minor=s.factor(nums[0]*D-nums[1]**2-nums[2]**2)
P=s.Poly(-minor/12,A,B)
assert all(i%2==j%2==0 for(i,j),_ in P.terms())
positive=s.Poly(sum(coefficient*U**(i//2)*V**(j//2)for(i,j),coefficient in P.terms()),V)
constant=s.factor(positive.nth(0))
expected=U*(U-36)**2*(U+12)**6*(2*U**3+51*U*U+1032*U+3120)
assert s.expand(constant-expected)==0
coefficient_checks=[]
for power in range(1,positive.degree()+1):
 polynomial=s.Poly(positive.nth(power),U)
 assert polynomial.eval(0)>0 and polynomial.count_roots(0,s.oo)==0
 coefficient_checks.append({'powerOfBSquared':power,'polynomialInASquared':str(polynomial.as_expr()),'atZero':str(polynomial.eval(0)),'nonnegativeHalfLineRootCount':0})
# Pure residual states need v=w²+t². The minor is strictly negative forB!=0;
# atB=0 it vanishes only atA=0,+6,-6. a0 was inconsistentinthischart.
survivors=[]
for av in [6,-6]:
 subs={A:av,B:0};vval,wval,tval=[s.simplify(num.subs(subs)/D.subs(subs))for num in nums]
 assert vval==4 and wval==-s.Rational(av,3) and tval==0
 polynomial=s.Poly(s.expand((z**4+av*z*z+1)**2*(z**4+wval*z*z+1)),z)
 partition=sorted([m for factor,m in polynomial.sqf_list()[1]for _ in range(factor.degree())],reverse=True)
 assert partition==[2]*6
 survivors.append({'a':str(av),'beta':'1','gamma':str(wval),'polynomial':str(polynomial.as_expr()),'partition':partition})
survivors.append({'a':'0','beta':'0','gamma':'1','polynomial':'z^2*(z^4+1)^2','partition':[2]*6})
# The three nontrivial character spaces are one-dimensional.
characters=[]
for R in [z**4-1,z*(z*z-1),z*(z*z+1)]:
 poly=s.Poly(s.expand(Q**2*R),z);pc=[poly.nth(i)for i in range(13)]
 rho={}
 for aa in range(4):
  for cc in range(4):rho[aa,cc]=s.expand(sum(s.binomial(9,q)*pc[aa+q]*s.conjugate(pc[cc+q])/(s.binomial(12,aa+q)*s.binomial(12,cc+q))for q in range(10)))
 equations=[s.binomial(3,aa)*rho[aa,aa]-rho[0,0]for aa in range(1,4)]
 for aa in range(4):
  for cc in range(aa+1,4):equations +=[s.re(rho[aa,cc]),s.im(rho[aa,cc])]
 equations=[s.Poly(s.expand(e),A,B).clear_denoms()[1].primitive()[1].as_expr()for e in equations if e!=0]
 basis=s.groebner(equations,A,B)
 assert list(basis)==[1]
 characters.append({'R':str(R),'equations':list(map(str,equations)),'exactGroebnerBasis':['1']})
output={
 'verifiedExactly':True,'scope':'No pure symmetric12-qubit3AC state of full partition(2,2,2,2,1,1,1,1) has a proper D2 rotational subgroup.',
 'globalExclusionProved':False,
 'normalFormCompleteness':'A D2 orbit on the sphere has size4 except the three size2 antipodal-axis orbits. Four doubled points are either two size2 orbits (excluded by the antipodal theorem) or one size4 orbit {u,-u,1/u,-1/u}. The latter has Q=z^4+a*z^2+1 with arbitrary complex a. Q² has trivial character. All residual quartic character spaces are span{z4+1,z2}, span{z4-1}, span{z(z²-1)}, span{z(z²+1)}.',
 'dependency':'antipodal-obstruction-verification.json',
 'trivialCharacter':{
   'variables':'a=A+iB; beta!=0 normalized to 1; gamma=X+iY. Both A,B and X,Y are arbitrary real.',
   'necessaryMomentEquations':[str(F),str(H),str(G)],'densityMatrixConstraint':[[str(v)for v in row]for row in C.tolist()],
   'recoveryDenominator':str(D),'denominatorStrictlyNegativeExceptAEqualsBEqualsZero':True,
   'densityRecoveryNumerators':list(map(str,nums)),
   'pureStateMinorNumerator':'-12*P(A^2,B^2)',
   'PAtBSquaredZero':str(constant),'positiveCoefficientChecks':coefficient_checks,
   'conclusion':'For B!=0 the residual2x2 density determinant is strictly negative. For B=0 pure states require A=±6 (A=0 is inconsistent in the beta=1 chart); all survivors have partition 2^6.',
   'vanishingBetaBranch':'beta=0 forces A=B=0 and gamma nonzero, giving z²(z⁴+1)². It is retained and has partition 2^6.'},
 'nontrivialCharacters':characters,'survivingDegenerateStates':survivors,
 'remainingPossibleProperGroupsAfterC3C4Restrictions':['C1','C2'],
 'limitation':'This excludes D2, not arbitrary asymmetric or solely C2-symmetric constellations.'}
(HERE/'analytic-d2-obstruction.json').write_text(json.dumps(output,indent=2)+'\n')
print(json.dumps({'verifiedExactly':True,'excludedSubgroup':'D2','onlySurvivingPartition':[2]*6,'remainingPossibleProperGroups':['C1','C2'],'globalExclusionProved':False}))
