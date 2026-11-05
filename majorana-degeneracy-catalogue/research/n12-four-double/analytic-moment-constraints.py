#!/usr/bin/env python3
"""Exact necessary constraints at double Majorana points of a 12-qubit 3-AC state.

These strengthen the search domain but do not exclude (2^4,1^4) globally.
Requires SymPy; output is analytic-moment-constraints.json in this directory.
"""
import json,math
from pathlib import Path
import sympy as s
HERE=Path(__file__).resolve().parent
r,q=s.symbols('r q', real=True)
# Rotate any double point to the root z=0, so c0=c1=0.
# Three-AC fixes the first three r moments to those of uniform {0,...,12}.
mean=lambda f:s.simplify(sum(f.subs(r,j) for j in range(13))/13)
L=s.expand((r-3)*(r-9)*(r-10))
U=s.expand(-(r-8)*(r-9)*(r-12))
assert L.subs(r,2)==-56 and mean(L)==-20
assert all(L.subs(r,j)>=0 for j in range(3,13))
assert U.subs(r,2)==420 and mean(U)==190
assert all(U.subs(r,j)>=0 for j in range(3,13))
lower=s.Rational(5,14);upper=s.Rational(19,42)
assert -mean(L)/56==lower and mean(U)/420==upper

def endpoint_weights(support):
    weights=s.symbols('w:'+str(len(support)))
    solution=s.solve([sum(w*j**d for w,j in zip(weights,support))-mean(r**d)
                      for d in range(4)],weights,dict=True)
    assert len(solution)==1
    values=[solution[0][w] for w in weights]
    assert all(v>0 for v in values)
    return {str(j):str(v) for j,v in zip(support,values)}
low_support=[2,3,9,10];high_support=[2,8,9,12]
low_weights=endpoint_weights(low_support);high_weights=endpoint_weights(high_support)
# At the lower endpoint only adjacent pairs 2->3 and 9->10 remain.
# J+ and {Jz,J+} give A+B=0 and -7A+7B=0. Both A,B are nonzero.
assert [2*j+1-12 for j in [2,9]]==[-7,7]
assert [s.sqrt((12-j)*(j+1)) for j in [2,9]]==[s.sqrt(30)]*2
assert s.Matrix([[1,1],[-7,7]]).det()==14
# At the upper endpoint only 8->9 contributes to J+ and is nonzero.
assert s.sqrt((12-8)*(8+1))==6
# For two double roots at0,t>0, put q=t². The constrained coefficient
# vector is orthogonal to g(t),g'(t), where P(z)=z²g(z).
A=sum(s.binomial(12,j)*q**(j-2) for j in range(2,13))
D=s.diff(A,q)+q*s.diff(A,q,2)
delta=s.expand(A*D-q*s.diff(A,q)**2)
# Maximum |c2|² in this vanishing subspace is 1-66D/delta.
separation=s.Poly(s.expand(9*delta-924*D),q).primitive()[1]
assert separation.degree()==18
coeffs=separation.all_coeffs();signs=[s.sign(v) for v in coeffs if v]
assert sum(a!=b for a,b in zip(signs,signs[1:]))==1
assert separation.eval(0)<0 and separation.LC()>0
interval=separation.refine_root(s.Rational(3237,10000),s.Rational(3238,10000),eps=s.Rational(1,10**30))
assert separation.count_roots(*interval)==1
root=s.nsolve(separation.as_expr(),q,s.Rational(324,1000),prec=60)
angle=2*math.atan(math.sqrt(float(root)))
output={
 'verifiedExactly':True,
 'scope':'Necessary conditions for every complex pure symmetric twelve-qubit third-order anticoherent state at any double Majorana point. These do not globally exclude partition (2,2,2,2,1,1,1,1).',
 'globalExclusionProved':False,
 'doublePointProbability':{
   'definition':'After rotating a double root to z=0, c0=c1=0 and p2=|c2|².',
   'strictLower':str(lower),'strictUpper':str(upper),
   'lowerPolynomial':str(L),'lowerUniformExpectation':str(mean(L)),
   'upperPolynomial':str(U),'upperUniformExpectation':str(mean(U)),
   'endpointSupports':{'lower':low_weights,'upper':high_weights},
   'lowerEndpointExcluded':'Only c2*c3bar and c9*c10bar contribute to J+. Their nonzero magnitudes and the independent {Jz,J+} equation force an impossible simultaneous A+B=0 and -7A+7B=0.',
   'upperEndpointExcluded':'The sole adjacent support pair8,9 gives a nonzero J+ expectation.'},
 'doublePointSeparation':{
   'coordinate':'q=tan(theta/2)^2, theta is the angular separation of two distinct double points.',
   'coefficientKernelA':str(A),'coefficientKernelD':str(D),'gramDeterminant':str(delta),
   'maximumP2':'1-66*D/(A*D-q*(dA/dq)^2)',
   'positivityJustification':'A*D-q*(A prime)^2 is a Gram determinant of two independent evaluation/derivative vectors for every q>0.',
   'definingPolynomial':str(separation.as_expr()),'uniquePositiveRoot':True,
   'isolatingInterval':list(map(str,interval)),
   'qLowerBoundApproximation':str(root),
   'angularThresholdRadiansApproximation':angle,'angularThresholdDegreesApproximation':angle*180/math.pi,
   'exactAngularThreshold':'2*atan(sqrt(q_star)), where q_star is the unique positive root of definingPolynomial.',
   'numericalScope':'Floating angles approximate the exact threshold; they are not asserted to be outward-rounded interval bounds. Antipodal pairs satisfy the separation bound automatically.',
   'conclusion':'Every pair of distinct double points has angular separation strictly greater than this bound.'},
 'interpretation':'The repeated points cannot themselves approach one another in an exact3AC state. Simple points may still coalesce into additional double points; the double-octahedral boundary therefore remains possible.'}
(HERE/'analytic-moment-constraints.json').write_text(json.dumps(output,indent=2)+'\n')
print(json.dumps({'verifiedExactly':True,'strictP2Bounds':[str(lower),str(upper)],'minimumDoublePointAngleDegrees':angle*180/math.pi,'globalExclusionProved':False}))
