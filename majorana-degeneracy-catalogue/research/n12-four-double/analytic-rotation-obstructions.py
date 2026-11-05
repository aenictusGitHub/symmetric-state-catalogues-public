#!/usr/bin/env python3
"""Exact C3 and C4 exclusions for twelve-qubit partition (2^4,1^4).

Does not assume all constellations have either symmetry, and does not settle
arbitrary C1, C2 or D2 configurations. Requires SymPy only.
"""
from pathlib import Path
import json
import sympy as s
HERE=Path(__file__).resolve().parent

def moments(support):
    w=s.symbols('w:'+str(len(support)))
    equations=[sum(v*r**m for v,r in zip(w,support))-sum(s.Integer(j)**m for j in range(13))/13 for m in range(4)]
    solutions=s.solve(equations,w,dict=True);assert len(solutions)==1
    return [solutions[0][v] for v in w]

c3=[]
for support in [[2,5,8,11],[1,4,7,10]]:
    probabilities=moments(support)
    assert s.Rational(-1,27) in probabilities
    c3.append({'support':support,'forcedProbabilities':list(map(str,probabilities))})
c4prob=moments([0,4,8,12]);assert c4prob==[s.Rational(5,32),s.Rational(11,32),s.Rational(11,32),s.Rational(5,32)]
ratio=s.sqrt(s.binomial(12,4)*c4prob[1]/c4prob[0]);assert ratio==33
r,c,h=s.symbols('r c h',real=True)
# After a phase rotation set u=r>0. From |u²v|=1, v=r^-2 exp(i phi).
# A=-(2u+v), B=u²+2uv are the monic cubic's middle coefficients.
absA2=4*r*r+r**-4+4*c/r
absB2=r**4+4/r**2+4*r*c
factor=(r*r-1)*(r**6-3*r**4-3*r*r+1+4*r**3*c)/r**4
assert s.simplify(absB2-absA2-factor)==0
cvalue=-(r**6-3*r**4-3*r*r+1)/(4*r**3)
assert s.simplify(cvalue-((r+1/r)*(6-(r+1/r)**2)/4))==0
assert s.simplify(absA2.subs(c,cvalue)-3*((r+1/r)**2-1))==0
assert s.expand(h**3-6*h-4-(h+2)*(h*h-2*h-2))==0
maximum=s.simplify(3*((1+s.sqrt(3))**2-1));assert maximum==9+6*s.sqrt(3)
assert maximum<25<ratio**2
out={
 'verifiedExactly':True,'scope':'Pure symmetric12-qubit states with full partition(2,2,2,2,1,1,1,1).',
 'globalExclusionProved':False,
 'C3':{'excluded':True,'orbitCounting':'Four doubled points require one three-point orbit plus one polar double point. Four simple points require one three-point orbit plus the other polar simple point. The resulting Dicke support is one of the two listed sets.',
       'momentSystems':c3,'contradiction':'Uniform moments through order3 force a negative probability.'},
 'C4':{'excluded':True,'orbitCounting':'Both the doubled set and the simple set must be single four-point orbits; no polar point is possible. Thus P(z)=(z^4-u)^2(z^4-v), with u,v nonzero and distinct.',
       'support':[0,4,8,12],'forcedProbabilities':list(map(str,c4prob)),
       'requiredMonicMiddleCoefficientMagnitude':'33',
       'doubleRootCubicLemma':'For (w-u)^2(w-v), |u²v|=1 and equal magnitudes of the two middle coefficients imply each squared magnitude is at most9+6sqrt(3).',
       'lemmaDerivation':'Rotate u to r>0 and write v=r^-2 exp(i phi). Equal magnitudes factor as (r²-1)[r6-3r4-3r2+1+4r3cos(phi)]=0. If r=1 the magnitude is at most3. Otherwise h=r+1/r>=2 obeys cos(phi)=h(6-h²)/4, hence h<=1+sqrt3, and the squared middle magnitude is3(h²-1)<=9+6sqrt3<25.',
       'contradiction':'The required magnitude33 exceeds the universal bound, which is less than5.'},
 'ordersAtLeast5':{'excluded':True,'argument':'A nonaxial orbit under an order-m proper rotation has exactlym points. With onlyfour doubled points and at mosttwo axial points, m>=5 is impossible.'},
 'remainingProperRotationGroups':['C1','C2','D2'],
 'groupClassificationReason':'Finite subgroups of SO(3) are cyclic, dihedral or polyhedral. If no element has order>=3, only C1,C2,D2 remain. Continuous rotation symmetry is impossible for eight distinct points.',
 'limitation':'The target partition remains unsettled for general geometries in the three remaining proper-group classes. No numerical failure is used as a global exclusion.'}
(HERE/'analytic-rotation-obstructions.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps({'verifiedExactly':True,'excludedOrders':['3','4','at least5'],'remainingProperGroups':out['remainingProperRotationGroups'],'globalExclusionProved':False}))
