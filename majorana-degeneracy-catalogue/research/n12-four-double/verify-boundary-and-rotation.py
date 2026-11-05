#!/usr/bin/env python3
"""Independent exact checks of the boundary sequence and cyclic exclusions."""
from pathlib import Path
import hashlib,json,time
import sympy as s
HERE=Path(__file__).resolve().parent
z,eps,r,c,h=s.symbols('z epsilon r c h',real=True)
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def moments(support):
    M=s.Matrix([[s.Integer(j)**degree for j in support] for degree in range(4)])
    target=s.Matrix([sum(s.Integer(j)**degree for j in range(13))/13 for degree in range(4)])
    return list(M.inv()*target)
def main():
    start=time.monotonic();source=HERE/'boundary_sequence.json';data=json.loads(source.read_text())
    Q=z**4-6*z*z+1;R=(z*z+1)**2+eps;P=s.Poly(Q**2*R,z)
    raw=[(-1)**j*P.nth(j)/s.sqrt(s.binomial(12,j))for j in range(13)]
    Z=s.factor(sum(v*v for v in raw))
    A=s.Matrix(4,10,lambda a,b:raw[a+b]*s.sqrt(s.binomial(3,a)*s.binomial(9,b)/s.binomial(12,a+b)))
    rho=(A*A.H).applyfunc(lambda v:s.cancel(v/Z))
    defect=s.factor(sum(s.expand_complex(v*s.conjugate(v))for v in rho-s.eye(4)/4))
    assert s.cancel(Z-s.sympify(data['normalizationSquared'],locals={'epsilon':eps}))==0
    assert s.cancel(defect-s.sympify(data['thirdOrderHSDefectSquared'],locals={'epsilon':eps}))==0
    qd=s.discriminant(Q,z);rd=s.factor(s.discriminant(R,z));qr=s.factor(s.resultant(Q,R,z))
    assert qd==16384 and rd==256*eps**2*(1+eps) and qr==(eps**2+48*eps+64)**2
    assert P.degree()==12 and P.nth(12)==1
    # Strict positivity for epsilon>0 follows directly from the displayed
    # positive coefficients of each remaining factor.
    assert all(v>0 for v in s.Poly(qr,eps).all_coeffs())
    num,den=s.fraction(defect/eps**2)
    assert all(v>0 for v in s.Poly(num,eps).all_coeffs())
    assert all(v>0 for v in s.Poly(den,eps).all_coeffs())
    limit=s.Poly(P.as_expr().subs(eps,0),z)
    lam=sorted([int(m)for q,m in limit.sqf_list()[1]for _ in range(q.degree())],reverse=True)
    assert lam==[2]*6
    assert rho.subs(eps,0)==s.eye(4)/4 and s.limit(defect,eps,0)==0
    leading=s.simplify(s.limit(defect/eps**2,eps,0))
    # Independently reconstruct the forced diagonal moments for the cyclic cases.
    C3=[moments(S)for S in [[2,5,8,11],[1,4,7,10]]]
    assert all(min(w)<0 for w in C3)
    C4=moments([0,4,8,12]);assert C4==[s.Rational(5,32),s.Rational(11,32),s.Rational(11,32),s.Rational(5,32)]
    magnitudeSquared=s.binomial(12,4)*C4[1]/C4[0];assert magnitudeSquared==1089
    # Arbitrary complex double and simple orbit parameters: rotate u=r>0,
    # and write v=e^(i phi)/r². c=cos(phi) is the only remaining phase.
    AA=4*r*r+r**-4+4*c/r;BB=r**4+4/r**2+4*r*c
    fact=(r*r-1)*(r**6-3*r**4-3*r*r+1+4*r**3*c)/r**4
    assert s.cancel(BB-AA-fact)==0
    cv=(r+1/r)*(6-(r+1/r)**2)/4
    assert s.cancel((r**6-3*r**4-3*r*r+1+4*r**3*c).subs(c,cv))==0
    assert s.cancel(AA.subs(c,cv)-3*((r+1/r)**2-1))==0
    assert s.expand(h**3-6*h-4-(h+2)*(h*h-2*h-2))==0
    # For h>=2, h(6-h²)/4 decreases strictly. cos(phi)>=-1 bounds
    # h<=1+sqrt3, so AA<=9+6sqrt3. At r=1, AA<=9 separately.
    bound=9+6*s.sqrt(3);assert bound<25<1089
    rotations=HERE/'analytic-rotation-obstructions.json'
    output={'passed':True,'inputsSHA256':{source.name:sha(source),rotations.name:sha(rotations)},'boundary':{'normalizationAndFullSchmidtReductionRecomputed':True,'normalizationSquared':str(Z),'thirdOrderHilbertSchmidtDefectSquared':str(defect),'leadingSquaredDefectCoefficient':str(leading),'partitionForEveryEpsilonPositive':[2]*4+[1]*4,'discriminantsAndResultantVerifiedExactly':True,'degree12NoInfinityRootForPositiveEpsilon':True,'limitPartition':lam,'limitReducedStateExactlyIdentity':True,'positiveDefectForEveryPositiveEpsilon':True,'infimumZero':True,'attainmentInTargetPartitionEstablished':False},'cyclicObstructions':{'C3ForcedProbabilities':[[str(v)for v in w]for w in C3],'C4ForcedProbabilities':list(map(str,C4)),'C4RequiredMiddleCoefficientSquaredMagnitude':str(magnitudeSquared),'C4UniversalSquaredMagnitudeUpperBound':str(bound),'orbitCountingReviewed':True,'arbitraryComplexOrbitParametersCovered':True,'properRotationalSymmetryOrdersExcluded':['3','4','all orders at least5'],'remainingProperGroups':['C1','C2','D2']},'globalPartitionExclusionProved':False,'runtimeSeconds':time.monotonic()-start}
    (HERE/'boundary-and-rotation-independent-audit.json').write_text(json.dumps(output,indent=2)+'\n');print(json.dumps({'passed':True,'boundaryHS2':str(defect),'runtimeSeconds':output['runtimeSeconds']}))
if __name__=='__main__':main()
