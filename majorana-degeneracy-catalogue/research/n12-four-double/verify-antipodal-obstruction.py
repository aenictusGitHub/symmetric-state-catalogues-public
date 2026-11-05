#!/usr/bin/env python3
"""Exact obstruction for four doubled points in two antipodal pairs.

SymPy only. Reconstructs all quadrupole equations and the three-body reduction
of the time-reversal-odd directions; imports no discovery code.
This is a restricted-family theorem, not a global exclusion for (2^4,1^4).
"""
from pathlib import Path
import json,time
import sympy as s
HERE=Path(__file__).resolve().parent
h,z,b,c,d,e,f,t=s.symbols('h z b c d e f t',real=True)
vars=[h,b,c,d,e,f]
def primitive(E):return s.factor(s.Poly(s.expand(E),*vars).clear_denoms()[1].primitive()[1].as_expr())
def interior_root_count(E):
    P=s.Poly(E,h,domain=s.QQ)
    for endpoint in [-1,1]:
        while P.eval(endpoint)==0:P=P.exquo(s.Poly(h-endpoint,h))
    return int(P.count_roots(-1,1))

def main():
    start=time.monotonic();Q=z**4+2*h*z**2+1
    R=(b+s.I*d)*z**4+(e+s.I*f)*z**3+c*z**2+(-e+s.I*f)*z+b-s.I*d
    P=s.Poly(s.expand(Q**2*R),z);p=[P.nth(i)for i in range(13)]
    weights=[s.expand(v*s.conjugate(v))/s.binomial(12,r)for r,v in enumerate(p)]
    zz=sum(((6-r)**2-14)*weights[r]for r in range(13))
    zp=sum(s.conjugate(p[r-1])*p[r]*s.Rational(r,s.binomial(12,r-1))*(13-2*r)for r in range(1,13))
    pp=sum(s.conjugate(p[r-2])*p[r]*s.Rational(r*(r-1),s.binomial(12,r-2))for r in range(2,13))
    E=list(map(primitive,[zz,s.re(zp),s.im(zp),s.re(pp),s.im(pp)]))
    stored=json.loads((HERE/'antipodal-reduction.json').read_text())
    assert all(s.expand(a-s.sympify(q,locals={str(v):v for v in vars}))==0 for a,q in zip(E,stored['conditions']))
    # Both e and f nonzero: t=f/e. Solve only a nonsingular real linear system.
    bc=s.solve([E[1].subs({e:1,f:t}),E[2].subs({e:1,f:t})],[b,c])
    D=(h*h-2*h+5)*(h*h+2*h+5);K=7*h**4+18*h*h+7
    A=h**6+2*h**5+9*h**4+12*h**3+51*h*h+34*h+115
    lhs=s.factor(E[4].subs({b:bc[b],c:bc[c],f:t*e}))
    rhs=-(t*t*D*K*e*e+16*(h*h+1)*(t*t*A+A.subs(h,-h))*d*d)/(t*D)
    assert s.cancel(lhs-rhs)==0
    # A(±h)>=115-2-12-34=67; D,K>0 on [-1,1].
    assert 115-2-12-34==67
    # Exactly one of e,f nonzero: remaining equation makes d=0.
    assert s.expand(E[1].subs(e,0)+4*d*f*(h**4+4*h**3+12*h*h+16*h+23))==0
    assert s.expand(E[2].subs(f,0)+4*d*e*(h**4-4*h**3+12*h*h-16*h+23))==0
    assert 23-4-16>0
    single_branches=[]
    for zero,active,off in [(e,f,2),(f,e,1)]:
        csol=s.solve(E[off].subs({zero:0,d:0,active:1}),c)[0]
        polys=[s.together(E[j].subs({zero:0,d:0,c:csol})).as_numer_denom()[0]for j in[0,3]]
        M=s.Matrix([[s.expand(q).coeff(b,2),s.expand(q).coeff(active,2)]for q in polys])
        determinant=s.factor(M.det());assert determinant!=0 and interior_root_count(determinant)==0
        single_branches.append({'zero':str(zero),'active':str(active),'reconstructedC':str(s.factor(csol)),'determinant':str(determinant),'interiorRealRootCount':0})
    # e=f=0, d!=0: the last equation fixes c; the diagonal equation is positive.
    csol=-2*b*h/(h*h+2)
    diagonal=s.factor(s.together(E[0].subs({e:0,f:0,c:csol})).as_numer_denom()[0])
    expected=(-16*h**8+1104*h**4+4320*h*h+4320)*b*b+(-16*h**8-48*h**6+1088*h**4+4416*h*h+4352)*d*d
    assert s.expand(diagonal-expected)==0
    assert 4320-16>0 and 4352-16-48>0
    # e=f=d=0, b!=0: normalize b=1 and eliminate c.
    even=[E[j].subs({e:0,f:0,d:0,b:1})for j in[0,3]]
    resultant=s.factor(s.resultant(*even,c));allowed=resultant/h**2
    assert interior_root_count(allowed)==0
    assert even[0].subs(h,0)==1080
    # b=d=e=f=0 permits only h=0, c!=0, the double octahedron.
    assert s.factor(E[0].subs({b:0,d:0,e:0,f:0}))==-2*c*c*h*h*(3*h*h+11)
    boundary=s.Poly(z*z*(z**4+1)**2,z);sqf=boundary.sqf_list()[1]
    lam=sorted([m for q,m in sqf for _ in range(q.degree())]+[12-boundary.degree()],reverse=True)
    lam=list(map(int,lam));assert lam==[2]*6
    cv=[(-1)**r*boundary.nth(r)/s.sqrt(s.binomial(12,r))for r in range(13)]
    norm=sum(v*v for v in cv);cv=[v/s.sqrt(norm)for v in cv]
    schmidt=s.Matrix(4,10,lambda a,b:cv[a+b]*s.sqrt(s.binomial(3,a)*s.binomial(9,b)/s.binomial(12,a+b)))
    assert s.simplify(schmidt*schmidt.H)==s.eye(4)/4
    # Arbitrary complex R: the TR-odd directions must vanish under rho3.
    RB=[z**4+1,z**2,s.I*(z**4-1),z**3-z,s.I*(z**3+z)]
    pp=[[s.Poly(s.expand(Q**2*q),z).nth(r)for r in range(13)]for q in RB]
    pairs=[(i,j)for i in range(5)for j in range(i+1,5)]
    rowdes=[('diag',a,a)for a in range(4)]+[(kind,a,c)for a in range(4)for c in range(a+1,4)for kind in['re','im']]
    M=[]
    for kind,a,cidx in rowdes:
        line=[]
        for i,j in pairs:
            v=s.I*(-1)**(a+cidx)*sum(s.binomial(9,q)*(pp[i][a+q]*s.conjugate(pp[j][cidx+q])-pp[j][a+q]*s.conjugate(pp[i][cidx+q]))/(s.binomial(12,a+q)*s.binomial(12,cidx+q))for q in range(10))
            line.append(s.expand(s.im(v)if kind=='im'else s.re(v)))
        scale=s.ilcm(*[s.denom(q)for v in line for q in s.Poly(v,h).all_coeffs()]);M.append([s.expand(v*scale)for v in line])
    M=s.Matrix(M);odd=json.loads((HERE/'odd-restriction.json').read_text())
    assert M==s.Matrix([[s.sympify(v,locals={'h':h})for v in row]for row in odd['matrix']])
    determinant=s.factor(M[odd['selectedRows'],:].det(method='domain-ge'))
    assert s.expand(determinant-s.sympify(odd['determinant'],locals={'h':h}))==0
    factors=[]
    for q,m in s.factor_list(determinant,h)[1]:
        count=int(s.Poly(q,h).count_roots(-1,1));assert count==0
        factors.append({'polynomial':str(q),'multiplicity':int(m),'closedIntervalRealRootCount':count})
    output={'passed':True,'scope':'No 3-AC pure symmetric 12-qubit state has partition (2,2,2,2,1,1,1,1) with its four doubled points forming two antipodal pairs. Simple points may be arbitrary complex points. This does not exclude non-antipodal doubled-point configurations.','normalForm':'P(z)=(z^4+2h*z^2+1)^2 R(z), -1<h<1; arbitrary complex quartic R, including lower degree and pole roots.','antipodalResidualForcedByOddReductionRank':10,'oddRestrictionDeterminant':str(determinant),'oddRestrictionFactors':factors,'bothEFNonzeroContradiction':'e² times a strictly positive coefficient plus d² times a strictly positive coefficient is zero; e is nonzero.','positiveAUniformLowerBound':67,'singleNonzeroBranches':single_branches,'dNonzeroDiagonalLowerBounds':{'bSquared':4304,'dSquared':4288},'realEvenResultant':str(resultant),'soleSurvivingProjectivePolynomial':'z^2*(z^4+1)^2','survivingPartition':lam,'survivingOrder3VerifiedExactly':True,'targetMultiplicityGuard':'The survivor has six double roots, including zero and infinity, so it is outside (2^4,1^4). h=±1 would merge roots of Q and is also excluded.','globalComplexExclusionProved':False,'runtimeSeconds':time.monotonic()-start}
    (HERE/'antipodal-obstruction-verification.json').write_text(json.dumps(output,indent=2)+'\n');print(json.dumps({'passed':True,'rank':10,'onlyPartition':lam,'runtimeSeconds':output['runtimeSeconds']}))
if __name__=='__main__':main()
