#!/usr/bin/env python3
"""Independent exact audit of the balanced cyclic AC2 exclusion.
This certificate first forces one dimension in each cyclic sector, then excludes it;
it does not bound arbitrary AC subspaces more strongly than the global bounds.
"""
from fractions import Fraction as F
from pathlib import Path
import json
import sympy as S

base=Path(__file__).resolve().parent
N,m=S.symbols('N m',positive=True,real=True)
j=N/2;c=N*(N+2)/12
th=S.pi/8;cs=S.cos(th);sn=S.sin(th)
f=[S.sin((a+1)*th) for a in range(7)]
window=[]
for a in range(7):
    left=f[a-1] if a else S.Integer(0)
    right=f[a+1] if a<6 else S.Integer(0)
    assert S.simplify(left+right-2*cs*f[a])==0
    assert S.simplify(right-left-2*sn*S.cos((a+1)*th))==0
    raw=(1+m*m/c)*f[a]+((j*(j+1)-m*m+m)*right+(j*(j+1)-m*m-m)*left)/(2*c)
    expected=(1+3*cs+(1-cs)*m*m/c)*f[a]+m/c*sn*S.cos((a+1)*th)
    difference_num=S.together(raw-expected).as_numer_denom()[0]
    assert all(S.simplify(v)==0 for v in S.Poly(difference_num,m,N).coeffs())
    # This squared inequality implies sin(theta)*|cot((a+1)theta)| <= 1.
    trig_margin=S.simplify(f[a]**2-sn**2*S.cos((a+1)*th)**2)
    assert trig_margin.is_nonnegative
    window.append({'a':a,'trigonometricMargin':str(trig_margin),'coefficientIdentity':True})
assert F(49,36)**2<2
# cos(pi/8)>11/12 follows by squaring and sqrt(2)>49/36.
Q=lambda z:z*z-F(8,5)*z
lat=[1-F(2*l,5) for l in range(5)]
qex=[Q(z)+F(3,5) for z in lat]
assert qex==[0,0,F(8,25),F(24,25),F(48,25)]
records=[]
for s_max,threshold in [(5,420),(1,309)]:
    leading=[F(12)+F(2*l*s_max,5) for l in range(5)]
    slopes=[F(8,5)-2*z for z in lat]
    first=sum(slopes[l]*leading[l] for l in range(1,5))
    second=sum(leading[l]**2 for l in range(1,5))
    B=F(15,4)-F(6,threshold+2)
    U=-F(3,5)+(F(16,5)+first/threshold+second/threshold**2)/B
    gap=F(1,3)-U
    assert B>0 and gap>0
    assert leading[0]/threshold<=F(2,5) # Band zero contributes no positive error.
    target=threshold//5+1
    assert target>=s_max+6 # All seven indices in every band exist.
    records.append({'deficitMaximum':s_max,'thresholdN':threshold,
                    'linearErrorCoefficient':str(first),'quadraticErrorCoefficient':str(second),
                    'B_at_threshold':str(B),'U_at_threshold':str(U),
                    'oneThirdMinusU':str(gap),'positiveGap':True})
assert records[0]['linearErrorCoefficient']=='584/5'
assert records[0]['quadraticErrorCoefficient']=='1176'
assert records[0]['oneThirdMinusU']=='872/1655325'
assert records[1]['linearErrorCoefficient']=='424/5'
assert records[1]['quadraticErrorCoefficient']=='3384/5'
assert records[1]['oneThirdMinusU']=='93358/738545535'
# A C_k character sector has at most five equally spaced magnetic nodes.
# Removing its top/bottom coefficient in a putative dimension-two subspace
# gives the following two variance bounds by E[(m-A)(B-m)]>=0.
k,d,zeta=S.symbols('k d zeta',positive=True,real=True)
assert S.expand(-(d-2*k)*(d+k)-(2*k*k+k*d-d*d))==0
assert S.expand(-(d-k)*(d+2*k)-(2*k*k-k*d-d*d))==0
sector_checks=[]
for deficit in range(1,6):
    gap=S.expand(((5*k-deficit)*(5*k-deficit+2)-24*k*k)/12)
    expected=(k*k-10*(deficit-1)*k+deficit*(deficit-2))/12
    assert S.expand(gap-expected)==0
    shifted=S.Poly(gap.subs(k,zeta+85).expand(),zeta)
    assert all(c>0 for c in shifted.coeffs())
    sector_checks.append({'deficit':deficit,'gapAfterK85Substitution':str(shifted.as_expr()),'positiveForAllK85':True})
assert S.simplify(((5*k-1)*(5*k+1)-24*k*k)/12-(k*k-1)/12)==0
report={'passed':True,'arithmetic':'Exact rational and symbolic trigonometric identities',
        'scope':'All second-order AC subspaces of dimension k=floor(N/5)+1 invariant under rotation through 2pi/k about any axis',
        'allNConclusion':'No such code for N>=420; on N=5k-1 already no such code for k>=62',
        'unrestrictedExclusion':False,
        'sectorRankLemma':{'description':'A dimension-two AC2 subspace on five step-k magnetic nodes would have c<=2k^2. In the stated ranges c>2k^2, so each cyclic sector has rank at most one. Dimension k then forces one dimension in every sector.','checks':sector_checks,'subsequenceGap':'(k^2-1)/12>0 for k>1'},
        'windowIdentities':window,'nodalExcesses':list(map(str,qex)),
        'thresholds':records,
        'monotonicity':'U(N)=-3/5+A(N)/B(N), where A is positive decreasing and B is positive increasing; hence U decreases above either stated threshold.',
        'projectorPremise':'AC2 makes C, Jx C, Jy C, Jz C mutually orthogonal; E=P+sum Ja P Ja/c is a projection, so every diagonal entry is <=1. Reindexing its weighted seven-residue sum yields the independently checked coefficient identity.'}
(base/'cyclic-exclusion-independent-audit.json').write_text(json.dumps(report,indent=2)+'\n')
print(json.dumps({'passed':True,'unrestrictedExclusion':False,'thresholdGaps':[r['oneThirdMinusU'] for r in records]}))
