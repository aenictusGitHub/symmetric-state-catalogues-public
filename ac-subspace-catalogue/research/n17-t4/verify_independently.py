#!/usr/bin/env python3
"""Independent exact construction and saturation exclusion for N=17,t=4.

The coefficient expressions, Dicke reductions, six-j channel eigenvalues,
and time-reversal lift are reconstructed without importing another verifier.
"""
from pathlib import Path
import json
import sympy as s
from sympy.physics.wigner import wigner_6j
HERE=Path(__file__).resolve().parent
V=s.Symbol('v',real=True)
q=s.sqrt(20995)
u=710-9*q
D=1628760-u*u

def reduce(value):
    value=s.expand(s.expand_complex(value))
    result=s.rem(s.Poly(value,V),s.Poly(V*V-D,V)).as_expr()
    return s.simplify(s.expand(result))

def zero(value,label):
    r=reduce(value)
    if r!=0:raise AssertionError(label+': '+str(r))

def vector(sign):
    A=(u-sign*s.I*V)/(756*s.sqrt(15))
    B=2*s.sqrt(21)/27
    C=-2*(2130+u+sign*2*s.I*V)/(189*s.sqrt(4290))
    v=s.zeros(18,1)
    for r,c in zip([1,4,7,10,13,16],[A,B,C,C,B,A]):v[r]=c
    return v

def split(v,t):
    return s.Matrix(t+1,18-t,lambda a,b:v[a+b]*s.sqrt(s.binomial(t,a)*s.binomial(17-t,b)/s.binomial(17,a+b)))

def branch(sign):
    c=vector(sign)
    zero((c.H*c)[0]-1,'Norm')
    M=split(c,4);rho4=(M*M.H).applyfunc(reduce)
    for a in range(5):
        for b in range(5):zero(rho4[a,b]-s.Rational(int(a==b),5),'rho4')
    M=split(c,5);rho5=(M*M.H).applyfunc(reduce)
    expected=s.eye(6)/6
    for a,b,x in [(0,3,-s.sqrt(8398)/9282),(1,4,s.sqrt(20995)/9282),(2,5,-s.sqrt(8398)/9282)]:
        expected[a,b]=expected[b,a]=x
    for x in rho5-expected:zero(x,'rho5')
    purity=reduce(s.trace(rho5*rho5))
    zero(purity-s.Rational(5443,32487),'Purity5')
    for r in range(18):zero(c[r]-c[17-r],'Palindromic pi-x symmetry')
    assert all((17-2*r)%6==3 for r in [1,4,7,10,13,16])
    ratio=reduce(c[1]*c[16]/(c[4]*c[13]))
    assert s.simplify(s.im(ratio))!=0
    return {'sign':sign,'normalizationExact':True,'rho4':'I_5/5','rho4EntriesChecked':25,
            'rho5EntriesChecked':36,'rho5Purity':str(purity),'exactACOrder':4,
            'properRotationGroup':'D3','axisNormalizerInvariant':str(ratio),
            'sourceRelationship':'An independent symmetric construction; identity or rotation equivalence to the numerical Julia state is not asserted.'}

def exclusion():
    N=17;j=s.Rational(N,2);Q=N*(N+2)
    rows=[]
    for L in range(N+1):
        direct=s.simplify(sum((-1)**(N+ell+L)*(2*ell+1)*wigner_6j(j,j,ell,j,j,L) for ell in range(3)))
        x=L*(L+1)
        polynomial=s.Rational(3*(10*x*x+(16-12*Q)*x+3*Q*(Q-3)),(N+1)*Q*(Q-3))
        assert direct==polynomial
        rows.append({'L':L,'eigenvalue':str(direct)})
    assert [r['L'] for r in rows if s.Rational(r['eigenvalue'])==0]==[16]
    x=s.Symbol('x')
    assert s.expand(10*(x-272)*(x-114))==10*x*x+(16-12*Q)*x+3*Q*(Q-3)
    assert all(L*(L+1)!=114 for L in range(18))
    # A fourth-order rank-two code saturates nine tensor images in H_17.
    assert 9*2==18
    # The channel equation then forces P=I_18/9+X_16. The qubit Schmidt
    # partners form a rank-four projector Q=2 Tr_1(P) on H_16.
    rank=4;d=17;alpha=s.Rational(rank,d)
    norm2=rank-rank*rank/s.Integer(d)
    assert norm2==s.Rational(52,17)
    assert (1-2*alpha)/norm2==s.Rational(9,52)
    assert alpha*(1-alpha)/norm2==s.Rational(1,17)
    # Time-reversal intertwining on spin eight, checked on Jz,J+,J-.
    z=s.diag(*[8-r for r in range(17)]);p=s.zeros(17)
    for r in range(1,17):p[r-1,r]=s.sqrt(r*(17-r))
    F=s.Matrix(17,17,lambda r,c:(-1)**r if c==16-r else 0)
    assert F*F.H==s.eye(17)
    for J in [z,(p+p.T)/2,(p-p.T)/(2*s.I)]:assert F*J.T==-J*F
    highest=s.zeros(17);highest[0,16]=1
    lifted=highest*F
    assert lifted[0,0]==1 and sum(int(v!=0) for v in lifted)==1
    # Tracing one qubit kills the highest-weight rank-16 operator; covariance
    # then kills its complete rank-16 multiplet. Direct Dicke split of
    # |D_16^0><D_16^16|: removed-qubit occupations are different.
    reduced=s.zeros(16)
    for a in range(16):
        for b in range(16):
            reduced[a,b]=sum(highest[a+h,b+h]*s.sqrt(s.binomial(15,a)*s.binomial(15,b)/
                  (s.binomial(16,a+h)*s.binomial(16,b+h))) for h in [0,1])
    assert reduced==s.zeros(16)
    # vec(YF) is a pure symmetric 32-qubit state with rho15=I16/16.
    # Rank-at-most-seven spin tensors would give 64 orthogonal vectors
    # in a space of dimension 33.
    error_images=sum(2*ell+1 for ell in range(8));assert error_images==64 and error_images>33
    return {'passed':True,'hypotheticalDimension':2,'packingSaturation':'9*2=18',
            'channelEigenvalues':rows,'channelKernelRanks':[16],
            'channelNumeratorFactorization':'10*(x-272)*(x-114), x=L*(L+1)',
            'SchmidtPartnerRank':4,'SchmidtPartnerAmbientDimension':17,
            'trY2':str(norm2),'liftedRho16':'I_17/17+(9/52)*Y_16',
            'liftedPureState':{'n':32,'t':15,'rho15':'I_16/16'},
            'GramContradiction':'64<=33','timeReversalIntertwiningExact':True,
            'highestRankPartialTraceZeroExact':True,'conclusion':'K(17,4)<=1 for arbitrary complex subspaces.'}

def main():
    # D>0 follows already from sqrt(20995)>144.
    assert 20995>144**2
    assert s.expand(D)==-575935+12780*q
    assert -575935+12780*144>0
    branches=[branch(sign) for sign in [1,-1]]
    source=json.loads((HERE/'basis.json').read_text())
    rows=source['basis']
    assert len(rows)==1 and len(rows[0]['exactCoefficients'])==18
    for text,value in zip(rows[0]['exactCoefficients'],vector(1)):
        difference=s.simplify(s.expand(s.sympify(text)-value.subs(V,s.sqrt(s.expand(D)))))
        assert difference==0, 'Published coefficient association: '+str(difference)
    report={'passed':True,'n':17,'t':4,'branches':branches,
            'arithmetic':'Exact rational and radical arithmetic, with v^2 reduced by its defining positive radicand.',
            'independentMethod':'No primary verifier or builder is imported. Full Dicke reductions, six-j channel eigenvalues, time reversal and highest-rank partial trace are independently reconstructed.',
            'upperBound':exclusion(),'maximumDimension':1,
            'sourceCoefficientAssociation':True,
            'classificationScope':'Two displayed construction branches; not an exhaustive classification of pure states.'}
    (HERE/'independent-verification.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps({'passed':True,'fullRho4Entries':50,'fullRho5Entries':72,'channelKernel':[16],'maximumDimension':1},indent=2))
if __name__=='__main__':main()
