#!/usr/bin/env python3
"""Verify the analytical family and certify a rational point near the Julia plane."""
from pathlib import Path
from fractions import Fraction as F
import math,json,hashlib,importlib.util
import sympy as s
import mpmath as m
HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[1]

def sqrt_interval(x,digits=55):
    x=F(x);scale=10**digits;z=math.isqrt((x.numerator*scale*scale)//x.denominator)
    return F(z,scale),F(z+1,scale)
def plus(a,b):return a[0]+b[0],a[1]+b[1]
def mul(a,b):
    xs=[x*y for x in a for y in b];return min(xs),max(xs)
def square(a):
    assert a[0]>=0
    return a[0]**2,a[1]**2

def verify():
    p=s.Symbol('p',real=True);q=s.Rational(241,768)-p
    P=[p-s.Rational(65,768),s.Rational(799,768)-4*p,6*p-s.Rational(731,768),s.Rational(255,256)-4*p,p]
    Q=[q+s.Rational(65,768),s.Rational(505,768)-4*q,6*q-s.Rational(341,768),s.Rational(539,768)-4*q,q]
    X=-(1024*p-255)*(3072*p-425)/s.Integer(2895532523520);Y=-p*(768*p-241)/s.Integer(13056);Z=(25690112*p*p-9670144*p+770185)/s.Integer(68941250560)
    A=s.Matrix([-1820,3,-s.Rational(14,33),1,0]);B=s.Matrix([-3,s.Rational(2,455),-s.Rational(1,2145),0,1]);K=s.Matrix([[1,1001,3003,91,0],[0,1092,10296,1092,0],[0,91,3003,1001,1]])
    assert K*A==s.zeros(3,1) and K*B==s.zeros(3,1)
    for i in range(5):assert s.factor(A[i]**2*X+B[i]**2*Y+2*A[i]*B[i]*Z-P[i]*Q[i]/s.binomial(17,4*i)/s.binomial(17,4*i+1))==0
    for off,probs in enumerate([P,Q]):
        for ell in range(4):assert s.factor(sum(probs[i]*(4*i+off)**ell for i in range(5))-sum(s.Rational(j**ell,18) for j in range(18)))==0
    # Exact rational enclosure proves that the source-near parameter choice lies in the family domain.
    pp=F(18370114016312017,10**17);rr=F(47980997065595143,10**17)
    rational=lambda z:F(str(s.factor(z.subs(p,s.Rational(pp.numerator,pp.denominator)))))
    ps=list(map(rational,P));qs=list(map(rational,Q));xx,yy,zz=map(rational,[X,Y,Z]);det=xx*yy-zz*zz
    assert all(v>0 for v in ps+qs) and xx>0 and det>0
    hs=[14,2002,2002,14];L2=[F(hs[i]**2)*ps[i+1]*qs[i]/math.comb(17,4*i+4)/math.comb(17,4*i+1) for i in range(4)]
    cosine=(1-rr*rr)/(1+rr*rr)
    root=sqrt_interval(L2[0]*L2[1]);R=plus((L2[0]+L2[1],)*2,mul(root,(2*cosine,)*2));d=mul(plus(R,(L2[2]-L2[3],)*2),(F(1,2),)*2)
    ds=square(d);delta=plus(mul((L2[2],)*2,R),(-ds[1],-ds[0]));assert R[0]>0 and delta[0]>0
    # Associate the downloadable expanded expressions with the already-proved
    # recurrence. Importing the constructor here is only a serialization check;
    # all family and domain identities above were independently reconstructed.
    spec=importlib.util.spec_from_file_location('n17_radical_export',HERE/'verify.py')
    constructor=importlib.util.module_from_spec(spec);spec.loader.exec_module(constructor)
    exact=constructor.construction(s.Rational(pp.numerator,pp.denominator),s.Rational(rr.numerator,rr.denominator),1)
    export_raw=(HERE/'source-reconstruction.json').read_bytes();export=json.loads(export_raw)
    assert export['exactBasis']==[[str(v) for v in row] for row in exact['rows']]
    assert export['parameters']=={'p':str(pp),'tanHalfBend':str(rr),'triangleBranch':1}
    # 80-digit rendering of this exact radical member, used only for numerical source identification.
    m.mp.dps=80;mp=lambda f:m.mpf(f.numerator)/f.denominator
    ps,qs=list(map(mp,ps)),list(map(mp,qs));xx,yy,zz=map(mp,[xx,yy,zz]);LL=list(map(mp,L2));r=mp(rr)
    z=m.sqrt(xx);w=(zz+1j*m.sqrt(xx*yy-zz*zz))/z
    xs=[complex(0)]*5
    xs=[m.mpf(str(A[i]))*z+m.mpf(str(B[i]))*w for i in range(5)]
    eta=(1-r*r+2j*r)/(1+r*r);S=m.sqrt(LL[0])+m.sqrt(LL[1])*eta;Rval=abs(S)**2;dv=(Rval+LL[2]-LL[3])/2
    edges=[m.sqrt(LL[0]),m.sqrt(LL[1])*eta,-S*(dv+1j*m.sqrt(LL[2]*Rval-dv*dv))/Rval];edges.append(-S-edges[2]);aa=[m.sqrt(ps[0])]
    for i in range(4):aa.append(aa[-1]*edges[i]/hs[i]/xs[i])
    bb=[m.conj(x/a) for x,a in zip(xs,aa)]
    rows=[]
    for off,ar in enumerate([aa,bb]):
        row=[m.mpc(0)]*18
        for i in range(5):row[4*i+off]=m.sqrt(math.comb(17,4*i+off))*ar[i]
        rows.append(row)
    src=next(v for v in json.loads((ROOT/'research/acliege/numerical-subspaces.json').read_text()) if v['n']==17 and v['t']==3)['sourceBasis']
    beta=m.mpf('0.11513821260197721');phases=[];errors=[]
    for row,source in zip(rows,src):
        rot=[z*m.exp(-1j*beta*(m.mpf(17)/2-r)) for r,z in enumerate(row)];ss=[m.mpc(str(re),str(im)) for re,im in source];phase=m.arg(sum(m.conj(a)*b for a,b in zip(rot,ss)));phases.append(str(phase));errors.append(float(max(abs(a*m.exp(1j*phase)-b) for a,b in zip(rot,ss))))
    assert max(errors)<5e-10
    output={'passed':True,'familyMomentIdentitiesExact':8,'familyCrossProductModulusIdentitiesExact':5,'crossProductKernelDimension':2,'familyParameterCount':2,'sourceParameters':{'p':str(pp),'tanHalfBend':str(rr),'triangleBranch':1},'sourceParameterDomainCertified':True,'sourceRationalDomainCertificate':{'R':list(map(str,R)),'triangleDiscriminant':list(map(str,delta)),'phaseDeterminant':str(det)},'sourceComparison':{'commonZRotation':str(beta),'basisGlobalPhases':phases,'maximumDickeCoefficientErrors':errors,'precisionDigits':80,'exactEqualityToDecimalSourceAsserted':False},'scope':'The simple catalogue representative uses p=5/24 and tanHalfBend=1. The exact rational parameter point recorded here is another member of the same explicitly constructed family and reconstructs the numerical Julia plane within its numerical accuracy.'}
    output['sourceExpandedRadicalsAssociatedExactly']=True
    output['sourceReconstructionSHA256']=hashlib.sha256(export_raw).hexdigest()
    (HERE/'family-verification.json').write_text(json.dumps(output,indent=2)+'\n')
    return output
if __name__=='__main__':print(json.dumps(verify(),indent=2))
