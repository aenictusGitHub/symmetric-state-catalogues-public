#!/usr/bin/env python3
"""Exact radical construction of a two-dimensional third-order AC subspace.
All ordered three-qubit cross reductions are reconstructed independently.
"""
from pathlib import Path
import hashlib,json,re
import numpy as np
import sympy as s
ROOT=Path(__file__).resolve().parent

def clean(x):return s.simplify(s.expand_complex(x))
def zero(x,label):
    if clean(x)!=0:raise AssertionError(label+': '+str(x))
def basis():
    z=(545-192*s.I*s.sqrt(10))/(7*s.sqrt(13585))
    f=-(55+4*s.I*s.sqrt(10))/(7*s.sqrt(65))
    g=(19+8*s.I*s.sqrt(10))/s.sqrt(1001)
    u=s.zeros(17,1);v=s.zeros(17,1)
    for r,c in {0:s.sqrt(38)/16,4:z*s.sqrt(10)/8,8:z*s.Rational(5,8),12:z*s.sqrt(10)/8,16:z*z*s.sqrt(38)/16}.items():u[r]=c
    for r,c in {2:s.sqrt(5)/4,6:f*s.sqrt(3)/4,10:g*s.sqrt(3)/4,14:f*g*s.sqrt(5)/4}.items():v[r]=c
    return [u,v],{'z':z,'f':f,'g':g}
def split(v,t):
    return s.Matrix(t+1,17-t,lambda a,b:v[a+b]*s.sqrt(s.binomial(t,a)*s.binomial(16-t,b)/s.binomial(16,a+b)))
def source_comparison(vectors):
    path=ROOT.parents[1]/'references/ACLiege.jl';text=path.read_text()
    block=text.split('function Π_16_2_3()',1)[1].split('\nend',1)[0]
    rows=[]
    for line in block.splitlines():
        if '[' not in line or 'im' not in line:continue
        pairs=re.findall(r'([-+]?\d+\.\d+)\s*([+-])\s*(\d+\.\d+)im',line)
        rows.append([complex(float(a),float(b)*(-1 if sign=='-' else 1))for a,sign,b in pairs])
    assert list(map(len,rows))==[5,4]
    supports=[[0,4,8,12,16],[2,6,10,14]]
    exact=[[complex(s.N(v[r],30))for r in support]for v,support in zip(vectors,supports)]
    beta=float(np.angle(rows[0][1]/rows[0][0]/(exact[0][1]/exact[0][0]))/4)
    errors=[];phases=[]
    for a,b,support in zip(rows,exact,supports):
        a=np.array(a);b=np.array(b)*np.exp(1j*np.array(support)*beta)
        phase=np.vdot(b,a);phase/=abs(phase)
        errors.append(float(np.max(abs(a-phase*b))));phases.append(float(np.angle(phase)+8*beta))
    assert max(errors)<1e-8
    return {'source':'references/ACLiege.jl, Π_16_2_3 at line 901','sourceSha256':hashlib.sha256(path.read_bytes()).hexdigest(),'convention':'source_a approximately exp(i*chi_a) exp(-i*beta*J_z) u_a, J_z|D_16^r>=(8-r)|D_16^r>','betaRadians':beta,'chiRadians':phases,'maximumCoefficientErrors':errors,'scope':'Approximate comparison with the supplied numerical vectors. The exact proof uses only the radical identities.'}
def verify():
    vectors,parameters=basis();saved=json.loads((ROOT/'basis.json').read_text())
    for name,value in parameters.items():zero(value*s.conjugate(value)-1,'Unit phase '+name)
    for v,record in zip(vectors,saved['basis']):
        for x,y in zip(v,record['exactCoefficients']):zero(x-s.sympify(y),'Saved exact coefficients')
    gram=s.Matrix.hstack(*vectors).H*s.Matrix.hstack(*vectors)
    for x in gram-s.eye(2):zero(x,'Orthonormality')
    mats=[split(v,3)for v in vectors];count=0
    for a,A in enumerate(mats):
        for b,B in enumerate(mats):
            target=s.eye(4)/4 if a==b else s.zeros(4)
            for x in A*B.H-target:zero(x,'Full ordered cross reduction');count+=1
    # Check the four explicit cancellation identities printed in proof.tex,
    # and their association with every entry of the displayed cross matrix.
    z=parameters['z'];F=s.conjugate(parameters['f']);G=s.conjugate(parameters['g'])
    E1=7*s.sqrt(19)+z*(s.sqrt(91)*F*G+15*s.sqrt(11)*F+9*s.sqrt(35)*G)
    E2=s.sqrt(91)*F*G+s.sqrt(11)*F+2*s.sqrt(35)*G
    E3=2*s.sqrt(35)*F+s.sqrt(11)*G+s.sqrt(91)
    E4=7*s.sqrt(19)*F*G*z+9*s.sqrt(35)*F+15*s.sqrt(11)*G+s.sqrt(91)
    for e in [E1,E2,E3,E4]:zero(e,'Displayed phase cancellation')
    displayed=s.zeros(4)
    displayed[0,2]=E1/896;displayed[1,3]=3*z*E2/448
    displayed[2,0]=3*z*E3/448;displayed[3,1]=z*E4/896
    for x in mats[0]*mats[1].H-displayed:zero(x,'Displayed cross matrix')
    purities=[]
    for v in vectors:
        rho=split(v,4)*split(v,4).H;purity=clean(s.trace(rho*rho));assert purity>s.Rational(1,5);purities.append(str(purity))
    assert purities==['3457/16562','3433/16562']
    report={'passed':True,'n':16,'t':3,'dimension':2,'arithmetic':'Exact rational and radical arithmetic; numerical source comparison is separate.','orthonormalityEntriesChecked':4,'fullOrderedCrossReductionEntriesChecked':count,'displayedPhaseCancellationIdentitiesChecked':4,'crossReductionIdentity':'R_3(|u_a><u_b|)=delta_ab*I_4/4','individualExactACOrders':[3,3],'fourthReductionPurities':purities,'upperBound':'6*k <= N = 16','maximumDimension':2,'sourceComparison':source_comparison(vectors),'basisSha256':hashlib.sha256((ROOT/'basis.json').read_bytes()).hexdigest()}
    return report
if __name__=='__main__':
    report=verify();(ROOT/'verification.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report,indent=2))
