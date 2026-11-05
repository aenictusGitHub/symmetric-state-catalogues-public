#!/usr/bin/env python3
"""Exact radical construction of an optimal N=17, t=3 AC plane.
No numerical solve is used in the existence proof.
"""
from pathlib import Path
import json, math
import numpy as np
import sympy as s
HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[1]

def construction(p=s.Rational(5,24),bend=s.S.One,branch=1):
    p=s.sympify(p);bend=s.sympify(bend);q=s.Rational(241,768)-p
    P=[p-s.Rational(65,768),s.Rational(799,768)-4*p,6*p-s.Rational(731,768),s.Rational(255,256)-4*p,p]
    Q=[q+s.Rational(65,768),s.Rational(505,768)-4*q,6*q-s.Rational(341,768),s.Rational(539,768)-4*q,q]
    X=s.factor(-(1024*p-255)*(3072*p-425)/s.Integer(2895532523520))
    Y=s.factor(-p*(768*p-241)/s.Integer(13056))
    Z=s.factor((25690112*p*p-9670144*p+770185)/s.Integer(68941250560))
    det=s.factor(X*Y-Z*Z)
    A=[-1820,3,-s.Rational(14,33),1,0];B=[-3,s.Rational(2,455),-s.Rational(1,2145),0,1]
    z=s.sqrt(X);w=(Z+s.I*s.sqrt(det))/z
    x=[s.factor(aa*z+bb*w) for aa,bb in zip(A,B)]
    h=[14,2002,2002,14]
    L2=[s.factor(h[i]**2*P[i+1]*Q[i]/s.binomial(17,4*i+4)/s.binomial(17,4*i+1)) for i in range(4)]
    eta=(1-bend*bend+2*s.I*bend)/(1+bend*bend)
    S=s.sqrt(L2[0])+s.sqrt(L2[1])*eta
    R=s.factor(L2[0]+L2[1]+2*s.sqrt(L2[0]*L2[1])*(1-bend*bend)/(1+bend*bend))
    d=(R+L2[2]-L2[3])/2;delta=s.factor(L2[2]*R-d*d)
    e=[s.sqrt(L2[0]),s.sqrt(L2[1])*eta,-S*(d+branch*s.I*s.sqrt(delta))/R]
    e.append(-S-e[2])
    a=[s.sqrt(P[0])]
    for i in range(4):a.append(s.factor(a[-1]*e[i]/h[i]/x[i]))
    b=[s.conjugate(xi/ai) for xi,ai in zip(x,a)]
    rows=[]
    for offset,ar in [(0,a),(1,b)]:
        row=[s.S.Zero]*18
        for i in range(5):row[4*i+offset]=s.factor(s.sqrt(s.binomial(17,4*i+offset))*ar[i])
        rows.append(row)
    return locals()

def zero(v):return s.simplify(s.expand_complex(v))==0

def split(row,t):
    return s.Matrix(t+1,18-t,lambda a,b:row[a+b]*s.sqrt(s.binomial(t,a)*s.binomial(17-t,b)/s.binomial(17,a+b)))

def verify():
    data=construction();P,Q,x,e,L2=data['P'],data['Q'],data['x'],data['e'],data['L2']
    assert all(v>0 for v in P+Q) and data['X']>0 and data['det']>0 and data['delta']>0
    for i in range(5):assert zero(x[i]*s.conjugate(x[i])-P[i]*Q[i]/s.binomial(17,4*i)/s.binomial(17,4*i+1))
    assert zero(sum(e)) and all(zero(ei*s.conjugate(ei)-ell) for ei,ell in zip(e,L2))
    K=s.Matrix([[1,1001,3003,91,0],[0,1092,10296,1092,0],[0,91,3003,1001,1]])
    assert all(zero(v) for v in K*s.Matrix(x))
    rows=data['rows'];U=s.Matrix.hstack(*map(s.Matrix,rows))
    assert all(zero(v) for v in U.H*U-s.eye(2))
    M=[split(row,3) for row in rows];count=0
    for i in range(2):
        for j in range(2):
            residual=M[i]*M[j].H-(s.eye(4)/4 if i==j else s.zeros(4))
            for val in residual:assert zero(val);count+=1
    purities=[]
    for row in rows:
        M4=split(row,4);R4=M4*M4.H
        purity=s.simplify(s.expand_complex(s.trace(R4*R4)))
        assert purity>s.Rational(1,5);purities.append(str(purity))
    result={'passed':True,'n':17,'t':3,'dimension':2,'orthonormalityExact':True,'complexCrossEntriesExactlyVerified':count,'wholeSpanAC':True,'P':list(map(str,P)),'Q':list(map(str,Q)),'parameters':{'p':'5/24','tanHalfBend':'1','triangleBranch':1},'phaseDeterminant':str(data['det']),'triangleDiscriminant':str(data['delta']),'rho4Purities':purities,'displayedBasisVectorsExactACOrder':3,'upperBound':'6k <= 17','maximumDimension':2,'method':'Direct exact normalized-Dicke splitting, including all four ordered cross reductions; phase and quadrilateral identities are checked separately.'}
    (HERE/'verification.json').write_text(json.dumps(result,indent=2)+'\n')
    return result

def save_basis():
    d=construction();basis=[]
    for j,row in enumerate(d['rows'],1):
        basis.append({'id':f'radical-n17-k2-t3-basis-{j}','n':17,'label':f'Analytical third-order AC plane · basis vector {j}','exactCoefficients':list(map(str,row)),'coefficients':[[float(s.re(c).evalf()),float(s.im(c).evalf())] for c in row],'support':[i for i,c in enumerate(row) if c!=0],'stateLatex':(r'|u_1\rangle=\sum_{i=0}^{4}\sqrt{\binom{17}{4i}}\,a_i|D_{17}^{(4i)}\rangle' if j==1 else r'|u_2\rangle=\sum_{i=0}^{4}\sqrt{\binom{17}{4i+1}}\,b_i|D_{17}^{(4i+1)}\rangle'),'coefficientDefinitionLatex':r'a_0=\sqrt{P_0},\qquad a_{i+1}=\frac{a_i e_i}{h_i x_i},\qquad b_i=\overline{\frac{x_i}{a_i}}.\quad\text{All constants are defined in the linked analytical construction.}'})
    out={'id':'radical-n17-k2-t3','n':17,'t':3,'dimension':2,'label':'Analytical two-parameter family · simple radical representative','evidence':'exact','basis':basis,'parameters':{'p':'5/24','tanHalfBend':'1','triangleBranch':1},'derivationLatex':r'\begin{aligned}\boldsymbol P&=\frac1{768}(95,159,229,125,160),\\\boldsymbol Q&=\frac1{768}(146,181,145,215,81),\\R_3(|u_a\rangle\langle u_b|)&=\delta_{ab}I_4/4,\qquad K(17,3)=2.\end{aligned}','provenance':'Exact radical member of an analytical family on the support of the reported Julia plane Π_17_2_3. The displayed simple parameter choice is not claimed to equal that numerical plane. A nearby exact family member reconstructs the source within its stated numerical accuracy; see the linked derivation and source comparison.','source':{'file':'ACLiege.jl','function':'Π_17_2_3','line':939},'construction':{key:[str(x) for x in d[key]] if isinstance(d[key],list) else str(d[key]) for key in ['P','Q','X','Y','Z','det','A','B','x','h','L2','R','d','delta','e','a','b']},'family':{'parameters':['p','tanHalfBend'],'q':'241/768-p','admissibility':'P_i,Q_i>0, X>0, XY-Z^2>0, R>0 and L_2^2 R-d^2>=0. The formula uses either sign of the triangle square root.','sourceExactParameters':{'p':'18370114016312017/100000000000000000','tanHalfBend':'47980997065595143/100000000000000000','triangleBranch':1},'sourceParameterCertificate':'family-verification.json'}}
    (HERE/'basis.json').write_text(json.dumps(out,indent=2)+'\n')
    return out

if __name__=='__main__':
    save_basis();print(json.dumps(verify(),indent=2))
