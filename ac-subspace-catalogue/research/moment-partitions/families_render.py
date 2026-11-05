#!/usr/bin/env python3
"""Numerical renderings of the rigorously certified radical prescription.
These direct spin checks supplement the all-N exact polynomial certificate.
"""
from pathlib import Path
from fractions import Fraction as Q
import numpy as np
import json

def data(n):
    if n<17:raise ValueError('N>=17 required')
    q,rho=divmod(n-1,8);rho+=1
    h=2*q+1 if rho<8 else 2*q+2
    offset=(rho-1)//2 if rho<8 else 2
    j=Q(n,2);v=j*(j+1)/(3*h*h)
    rows=[]
    for i in range(q):
        x=(j-offset-2*i)/h
        w=[(v+x*x)/4-Q(35,36)*x+Q(5,6),-(v+x*x)/4+Q(5,12)*x+Q(1,2),
           -(v+x*x)/4+Q(13,12)*x-Q(1,2),(v+x*x)/4-Q(19,36)*x+Q(1,6)]
        if rho<8:w=[a+(x-Q(3,2))/h*b for a,b in zip(w,[-Q(1,3),1,-1,Q(1,3)])]
        r=[offset+2*i+l*h for l in range(4)]
        assert all(a>0 for a in w)
        assert sum(w)==1
        assert sum(a*(j-b)for a,b in zip(w,r))==0
        assert sum(a*(j-b)**2 for a,b in zip(w,r))==j*(j+1)/3
        rows.append({'indices':r,'weights':[str(a)for a in w]})
    edges=[]
    for i in range(q-1):
        rr=rows[i]['indices']
        ww=[Q(a)for a in rows[i]['weights']]
        ww_next=[Q(a)for a in rows[i+1]['weights']]
        squares=[ww[l]*ww_next[l]*(r+1)*(r+2)*(n-r)*(n-r-1)for l,r in enumerate(rr)]
        edges.append({'fromRow':i,'toRow':i+1,'squaredLengths':[str(a)for a in squares]})
    return {'n':n,'k':q,'t':2,'h':h,'offset':offset,'rows':rows,
        'exactBasisDefinition':{
          'convention':'|u_i> = sum_l sqrt(w_i,l) eta_i,l |D_N^(r_i,l)>',
          'initialPhases':'eta_0,l = 1',
          'edgeLengths':'(a,b,c,d) = componentwise positive square roots of squaredLengths',
          'diagonal':'D=(max(abs(a-b),abs(c-d))+min(a+b,c+d))/2',
          'pair':'Z(u,v)=(C_u+i sqrt(1-C_u^2), C_v-i sqrt(1-C_v^2)), with C_u=(D^2+u^2-v^2)/(2Du), C_v=(D^2+v^2-u^2)/(2Dv)',
          'increments':'(z_0,z_1)=Z(a,b), (z_2,z_3)=-Z(c,d)',
          'recursion':'eta_(i+1),l=eta_i,l z_l',
          'scope':'These finite radical recurrences define exact coefficients; decimal rendering is not the certificate.'},
        'phaseEdges':edges}

def basis(n):
    record=data(n);q=record['k']
    r=np.array([a['indices']for a in record['rows']])
    w=np.array([[float(Q(a))for a in b['weights']]for b in record['rows']])
    phases=np.ones((q,4),dtype=complex);margins=[]
    for i in range(q-1):
        rr=r[i]
        L=np.sqrt(w[i]*w[i+1]*(rr+1)*(rr+2)*(n-rr)*(n-rr-1))
        a,b,c,d=L;D=(max(abs(a-b),abs(c-d))+min(a+b,c+d))/2
        def pair(u,v):
            cu=(D*D+u*u-v*v)/(2*D*u);cv=(D*D+v*v-u*u)/(2*D*v)
            return cu+1j*np.sqrt(1-cu*cu),cv-1j*np.sqrt(1-cv*cv)
        z0,z1=pair(a,b);z2,z3=pair(c,d);zz=np.array([z0,z1,-z2,-z3])
        assert abs(L@zz)<1e-7
        phases[i+1]=phases[i]*zz
        margins.append(float(sum(L)-2*max(L)))
    B=np.zeros((n+1,q),complex)
    for i in range(q):B[r[i],i]=np.sqrt(w[i])*phases[i]
    return B,record,margins

def check(n):
    B,record,margins=basis(n);q=record['k'];j=n/2
    z=j-np.arange(n+1);Jp=np.diag(np.sqrt(np.arange(1,n+1)*(n+1-np.arange(1,n+1))),1)
    A=[B.conj().T@B-np.eye(q), B.conj().T@(z[:,None]*B),
       B.conj().T@((z*z)[:,None]*B)-j*(j+1)/3*np.eye(q),B.conj().T@Jp@B,
       B.conj().T@(Jp@Jp)@B,B.conj().T@((z[:,None]+z[None,:])*Jp)@B]
    residual=max(float(np.max(abs(a)))for a in A)
    assert residual<1e-8
    return {'n':n,'k':q,'maxAbsoluteCompressionResidual':residual,'minimumQuadrilateralMargin':min(margins)}

if __name__=='__main__':
    records=[check(n)for n in range(17,161)]
    directory=Path(__file__).resolve().parent
    (directory/'families-render-checks.json').write_text(json.dumps({'scope':'Numerical renderings only; theorem is proved by the exact polynomial certificate','passed':True,'records':records},indent=2)+'\n')
    (directory/'families-example-weights.json').write_text(json.dumps({'examples':[data(n)for n in [25,33,41,49,57,65,73,81]]},indent=2)+'\n')
    print('Direct full spin compressions passed N=17,...,160.')
