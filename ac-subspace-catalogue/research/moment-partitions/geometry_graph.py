#!/usr/bin/env python3
"""Exact radical construction for two diagonal scalar compressions.
Uses elementary Fan-Pall interlacing; no floating point search is needed.
"""
from fractions import Fraction as F
import sympy as s

def common_spectrum(a,b,k):
    a=sorted(a);b=sorted(b)
    if k>min(len(a),len(b)):return None
    lam=[max(a[i],b[i]) for i in range(k)]
    if all(lam[i]<=min(a[len(a)-k+i],b[len(b)-k+i]) for i in range(k)):return lam
    return None

def codimension_one(a,b):
    """Q.T Q=I, Q.T diag(a) Q=diag(b), for interlacing distinct spectra."""
    a=list(map(s.Rational,a));b=list(map(s.Rational,b));m=len(a)
    assert len(b)==m-1 and all(a[i]<=b[i]<=a[i+1] for i in range(m-1))
    assert len(set(a))==m and len(set(b))==m-1
    shared=set(a)&set(b);active=[i for i,x in enumerate(a) if x not in shared]
    ar=[a[i] for i in active];br=[x for x in b if x not in shared]
    weights=[s.prod(x-y for y in br)/s.prod(x-z for z in ar if z!=x) for x in ar]
    assert all(x>0 for x in weights)
    Q=s.zeros(m,m-1)
    for j,y in enumerate(b):
        if y in shared:Q[a.index(y),j]=1;continue
        norm=s.sqrt(sum(w/(x-y)**2 for w,x in zip(weights,ar)))
        for i,x,w in zip(active,ar,weights):Q[i,j]=s.sqrt(w)/(x-y)/norm
    return Q

def compress(a,lam):
    """Radical isometry realizing a rational interlacing compression spectrum."""
    a=list(map(s.Rational,a));lam=list(map(s.Rational,lam));m=len(a);k=len(lam)
    assert all(a[i]<=lam[i]<=a[m-k+i] for i in range(k))
    if m==k:
        assert a==lam;return s.eye(k)
    r=m-1-k
    mid=[max(a[i],lam[i-r]) if i>=r else a[i] for i in range(m-1)]
    return codimension_one(a,mid)*compress(mid,lam)

def prescription(n,residue=0,k=None):
    levels=[F(n,2)-r for r in range(residue,n+1,3)]
    c=F(n*(n+2),12)
    eps=F(1,100*(n+1)*(c+1)) if 0 in levels else F(0)
    first=[x+eps*(x*x-c) for x in levels]
    assert all(x for x in first)
    pos=sorted([i for i,x in enumerate(first) if x>0],key=lambda i:(levels[i]**2-c)/first[i])
    neg=sorted([i for i,x in enumerate(first) if x<0],key=lambda i:(levels[i]**2-c)/first[i])
    a=[(levels[i]**2-c)/first[i] for i in pos]
    b=[(levels[i]**2-c)/first[i] for i in neg]
    if k is None:k=max([z for z in range(1,min(len(a),len(b))+1) if common_spectrum(a,b,z) is not None],default=0)
    lam=common_spectrum(a,b,k)
    if lam is None:return None
    return dict(n=n,residue=residue,k=k,levels=levels,c=c,eps=eps,first=first,pos=pos,neg=neg,a=a,b=b,lam=lam)

def construct(n,residue=0,k=None):
    p=prescription(n,residue,k);assert p is not None
    U=compress(p['a'],p['lam']);V=compress(p['b'],p['lam']);Fmat=s.zeros(len(p['levels']),p['k'])
    for i,ix in enumerate(p['pos']):Fmat[ix,:]=U[i,:]/s.sqrt(s.Rational(p['first'][ix]))
    for i,ix in enumerate(p['neg']):Fmat[ix,:]=V[i,:]/s.sqrt(-s.Rational(p['first'][ix]))
    # Gram-Schmidt uses only field operations and square roots.
    columns=[]
    for col in range(p['k']):
        v=Fmat[:,col]
        for u in columns:v=v-u*(u.T*v)[0]
        norm=s.sqrt(s.simplify((v.T*v)[0]));columns.append(v/norm)
    return p,s.Matrix.hstack(*columns)

if __name__=='__main__':
    for n in (12,13,21,22,30,31,35,39,40,42,48,49,120):
        print(n,[prescription(n,r)['k'] for r in range(3)])
