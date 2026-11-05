#!/usr/bin/env python3
"""Exact nested AC2 construction in the reflection-pair AC1 parent.
This certifies a subfamily, not optimality or a classification of AC1 parents.
"""
from pathlib import Path
import json
import sympy as s
ROOT=Path(__file__).resolve().parent

def simp(x): return s.simplify(s.expand_complex(x))

def v(n,x):
    x=s.Rational(x);j=s.Rational(n,2);u=s.zeros(n+1,1)
    u[int(j-x)]=s.sqrt((x+1)/(2*x+1))
    u[int(j+x+1)]=s.sqrt(x/(2*x+1))
    return u

def pair(n,x,y):
    j=s.Rational(n,2);c=j*(j+1)/3;a=x*(x+1);b=y*(y+1)
    assert a<c<b
    return s.sqrt((b-c)/(b-a))*v(n,x)+s.sqrt((c-a)/(b-a))*v(n,y)

def verify(n,columns):
    j=s.Rational(n,2);c=j*(j+1)/3;V=s.Matrix.hstack(*columns);k=len(columns)
    Z=s.diag(*[j-r for r in range(n+1)]);R=s.zeros(n+1)
    for r in range(1,n+1):R[r-1,r]=s.sqrt(r*(n-r+1))
    X=(R+R.H)/2;Y=(R-R.H)/(2*s.I)
    assert (V.H*V-s.eye(k)).applyfunc(simp)==s.zeros(k)
    first=[];second=[]
    for name,A in zip(('x','y','z'),(X,Y,Z)):
        assert (V.H*A*V).applyfunc(simp)==s.zeros(k)
        first.append(name)
    for i,A in enumerate((X,Y,Z)):
        for l,B in enumerate((X,Y,Z)):
            assert (V.H*A*B*V-c*int(i==l)*s.eye(k)).applyfunc(simp)==s.zeros(k)
            second.append([i,l])
    def split(u):
        return s.Matrix(3,n-1,lambda a,b:u[a+b]*s.sqrt(s.binomial(2,a)*s.binomial(n-2,b)/s.binomial(n,a+b)))
    M=[split(u) for u in columns]
    for i,A in enumerate(M):
        for l,B in enumerate(M):assert (A*B.H-int(i==l)*s.eye(3)/3).applyfunc(simp)==s.zeros(3)
    return {'orthonormal':True,'firstMomentMatrices':first,'secondMomentProducts':second,'allTwoQubitCrossReductions':True,'crossReductionEntries':9*k*k}

def main():
    n=18;columns=[v(n,5),pair(n,s.Integer(2),s.Integer(8))]
    checks=verify(n,columns)
    record={'n':18,'dimension':2,'order':2,'maximal':False,'parent':'Canonical seven-dimensional first-order reflection-pair subspace',
        'magneticModules':[2,5,8],'moduleSecondMoments':[6,30,72],'targetSecondMoment':30,
        'basisInParent':['v_5','sqrt(7/11)*v_2+sqrt(4/11)*v_8'],
        'basis':[{'exactCoefficients':[str(simp(c)) for c in u],'support':[r for r,c in enumerate(u) if c!=0]} for u in columns],
        'verification':checks,
        'generalRecipe':{'assumptions':'Choose positive x values spaced at least three apart; take one v_x from each reflected AC1 module. All magnetic support positions must be separated by at least three.',
            'module':'v_x=sqrt((x+1)/(2*x+1))|x>+sqrt(x/(2*x+1))|-x-1>',
            'target':'c=j*(j+1)/3','moduleMoment':'h_x=x*(x+1)',
            'pair':'sqrt((h_y-c)/(h_y-h_x))*v_x+sqrt((c-h_x)/(h_y-h_x))*v_y for h_x<c<h_y',
            'dimension':'number(h=c)+min(number(h<c),number(h>c))',
            'scope':'Exact maximal dimension within this one-vector-per-module subspace; not an unrestricted maximum.'}}
    (ROOT/'nested-examples.json').write_text(json.dumps(record,indent=2)+'\n')
    print(json.dumps({'n':n,'dimension':len(columns),'verified':True,'crossReductionEntries':checks['crossReductionEntries']}))
if __name__=='__main__':main()
