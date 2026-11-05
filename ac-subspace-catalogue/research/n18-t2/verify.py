#!/usr/bin/env python3
"""Exact whole-subspace verification of the analytical N=18, k=4, t=2 basis."""
from pathlib import Path
import json,hashlib
import sympy as s
ROOT=Path(__file__).resolve().parent

def split(c,t):
    return s.Matrix(t+1,19-t,lambda a,b:c[a+b]*s.sqrt(s.binomial(t,a)*s.binomial(18-t,b)/s.binomial(18,a+b)))

def zero(x):return s.simplify(s.expand_complex(x))==0

def verify():
    raw=(ROOT/'basis.json').read_bytes();data=json.loads(raw)
    assert (data['n'],data['t'],data['dimension'])==(18,2,4)
    rows=[[s.sympify(z) for z in b['exactCoefficients']] for b in data['basis']]
    B=s.Matrix.hstack(*map(s.Matrix,rows))
    assert all(zero(z) for z in B.H*B-s.eye(4))
    matrices=[split(c,2) for c in rows];count=0
    for i in range(4):
        for j in range(4):
            D=matrices[i]*matrices[j].H-(s.eye(3)/3 if i==j else s.zeros(3))
            for z in D:assert zero(z);count+=1
    assert count==144
    assert rows[0]==rows[0][::-1]
    assert all(zero(a+b) for a,b in zip(rows[1],rows[1][::-1]))
    assert rows[3]==rows[2][::-1]
    weight_lists=[];third_purities=[]
    for c in rows:
        weight_lists.append([str(s.simplify(s.expand_complex(z*s.conjugate(z)))) for z in c])
        M=split(c,3);rho=(M*M.H).applyfunc(lambda z:s.simplify(s.expand_complex(z)))
        third_purities.append(str(s.simplify(s.trace(rho*rho))))
    # Exact inertia upper bound specialized from the general auxiliary-spin theorem.
    j=s.Integer(9);a=s.Rational(3,2)
    spins=[j-a+i for i in range(4)]
    x=[(L*(L+1)-j*(j+1)-a*(a+1))/2 for L in spins]
    ev=[s.factor(v*v-7*v-s.Rational(225,2)) for v in x]
    assert ev==[s.Rational(435,2),-s.Rational(99,4),-s.Rational(249,2),-s.Rational(99,4)]
    positive=sum(int(2*L+1) for L,e in zip(spins,ev) if e>0)
    negative=sum(int(2*L+1) for L,e in zip(spins,ev) if e<0)
    assert (positive,negative)==(16,60)
    return {'passed':True,'n':18,'dimension':4,'order':2,'maximumDimension':4,'basisSHA256':hashlib.sha256(raw).hexdigest(),'orthonormal':True,'complexCrossEntriesExactlyVerified':count,'fullCrossReductionEntriesChecked':count,'allOrderedBasisPairs':16,'basisVectorThirdReductionPurities':third_purities,'squaredDickeAmplitudes':weight_lists,'reflectionStructureExactlyVerified':True,'inertiaUpperBound':{'auxiliarySpin':'3/2','witness':'X^2-7X-225I/2','eigenvalues':list(map(str,ev)),'positiveDimension':positive,'negativeDimension':negative,'conclusion':'4k<=16, hence k<=4'},'arithmetic':'Exact rational and radical arithmetic; all 144 complex cross entries reconstructed from normalized Dicke splitting. No numerical residual establishes existence.'}

if __name__=='__main__':
    report=verify();(ROOT/'verification.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps({k:v for k,v in report.items() if k!='squaredDickeAmplitudes'},indent=2))
