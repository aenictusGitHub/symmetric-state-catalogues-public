#!/usr/bin/env python3
"""Exact obstruction to embedding the N=18 optimal AC2 code in a largest AC1 code."""
from pathlib import Path
import json, hashlib, sympy as s
ROOT=Path(__file__).resolve().parent
DATA=ROOT.parent/'n18-t2/basis.json'
simp=lambda z:s.simplify(s.expand_complex(z))

def main():
    data=json.loads(DATA.read_text())
    basis_sha=hashlib.sha256(DATA.read_bytes()).hexdigest()
    source_verification=json.loads((DATA.parent/'verification.json').read_text())
    assert source_verification['passed'] and source_verification['basisSHA256']==basis_sha
    B=s.Matrix.hstack(*[s.Matrix([s.sympify(v) for v in b['exactCoefficients']]) for b in data['basis']])
    P=B*B.H
    def pp(a,b):return P[a,b] if 0<=a<=18 and 0<=b<=18 else s.S.Zero
    def rr(a,b):
        term=(1+s.Rational((9-a)*(9-b),30))*pp(a,b)
        if a<18 and b<18:term+=s.sqrt((a+1)*(18-a)*(b+1)*(18-b))*pp(a+1,b+1)/60
        if a>0 and b>0:term+=s.sqrt(a*(19-a)*b*(19-b))*pp(a-1,b-1)/60
        return simp(int(a==b)-term)
    # R=I-P-(Jx P Jx+Jy P Jy+Jz P Jz)/30, the complement of four orthogonal error spaces.
    row=[rr(1,b) for b in range(19)]
    weights=[simp(z*s.conjugate(z)) for z in row]
    diag=s.factor(sum((9-r)*z for r,z in enumerate(weights)))
    assert diag==s.Rational(748,405)
    assert row[1]==s.Rational(1496,3645)
    # Three normalized columns of R give a particularly simple basis of its rank-three range.
    Q=s.Matrix.hstack(*[s.Matrix([rr(a,b) for a in range(19)])/s.sqrt(rr(b,b)) for b in [1,0,2]])
    assert all(simp(z)==0 for z in Q.H*Q-s.eye(3))
    assert all(simp(z)==0 for z in B.H*Q)
    Jz=s.diag(*range(9,-10,-1));Jp=s.zeros(19)
    for a in range(18):Jp[a,a+1]=s.sqrt((a+1)*(18-a))
    checks=[B.H*B-s.eye(4),B.H*Jz*B,B.H*Jp*B,
            B.H*Jz**2*B-30*s.eye(4),B.H*Jp**2*B,
            B.H*(Jz*Jp+Jp*Jz)*B]
    assert all(simp(z)==0 for check in checks for z in check)
    R=s.Matrix(19,19,rr)
    assert all(simp(z)==0 for z in R-Q*Q.H)
    assert simp(s.trace(R))==3
    rz=(Q.H*Jz*Q).applyfunc(simp);rp=(Q.H*Jp*Q).applyfunc(simp)
    print('compressed Jz',rz); print('compressed J+',rp)
    assert rz==s.diag(s.Rational(9,2),0,-s.Rational(9,2))
    target=s.zeros(3);target[1,0]=9/s.sqrt(2);target[2,1]=-9*s.I/s.sqrt(2)
    assert all(simp(z)==0 for z in rp-target)
    report={'passed':True,'source':'../n18-t2/basis.json','basisSHA256':basis_sha,'sourceExactCertificateHashMatches':True,'inputAC2SpinConditionsIndependentlyRechecked':True,'complementProjectorExactlyVerified':True,'n':18,'AC2dimension':4,'maximumAC1dimension':7,'errorComplementDimension':3,'nonzeroComplementSpinEntry':'748/405','complementProjectorEntry':'1496/3645','compressedJz':[[str(z) for z in row] for row in rz.tolist()],'compressedJplus':[[str(z) for z in row] for row in rp.tolist()],'largestAC1ParentDimensionOfThisCode':5,'notContainedInAnyMaximumDimensionAC1':True,'scope':'This exact AC2 code; not all N=18 AC2 codes. The general extension criterion applies to every AC2 subspace.'}
    (ROOT/'extension-obstruction.json').write_text(json.dumps(report,indent=2)+'\n')
if __name__=='__main__':main()
