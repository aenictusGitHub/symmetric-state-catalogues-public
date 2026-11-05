#!/usr/bin/env python3
"""Independent exact audit of the all-N auxiliary-spin dimension bound.

Requires SymPy only. The universal proof uses a rank-3/2 auxiliary spin and a
quadratic Hermitian form. Exact finite spin blocks cross-check the angular-
momentum eigenvalues and their multiplicities; no optimizer is invoked.
"""
import hashlib,json,time
from pathlib import Path
import sympy as s
BASE=Path(__file__).resolve().parent

def verify():
    start=time.time();n=s.symbols('n',integer=True,positive=True);y=s.symbols('y',nonnegative=True)
    spin=n/2;aux=s.Rational(3,2);Q=n*(n+2)
    total=[spin-aux+i for i in range(4)]
    xs=[s.factor((J*(J+1)-spin*(spin+1)-aux*(aux+1))/2)for J in total]
    expected_x=[-3*(n+2)/4,-(n+8)/4,(n-6)/4,3*n/4]
    assert all(s.simplify(a-b)==0 for a,b in zip(xs,expected_x))
    ws=[s.factor(x*x-(n-4)*x/2-5*Q/16)for x in xs]
    expected_w=[(n+2)*(5*n-3)/8,-n*(n-7)/8,-(3*n*n+n+6)/8,-n*(n-7)/8]
    assert all(s.factor(a-b)==0 for a,b in zip(ws,expected_w))
    dims=[s.expand(2*J+1)for J in total];assert dims==[n-2,n,n+2,n+4]
    # These identities establish strict inertia for every n>=8, not only
    # for the finite spin matrices subsequently checked.
    for p in [ws[0],-ws[1],-ws[2],-ws[3]]:
        assert all(c>0 for c in s.Poly(p.subs(n,y+8),y).all_coeffs())
    assert s.expand(sum(dims))==4*(n+1)
    assert s.simplify(spin*(spin+1)*aux*(aux+1)/3-5*Q/16)==0
    # Independent operator-algebra check: antisymmetric parts in X^2 supply
    # -X/2, leaving only rank-one and traceless rank-two input-spin tensors.
    for a in range(3):
        for b in range(3):
            assert sum(s.LeviCivita(i,j,a)*s.LeviCivita(i,j,b)for i in range(3)for j in range(3))==2*int(a==b)
    Sz=s.diag(s.Rational(3,2),s.Rational(1,2),-s.Rational(1,2),-s.Rational(3,2));Sp=s.zeros(4)
    for a in range(1,4):Sp[a-1,a]=s.sqrt(a*(4-a))
    Sx=(Sp+Sp.T)/2;Sy=(Sp-Sp.T)/(2*s.I);S=[Sx,Sy,Sz]
    assert sum((A*A for A in S),s.zeros(4))==s.Rational(15,4)*s.eye(4)
    for a in range(3):
        for b in range(3):assert S[a]*S[b]-S[b]*S[a]==s.I*sum((s.LeviCivita(a,b,c)*S[c]for c in range(3)),s.zeros(4))
    # For any 2-AC projector, P J_a P=0 and P J_a J_b P=c delta_ab P.
    # The following exact replacement is the whole compression of W.
    compression=spin*(spin+1)/3*sum((A*A for A in S),s.zeros(4))-5*Q/16*s.eye(4)
    assert compression.applyfunc(s.simplify)==s.zeros(4)
    # Exact physical spin matrices in each conserved total-magnetic block.
    lam=s.symbols('lambda');block_count=0;controls=[]
    for N in range(3,19):
        j=s.Rational(N,2);roots=[x.subs(n,N)for x in xs];Js=[J.subs(n,N)for J in total];counts=[0]*4
        for rsum in range(N+4):
            states=[(r,rsum-r)for r in range(N+1)if 0<=rsum-r<4];ix={v:i for i,v in enumerate(states)};B=s.zeros(len(states))
            for col,(r,a)in enumerate(states):
                B[col,col]=(j-r)*(aux-a)
                if r>0 and a<3:B[ix[r-1,a+1],col]=s.sqrt(r*(N-r+1)*(a+1)*(3-a))/2
                if r<N and a>0:B[ix[r+1,a-1],col]=s.sqrt((r+1)*(N-r)*a*(4-a))/2
            M=j+aux-rsum
            active=[i for i,J in enumerate(Js)if J>=abs(M)]
            assert len(active)==len(states)
            expected=s.prod(lam-roots[i]for i in active)
            actual=B.charpoly(lam).as_expr()
            assert s.Poly(s.expand(actual-expected),lam).is_zero
            for i in active:counts[i]+=1
            block_count+=1
        assert counts==[int(v.subs(n,N))for v in dims]
        signs=[s.sign(w.subs(n,N))for w in ws]
        if N>=8:assert signs==[1,-1,-1,-1]
        controls.append({'n':N,'sectorDimensions':counts,'WSigns':list(map(int,signs))})
    # The known N=7 dimension-two construction is a boundary control: W has
    # zero sectors there, so the strict n>=8 argument correctly does not apply.
    assert[ws[i].subs(n,7)for i in [1,3]]==[0,0]
    catalogue=BASE.parent.parent/'data'/'catalogue.json';raw=catalogue.read_bytes();data=json.loads(raw)
    contradictions=[];improvements=[];tight=[]
    for cell in data['entries']:
        N=cell['n'];t=cell['t']
        if N<8 or t<2:continue
        upper=(N-2)//4
        if cell['lower']>upper:contradictions.append({'n':N,'t':t,'lower':cell['lower'],'upper':upper})
        if upper<cell['upper']:improvements.append({'n':N,'t':t,'lower':cell['lower'],'oldUpper':cell['upper'],'newUpper':upper,'newlySettled':upper==cell['lower']})
        if t==2 and cell['lower']==upper:tight.append({'n':N,'dimension':upper})
    assert not contradictions
    for N,k in [(13,2),(16,3),(17,3)]:
        cell=next(c for c in data['entries']if(c['n'],c['t'])==(N,2));assert cell['lower']==k and(N-2)//4==k
    return {'passed':True,'domain':'Every integer N>=8; every AC order t>=2','bound':'4k<=N-2',
      'proof':'Compression of W vanishes on the 4k-dimensional code tensor auxiliary-spin space. W has no zero eigenvalues and positive index N-2. Projection onto its positive spectral space is injective on a totally isotropic subspace.',
      'auxiliarySpin':'3/2','polynomial':'W=X^2-(N-4)X/2-5N(N+2)I/16','symbolicEigenvalues':list(map(str,ws)),
      'symbolicSectorDimensions':list(map(str,dims)),'universalStrictSignCheck':True,'wholeSubspaceCompressionIdentityVerified':True,
      'finiteExactTotalMagneticBlocksChecked':block_count,'finiteSpinControls':controls,
      'n13MaximumDimension':2,'newlySettledExactExamples':[{'n':13,'t':2,'dimension':2},{'n':16,'t':2,'dimension':3},{'n':17,'t':2,'dimension':3}],
      'catalogueSha256':hashlib.sha256(raw).hexdigest(),'exactWitnessContradictions':contradictions,'upperImprovementsStillToApply':improvements,'tightCatalogueSecondOrderCells':tight,'runtimeSeconds':time.time()-start}

if __name__=='__main__':
    report=verify();(BASE/'independent-verification.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps({k:v for k,v in report.items()if k!='finiteSpinControls'},indent=2))
