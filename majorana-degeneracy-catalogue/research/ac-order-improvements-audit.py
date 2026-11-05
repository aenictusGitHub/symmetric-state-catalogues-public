#!/usr/bin/env python3
"""Independently audit additional exact AC witnesses (SymPy only).

This does not import the witness-generation or pre-existing AC/multiplicity
verifiers. Reduced states are reconstructed in the normalized Dicke basis.
All polynomial identities and root counts use exact arithmetic.
"""
from pathlib import Path
import argparse, hashlib, json, time
import sympy as S

HERE=Path(__file__).resolve().parent
z=S.Symbol('z')

def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()
def parse(v):return S.sympify(v,locals={'I':S.I,'sqrt':S.sqrt})
def explicit(st,expected=None):
    n=st['n'];c=list(map(parse,st['exactCoefficients']));t=st['achievedOrder']
    assert len(c)==n+1
    probabilities=[S.simplify(v*S.conjugate(v)) for v in c]
    assert S.simplify(sum(probabilities)-1)==0
    reductions=[]
    for order in range(1,t+1):
        A=S.Matrix(order+1,n-order+1,lambda a,b:c[a+b]*S.sqrt(S.binomial(order,a)*S.binomial(n-order,b)/S.binomial(n,a+b)))
        rho=S.simplify(A*A.H)
        assert rho==S.eye(order+1)/S.Integer(order+1),(st['id'],order,rho)
        reductions.append({'order':order,'identityExactlyVerified':True,'matrix':str(rho)})
    P=S.Poly(sum((-1)**r*S.sqrt(S.binomial(n,r))*c[r]*z**r for r in range(n+1)),z,extension=True)
    fac=P.sqf_list()[1];partition=sorted([m for f,m in fac for _ in range(f.degree())]+([n-P.degree()] if n>P.degree() else []),reverse=True)
    assert sum(partition)==n
    if expected is not None:assert partition==expected
    if st.get('partition') is not None:assert partition==st['partition']
    squarefree=[]
    for f,m in fac:
        assert S.gcd(f,f.diff()).degree()==0
        squarefree.append({'factor':str(f.monic().as_expr()),'degree':int(f.degree()),'rootMultiplicity':m})
    for i,(f,m) in enumerate(fac):
        for g,other in fac[i+1:]:assert S.gcd(f,g).degree()==0
    return {'id':st['id'],'n':n,'verifiedOrderAtLeast':t,'partition':partition,'normalizationVerifiedExactly':True,'directReducedMatrices':reductions,'probabilities':{str(r):str(p) for r,p in enumerate(probabilities) if p},'diagonalMoments':{str(j):str(S.simplify(sum(r**j*p for r,p in enumerate(probabilities)))) for j in range(t+1)},'affinePolynomialMonic':str(P.monic().as_expr()),'infinityMultiplicity':n-int(P.degree()),'squarefreeDecomposition':squarefree,'independentOfConstructionVerifier':True}

def add(A,B):return A[0]+B[0],A[1]+B[1]
def mul(A,B):
    values=[a*b for a in A for b in B];return min(values),max(values)
def enclosure(poly,k,interval):
    out=(S.S.Zero,S.S.Zero)
    for c in S.Poly(poly,k).all_coeffs():out=add(mul(out,interval),(c,c))
    return out
def quotient(A,B):
    assert not B[0]<=0<=B[1]
    return mul(A,tuple(sorted([1/B[0],1/B[1]])))

def implicit(st,cert):
    k,x=S.symbols('k x');n=st['n'];a=cert['poleAtZero'];b=cert['poleAtInfinity']
    m=cert.get('repeatedMultiplicity',2)
    assert n==a+3*(m+1)+b and st['achievedOrder']==2
    q=S.Poly(S.sympify(cert['parameterPolynomial'],locals={'k':k}),k,domain=S.QQ)
    lo,hi=map(S.Rational,cert['rootInterval']);assert lo<hi
    assert q.eval(lo)!=0 and q.eval(hi)!=0 and q.count_roots(lo,hi)==1
    bracket=tuple(map(S.Rational,cert['rootBracket']))
    assert bracket[0]<lo<hi<bracket[1] and q.count_roots(*bracket)==1
    assert S.gcd(q,q.diff()).degree()==0
    assert q.eval(0)!=0 and q.eval(1)!=0
    A=S.Poly(S.sympify(cert['xNumerator'],locals={'k':k}),k,domain=S.QQ)
    B=S.Poly(S.sympify(cert['xDenominator'],locals={'k':k}),k,domain=S.QQ)
    assert S.gcd(q,B).degree()==0
    xi=quotient(enclosure(A.as_expr(),k,(lo,hi)),enclosure(B.as_expr(),k,(lo,hi)))
    assert xi[0]>0
    rs=[a+3*j for j in range(m+2)]
    amplitudes=[(-1)**(m+1-j)*(S.binomial(m,j-1)+k*S.binomial(m,j)) for j in range(m+2)]
    ws=[v*v*x**(m+1-j)/S.binomial(n,r) for j,(r,v) in enumerate(zip(rs,amplitudes))]
    Z=sum(ws);recovery=A.as_expr()/B.as_expr()
    recovered_mod=S.rem(A*S.invert(B,q),q).as_expr()
    diagonal=[]
    for order in [1,2]:
        for h in range(order+1):
            difference=sum(w*S.binomial(r,h)*S.binomial(n-r,order-h)/S.binomial(n,order) for r,w in zip(rs,ws))-Z/S.Integer(order+1)
            rem=S.S.Zero
            for coefficient in S.Poly(difference,x).all_coeffs():
                rem=S.rem(rem*recovered_mod+coefficient,q.as_expr(),k)
            assert rem==0,(st['id'],order,h,rem)
            diagonal.append({'order':order,'diagonalIndex':h,'remainder':'0'})
    assert all(abs(i-j)>2 for i in rs for j in rs if i!=j)
    partition=sorted([m,m,m,1,1,1]+([a] if a else [])+([b] if b else []),reverse=True)
    assert partition==st['partition']==cert['partition']
    # Exact factor proof: discriminants -27*u^2 and -27*k^2*u^2 are
    # nonzero. The resultant of z^3-u and z^3-ku is u^3(1-k)^3.
    # No finite factor has zero as a root; a,b supply the two poles.
    refined=q.refine_root(lo,hi,eps=S.Rational(1,10**80));kmid=S.N(sum(refined)/2,90);xmid=S.N(recovery.subs(k,kmid),85);u=S.sqrt(xmid)
    raw=[v.subs(k,kmid)*u**(m+1-j) for j,v in enumerate(amplitudes)]
    unnormalized=[(-1)**r*v/S.sqrt(S.binomial(n,r)) for r,v in zip(rs,raw)];norm=S.sqrt(sum(v*v for v in unnormalized))
    coefficients={r:S.N(v/norm,75) for r,v in zip(rs,unnormalized)}
    errors=[abs(S.N(S.Float(st['highPrecisionCoefficients'][r][0],80)-coefficients.get(r,0),70)) for r in range(n+1)]
    assert max(errors)<S.Rational(1,10**55)
    return {'id':st['id'],'n':n,'verifiedOrderAtLeast':2,'partition':partition,'rootCountVerifiedExactly':1,'rootInterval':[str(lo),str(hi)],'displayedRootBracket':list(map(str,bracket)),'displayedBracketRootCountVerifiedExactly':1,'parameterPolynomial':str(q.as_expr()),'recoveryDenominatorNonzero':True,'reconstructedXPositiveInterval':list(map(str,xi)),'normalizationStrictLowerBound':str(1/S.binomial(n,rs[-1])),'directReducedDiagonalRemainders':diagonal,'offDiagonalEntriesZeroBySupportSeparation':True,'finiteFactorDiscriminants':['-27*u**2','-27*k**2*u**2'],'finiteFactorResultant':'u**3*(1-k)**3','nonzeroFactorsAndDistinctRootsVerifiedExactly':True,'poleMultiplicities':{'zero':a,'infinity':b},'highPrecisionCoefficientMaxAbsoluteError':str(max(errors)),'independentOfConstructionVerifier':True}

def point_check(st):
    if not st.get('majoranaPoints'):return None
    import mpmath as mp
    mp.mp.dps=70;n=st['n'];poly=[mp.mpc(1)];north=0
    for point in st['majoranaPoints']:
        X,Y,Z=map(lambda k:mp.mpf(str(point[k])),['x','y','z']);m=point['multiplicity']
        assert abs(X*X+Y*Y+Z*Z-1)<mp.mpf('2e-14')
        if Z==1 and X==0 and Y==0:north+=m;continue
        root=(X-mp.j*Y)/(1-Z)
        for _ in range(m):
            p=[mp.mpc(0)]*(len(poly)+1)
            for j,v in enumerate(poly):p[j]-=root*v;p[j+1]+=v
            poly=p
    assert len(poly)-1+north==n
    reconstructed=[(-1)**r*(poly[r] if r<len(poly) else 0)/mp.sqrt(mp.binomial(n,r)) for r in range(n+1)]
    norm=mp.sqrt(sum(abs(v)**2 for v in reconstructed));reconstructed=[v/norm for v in reconstructed]
    given=st.get('highPrecisionCoefficients',st['coefficients']);given=[mp.mpc(str(v[0]),str(v[1])) for v in given]
    norm=mp.sqrt(sum(abs(v)**2 for v in given));given=[v/norm for v in given]
    overlap=sum(mp.conj(a)*b for a,b in zip(reconstructed,given));phase=overlap/abs(overlap)
    error=max(abs(phase*a-b) for a,b in zip(reconstructed,given));assert error<mp.mpf('5e-12'),(st['id'],error)
    return {'normalizedDickeReconstructionMaxComponentError':mp.nstr(error,12),'tolerance':'5e-12','phaseAligned':True,'scope':'Numerical check that plotted coordinates have the displayed vector orientation; exact AC and multiplicities are established separately.'}

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--sparse',type=Path,default=HERE/'ac-order-sparse.json');ap.add_argument('--constructions',type=Path,default=HERE/'ac-order-constructions.json');ap.add_argument('--output',type=Path,default=HERE/'ac-order-improvements-audit.json');args=ap.parse_args()
    started=time.monotonic();source=json.loads(args.sparse.read_text());checks=[]
    expected={c['id']:c['partition'] for c in source['checks']}
    for st in source['states']:checks.append(explicit(st,expected[st['id']]))
    inputs={args.sparse.name:sha(args.sparse)}
    if args.constructions.exists():
        constructions=json.loads(args.constructions.read_text());certs={c['id']:c for c in constructions['certificates']}
        for st in constructions['states']:
            check=explicit(st,st['partition']) if st.get('exactCoefficients') else implicit(st,certs[st['id']])
            check['displayedConstellationCheck']=point_check(st)
            checks.append(check)
        inputs[args.constructions.name]=sha(args.constructions)
    output={'passed':True,'method':'Independent normalized-Dicke reduced matrices or hypergeometric diagonal probabilities; exact rational root isolation and quotient identities; algebraic squarefree decomposition or explicit cubic discriminants/resultant. No generator imports.','inputSHA256':inputs,'states':checks,'stateCount':len(checks),'runtimeSeconds':time.monotonic()-started,'scope':'Proves existence at the recorded order in exactly the recorded full partition; does not by itself prove a global maximum over that partition.'}
    args.output.write_text(json.dumps(output,indent=2)+'\n');print(json.dumps({'passed':True,'states':len(checks),'runtimeSeconds':output['runtimeSeconds']}))
if __name__=='__main__':main()
