#!/usr/bin/env python3
"""Independent direct-reduction audit of the D2 obstruction (SymPy only)."""
from pathlib import Path
import hashlib,json,runpy,time
import sympy as s
HERE=Path(__file__).resolve().parent
z,A,B,X,Y,v,U=s.symbols('z A B X Y v U',real=True)
def expression(text):return s.sympify(text,locals={str(q):q for q in[z,A,B,X,Y,v,U]})
def necessary(Q,R):
    P=s.Poly(s.expand(Q**2*R),z);p=[P.nth(r)for r in range(13)]
    norm=s.expand(sum(q*s.conjugate(q)/s.binomial(12,r)for r,q in enumerate(p)))
    r00=sum(s.binomial(9,q)*p[q]*s.conjugate(p[q])/s.binomial(12,q)**2 for q in range(10))
    r02=sum(s.binomial(9,q)*p[q]*s.conjugate(p[q+2])/(s.binomial(12,q)*s.binomial(12,q+2))for q in range(10))
    return [s.Poly(s.expand(q),A,B).clear_denoms()[1].primitive()[1].as_expr()for q in [r00-norm/4,s.re(r02),s.im(r02)]]
def main():
    start=time.monotonic()
    # This is our own independent Schmidt-reduction derivation, not the
    # obstruction author's collective-moment implementation.
    runpy.run_path(str(HERE/'independent-d2-reduction.py'))
    own=json.loads((HERE/'independent-d2-reduction.json').read_text());source=HERE/'analytic-d2-obstruction.json';other=json.loads(source.read_text())
    det=expression(own['linearDeterminant']);assert all(c>0 for c in s.Poly(det,A,B).coeffs())
    assert s.expand(det+expression(other['trivialCharacter']['recoveryDenominator']))==0
    c,d=s.symbols('c d',real=True);locals={'A':A,'B':B,'c':c,'d':d,'v':v}
    rec=[s.sympify(q,locals=locals)for q in own['linearSolution']]
    defect=s.sympify(own['rankOneDefect'],locals=locals)
    assert s.cancel(defect-(rec[2]-rec[0]**2-rec[1]**2))==0
    P=s.Poly(s.cancel(-defect*det**2/12),A,B)
    checks=[]
    for (power,),coef in s.Poly(P.as_expr(),B).terms():
        assert power%2==0
        terms=s.Poly(coef,A).terms();assert all(m[0]%2==0 for m,q in terms)
        q=s.Poly(sum(value*U**(m[0]//2)for m,value in terms),U)
        if power==0:
            assert s.expand(q.as_expr()-U*(U-36)**2*(U+12)**6*(2*U**3+51*U**2+1032*U+3120))==0
        else:
            assert q.eval(0)>0 and q.count_roots(0,s.oo)==0
            recorded=next(r for r in other['trivialCharacter']['positiveCoefficientChecks']if r['powerOfBSquared']==power//2)
            assert s.expand(q.as_expr()-expression(recorded['polynomialInASquared']))==0
            checks.append({'powerOfBSquared':power//2,'positiveAtZero':str(q.eval(0)),'nonnegativeHalfLineRoots':0})
    # The zero-denominator branch A=B=0 cannot solve beta=1 conditions.
    eq=[s.sympify(q,locals=locals)for q in own['equations']]
    assert any(q.subs({A:0,B:0})!=0 and not q.subs({A:0,B:0}).free_symbols for q in eq)
    survivors=[]
    Q=z**4+(A+s.I*B)*z*z+1
    for av in [-6,6]:
        recovered=[s.simplify(q.subs({A:av,B:0}))for q in rec]
        assert recovered==[s.Rational(-av,3),0,4]
        p=s.Poly((Q**2*(z**4+1+(recovered[0]+s.I*recovered[1])*z*z)).subs({A:av,B:0}),z)
        lam=sorted([int(m)for q,m in p.sqf_list()[1]for _ in range(q.degree())],reverse=True);assert lam==[2]*6
        survivors.append({'a':av,'gamma':str(recovered[0]),'partition':lam})
    # beta=0 is retained; two off-diagonal equations force A=B=0.
    beta0=necessary(Q,z*z)
    assert s.groebner(beta0,A,B).polys==s.groebner([A,B],A,B).polys
    # All nontrivial residual characters are one-dimensional and are tested
    # directly, allowing arbitrary complex doubled-orbit parameter A+iB.
    characters=[]
    for R in [z**4-1,z*(z*z-1),z*(z*z+1)]:
        equations=necessary(Q,R);G=s.groebner(equations,A,B)
        assert len(G.polys)==1 and G.polys[0].as_expr()==1
        characters.append({'residual':str(R),'necessaryEquations':list(map(str,equations)),'groebnerBasis':['1']})
    out={'passed':True,'sourceSHA256':hashlib.sha256(source.read_bytes()).hexdigest(),'scope':'Independent audit that no target-partition 3-AC state can have a D2 proper rotational subgroup. The general C1/C2 cases remain unresolved.','normalFormAndAllFourCharactersReviewed':True,'directSchmidtReductionRecomputed':True,'recoveryDeterminantMatchesUpToSign':True,'rankOnePositivityCertificates':checks,'betaZeroBranchRetained':True,'survivingBetaNonzeroCases':survivors,'nontrivialCharacters':characters,'dependency':'The doubled set consisting of two size2 orbits is excluded by the separately verified antipodal obstruction.','globalComplexPartitionExclusionProved':False,'runtimeSeconds':time.monotonic()-start}
    (HERE/'d2-independent-audit.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps({'D2AuditPassed':True,'runtimeSeconds':out['runtimeSeconds']}))
if __name__=='__main__':main()
