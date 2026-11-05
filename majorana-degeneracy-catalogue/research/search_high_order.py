#!/usr/bin/env python3
"""Reproducible fixed-partition searches; numerical AC is never an exact proof.

P(z)=product_m F_m(z)^m, where F_m has one root for each point of multiplicity m.
The factors are monic and optimized in complex coefficients. Exact Gaussian-
rational gcd checks certify the partition of a serialized candidate polynomial,
while its AC residual remains a numerical statement. Root collisions are reported
explicitly: approaching a more degenerate AC state does not establish existence
inside the requested partition.
"""
import argparse
import json
import math
import time
from collections import Counter
from pathlib import Path

import numpy as np
from scipy.optimize import least_squares


def power(poly, exponent):
    answer=np.array([1.0+0j])
    for _ in range(exponent): answer=np.convolve(answer,poly)
    return answer


class FactorModel:
    def __init__(self, partition, order, allowed=None):
        self.partition=sorted(partition,reverse=True)
        self.n=sum(partition); self.t=order
        self.groups=sorted(Counter(partition).items(),reverse=True)
        self.parameters=[]
        for group,(multiplicity,degree) in enumerate(self.groups):
            for k in (allowed.get(multiplicity,range(degree)) if allowed else range(degree)):
                self.parameters.append((group,k))
        self.dimension=len(self.parameters)
        self.sqrtbinom=np.sqrt([math.comb(self.n,r) for r in range(self.n+1)])
        self.indices=np.array([[a+b for b in range(self.n-order+1)] for a in range(order+1)])
        self.weights=np.array([[math.sqrt(math.comb(order,a)*math.comb(self.n-order,b)/math.comb(self.n,a+b))
                              for b in range(self.n-order+1)] for a in range(order+1)])
        self.cached_x=None

    def factors(self,x):
        result=[np.r_[np.zeros(degree,dtype=complex),1+0j] for _,degree in self.groups]
        values=x[:self.dimension]+1j*x[self.dimension:]
        for value,(g,k) in zip(values,self.parameters): result[g][k]=value
        return result

    def encode(self,factors):
        values=np.array([factors[g][k] for g,k in self.parameters])
        return np.r_[values.real,values.imag]

    def evaluate(self,x):
        if self.cached_x is not None and np.array_equal(x,self.cached_x): return
        factors=self.factors(x); powers=[power(f,m) for f,(m,_) in zip(factors,self.groups)]
        polynomial=np.array([1+0j])
        for p in powers: polynomial=np.convolve(polynomial,p)
        derivatives=[]
        for group,k in self.parameters:
            m,_=self.groups[group]
            derivative=power(factors[group],m-1)*m
            for g,p in enumerate(powers):
                if g!=group: derivative=np.convolve(derivative,p)
            derivative=np.pad(derivative,(k,self.n+1-len(derivative)-k))
            derivatives.append(derivative)
        derivatives=np.array(derivatives).T
        derivatives=np.concatenate([derivatives,1j*derivatives],axis=1)
        z=polynomial/self.sqrtbinom
        dz=derivatives/self.sqrtbinom[:,None]
        norm=np.linalg.norm(z); c=z/norm
        dc=dz/norm-c[:,None]*(np.real(np.conj(c)@dz)/norm)[None,:]
        A=c[self.indices]*self.weights
        dA=dc[self.indices,:]*self.weights[:,:,None]
        rho=A@A.conj().T
        drho=np.einsum('akp,bk->abp',dA,A.conj())+np.einsum('ak,bkp->abp',A,dA.conj())
        residual=[];jac=[]
        for a in range(self.t+1):
            residual.append(rho[a,a].real-1/(self.t+1));jac.append(drho[a,a].real)
            for b in range(a):
                residual.extend([math.sqrt(2)*rho[a,b].real,math.sqrt(2)*rho[a,b].imag])
                jac.extend([math.sqrt(2)*drho[a,b].real,math.sqrt(2)*drho[a,b].imag])
        self.cached_x=x.copy();self.residual=np.array(residual);self.jacobian=np.array(jac)
        self.coefficients=c;self.polynomial=polynomial

    def fun(self,x): self.evaluate(x);return self.residual
    def jac(self,x): self.evaluate(x);return self.jacobian


def separation(factors):
    roots=[root for factor in factors for root in np.roots(factor[::-1])]
    minimum=min((2*abs(a-b)/math.sqrt((1+abs(a)**2)*(1+abs(b)**2))
                 for i,a in enumerate(roots) for b in roots[:i]),default=2.0)
    return float(minimum),roots


def strings(values): return [[format(float(z.real),'.17g'),format(float(z.imag),'.17g')] for z in values]


def certify_partition(factors,groups):
    import sympy as s
    z=s.Symbol('z');polynomials=[]
    for factor in factors:
        coefficients=[s.Rational(a)+s.I*s.Rational(b) for a,b in strings(factor)]
        polynomials.append(s.Poly(sum(a*z**k for k,a in enumerate(coefficients)),z,extension=s.I))
    self_gcd=[s.degree(s.gcd(p,p.diff())) for p in polynomials]
    pair_gcd=[s.degree(s.gcd(a,b)) for i,a in enumerate(polynomials) for b in polynomials[:i]]
    passed=all(v==0 for v in self_gcd+pair_gcd)
    return {'passed':passed,'method':'Exact gcds over Gaussian rationals for the serialized monic factors.',
            'selfGcdDegrees':list(map(int,self_gcd)),'pairwiseGcdDegrees':list(map(int,pair_gcd)),
            'partition':sorted([m for m,d in groups for _ in range(d)],reverse=True) if passed else None,
            'scope':'Exact multiplicities of the serialized factorized polynomial; no exact AC claim.'}


def initial_factors(model,rng,start,case):
    if case=='selected' and start<12:
        # The known order2 tetrahedral construction and an order3 boundary
        # double-octahedron, with perturbations of several magnitudes.
        if start%2==0: base=[np.array([1,0,2j*math.sqrt(3),0,1]),np.array([1,0,-2j*math.sqrt(3),0,1])]
        else: base=[np.array([1,0,-6,0,1],dtype=complex),np.array([1,0,2,0,1],dtype=complex)]
        strength=[0,0.001,0.03,0.15,0.6,1.5][start//2]
        return [p+np.r_[strength*(rng.normal(size=len(p)-1)+1j*rng.normal(size=len(p)-1)),0] for p in base]
    factors=[]
    for m,d in model.groups:
        roots=np.exp(rng.normal(scale=.65,size=d)+1j*rng.uniform(0,2*np.pi,size=d))
        p=np.poly(roots)[::-1].astype(complex)
        allowed={k for g,k in model.parameters if model.groups[g][0]==m}
        for k in range(d):
            if k not in allowed:p[k]=0
        factors.append(p)
    return factors


def run_case(case,partition,order,starts,seed,allowed=None,max_nfev=1200):
    rng=np.random.default_rng(seed);model=FactorModel(partition,order,allowed)
    runs=[];start_time=time.time()
    for start in range(starts):
        factors=initial_factors(model,rng,start,case)
        x0=model.encode(factors)
        # Keep the finite affine chart bounded; this is a search restriction,
        # never used as a nonexistence argument.
        x0=np.clip(x0,-99,99)
        tick=time.time()
        result=least_squares(model.fun,x0,jac=model.jac,bounds=(-100,100),max_nfev=max_nfev,
                             ftol=2e-13,xtol=2e-13,gtol=2e-13)
        model.evaluate(result.x);factors=model.factors(result.x)
        residual=float(np.linalg.norm(model.residual));gap,_=separation(factors)
        run={'start':start,'HSResidual':residual,'minimumDistinctRootChordSeparation':gap,
             'provisionalCandidateForRefinement':residual<1e-10 and gap>1e-3,
             'nearCollision':gap<=1e-3,'evaluations':result.nfev,'termination':result.message,
             'maximumFactorParameter':float(np.max(np.abs(result.x))),
             'seconds':round(time.time()-tick,3),
             'factors':[{'multiplicity':m,'degree':d,'coefficients':strings(p)} for p,(m,d) in zip(factors,model.groups)],
             'normalizedDickeCoefficients':strings(model.coefficients)}
        if run['provisionalCandidateForRefinement']:
            run['partitionCertificate']=certify_partition(factors,model.groups)
        runs.append(run)
        print(case,start,'HS',format(residual,'.4g'),'gap',format(gap,'.4g'),
              'provisional',run['provisionalCandidateForRefinement'],flush=True)
    return {'case':case,'n':sum(partition),'partition':partition,'targetOrder':order,'seed':seed,
            'parameterization':'Product of monic complex factors with the requested multiplicities.',
            'allowedFactorIndices':allowed,'maximumEvaluationsPerStart':max_nfev,
            'starts':starts,'runtimeSeconds':round(time.time()-start_time,3),
            'bestHSResidual':min(r['HSResidual'] for r in runs),
            'provisionalCandidatesForRefinement':sum(r['provisionalCandidateForRefinement'] for r in runs),
            'runs':runs}


def ansatz_certificate():
    import sympy as s
    a,b=s.symbols('a b');A=2*a+b;B=a*a+2*a*b+3;C=a*a*b+4*a+2*b
    F=s.expand(33*A+3*A*B+B*C)
    G=s.expand(8712+12*A*A-8*B*B-3*C*C)
    resultant=s.factor(s.resultant(F,G,b))
    special_gcds={str(value):str(s.Poly(F.subs(a,value),b).gcd(s.Poly(G.subs(a,value),b)).monic().as_expr())
                  for value in (0,6,-6)}
    return {'ansatz':'(z^4+a*z^2+1)^2*(z^4+b*z^2+1), a,b real',
            'constraints':[str(F),str(G)],'resultantEliminatingB':str(resultant),
            'exactCommonFactorInBAtPossibleRealA':special_gcds,
            'realSolutions':[['6','-2'],['-6','2']],
            'interpretation':'The a=0 resultant root is extraneous; the remaining polynomial in a^2 has positive coefficients. Both real solutions make the nominal simple quartic have repeated roots, yielding partition(2,2,2,2,2,2). This excludes only this shared-D2 ansatz, not the full requested partition.'}


def refine_candidate(case,run,max_steps=32):
    """High-precision residual/Jacobian with an SVD pseudoinverse; no existence claim."""
    import mpmath as mp
    mp.mp.dps=85
    model=FactorModel(case['partition'],case['targetOrder'],case['allowedFactorIndices'])
    # JSON integer mapping keys may have become strings when loading a report.
    if case['allowedFactorIndices']:
        allowed={int(k):v for k,v in case['allowedFactorIndices'].items()}
        model=FactorModel(case['partition'],case['targetOrder'],allowed)
    initial=[[mp.mpc(a,b) for a,b in f['coefficients']] for f in run['factors']]
    def conv(a,b):
        out=[mp.mpc(0)]*(len(a)+len(b)-1)
        for i,u in enumerate(a):
            for j,v in enumerate(b):out[i+j]+=u*v
        return out
    def powpoly(a,m):
        out=[mp.mpc(1)]
        for _ in range(m):out=conv(out,a)
        return out
    values=[initial[g][k] for g,k in model.parameters]
    x=[mp.re(v) for v in values]+[mp.im(v) for v in values]
    n,t,d=model.n,model.t,model.dimension
    choose=[mp.sqrt(math.comb(n,r)) for r in range(n+1)]
    weights=[[mp.sqrt(mp.mpf(math.comb(t,a))*math.comb(n-t,k)/math.comb(n,a+k))
              for k in range(n-t+1)] for a in range(t+1)]
    def evaluate(x,with_jac=True):
        factors=[[mp.mpc(0)]*degree+[mp.mpc(1)] for m,degree in model.groups]
        for j,(group,k) in enumerate(model.parameters):factors[group][k]=mp.mpc(x[j],x[j+d])
        powers=[powpoly(p,m) for p,(m,_) in zip(factors,model.groups)]
        polynomial=[mp.mpc(1)]
        for p in powers:polynomial=conv(polynomial,p)
        z=[polynomial[r]/choose[r] for r in range(n+1)]
        norm=mp.sqrt(mp.fsum(abs(v)**2 for v in z));c=[v/norm for v in z]
        dc=None
        if with_jac:
            dz=[]
            for group,k in model.parameters:
                m,_=model.groups[group];p=[v*m for v in powpoly(factors[group],m-1)]
                for h,power_h in enumerate(powers):
                    if h!=group:p=conv(p,power_h)
                p=[mp.mpc(0)]*k+p+[mp.mpc(0)]*(n+1-k-len(p))
                dz.append([p[r]/choose[r] for r in range(n+1)])
            dz += [[mp.j*v for v in p] for p in dz]
            dc=[]
            for p in dz:
                projection=mp.re(mp.fsum(mp.conj(c[r])*p[r] for r in range(n+1)))/norm
                dc.append([p[r]/norm-c[r]*projection for r in range(n+1)])
        F=[];J=[]
        for a in range(t+1):
            for b in range(a+1):
                rho=mp.fsum(weights[a][k]*weights[b][k]*c[a+k]*mp.conj(c[b+k]) for k in range(n-t+1))
                derivative=[]
                if with_jac:
                    derivative=[mp.fsum(weights[a][k]*weights[b][k]*(p[a+k]*mp.conj(c[b+k])+c[a+k]*mp.conj(p[b+k]))
                                        for k in range(n-t+1)) for p in dc]
                if a==b:
                    F.append(mp.re(rho)-mp.mpf(1)/(t+1))
                    if with_jac:J.append([mp.re(v) for v in derivative])
                else:
                    F.extend([mp.sqrt(2)*mp.re(rho),mp.sqrt(2)*mp.im(rho)])
                    if with_jac:J.extend([[mp.sqrt(2)*mp.re(v) for v in derivative],[mp.sqrt(2)*mp.im(v) for v in derivative]])
        return mp.matrix(F),mp.matrix(J) if with_jac else None,factors,c
    history=[];started=time.time()
    for step in range(max_steps):
        F,J,factors,c=evaluate(x)
        residual=mp.norm(F)
        gap,_=separation([np.array([complex(v) for v in p]) for p in factors])
        history.append({'step':step,'HSResidual':mp.nstr(residual,35),'doublePrecisionRootGap':gap})
        print('refine',case['case'],run['start'],step,mp.nstr(residual,6),'gap',gap,flush=True)
        if residual<mp.mpf('1e-65'):break
        # Fix rotational gauge without sending an almost-antipodal root to
        # infinity: constrain the Newton correction to be orthogonal to the
        # three infinitesimal SU(2) actions on the monic factor coefficients.
        tangent=[]
        for generator in range(3):
            values=[]
            for group,k in model.parameters:
                p=factors[group];degree=len(p)-1
                left=(degree-k+1)*p[k-1] if k else mp.mpc(0)
                right=(k+1)*p[k+1]
                value=(left-right-p[k]*p[-2]) if generator==0 else \
                      -mp.j*(left+right-p[k]*p[-2]) if generator==1 else mp.j*(degree-k)*p[k]
                values.append(value)
            vector=[mp.re(v) for v in values]+[mp.im(v) for v in values]
            for previous in tangent:
                projection=mp.fdot(vector,previous)
                vector=[v-projection*w for v,w in zip(vector,previous)]
            norm=mp.sqrt(mp.fdot(vector,vector))
            if norm>mp.mpf('1e-60'):tangent.append([v/norm for v in vector])
        augmentedJ=mp.matrix([[J[i,j] for j in range(2*d)] for i in range(J.rows)]+tangent)
        augmentedF=list(F)+[mp.mpf(0)]*len(tangent)
        U,S,V=mp.svd_r(augmentedJ,full_matrices=False)
        retained=[i for i,s in enumerate(S) if s>S[0]*mp.mpf('1e-65')]
        projection=[mp.fsum(U[k,i]*augmentedF[k] for k in range(len(augmentedF)))/S[i] for i in retained]
        delta=[-mp.fsum(V[i,j]*value for i,value in zip(retained,projection)) for j in range(2*d)]
        accepted=False
        for power in range(12):
            trial=[v+dv/(2**power) for v,dv in zip(x,delta)]
            nextF,_,_,_=evaluate(trial,False)
            if mp.norm(nextF)<residual:x=trial;accepted=True;break
        if not accepted:break
    F,_,factors,c=evaluate(x,False)
    exact_roots=[];root_error=None
    try:
        for p in factors:
            roots=mp.polyroots(list(reversed(p)),maxsteps=600,extraprec=100,error=False)
            exact_roots.extend(roots)
        gap=min(2*abs(a-b)/mp.sqrt((1+abs(a)**2)*(1+abs(b)**2)) for i,a in enumerate(exact_roots) for b in exact_roots[:i])
    except Exception as error:gap=None;root_error=repr(error)
    final_residual=mp.norm(F)
    return {'case':case['case'],'start':run['start'],'workingPrecision':85,'history':history,
            'rotationGauge':'Newton corrections orthogonal to infinitesimal SU(2) directions in the original finite chart.',
            'finalHSResidual':mp.nstr(final_residual,45),'finalRootGap':mp.nstr(gap,35) if gap is not None else None,
            'rootSolverError':root_error,'runtimeSeconds':round(time.time()-started,3),
            'factors':[{'multiplicity':m,'degree':degree,'coefficients':[[mp.nstr(mp.re(v),75),mp.nstr(mp.im(v),75)] for v in p]}
                       for p,(m,degree) in zip(factors,model.groups)],
            'normalizedDickeCoefficients':[[mp.nstr(mp.re(v),75),mp.nstr(mp.im(v),75)] for v in c],
            'numericalOnly':True,'strongSeparatedCandidate':bool(gap is not None and gap>mp.mpf('.001') and final_residual<mp.mpf('1e-60'))}


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,default=Path(__file__).with_name('high-order-search.json'))
    parser.add_argument('--selected-starts',type=int,default=32)
    parser.add_argument('--three-double-starts',type=int,default=24)
    parser.add_argument('--c3-starts',type=int,default=32)
    parser.add_argument('--max-nfev',type=int,default=1200)
    parser.add_argument('--refine-existing',action='store_true')
    parser.add_argument('--refinement-steps',type=int,default=32)
    args=parser.parse_args()
    if args.refine_existing:
        output=json.loads(args.output.read_text());output['highPrecisionRefinements']=[]
        output['restrictedAnsatzCertificate']=ansatz_certificate()
        for case in output['cases']:
            if 'wellSeparatedCandidates' in case:
                case['provisionalCandidatesForRefinement']=case.pop('wellSeparatedCandidates')
            for old_run in case['runs']:
                if 'wellSeparatedNumericalCandidate' in old_run:
                    old_run['provisionalCandidateForRefinement']=old_run.pop('wellSeparatedNumericalCandidate')
            if case['case'] not in ('selected','three-double-general'):continue
            eligible=[r for r in case['runs'] if r['HSResidual']<1e-8 and r['minimumDistinctRootChordSeparation']>1e-5]
            if eligible:
                run=max(eligible,key=lambda r:r['minimumDistinctRootChordSeparation'])
                output['highPrecisionRefinements'].append(refine_candidate(case,run,args.refinement_steps))
                args.output.write_text(json.dumps(output,indent=2)+'\n')
        output['outcome']={
            'totalSearchStarts':sum(c['starts'] for c in output['cases']),
            'highPrecisionDiagnostics':len(output['highPrecisionRefinements']),
            'strongSeparatedNumericalCandidates':sum(r['strongSeparatedCandidate'] for r in output['highPrecisionRefinements']),
            'newExactACWitnesses':0,'selectedRigorousRangeUnchanged':[2,3],
            'conclusion':'No pattern-preserving third-order AC witness was established. The smallest selected-row residuals occur near additional root collisions; higher-precision diagnostics did not converge to a separated high-accuracy state. The selected row remains2–3.',
            'limitations':'This finite search and the restricted-ansatz exclusion are not global nonexistence proofs. Provisional double-precision candidates are not catalogue lower-bound improvements.'}
        args.output.write_text(json.dumps(output,indent=2)+'\n')
        return
    model=FactorModel([2]*4+[1]*4,3)
    rng=np.random.default_rng(571);x=model.encode(initial_factors(model,rng,14,'selected'))
    direction=rng.normal(size=len(x));direction/=np.linalg.norm(direction)
    analytic=model.jac(x)@direction;h=1e-5
    numeric=(model.fun(x+h*direction)-model.fun(x-h*direction))/(2*h)
    jac_error=float(np.max(np.abs(analytic-numeric)))
    assert jac_error<1e-7,jac_error
    output={'scientificStatus':'Exploratory numerical searches, not exact AC existence or global nonexistence proofs.',
            'selectedRowRigorousRange':[2,3],'selectedRowPartition':[2,2,2,2,1,1,1,1],
            'analyticJacobianDirectionalCheckMaximumError':jac_error,'restrictedAnsatzCertificate':ansatz_certificate(),'cases':[]}
    plans=[('selected',[2]*4+[1]*4,args.selected_starts,None,173),
           ('three-double-c3',[2]*3+[1]*6,args.c3_starts,{2:[0],1:[0,3]},1203),
           ('three-double-general',[2]*3+[1]*6,args.three_double_starts,None,271)]
    for case,partition,starts,allowed,seed in plans:
        if starts:
            output['cases'].append(run_case(case,partition,3,starts,seed,allowed,args.max_nfev))
            args.output.write_text(json.dumps(output,indent=2)+'\n')
    print('DONE',[(c['case'],c['provisionalCandidatesForRefinement'],c['bestHSResidual']) for c in output['cases']],flush=True)


if __name__=='__main__':main()
