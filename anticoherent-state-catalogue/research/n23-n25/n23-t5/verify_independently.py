"""Verify a rational Banach certificate for an at least 5-AC 23-qubit state.

Standard library only. The reduced-state equations are reconstructed directly;
the square Jacobian is reconstructed by exact central finite differences.
No numerical library, numerical tolerance, or discovery import is used.
"""
from pathlib import Path
from fractions import Fraction as Q
from math import comb
import json,hashlib,time

here=Path(__file__).resolve().parent
start=time.monotonic()
raw=(here/'certificate.json').read_bytes();cert=json.loads(raw)
n,t=cert['n'],cert['t'];assert(n,t)==(23,5)
d=n+1;m=(t+1)**2
x0=list(map(Q,cert['center']));scales=cert['scaleFactors']
assert len(x0)==2*d and len(scales)==d
assert all(isinstance(x,int) and x>0 and x&(x-1)==0 for x in scales)
active=cert['activeVariables'];fixed=cert['fixedVariables']
assert len(active)==m and len(fixed)==2*d-m
assert sorted(active+fixed)==list(range(2*d))

def reduction_equations(x):
    # Division by sqrt(C(t,a)C(t,b)) removes all irrational coefficients.
    # c_r = sqrt(C(n,r))/d_r * (x_r+i*x_(d+r)).
    diagonal=[];upper=[]
    for a in range(t+1):
        diagonal.append(sum(Q(comb(n-t,l),scales[a+l]**2)*
            (x[a+l]**2+x[d+a+l]**2) for l in range(n-t+1))
            -Q(1,(t+1)*comb(t,a)))
        for b in range(a+1,t+1):
            real=Q(0);imag=Q(0)
            for l in range(n-t+1):
                i,j=a+l,b+l;w=Q(comb(n-t,l),scales[i]*scales[j])
                real+=w*(x[i]*x[j]+x[d+i]*x[d+j])
                imag+=w*(x[d+i]*x[j]-x[i]*x[d+j])
            upper.extend([real,imag])
    return diagonal+upper

f0=reduction_equations(x0)
assert len(f0)==m
# Quadratic central differences are exact, including for nonzero fixed values.
columns=[]
for k in active:
    plus=x0.copy();minus=x0.copy();plus[k]+=1;minus[k]-=1
    fp,fm=reduction_equations(plus),reduction_equations(minus)
    columns.append([(p-q)/2 for p,q in zip(fp,fm)])
bits=cert['inverseDyadicBits'];assert isinstance(bits,int) and bits>0
A=[[Q(int(x),1<<bits) for x in row] for row in cert['inverseNumerators']]
assert len(A)==m and all(len(row)==m for row in A)
eta=max(abs(sum(a*f for a,f in zip(row,f0))) for row in A)
beta=max(sum(abs(Q(i==j)-sum(a*x for a,x in zip(A[i],columns[j])))
    for j in range(m)) for i in range(m))
anorm=max(sum(abs(x) for x in row) for row in A)
# Each equation has two quadratic monomials per l. Each monomial contributes
# at most twice its absolute coefficient to the sum of absolute Hessian entries.
# Scales are >=1, so H <= 4*sum_l C(n-t,l) = 2**(n-t+2).
hbound=Q(1<<(n-t+2))
radius=Q(cert['radius']);assert radius>0
q=beta+anorm*hbound*radius
assert beta<1 # A J0 is invertible, hence A is invertible.
assert 0<=q<1
assert eta+q*radius<radius
assert eta<Q(4,10**66) and beta<Q(1,10**76) and anorm<1600 and q<Q(1,10**55)
for name,value in [('residual',eta),('inverseDefect',beta),('inverseInfinityNorm',anorm)]:
    assert value==Q(cert['bounds'][name])
report={'passed':True,'n':n,'t':t,'equationCount':m,'variableCount':2*d,
    'activeVariableCount':m,'fixedVariableCount':len(fixed),
    'certificateSHA256':hashlib.sha256(raw).hexdigest(),
    'method':'Reconstructed rational reduced-state equations; exact central-difference Jacobian; exact residual and approximate-inverse defect; universal integer Hessian bound; Banach contraction and strict self-mapping.',
    'proofBounds':{'eta':str(eta),'beta':str(beta),'inverseInfinityNorm':str(anorm),
        'universalHessianBound':str(hbound),'radius':str(radius),'contractionFactor':str(q),
        'selfMappingMargin':str(radius-eta-q*radius)},
    'conclusion':'A unique zero exists in the stated active-coordinate ball with the twelve other coordinates fixed. The exactly defined state has rho_5=I_6/6 and is normalized.',
    'coefficientEvidence':'Exact implicit algebraic definition by a rational polynomial system and isolating box; saved floating components are approximations, not closed-form radicals.',
    'elapsedSeconds':time.monotonic()-start}
(here/'independent-verification.json').write_text(json.dumps(report,indent=2)+'\n')
print(json.dumps({k:v for k,v in report.items() if k!='proofBounds'},indent=2))
print('Verified conservative bounds:',{'eta<4e-66':eta<Q(4,10**66),'beta<1e-76':beta<Q(1,10**76),'inverseNorm<1600':anorm<1600,'contraction<1e-55':q<Q(1,10**55)},flush=True)
