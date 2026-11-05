#!/usr/bin/env python3
"""Build exact algebraic certificates from numerical parity candidates."""
from pathlib import Path
from math import comb,log2,sqrt
from fractions import Fraction as F
import json,argparse
import numpy as np
import mpmath as mp
from scipy.linalg import qr
from scipy.optimize import least_squares
from verify_parity_certificate import reconstruct,exact_bounds,verify
HERE=Path(__file__).resolve().parent

def build(n,k,parity=0):
    raw=json.loads((HERE/f'search-parity-n{n}-k{k}-p{parity}.json').read_text())
    B=np.array(raw['basis']);B=B[:,:,0]+1j*B[:,:,1]
    indices=list(range(parity,n+1,2));d=len(indices)
    exponents=[max(0,round(log2(sqrt(comb(n,r)))))for r in indices]
    weights=np.array([sqrt(comb(n,r))/2**e for r,e in zip(indices,exponents)])
    X=B[indices]/weights[:,None];v=np.r_[X.real.ravel(),X.imag.ravel()]
    indices,size,L,eq=reconstruct(n,k,parity,exponents);m=len(eq)
    terms=[(np.array([t[0]for t in e['terms']]),np.array([t[1]for t in e['terms']]),np.array([t[2]/L for t in e['terms']]),e['constant']/L)for e in eq]
    def fun_all(v):return np.array([c+np.sum(w*v[a]*v[b])for a,b,w,c in terms])
    def jac_all(v):
        J=np.zeros((m,size))
        for row,(a,b,w,c)in enumerate(terms):
            np.add.at(J[row],a,w*v[b]);np.add.at(J[row],b,w*v[a])
        return J
    J=jac_all(v);_,R,piv=qr(J,pivoting=True,mode='economic');unknown=piv[:m].tolist()
    assert np.linalg.matrix_rank(J)==m,(n,k,np.linalg.matrix_rank(J),m)
    den=10**16;v=np.array([round(z*den)/den for z in v])
    def full(x):
        result=v.copy();result[unknown]=x;return result
    sol=least_squares(lambda x:fun_all(full(x)),v[unknown],jac=lambda x:jac_all(full(x))[:,unknown],max_nfev=100,ftol=1e-14,gtol=1e-14,xtol=1e-14)
    v=full(sol.x);center=[str(round(z*den))for z in v];v=np.array(list(map(int,center)),float)/den
    inverse=np.linalg.inv(jac_all(v)[:,unknown]);ad0=10**16
    inversenums=[[str(round(z*ad0))for z in row]for row in inverse]
    # A fixed rational approximate inverse suffices for high-precision
    # preconditioned iteration; exact contraction below remains the proof.
    mp.mp.dps=100
    vm=[mp.mpf(a)/den for a in center]
    Am=[[mp.mpf(a)/ad0 for a in row]for row in inversenums]
    for step in range(9):
        fm=[(e['constant']+sum(mp.mpf(w)*vm[a]*vm[b]for a,b,w in e['terms']))/L for e in eq]
        correction=[sum(a*f for a,f in zip(row,fm))for row in Am]
        for slot,change in zip(unknown,correction):vm[slot]-=change
        if max(abs(z)for z in correction)<mp.mpf('1e-90'):break
    den=10**80;center=[str(int(mp.nint(z*den)))for z in vm]
    v=np.array([float(z)for z in vm])
    certificate={'n':n,'dimension':k,'order':2,'parity':parity,'scaleExponents':exponents,'systemScale':str(L),
        'coordinateConvention':'First d*k real coordinates, then d*k imaginary coordinates, each in row-major parity-index then basis-column order. Dicke amplitude at r equals sqrt(binomial(N,r))/2^e_r times its complex coordinate.',
        'unknownCoordinateIndices':unknown,'centerDenominator':str(den),'fullCenterNumerators':center,
        'inverseDenominator':str(ad0*L),'inverseNumerators':inversenums,'radius':'1/'+str(10**50),
        'fixedCoordinates':'Every coordinate outside unknownCoordinateIndices is fixed exactly to its rational center.',
        'equationConvention':'Three Hermitian upper-triangle systems for Gram-I, Jz/j and (Jz^2-cI)/j^2, followed by all ordered real and imaginary J+^2/j^2 entries, using the common integer systemScale.'}
    bounds,equations=exact_bounds(certificate);certificate['bounds']={key:str(val)for key,val in bounds.items()}
    dest=HERE/f'n{n}-k{k}';dest.mkdir(exist_ok=True)
    (dest/'certificate.json').write_text(json.dumps(certificate,indent=2)+'\n')
    (dest/'system.json').write_text(json.dumps({'n':n,'dimension':k,'equations':equations},indent=2)+'\n')
    coord=np.array(list(map(int,center)),float)/den;Z=coord[:d*k].reshape(d,k)+1j*coord[d*k:].reshape(d,k);approx=np.zeros((n+1,k),complex);approx[indices]=weights[:,None]*Z
    known=set(unknown);position={slot:i+1 for i,slot in enumerate(unknown)}
    def atom(slot):return 'x'+str(position[slot])if slot in known else '('+str(F(int(center[slot]),den))+')'
    records=[]
    for col in range(k):
        exact=['0']*(n+1)
        for row,r in enumerate(indices):
            a=row*k+col;b=d*k+a
            exact[r]=f'sqrt({comb(n,r)})/({2**exponents[row]})*({atom(a)}+I*{atom(b)})'
        records.append({'exactCoefficients':exact,'coefficients':[[float(z.real),float(z.imag)]for z in approx[:,col]],'support':indices})
    intervals={'x'+str(i+1):[str(F(int(center[slot]),den)-F(certificate['radius'])),str(F(int(center[slot]),den)+F(certificate['radius']))]for i,slot in enumerate(unknown)}
    basis={'id':f'algebraic-parity-n{n}-k{k}-t2','n':n,'t':2,'dimension':k,'basis':records,
        'label':f'Exact algebraic parity construction at N = {n}',
        'exactnessKind':'isolated-algebraic-root',
        'provenance':'Independent nested construction in a variable first-order AC parent within one Dicke-parity sector. Exact rational contraction certifies the full second-order subspace.',
        'derivationLatex':r'c_{r,a}=\frac{\sqrt{\binom Nr}}{2^{e_r}}(x_{r,a}+i y_{r,a}),\qquad r\equiv0\pmod2,\qquad c_{r,a}=0\text{ otherwise}.',
        'algebraicDefinition':{'type':'Unique zero of a rational quadratic system in a certified rational box, after fixing the other coordinates.','certificate':'certificate.json','system':'system.json','variables':['x'+str(i+1)for i in range(m)],'parameterIntervals':intervals,'unknownCoordinateIndices':unknown,'fullCenterNumerators':center,'centerDenominator':str(den),'radius':certificate['radius'],'scaleExponents':exponents,'parityIndices':indices},
        'scope':'Exact algebraic existence; short radical formulas and unrestricted optimality are not asserted.'}
    (dest/'basis.json').write_text(json.dumps(basis,indent=2)+'\n')
    report=verify(dest/'certificate.json');(dest/'verification.json').write_text(json.dumps(report,indent=2)+'\n')
    print(n,k,'proved',m,'unknowns','residual',np.linalg.norm(fun_all(v)),{key:float(val)for key,val in bounds.items()},flush=True)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--n',type=int,required=True);p.add_argument('--k',type=int,required=True);p.add_argument('--parity',type=int,default=0);a=p.parse_args();build(a.n,a.k,a.parity)
