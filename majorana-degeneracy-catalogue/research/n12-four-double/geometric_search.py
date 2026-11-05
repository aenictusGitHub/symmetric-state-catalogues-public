#!/usr/bin/env python3
"""Separated, arbitrary-complex Majorana searches for N12 partition2^4 1^4.

Eight normalized physical spinors carry weights(2,2,2,2,1,1,1,1). The first
double star is fixed at north and the second on a meridian, removing rotations.
Pairwise chord-distance inequalities forbid the root-collision boundary.
Numerical failures do not prove nonexistence; numerical successes are not proofs.
"""
import argparse
import itertools
import json
import math
import time
from pathlib import Path

import numpy as np
from scipy.optimize import differential_evolution, minimize

HERE=Path(__file__).resolve().parent
WEIGHTS=np.array([2,2,2,2,1,1,1,1])
PAIRS=list(itertools.combinations(range(8),2))


def normalized(v):return np.asarray(v,dtype=float)/np.linalg.norm(v)


def gauge(points):
    p=np.asarray(points,dtype=float).copy();p/=np.linalg.norm(p,axis=1)[:,None]
    # Relabel double stars to avoid an antipodal meridian reference.
    j=max(range(1,4),key=lambda k:np.linalg.norm(np.cross(p[0],p[k])))
    p[[1,j]]=p[[j,1]]
    ez=p[0];ex=normalized(p[1]-np.dot(p[1],ez)*ez);ey=np.cross(ez,ex)
    p=p@np.stack([ex,ey,ez],axis=1)
    theta=np.arccos(np.clip(p[:,2],-1,1));phi=np.arctan2(p[:,1],p[:,0])
    return np.r_[theta[1],np.column_stack([theta[2:],phi[2:]]).ravel()]


class Model:
    def __init__(self,double_gap=0):
        self.double_gap=double_gap
        self.parameters=[(1,0)]+[(i,j) for i in range(2,8) for j in (0,1)]
        self.binomial=np.sqrt([math.comb(12,r) for r in range(13)])
        self.indices=np.array([[a+k for k in range(10)] for a in range(4)])
        self.rweights=np.array([[math.sqrt(math.comb(3,a)*math.comb(9,k)/math.comb(12,a+k)) for k in range(10)] for a in range(4)])
        self.last=None;self.has_jac=False

    def evaluate(self,x,jac=False):
        if self.last is not None and np.array_equal(self.last,x) and (not jac or self.has_jac):return
        angles=np.zeros((8,2))
        for value,(i,j) in zip(x,self.parameters):angles[i,j]=value
        th,ph=angles.T
        self.points=np.column_stack([np.sin(th)*np.cos(ph),np.sin(th)*np.sin(ph),np.cos(th)])
        spinors=np.column_stack([np.cos(th/2),np.exp(1j*ph)*np.sin(th/2)])
        factors=[np.convolve(p,p) if w==2 else p for p,w in zip(spinors,WEIGHTS)]
        prefix=[np.array([1+0j])]
        for f in factors:prefix.append(np.convolve(prefix[-1],f))
        suffix=[None]*9;suffix[8]=np.array([1+0j])
        for i in range(7,-1,-1):suffix[i]=np.convolve(factors[i],suffix[i+1])
        z=prefix[-1]/self.binomial;norm=np.linalg.norm(z);c=z/norm
        A=c[self.indices]*self.rweights;rho=A@A.conj().T
        F=[]
        for a in range(4):
            F.append(rho[a,a].real-.25)
            for b in range(a):F.extend([math.sqrt(2)*rho[a,b].real,math.sqrt(2)*rho[a,b].imag])
        self.F=np.array(F);self.coefficients=c
        self.constraints=np.array([np.sum((self.points[i]-self.points[j])**2) for i,j in PAIRS])
        if jac:
            dc=[];dp=np.zeros((8,3,13))
            for column,(i,kind) in enumerate(self.parameters):
                if kind==0:
                    derivative=np.array([-.5*np.sin(th[i]/2),.5*np.exp(1j*ph[i])*np.cos(th[i]/2)])
                    dp[i,:,column]=[np.cos(th[i])*np.cos(ph[i]),np.cos(th[i])*np.sin(ph[i]),-np.sin(th[i])]
                else:
                    derivative=np.array([0,1j*spinors[i,1]])
                    dp[i,:,column]=[-np.sin(th[i])*np.sin(ph[i]),np.sin(th[i])*np.cos(ph[i]),0]
                derivative=np.convolve(spinors[i],derivative)*2 if WEIGHTS[i]==2 else derivative
                dz=np.convolve(np.convolve(prefix[i],derivative),suffix[i+1])/self.binomial
                dc.append(dz/norm-c*np.real(np.vdot(c,dz))/norm)
            dc=np.array(dc).T;dA=dc[self.indices,:]*self.rweights[:,:,None]
            drho=np.einsum('akp,bk->abp',dA,A.conj())+np.einsum('ak,bkp->abp',A,dA.conj())
            J=[]
            for a in range(4):
                J.append(drho[a,a].real)
                for b in range(a):J.extend([math.sqrt(2)*drho[a,b].real,math.sqrt(2)*drho[a,b].imag])
            self.J=np.array(J)
            self.CJ=np.array([2*(self.points[i]-self.points[j])@(dp[i]-dp[j]) for i,j in PAIRS])
        self.last=x.copy();self.has_jac=jac

    def fun(self,x):self.evaluate(x);return float(self.F@self.F)
    def grad(self,x):self.evaluate(x,True);return 2*self.J.T@self.F
    def thresholds(self,gap):return np.array([max(gap,self.double_gap)**2 if i<4 and j<4 else gap*gap for i,j in PAIRS])
    def con(self,x,gap):self.evaluate(x);return self.constraints-self.thresholds(gap)
    def conjac(self,x,gap):self.evaluate(x,True);return self.CJ


def random_points(rng,gap):
    for attempt in range(10000):
        if gap>.75:
            cube=np.array(list(itertools.product([-1.,1.],repeat=3)))/math.sqrt(3)
            p=cube+.12*rng.normal(size=(8,3))
            # Random labels avoid assigning the same tetrahedron to double stars.
            p=p[rng.permutation(8)]
        else:p=rng.normal(size=(8,3))
        p/=np.linalg.norm(p,axis=1)[:,None]
        if min(np.linalg.norm(p[i]-p[j]) for i,j in PAIRS)>gap:return p
    raise RuntimeError('Could not initialize separated points')


def split_octahedron(rng,gap,adjacent):
    octa=np.array([[1,0,0],[-1,0,0],[0,1,0],[0,-1,0],[0,0,1],[0,0,-1]],dtype=float)
    omitted=[0,2] if adjacent else [0,1]
    points=[p for i,p in enumerate(octa) if i not in omitted]
    opening=math.asin(min(.8,.55*gap))
    for i in omitted:
        p=octa[i];v=rng.normal(size=3);v=normalized(v-np.dot(v,p)*p)
        points.extend([math.cos(opening)*p+math.sin(opening)*v,math.cos(opening)*p-math.sin(opening)*v])
    return np.array(points)


def summarize(model,x,kind,gap,seed,seconds,details):
    model.evaluate(x,True);c=model.coefficients
    time_reversed=np.array([(-1)**r*np.conj(c[12-r]) for r in range(13)])
    actual_gap=float(math.sqrt(max(0,min(model.constraints))))
    double_gap=min(np.linalg.norm(model.points[i]-model.points[j]) for i,j in PAIRS if i<4 and j<4)
    residual=float(np.linalg.norm(model.F))
    feasible=bool(actual_gap>=gap-1e-7 and double_gap>=model.double_gap-1e-7)
    return {'kind':kind,'seed':seed,'enforcedSeparation':gap,'actualMinimumSeparation':actual_gap,
            'enforcedDoublePointSeparation':model.double_gap,'actualDoublePointSeparation':float(double_gap),
            'HSResidual':residual,'feasible':feasible,
            'strongNumericalCandidate':feasible and residual<1e-11,
            'timeReversalInfidelity':float(max(0,1-abs(np.vdot(c,time_reversed))**2)),
            'variables':x.tolist(),'weightedPoints':[{'x':float(p[0]),'y':float(p[1]),'z':float(p[2]),'multiplicity':int(w)} for p,w in zip(model.points,WEIGHTS)],
            'normalizedDickeCoefficients':[[float(v.real),float(v.imag)] for v in c],
            'jacobianSingularValues':np.linalg.svd(model.J,compute_uv=False).tolist(),
            'seconds':round(seconds,3),**details}


def local_search(points,gap,kind,seed,maxiter=1200,double_gap=0):
    model=Model(double_gap);x=gauge(points);started=time.time()
    result=minimize(model.fun,x,jac=model.grad,method='SLSQP',
                    constraints=[{'type':'ineq','fun':lambda x:model.con(x,gap),'jac':lambda x:model.conjac(x,gap)}],
                    bounds=[(.005,math.pi-.005)]+[(None,None)]*12,
                    options={'maxiter':maxiter,'ftol':2e-15,'disp':False})
    return summarize(model,result.x,kind,gap,seed,time.time()-started,
                     {'iterations':result.nit,'optimizerSuccess':bool(result.success),'message':result.message})


def global_search(gap,seed,iterations,double_gap=0):
    model=Model(double_gap);started=time.time()
    def objective(x):
        value=model.fun(x);violation=np.maximum(0,model.thresholds(gap)-model.constraints)
        return value+10*np.dot(violation,violation)
    bounds=[(.05,math.pi-.05)]+[(.03,math.pi-.03),(-math.pi,math.pi)]*6
    result=differential_evolution(objective,bounds,seed=seed,popsize=12,maxiter=iterations,
                                  tol=1e-9,polish=False,updating='immediate',workers=1)
    model.evaluate(result.x);points=model.points.copy()
    answer=local_search(points,gap,'differential-evolution-plus-constrained-polish',seed,double_gap=double_gap)
    answer['globalIterations']=result.nit;answer['globalEvaluations']=result.nfev
    answer['seconds']=round(time.time()-started,3)
    return answer


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--output',type=Path,default=HERE/'geometric-search-results.json')
    p.add_argument('--random-starts',type=int,default=120)
    p.add_argument('--split-starts',type=int,default=40)
    p.add_argument('--global-runs',type=int,default=3)
    p.add_argument('--global-iterations',type=int,default=250)
    p.add_argument('--double-gap',type=float,default=0)
    a=p.parse_args();rng=np.random.default_rng(12024411)
    model=Model();x=gauge(random_points(rng,.3));direction=rng.normal(size=13);direction/=np.linalg.norm(direction)
    model.evaluate(x,True);analytic=model.J@direction;constraint_analytic=model.CJ@direction
    h=1e-5;model.evaluate(x+h*direction);plus=model.F.copy();cplus=model.constraints.copy()
    model.evaluate(x-h*direction);minus=model.F.copy();cminus=model.constraints.copy()
    derivative_errors={'AC':float(max(abs(analytic-(plus-minus)/(2*h)))),
                       'separation':float(max(abs(constraint_analytic-(cplus-cminus)/(2*h))))}
    assert max(derivative_errors.values())<1e-7
    output={'scientificStatus':'Arbitrary-complex numerical searches with enforced distinct-star separation; neither existence nor nonexistence is proved by an optimizer.',
            'n':12,'partition':[2,2,2,2,1,1,1,1],'targetOrder':3,'seed':12024411,'doublePointSeparationConstraint':a.double_gap,
            'parameterization':'Eight weighted normalized physical spinors, with three spatial rotation freedoms fixed. Thirteen real variables remain.',
            'derivativeChecks':derivative_errors,'runs':[]}
    def save(answer):
        output['runs'].append(answer)
        feasible=[r for r in output['runs'] if r['feasible']]
        output['summary']={'attempts':len(output['runs']),'feasible':len(feasible),
                           'strongNumericalCandidates':sum(r['strongNumericalCandidate'] for r in output['runs']),
                           'bestFeasibleHSResidual':min((r['HSResidual'] for r in feasible),default=None)}
        a.output.write_text(json.dumps(output,indent=2)+'\n')
        print(len(output['runs']),answer['kind'],'gap',answer['actualMinimumSeparation'],
              'HS',answer['HSResidual'],'TR',answer['timeReversalInfidelity'],flush=True)
    for i in range(a.random_starts):
        gap=[.15,.35,.6,.9][i%4]
        save(local_search(random_points(rng,gap),gap,'separated-random-multistart',12024411+i,double_gap=a.double_gap))
    for i in range(a.split_starts):
        gap=[.5,.2,.1,.05,.02][i%5]
        save(local_search(split_octahedron(rng,gap,bool((i//5)%2)),gap,
                          'split-octahedron-adjacent' if (i//5)%2 else 'split-octahedron-opposite',130000+i,double_gap=a.double_gap))
    for i in range(a.global_runs):save(global_search([.15,.35,.6][i%3],140000+i,a.global_iterations,a.double_gap))
    print('DONE',output['summary'],flush=True)


if __name__=='__main__':main()
