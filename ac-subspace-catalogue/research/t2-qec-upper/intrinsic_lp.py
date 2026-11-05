"""AC2 specialization of intrinsic SU(2) MacWilliams enumerator constraints.
Exact rational Racah transform, numerical feasibility only unless separately certified.
"""
from fractions import Fraction as F
from pathlib import Path
import numpy as np,json,argparse,time
from scipy.optimize import linprog
HERE=Path(__file__).resolve().parent

def racah(n,a,b):
    term=F(1);total=term
    for h in range(min(a,b)):
        term*=F((-a+h)*(a+1+h)*(-b+h)*(b+1+h), (h+1)**2*(-n+h)*(n+2+h))
        total+=term
    return F(2*b+1,n+1)*total

def matrix(n):
    # Detailed balance avoids about half the summations.
    out=np.zeros((n+1,n+1))
    for a in range(n+1):
        for b in range(a,n+1):
            f=racah(n,a,b);out[b,a]=float(f);out[a,b]=float(f*F(2*a+1,2*b+1))
    return out

def check(n,k,mat,shadow=True,nested=False):
    d=n+1;I=np.eye(d);Aeq=[np.ones(d),I[0],I[1],I[2],mat[1],mat[2]];beq=[1,k/d,0,0,0,0]
    Au=[I-k*mat,-mat];bu=[np.zeros(d),np.zeros(d)]
    if n%2 and k%2:
        Au.append(np.array([[(-1)**l for l in range(d)]]));bu.append(np.array([(k-1)/k]))
    if shadow:
        sh=mat*np.array([(-1)**l for l in range(d)])
        Au.extend([-sh,sh]);bu.extend([np.zeros(d),np.array([min(2*l+1,k*(k+(-1)**(n-l))/2)for l in range(d)])/k])
    if nested:
        c=n*(n+2)/12;r=d-4*k
        if r<0:return {'n':n,'k':k,'status':2,'reason':'Negative error complement rank'}
        fac=[np.ones(d),np.array([3-l*(l+1)/(2*c)for l in range(d)]),np.array([-4+l*(l+1)/(2*c)for l in range(d)])]
        fac[1][0]=3;fac[2][0]=r/k;dims=[k,3*k,r]
        for a in range(3):
            for b in range(a,3):
                f=fac[a]*fac[b];B=mat*f;S=B*np.array([(-1)**l for l in range(d)])
                Au.extend([-B,-S,S]);bu.extend([np.zeros(d),np.zeros(d),np.array([2*l+1 for l in range(d)])/k])
                if a==b:Au.append(np.diag(f)-dims[a]*B);bu.append(np.zeros(d))
    Au=np.vstack(Au);bu=np.concatenate(bu);Ae=np.array(Aeq);be=np.array(beq)
    result=linprog(np.zeros(d),A_ub=Au,b_ub=bu,A_eq=Ae,b_eq=be,bounds=(0,None),method='highs')
    rec={'n':n,'k':k,'status':result.status,'success':result.success,'message':result.message,'shadow':shadow,'nested':nested}
    if result.success:rec.update(point=result.x.tolist(),equalityResidual=float(np.max(abs(Ae@result.x-be))),minimumSlack=float(np.min(bu-Au@result.x)))
    return rec

def main():
    p=argparse.ArgumentParser();p.add_argument('n',nargs='+',type=int);args=p.parse_args();records=[]
    for n in args.n:
        start=time.monotonic();M=matrix(n);k=min((n-2)//4,(n+5)//5)
        for nested in [False,True]:
            r=check(n,k,M,nested=nested);records.append(r);print(n,k,nested,r['status'],round(time.monotonic()-start,2),flush=True)
    target=HERE/'intrinsic-lp-scan.json'
    previous=json.loads(target.read_text())['records'] if target.exists() else []
    merged={(r['n'],r['k'],r['nested']):r for r in previous+records}
    target.write_text(json.dumps({'scope':'Exact rational transform; floating-point LP feasibility, not existence of physical codes or a nonexistence certificate.','records':[merged[key]for key in sorted(merged)]},indent=2)+'\n')
if __name__=='__main__':main()
