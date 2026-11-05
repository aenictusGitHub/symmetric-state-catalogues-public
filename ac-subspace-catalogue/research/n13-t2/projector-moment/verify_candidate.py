#!/usr/bin/env python3
"""Check the saved numerical relaxation without an optimizer. Not an exact proof."""
from generate import *
def verify():
 cand=json.loads((HERE/'candidate.json').read_text());n=cand['n'];k=cand['dimension']
 path=ROOT/f'p3-n{n}-k{k}-maps.npz'
 if not path.exists():build(n,k)
 words,indices,full,schur,vmeta=build(n,k,True)
 meta=json.loads((ROOT/f'p3-n{n}-k{k}-metadata.json').read_text());M=np.load(path)
 Cs=[np.array(a['matrix'])for a in cand['blocks']]
 for B,C in zip(full,cand['blocks']):assert B['J']==C['J2']and B['aa']==C['pairSpins']
 x=np.zeros(len(vmeta));input_mats=[]
 for b in schur:
  H=b['H'][0];input_mats.append(H.T@Cs[b['ib']]@H)
 for i,(sb,a,b)in enumerate(vmeta):x[i]=input_mats[sb][a,b]
 y=M['direct']@x;theta=M['theta']@x;z=M['iz']@x;zt=M['iztheta']@x
 ident=np.array([float(a==b)for ib,a,b in meta['fullmeta']]);one=1/(k*k*(n+1));zero=1/k**3
 fm={tuple(v):i for i,v in enumerate(meta['fullmeta'])};tests=[]
 for ib,C in enumerate(input_mats):tests.append(('input-Schur',ib,float(np.linalg.eigvalsh(C)[0])))
 locals=[('PPT-PPP',theta)]
 for i in range(3):locals.extend([(f'EPP-{i}',z[i]-y),(f'PPT-EPP-{i}',zt[i]-theta)])
 for i,j in[(0,1),(0,2),(1,2)]:locals.extend([(f'EEP-{i}{j}',y-z[i]-z[j]+one*ident),(f'PPT-EEP-{i}{j}',theta-zt[i]-zt[j]+one*ident)])
 locals.extend([('EEE',z.sum(axis=0)-y+(zero-3*one)*ident),('PPT-EEE',zt.sum(axis=0)-theta+(zero-3*one)*ident)])
 for name,v in locals:
  for ib,b in enumerate(meta['full']):
   d=len(b['aa']);C=np.array([[v[fm[ib,min(a,c),max(a,c)]]for c in range(d)]for a in range(d)])
   tests.append((name,ib,float(np.linalg.eigvalsh(C)[0])))
 report={'numericalChecksPassed':True,'n':n,'dimension':k,'realMomentVariables':len(x),'PSDblocksChecked':len(tests),'trace':float(M['trace']@x),'maximumEqualityResidual':float(np.max(np.abs(M['A']@x))),'minimumEigenvalue':min(v[2]for v in tests),'interpretation':'Numerically consistent relaxation only. Eigenvalues within1e-10 ofzero are treated as roundoff for this check; this is not an exact certificate and is not a physical code.'}
 assert abs(report['trace']-1)<1e-9
 assert report['maximumEqualityResidual']<1e-10
 assert report['minimumEigenvalue']>-1e-10
 (ROOT/'candidate-verification.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report));return report
if __name__=='__main__':verify()
