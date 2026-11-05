#!/usr/bin/env python3
import os
os.environ.setdefault('OPENBLAS_NUM_THREADS','1')
os.environ.setdefault('OMP_NUM_THREADS','1')
import cvxpy as cp,numpy as np,json,sys,time
from pathlib import Path
from scipy.linalg import svd
HERE=Path(__file__).resolve().parent
ROOT=Path(os.environ.get('AC_P3_OUTPUT',str(HERE/'generated'))).resolve()
ROOT.mkdir(parents=True,exist_ok=True);n=int(sys.argv[1])if len(sys.argv)>1 else 13;k=int(sys.argv[2])if len(sys.argv)>2 else 3
M=np.load(ROOT/f'p3-n{n}-k{k}-maps.npz');meta=json.loads((ROOT/f'p3-n{n}-k{k}-metadata.json').read_text());A=M['A'];U,S,Vh=svd(A,full_matrices=False);rank=sum(S>S[0]*1e-10);K=Vh[rank:].T
print('rank',rank,'kernel',K.shape,'singularboundary',S[max(0,rank-3):rank+3],flush=True)
y=cp.Variable(K.shape[1]);x=K@y;delta=cp.Variable();cons=[M['trace']@x==1];blocks=[];labels=[]
vm={tuple(v):i for i,v in enumerate(meta['variables'])};fm={tuple(v):i for i,v in enumerate(meta['fullmeta'])}
for ib,b in enumerate(meta['schur']):
 d=b['dim'];C=cp.bmat([[x[vm[ib,min(a,c),max(a,c)]]for c in range(d)]for a in range(d)])
 blocks.append(C);labels.append(('schur',ib));cons.append(C>>0)
Y=M['direct'];T=M['theta'];Z=M['iz'];ZT=M['iztheta'];ident=np.array([float(a==b)for ib,a,b in meta['fullmeta']]);one=1/(k*k*(n+1));zero=1/(k**3)
localizers=[('PPT-PPP',T,0.)]
for ii in range(3):
 localizers.append((f'EPP-{ii}',Z[ii]-Y,0.))
 localizers.append((f'PPT-EPP-{ii}',ZT[ii]-T,0.))
for ii,jj in [(0,1),(0,2),(1,2)]:
 localizers.append((f'EEP-{ii}{jj}',Y-Z[ii]-Z[jj],one))
 localizers.append((f'PPT-EEP-{ii}{jj}',T-ZT[ii]-ZT[jj],one))
localizers.extend([('EEE',Z.sum(axis=0)-Y,zero-3*one),('PPT-EEE',ZT.sum(axis=0)-T,zero-3*one)])
for name,L,const in localizers:
 v=L@x+const*ident
 for ib,b in enumerate(meta['full']):
  d=len(b['aa']);C=cp.bmat([[v[fm[ib,min(a,c),max(a,c)]]for c in range(d)]for a in range(d)])
  blocks.append(C);labels.append((name,ib));cons.append(C>>0)
prob=cp.Problem(cp.Minimize(cp.sum_squares(y)*1e-5),cons)
t0=time.monotonic()
try:prob.solve(solver='CLARABEL',tol_gap_abs=1e-9,tol_feas=1e-9,tol_gap_rel=1e-9,max_iter=200,time_limit=90.,verbose=False)
except Exception as e:print('solver-error',repr(e),flush=True)
report={'n':n,'k':k,'variables':len(meta['variables']),'rank':int(rank),'kernelDimension':K.shape[1],'PSDblocks':len(blocks),'status':prob.status,'elapsedSeconds':time.monotonic()-t0}
if y.value is not None:
 xx=K@y.value;eigs=[np.linalg.eigvalsh(C.value)[0]for C in blocks]
 report.update(minEigenvalue=float(min(eigs)),minimumBlock=labels[int(np.argmin(eigs))],equalityResidual=float(np.max(np.abs(A@xx))),trace=float(M['trace']@xx))
 np.savez(ROOT/f'p3-n{n}-k{k}-solution.npz',x=xx,minimumEigenvalues=eigs,K=K)
print(json.dumps(report),flush=True);(ROOT/f'p3-n{n}-k{k}-report.json').write_text(json.dumps(report,indent=2)+'\n')
