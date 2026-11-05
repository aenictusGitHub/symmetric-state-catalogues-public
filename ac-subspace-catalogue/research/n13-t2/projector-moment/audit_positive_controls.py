#!/usr/bin/env python3
"""Populate moments from actual known codes and independently check every map."""
from generate import *
def audit(n):
 fixture=next(v for v in json.loads((HERE/'positive-control-bases.json').read_text())if v['n']==n);k=fixture['dimension'];d=n+1
 if not(ROOT/f'p3-n{n}-k{k}-maps.npz').exists():build(n,k)
 words,indices,full,schur,vmeta=build(n,k,True)
 B=np.array([[complex(*z)for z in row]for row in fixture['dickeColumns']]).T
 P=B@B.conj().T;F=np.fliplr(np.diag((-1.)**np.arange(d)));Pt=F@P.T@F.T
 def twirl(mats):
  out=[]
  for blk in full:
   C=np.zeros((len(blk['aa']),len(blk['aa'])),complex)
   for M,U in blk['U'].items():
    ws=np.array(words[M]);A=np.ones((len(ws),len(ws)),complex)
    for leg in range(3):A*=mats[leg][ws[:,leg,None],ws[None,:,leg]]
    C+=U.T@A@U/(blk['J']+1)
   out.append(C)
  return out
 C=twirl([P,P,P]);x=np.zeros(len(vmeta))
 for i,(sb,a,b)in enumerate(vmeta):
  info=schur[sb];H=info['H'][0];local=H.T@C[info['ib']]@H;x[i]=local[a,b].real
 maps=np.load(ROOT/f'p3-n{n}-k{k}-maps.npz');meta=json.loads((ROOT/f'p3-n{n}-k{k}-metadata.json').read_text())
 extract=lambda blocks:np.array([blocks[ib][a,b].real for ib,a,b in meta['fullmeta']])
 ct=twirl([Pt,P,P]);ci=[twirl([np.eye(d)if leg==ii else P for leg in range(3)])for ii in range(3)]
 cit=[twirl([np.eye(d),P,P]),twirl([Pt,np.eye(d),P]),twirl([Pt,P,np.eye(d)])]
 report={'n':n,'dimension':k,'source':fixture['source'],'equationResidual':float(np.max(np.abs(maps['A']@x))),'trace':float(maps['trace']@x),'directMomentMapError':float(np.max(np.abs(maps['direct']@x-extract(C)))),'partialTimeReversalMapError':float(np.max(np.abs(maps['theta']@x-extract(ct)))),'identityTimesPairMomentMapError':max(float(np.max(np.abs(maps['iz'][ii]@x-extract(ci[ii]))))for ii in range(3)),'partialTimeReversalIdentityPairMomentError':max(float(np.max(np.abs(maps['iztheta'][ii]@x-extract(cit[ii]))))for ii in range(3))}
 assert report['equationResidual']<1e-11
 assert abs(report['trace']-k**3)<1e-11
 assert max(report[q]for q in report if 'Error'in q)<1e-11
 (ROOT/f'n{n}-positive-control.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report));return report
if __name__=='__main__':
 for n in[7,14]:audit(n)
