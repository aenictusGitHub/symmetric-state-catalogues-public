"""Exact Gaussian-integer polynomial evaluations; all int64 operations have checked bounds."""
import os;os.environ['OPENBLAS_NUM_THREADS']='1'
from pathlib import Path
import json,numpy as np
from scipy import sparse
from math import comb
ROOT=Path(__file__).parent

def polarized(mons,points):
 x=np.asarray(points,dtype=np.int64);n=len(x);maximum=int(np.max(abs(x)));d=len(mons[0]);p=sum(i<12 for i in mons[0]);q=d-p;assert (1<<d)*(2*maximum)**d<2**60
 y=x[:,::-1,:].copy();y[:,:,1]*=-1;y*=((-1)**np.arange(12))[None,:,None];xb=x.copy();xb[:,:,1]*=-1;yb=y.copy();yb[:,:,1]*=-1;a=np.concatenate([x,xb],axis=1);b=np.concatenate([y,yb],axis=1);out=np.zeros(((p+1)*(q+1),2,len(mons),n),dtype=np.int64)
 for mask in range(1<<d):
  re=np.ones((len(mons),n),dtype=np.int64);im=np.zeros_like(re)
  for j in range(d):
   values=(b if mask>>j&1 else a)[:,[m[j]for m in mons],:];u=values[:,:,0].T;v=values[:,:,1].T;re,im=re*u-im*v,re*v+im*u
  i=(mask&((1<<p)-1)).bit_count();j=(mask>>p).bit_count();out[i*(q+1)+j,0]+=re;out[i*(q+1)+j,1]+=im
 return out

def evaluate(raw,ib,points,folder):
 B=raw[ib];path=folder/f'{ib}.npy'
 if path.exists():return
 mons=sorted(set(tuple(m)for c in B['covariants']for p in c for m,v in p));mi={m:i for i,m in enumerate(mons)};mv=polarized(mons,points);pols=mv.shape[0];ns=len(points);dim=B['dimension'];L=B['L2'];row=[];col=[];coef=[]
 for a,cov in enumerate(B['covariants']):
  for r,p in enumerate(cov):
   for m,c in p:row.append(a*(L+1)+r);col.append(mi[tuple(m)]);coef.append(int(c))
 count=max(abs(v).bit_length()for v in coef);bits=8;chunks=(count+bits-1)//bits;mm=mv.transpose(2,0,1,3).reshape(len(mons),-1);maximum=int(np.max(abs(mm)));vals=[]
 for k in range(chunks):
  cc=[(1 if v>=0 else -1)*((abs(v)>>(bits*k))&255)for v in coef];bound=[0]*(dim*(L+1))
  for r,v in zip(row,cc):bound[r]+=abs(v)
  assert max(bound)*maximum<2**62
  M=sparse.csr_matrix((np.array(cc,dtype=np.int64),(row,col)),shape=(dim*(L+1),len(mons)),dtype=np.int64);vv=M@mm;assert vv.dtype==np.int64;vals.append(vv)
 carry=np.zeros_like(vals[0]);lo=np.zeros_like(carry);mid=np.zeros_like(carry)
 for k in range(8):
  if k<chunks:carry+=vals[k]
  assert int(np.max(abs(carry)))<2**63-1
  digit=carry&255
  if k<4:lo|=digit<<(8*k)
  else:mid|=digit<<(8*(k-4))
  carry>>=8
 assert chunks<=8 and int(np.max(abs(carry)))<2**62
 # [point, covariant, spin component, logical component, real/imag, limb]
 out=np.stack([lo,mid,carry],axis=-1).reshape(dim,L+1,pols,2,ns,3).transpose(4,0,1,2,3,5).copy();np.save(path,out);print('EXACT',folder.name,ib,out.shape,flush=True)

if __name__=='__main__':
 import sys
 c=json.load(open(ROOT/'sos-certificate.json'));points=c['evaluationPoints'];which=sys.argv[1]if len(sys.argv)>1 else'source';raw=json.load(open(ROOT/('source-basis.json'if which=='source'else'holomorphic-basis.json')))['blocks'];folder=ROOT/(which+'-integer-evals');folder.mkdir(exist_ok=True)
 ids=json.load(open(ROOT/'covariants.json'))['sourceBlocks']if which=='source'else range(len(raw))
 for ib in ids:evaluate(raw,ib,points,folder)
