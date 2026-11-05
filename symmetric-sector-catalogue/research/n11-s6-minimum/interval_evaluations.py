"""Rigorous averaged SOS evaluations using directed-rounding FLINT/Arb."""
import os;os.environ['OPENBLAS_NUM_THREADS']='1'
from exact_field import *
from flint import arb,ctx
from fractions import Fraction as F
from math import comb,factorial
import ctypes as C,numpy as np,time
ROOT=Path(__file__).parent;ctx.prec=1024;PREC=256;start=time.time()
cert=json.load(open(ROOT/'sos-certificate.json'));definition=json.load(open(ROOT/'covariants.json'));raw=json.load(open(ROOT/'source-basis.json'))['blocks'];ker=json.load(open(ROOT/'source-kernels.json'))['blocks']
lo,hi=map(fmpq,DATA['rootInterval']);assert MOD(lo)*MOD(hi)<0
for _ in range(1150):
 mid=(lo+hi)/2
 if MOD(lo)*MOD(mid)<0:hi=mid
 else:lo=mid
ball=arb(lo).union(arb(hi))
def ev(z):
 out=arb(0)
 for c in reversed(K(z).p.coeffs()):out=out*ball+arb(c)
 return out
orig=arb(fmpq(DATA['rootInterval'][0])).union(arb(fmpq(DATA['rootInterval'][1])));d=arb(0)
for c in reversed(MOD.derivative().coeffs()):d=d*orig+arb(c)
assert not d.contains(0)and ev(U)>0 and ev(Q)>0
h=C.CDLL(str(ROOT/('arb_helper.dylib'if sys.platform=='darwin'else'arb_helper.so')));vp=C.c_void_p;lng=C.c_long;ip=C.POINTER(C.c_int64)
for name,args,result in [('newvec',[lng],vp),('delvec',[vp,lng],None),('newmat',[lng,lng],vp),('delmat',[vp],None),('matentry',[vp,lng,lng],vp),('vecentry',[vp,lng],vp),('setvalue',[vp,C.c_char_p,C.c_char_p,lng],C.c_int),('realvalue',[vp,lng],vp),('freestring',[vp],None),('eval_source',[vp,vp,ip,lng,lng,lng],None),('eval_block',[vp,C.POINTER(vp),ip,ip,vp,lng,vp,vp,vp,ip,lng,lng,lng,lng,ip,lng,vp,lng],None)]:
 fn=getattr(h,name);fn.argtypes=args;fn.restype=result

def put(ptr,z):
 if not isinstance(z,arb):z=arb(fmpq(str(z)))
 assert h.setvalue(ptr,z.str(110).encode(),b'0',PREC)==0

def get(ptr):
 p=h.realvalue(ptr,100);s=C.string_at(p).decode();h.freestring(p);return arb(s)
def ptr(a):return a.ctypes.data_as(ip)
sources=(vp*len(raw))();sourcepq=np.array([[int(b['kind'][0]),int(b['kind'][1])]for b in raw],dtype=np.int64);mixes={};arrays={}
for bn,ib in enumerate(definition['sourceBlocks']):
 B=raw[ib];d=B['dimension'];e=ker[ib]['dimension'];J=B['L2'];pol=(int(B['kind'][0])+1)*(int(B['kind'][1])+1);sources[ib]=h.newmat(e,(J+1)*pol);R=[[K(v)for v in col]for col in ker[ib]['kernelColumns']];sc=[F(v)for v in definition['sourceKernelScales'][bn]];W=[[F(v)for v in row]for row in definition['sourceMixing'][bn]];M=h.newmat(e,d)
 for a in range(e):
  for b in range(d):put(h.matentry(M,a,b),ev(sum((R[j][b]*K(W[j][a]/sc[j])for j in range(e)),K())))
 mixes[ib]=M;arrays[ib]=np.load(ROOT/'source-integer-evals'/f'{ib}.npy',mmap_mode='r')
print('INITIALIZED SOURCES',time.time()-start,flush=True)
def ff(a,b):return factorial(a)//factorial(a-b)if a>=b>=0 else 0
variables=[(ib,a,b)for ib,B in enumerate(definition['blocks'])for b in range(B['dimension'])for a in range(b+1)];pivs={v:i for i,v in enumerate(cert['pivotVariables'])};blocks=[];offset=0
for bi,B in enumerate(definition['blocks']):
 d=B['dimension'];L=B['L2'];p,q=map(int,B['kind']);pol=(p+1)*(q+1);ops=[];coefs=[]
 for ci,ch in enumerate(B['channels']):
  ib=ch['sourceBlock'];j=ch['sourceChannel'];J=raw[ib]['L2'];shift=(J+11-L)//2
  for a in range(J+1):
   for b in range(12):
    r=a+b-shift
    if not 0<=r<=L:continue
    v=sum((-1)**k*comb(shift,k)*ff(J-a,shift-k)*ff(a,k)*ff(11-b,k)*ff(b,shift-k)for k in range(shift+1))*comb(11,b)
    if v:ops.append([ci,ib,j,a,b,r]);coefs.append(F(ch['scale'])*v)
 ops=np.array(ops,dtype=np.int64);coef=h.newvec(len(coefs))
 for j,v in enumerate(coefs):put(h.vecentry(coef,j),v)
 mix=h.newmat(d,d)
 for a in range(d):
  for b in range(d):put(h.matentry(mix,a,b),F(B['mixing'][b][a]))
 M=[[F(1,6*comb(5,i2+j))if i+j2==i2+j else F(0)for i2 in range(p+1)for j2 in range(q+1)]for i in range(p+1)for j in range(q+1)];LL=[[F(a==b)for b in range(pol)]for a in range(pol)];diag=[]
 for a in range(pol):
  value=M[a][a]-sum(LL[a][j]**2*diag[j]for j in range(a));assert value>0;diag.append(value)
  for b in range(a+1,pol):LL[b][a]=(M[b][a]-sum(LL[a][j]*LL[b][j]*diag[j]for j in range(a)))/value
 assert all(sum(LL[a][j]*diag[j]*LL[b][j]for j in range(pol))==M[a][b]for a in range(pol)for b in range(pol))
 logical=h.newmat(pol,pol);weights=h.newvec((L+1)*pol)
 for a in range(pol):
  for b in range(pol):put(h.matentry(logical,a,b),LL[a][b])
 for r in range(L+1):
  for j in range(pol):put(h.vecentry(weights,r*pol+j),diag[j]/comb(L,r))
 pp=[]
 for b in range(d):
  for a in range(b+1):
   if offset in pivs:pp.append([pivs[offset],a,b])
   offset+=1
 pp=np.array(pp,dtype=np.int64).reshape(-1,3);blocks.append((B,ops,coef,mix,logical,weights,pp));print('INITIALIZED BLOCK',bi,flush=True)
rows=[];points=cert['evaluationPoints'];begin=int(os.environ.get('BEGIN_POINT','0'));stop=int(os.environ.get('STOP_POINT',len(points)));saved=ROOT/'interval-evaluations.json'
if saved.exists():rows=json.load(open(saved));begin=max(begin,len(rows))
for point in range(begin,stop):
 for ib in definition['sourceBlocks']:
  B=raw[ib];d=B['dimension'];pol=(int(B['kind'][0])+1)*(int(B['kind'][1])+1);digits=np.ascontiguousarray(arrays[ib][point]);h.eval_source(sources[ib],mixes[ib],ptr(digits),d,(B['L2']+1)*pol,PREC)
 answer=h.newvec(287);xx=np.asarray(points[point],dtype=np.int64)
 for B,ops,coef,mix,logical,weights,pp in blocks:
  h.eval_block(None,sources,ptr(sourcepq),ptr(ops),coef,len(ops),mix,logical,weights,ptr(xx),B['dimension'],B['L2'],int(B['kind'][0]),int(B['kind'][1]),ptr(pp),len(pp),answer,PREC)
 q0=sum(comb(11,i)*(int(a)**2+int(b)**2)for i,(a,b)in enumerate(points[point]));norm=arb(q0**5);row=[(get(h.vecentry(answer,j))/norm).str(90)for j in range(287)];h.delvec(answer,287);rows.append(row);saved.write_text(json.dumps(rows,separators=(',',':'))+'\n')
 print('POINT',point+1,'seconds',round(time.time()-start),flush=True)
