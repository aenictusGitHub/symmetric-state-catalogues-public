from pathlib import Path
import re,ast,cmath,json,unicodedata,math
import numpy as np
ROOT=Path(__file__).resolve().parent
src=ROOT.parents[1]/'references/ACLiege.jl';raw=src.read_text();clean='\n'.join(l.split('#')[0] for l in raw.splitlines())
def parse(expression,env):
 expression=expression.replace('−','-').replace('^','**');expression=re.sub(r'(?<=[0-9)])im\b','*im',expression)
 def ev(a):
  if isinstance(a,ast.Constant) and isinstance(a.value,(int,float)):return a.value
  if isinstance(a,ast.Name):return env[a.id]
  if isinstance(a,ast.List):return [ev(x) for x in a.elts]
  if isinstance(a,ast.UnaryOp):return -ev(a.operand) if isinstance(a.op,ast.USub) else ev(a.operand)
  if isinstance(a,ast.BinOp):
   x,y=ev(a.left),ev(a.right)
   if isinstance(a.op,ast.Add):return x+y
   if isinstance(a.op,ast.Sub):return x-y
   if isinstance(a.op,ast.Mult):return x*y
   if isinstance(a.op,ast.Div):return x/y
   if isinstance(a.op,ast.Pow):return x**y
  if isinstance(a,ast.Call) and isinstance(a.func,ast.Name) and len(a.args)==1:
   f={'sqrt':cmath.sqrt,'exp':cmath.exp,'atan':cmath.atan,'acos':cmath.acos,'normalize':lambda v:np.array(v)/np.linalg.norm(v)}[a.func.id]
   return f(ev(a.args[0]))
  raise ValueError(ast.dump(a))
 return ev(ast.parse(expression.strip(),mode='eval').body)
def red(a,b,t):
 n=len(a)-1;v=np.zeros((t+1,t+1),complex)
 for i in range(t+1):
  for j in range(t+1):
   v[i,j]=sum(a[i+l]*b[j+l].conjugate()*math.comb(n-t,l)*math.sqrt(math.comb(t,i)*math.comb(t,j)/(math.comb(n,i+l)*math.comb(n,j+l))) for l in range(n-t+1))
 return v
records=[]
for match in re.finditer(r'^function (Π_(\d+)_(\d+)_(\d+)\w*)\(.*?^end',clean,re.M|re.S):
 name,n,k,t=match.group(1),int(match.group(2)),int(match.group(3)),int(match.group(4));body=match.group(0)
 if not re.search(r'\d\.\d{5}',body):continue
 env={'im':1j};before,ret=body.split('return',1)
 for line in before.splitlines()[1:]:
  if '=' in line:
   key,exp=line.split('=',1);env[unicodedata.normalize('NFKC',key.strip())]=parse(exp,env)
 a=np.array(parse(ret.rsplit('end',1)[0],env),complex)
 if a.shape[1]!=n+1:
  records.append({'id':'acliege-'+name,'n':n,'t':t,'reportedDimension':k,'sourceFunction':name,'sourceLine':clean[:match.start()].count('\n')+1,'accepted':False,'shapeError':{'expected':n+1,'actual':a.shape[1]},'sourceBasis':[[[float(z.real),float(z.imag)] for z in v] for v in a]});print('INVALID',name,a.shape,flush=True);continue
 rawgram=np.linalg.norm(a.conj()@a.T-np.eye(len(a)));sing=np.linalg.svd(a,compute_uv=False);rank=int(sum(sing>1e-8))
 # Use SVD to obtain an orthonormal basis for the actual represented span.
 _,_,q=np.linalg.svd(a,full_matrices=False);basis=q[:rank]
 err=max(np.linalg.norm(red(x,y,t)-(np.eye(t+1)/(t+1) if i==j else 0)) for i,x in enumerate(basis) for j,y in enumerate(basis))
 rawerr=max(np.linalg.norm(red(x,y,t)-np.eye(t+1)*np.vdot(y,x)/(t+1)) for x in a for y in a)
 rec={'id':'acliege-'+name,'n':n,'t':t,'dimension':rank,'reportedDimension':k,'sourceFunction':name,'sourceLine':clean[:match.start()].count('\n')+1,'sourceVectorCount':len(a),'sourceGramError':rawgram,'crossReductionError':err,'sourceCrossReductionError':rawerr,'accepted':bool(rank==k and err<1e-7),'basis':[[[float(z.real),float(z.imag)] for z in v] for v in basis],'sourceBasis':[[[float(z.real),float(z.imag)] for z in v] for v in a],'singularValues':sing.tolist()}
 records.append(rec);print(name,'vectors',len(a),'rank',rank,'gram',rawgram,'cross',err,flush=True)
(ROOT/'numerical-subspaces.json').write_text(json.dumps(records,indent=2)+'\n')
