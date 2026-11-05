from survey import *
def check(n,root):
 s,V,M,P=get(n);k=s['dimension'];m=k-1
 D=mp.matrix([[sum(P[r][j]*mp.factorial(r)/mp.factorial(r-l)*root**(r-l)for r in range(l,n+1))for j in range(k)]for l in range(m)])
 Dnp=np.array([[complex(D[a,j])for j in range(k)]for a in range(m)]);Dnp/=np.linalg.norm(Dnp,axis=1)[:,None];_,_,vh=np.linalg.svd(Dnp);pivot=int(np.argmax(abs(vh[-1])));keep=[j for j in range(k)if j!=pivot]
 sol=mp.lu_solve(mp.matrix([[D[a,j]for j in keep]for a in range(m)]),mp.matrix([-D[a,pivot]for a in range(m)]));alpha=[mp.mpc(1)if j==pivot else sol[keep.index(j)]for j in range(k)];alpha=[x/mp.sqrt(sum(abs(y)**2 for y in alpha))for x in alpha]
 p=[sum(P[r][j]*alpha[j]for j in range(k))for r in range(n+1)];hi=max(r for r in range(n+1)if any(P[r][j]!=0 for j in range(k)));q=p[:hi+1]
 for _ in range(m):
  b=[q[-1]]
  for a in reversed(q[:-1]):b.append(a+root*b[-1])
  assert abs(b[-1])<mp.mpf('1e-80');q=list(reversed(b[:-1]))
 roots=[root]*m+mp.polyroots(q[::-1],maxsteps=4000,extraprec=100);return points_bary(roots,n,n-hi)[2]
if __name__=='__main__':
 out=[]
 for n in [10,11,13,14]:
  file=ROOT/f'n{n}-lowerfold.json'
  if not file.exists():continue
  data=json.loads(file.read_text());best=json.loads((ROOT/f'n{n}-result.json').read_text());ref=mp.mpf(best['distance'])
  for row in data['attempts']:
   d=check(n,mp.mpc(*map(str,row['root'])));out.append(dict(n=n,root=row['root'],distance=mp.nstr(d,70),differenceFromWronskian=mp.nstr(d-ref,20)))
  print(n,max(mp.mpf(r['differenceFromWronskian'])for r in out if r['n']==n),flush=True)
 (ROOT/'lowerfold-refined-check.json').write_text(json.dumps(out,indent=2))
