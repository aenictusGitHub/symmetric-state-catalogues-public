from search import *
z=sy.symbols('z');out=[]
for n in [5,7,8,9]:
 r=ROWS[n];k=len(r['basis']);support=[]
 for b in r['basis']:support.append({i:sy.simplify((-1)**i*sy.sqrt(sy.binomial(n,i))*sy.sympify(c))for i,c in enumerate(b['exactCoefficients'])if c!='0'})
 co={}
 for ix in itertools.product(*[list(x)for x in support]):
  if len(set(ix))!=k:continue
  v=sy.prod(support[j][ix[j]]for j in range(k))*sy.prod(ix[b]-ix[a]for a in range(k)for b in range(a+1,k));e=sum(ix)-k*(k-1)//2;co[e]=co.get(e,0)+v
 co={e:sy.simplify(sy.expand_complex(v))for e,v in co.items()};co={e:v for e,v in co.items()if v!=0};low,hi=min(co),max(co);lead=co[hi];co={e-low:sy.simplify(sy.expand_complex(v/lead))for e,v in co.items()};f=sy.Poly.from_dict({(e,):v for e,v in co.items()},z);W=sy.Poly(f.as_expr(),z,extension={5:[sy.sqrt(3)],7:[sy.sqrt(41)*sy.I],8:[sy.sqrt(21)],9:[sy.sqrt(35),sy.I]}[n]);g=W.gcd(W.diff()).degree()
 rec=dict(n=n,zeroFactorDegree=low,degree=f.degree(),coefficients={str(e):str(v)for e,v in co.items()},polynomial=str(f.as_expr()),squareFree=(g==0),gcdDegree=g);out.append(rec);print(n,rec,flush=True)
save('exact-wronskians.json',out)
