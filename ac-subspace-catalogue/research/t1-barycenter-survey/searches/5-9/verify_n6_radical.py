from search import *
r=json.loads((ROOT/'n6-radical-result.json').read_text());n=6;P=Problem(n);c=[mp.mpc(*v)for v in r['dickeCoefficients']];a=[mp.mpc(*v)for v in r['basisCoordinates']]
member=max(abs(sum(P.MV[i,j]*a[j]for j in range(P.k))-c[i])for i in range(n+1));gram=abs(sum(abs(v)**2 for v in c)-1)
A=mp.matrix([[c[i+j]*mp.sqrt(mp.binomial(n-1,j)/mp.binomial(n,i+j))for j in range(n)]for i in range(2)]);rho=A*A.H;ac=max(abs(rho[i,j]-mp.mpf(int(i==j))/2)for i in range(2)for j in range(2))
pol=[mp.mpc(1)];stars=[]
for rr in r['roots']:
 z=mp.mpc(*rr['root']);stars +=[xyz(z)]*rr['multiplicity']
 for _ in range(rr['multiplicity']):
  temp=[mp.mpc(0)]*(len(pol)+1)
  for j,v in enumerate(pol):temp[j]-=z*v;temp[j+1]+=v
  pol=temp
expected=[(-1)**j*mp.sqrt(mp.binomial(n,j))*c[j]for j in range(n+1)];pivot=max(range(n+1),key=lambda j:abs(expected[j]));scale=expected[pivot]/pol[pivot];reconstruct=max(abs(pol[j]*scale-expected[j])for j in range(n+1));B=[sum(v[j]for v in stars)/n for j in range(3)];distance=mp.sqrt(sum(v*v for v in B));derr=abs(distance-mp.mpf(r['distance']))
assert max(member,gram,ac,reconstruct,derr)<mp.mpf('1e-70')
r['verification']={'passed':True,'membershipMaxResidual':mp.nstr(member,15),'normalizationResidual':mp.nstr(gram,15),'oneQubitReductionMaxResidual':mp.nstr(ac,15),'starReconstructionMaxResidual':mp.nstr(reconstruct,15),'distanceResidual':mp.nstr(derr,15)}
r['stars']=[{'xyz':[mp.nstr(x,80)for x in xyz(mp.mpc(*p['root']))],'multiplicity':p['multiplicity']}for p in r['roots']];r['sphereConvention']='(2Re(z),-2Im(z),|z|^2-1)/(1+|z|^2); infinity is north';r['precisionDigits']=90
summary=json.loads((ROOT/'n6-summary.json').read_text());r['searchStarts']=summary['searchStarts'];r['searchMethods']={'BFGS':240,'DEProjectiveCharts':3,'DERepeatedRoot':6};r['numericalSearchMaxDistance']=summary['bestSearch']['distance']
save('n6-radical-result.json',r);print(json.dumps(r['verification'],indent=2))
