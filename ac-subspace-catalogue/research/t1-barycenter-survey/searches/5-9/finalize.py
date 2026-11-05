from search import *
records=[]
for n in range(5,10):
 if n==6:r=json.loads((ROOT/'n6-refined.json').read_text())
 elif n==8:
  P=Problem(n);z=((23+5*mp.sqrt(22))/mp.sqrt(21))**(mp.mpf(1)/5);r=P.mpstate(z,3);r['exactRepeatedRoot']='((23+5*sqrt(22))/sqrt(21))**(1/5)';r['exactWronskianPolynomial']='21*z**10-46*sqrt(21)*z**5-21'
 else:r=max(json.loads((ROOT/f'n{n}-wronskian.json').read_text())['candidates'],key=lambda r:mp.mpf(r['distance']))
 if n==5:r['exactRepeatedRoot']='-(3*sqrt(3)+2*sqrt(7))**(1/4)';r['exactWronskianPolynomial']='z**8-6*sqrt(3)*z**4-1'
 summary=json.loads((ROOT/f'n{n}-summary.json').read_text());r['searchStarts']=summary['searchStarts'];r['searchMethods']={'BFGS':240,'DEProjectiveCharts':r['dimension'],'DERepeatedRoot':6};r['numericalSearchMaxDistance']=summary['bestSearch']['distance'];r['precisionDigits']=90
 c=[mp.mpc(*v)for v in r['dickeCoefficients']];a=[mp.mpc(*v)for v in r['basisCoordinates']];P=Problem(n);member=max(abs(sum(P.MV[i,j]*a[j]for j in range(P.k))-c[i])for i in range(n+1));gram=abs(sum(abs(v)**2 for v in c)-1)
 A=mp.matrix([[c[i+j]*mp.sqrt(mp.binomial(n-1,j)/mp.binomial(n,i+j))for j in range(n)]for i in range(2)]);rho=A*A.H;ac=max(abs(rho[i,j]-mp.mpf(int(i==j))/2)for i in range(2)for j in range(2))
 pol=[mp.mpc(1)];stars=[];nf=0
 for rr in r['roots']:
  if rr['root']is None:stars +=[[mp.mpf(0),mp.mpf(0),mp.mpf(1)]]*rr['multiplicity'];continue
  z=mp.mpc(*rr['root']);stars +=[xyz(z)]*rr['multiplicity']
  for _ in range(rr['multiplicity']):
   temp=[mp.mpc(0)]*(len(pol)+1)
   for j,v in enumerate(pol):temp[j]-=z*v;temp[j+1]+=v
   pol=temp
 while len(pol)<n+1:pol.append(mp.mpc(0))
 expected=[(-1)**j*mp.sqrt(mp.binomial(n,j))*c[j]for j in range(n+1)];pivot=max(range(n+1),key=lambda j:abs(expected[j]));scale=expected[pivot]/pol[pivot];reconstruct=max(abs(pol[j]*scale-expected[j])for j in range(n+1));B=[sum(v[j]for v in stars)/n for j in range(3)];distance=mp.sqrt(sum(v*v for v in B));derr=abs(distance-mp.mpf(r['distance']))
 assert max(member,gram,ac,reconstruct,derr)<mp.mpf('1e-70'),(n,member,gram,ac,reconstruct,derr)
 r['verification']={'passed':True,'membershipMaxResidual':mp.nstr(member,15),'normalizationResidual':mp.nstr(gram,15),'oneQubitReductionMaxResidual':mp.nstr(ac,15),'starReconstructionMaxResidual':mp.nstr(reconstruct,15),'distanceResidual':mp.nstr(derr,15)}
 r['stars']=[{'xyz':[mp.nstr(x,80)for x in xyz(mp.mpc(*p['root']))]if p['root']is not None else ['0','0','1'],'multiplicity':p['multiplicity']}for p in r['roots']]
 r['sphereConvention']='(2Re(z),-2Im(z),|z|^2-1)/(1+|z|^2); infinity is north'
 save(f'n{n}-result.json',r);records.append(r)
save('results.json',{'records':records,'scope':'Numerical maxima within the selected specified maximal first-order AC subspaces. No global optimality proof.','totalSearchStarts':sum(r['searchStarts']for r in records)})
for r in records:print(r['n'],r['distance'],r['multiplicityPattern'],r['verification'])
