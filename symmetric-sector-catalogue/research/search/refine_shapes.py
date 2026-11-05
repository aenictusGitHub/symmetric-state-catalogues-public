import numpy as np,json,math
from scipy.optimize import least_squares,minimize
from optimize import ROOT,Sector
from verify_independent import pauli_sector

for n,k,support,signs in [(11,6,[0,3,6,9],[-1,1,1,1]),(12,8,[1,4,7,10],[1,1,1,1]),(12,6,[0,2,4,6,8,10,12],[1,1,-1,1,1,1,-1])]:
 obj=Sector(n,k)
 if n==12:
  raw=min(json.loads((ROOT/f'n{n}-s{k}-min-patterns.json').read_text()),key=lambda e:e['value'])
 else:raw=json.loads((ROOT/f'n{n}-s{k}-min-sparse.json').read_text())['best']
 c=np.sqrt(np.array(raw['weights']))
 if n==12 and k==6:
  basis=np.zeros((n+1,4));basis[0,0]=basis[2,1]=basis[8,2]=1/math.sqrt(2);basis[12,0]=basis[4,2]=-1/math.sqrt(2);basis[10,1]=1/math.sqrt(2);basis[6,3]=1
  x=np.array([c[0]*math.sqrt(2),c[2]*math.sqrt(2),c[4]*math.sqrt(2),c[6]])
 else:
  basis=np.zeros((n+1,4))
  for j,(r,sign) in enumerate(zip(support,signs)):basis[r,j]=sign
  x=c[support]
 x/=x[-1]
 def func(z):
  xx=np.r_[z,1.];cc=basis@xx;v,g=obj.value_grad(np.r_[cc,np.zeros(n+1)]);return v,(basis.T@g[:n+1])[:-1]
 opt=least_squares(lambda z:func(z)[1],x[:-1],xtol=1e-14,gtol=1e-13,ftol=1e-14,max_nfev=3000,diff_step=1e-5)
 cc=basis@np.r_[opt.x,1.];cc/=np.linalg.norm(cc)
 result={'n':n,'k':k,'mode':'min','value':float(pauli_sector(cc,k)),'coefficients':[[float(a),0.] for a in cc],'support':support,'weights':[float(a*a) for a in cc],'gradientNorm':float(np.linalg.norm(opt.fun)),'status':'numerical'}
 (ROOT/f'n{n}-s{k}-min-refined.json').write_text(json.dumps(result,indent=2)+'\n')
 print(n,k,result['value'],'grad',result['gradientNorm'],'coeff',cc[support],flush=True)
