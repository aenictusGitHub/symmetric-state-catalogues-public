#!/usr/bin/env python3
"""Match refined coefficients to nearby simple Majorana points numerically.
Repeated poles are held exactly; independent Rouché audit is separate.
"""
import json,math,time
from pathlib import Path
import mpmath as mp
HERE=Path(__file__).resolve().parent;mp.mp.dps=110

def evaluate(q,z):
 value=q[-1];deriv=mp.mpc(0)
 for a in reversed(q[:-1]):deriv=deriv*z+value;value=value*z+a
 return value,deriv

def compute(state):
 start=time.time();n=state['n'];c=[mp.mpc(*x) for x in state['refinedCoefficients']];original=[tuple(mp.mpf(x) for x in row) for row in state['originalPoints']];poles={x['pole']:x['multiplicity'] for x in state['repeatedPoles']};low=poles.get('south',0);high=poles.get('north',0)
 assert all(c[i]==0 for i in range(low)) and all(c[n-i]==0 for i in range(high))
 q=[(-1)**r*mp.sqrt(math.comb(n,r))*c[r] for r in range(low,n-high+1)];scale=max(map(abs,q));q=[x/scale for x in q];points=[];iterations=[];maxcorr=mp.mpf(0);maxshift=mp.mpf(0);root_errors=[]
 for index,(x,y,z) in enumerate(original):
  pole='north' if z>0 else 'south'
  if x==0 and y==0 and abs(z)==1 and pole in poles:
   pt=(mp.mpf(0),mp.mpf(0),mp.mpf(1 if pole=='north' else -1));steps=0;correction=mp.mpf(0)
  else:
   north=z>=0;w=(x+mp.j*y)/(1+z) if north else (x-mp.j*y)/(1-z);coeff=list(reversed(q)) if north else q;steps=0
   for step in range(24):
    value,deriv=evaluate(coeff,w)
    if deriv==0:raise ValueError(f'zero derivative {state["id"]} point{index}')
    delta=value/deriv;w-=delta;steps+=1
    if abs(delta)<mp.mpf('1e-75'):break
   value,deriv=evaluate(coeff,w);correction=abs(value/deriv)
   if correction>mp.mpf('1e-65'):root_errors.append({'pointIndex':index,'chartRootCorrection':mp.nstr(correction,20)})
   den=1+abs(w)**2
   pt=(2*mp.re(w)/den,2*mp.im(w)/den,(1-abs(w)**2)/den) if north else (2*mp.re(w)/den,-2*mp.im(w)/den,(abs(w)**2-1)/den)
  shift=mp.sqrt(sum((a-b)**2 for a,b in zip(pt,(x,y,z))));maxshift=max(maxshift,shift);maxcorr=max(maxcorr,correction);points.append(pt);iterations.append(steps)
 grouped=[];seen=set()
 for pt in points:
  pole='north' if pt[2]>0 else 'south'
  if pt[0]==0 and pt[1]==0 and abs(pt[2])==1 and pole in poles:
   if pole in seen:continue
   seen.add(pole);multiplicity=poles[pole]
  else:multiplicity=1
  grouped.append({'x':float(pt[0]),'y':float(pt[1]),'z':float(pt[2]),'multiplicity':multiplicity})
 partition=sorted([p['multiplicity'] for p in grouped],reverse=True);assert partition==state['partition']
 # Every source point retains its index. Check distinct refined locations stay far apart.
 distinct=[tuple(mp.mpf(p[k]) for k in ('x','y','z')) for p in grouped];separation=min([mp.sqrt(sum((a-b)**2 for a,b in zip(p,r))) for i,p in enumerate(distinct) for r in distinct[:i]]+[mp.inf])
 result={'id':state['id'],'n':n,'partition':partition,'refinedMajoranaPoints':grouped,'highPrecisionRefinedPoints':[[mp.nstr(a,75) for a in pt] for pt in points],'maxPointDisplacement':mp.nstr(maxshift,30),'minimumDistinctPointSeparation':mp.nstr(separation,25),'maximumFinalChartRootCorrection':mp.nstr(maxcorr,25),'maximumNewtonSteps':max(iterations,default=0),'numericalRootErrors':root_errors,'method':'110-digit Newton correction of each simple root from its original sphere point, using reciprocal charts; repeated poles are fixed exactly. This is numerical matching, with an independent Rouché check supplied separately.','runtimeSeconds':time.time()-start}
 print(state['id'],'maxshift',mp.nstr(maxshift,6),'rootcorr',mp.nstr(maxcorr,5),'steps',max(iterations),'errors',len(root_errors),'seconds',round(time.time()-start,2),flush=True);return result

def main():
 results=[];errors=[]
 for state in json.loads((HERE/'refined.json').read_text())['records']:
  try:results.append(compute(state))
  except Exception as e:errors.append({'id':state['id'],'error':repr(e)});print('ERROR',state['id'],repr(e),flush=True)
  (HERE/'refined-points.json').write_text(json.dumps({'records':results,'errors':errors},indent=2))
 print('DONE',len(results),'errors',len(errors),flush=True)
if __name__=='__main__':main()
