#!/usr/bin/env python3
"""Refine only reported AC orders; preserve structural zeros at repeated poles.
The output remains numerical, never an exact-state certificate.
"""
import argparse,json,math,time
from pathlib import Path
import numpy as np
from scipy.linalg import svd
import mpmath as mp
HERE=Path(__file__).resolve().parent;mp.mp.dps=100

def strings(c,digits=85):return [[mp.nstr(mp.re(x),digits),mp.nstr(mp.im(x),digits)] for x in c]

def setup(n,t,active):
 pairs=[];sqrt2=mp.sqrt(2)
 for a in range(t+1):
  for b in range(a,t+1):
   terms=[(a+q,b+q,mp.sqrt(mp.mpf(math.comb(t,a))*math.comb(t,b))*math.comb(n-t,q)/mp.sqrt(mp.mpf(math.comb(n,a+q))*math.comb(n,b+q))) for q in range(n-t+1)]
   pairs.append((a,b,terms))
 def residual(c):
  F=[];abs2=[mp.re(x)**2+mp.im(x)**2 for x in c]
  for a,b,terms in pairs:
   if a==b:F.append(mp.fdot([w for r,s,w in terms],[abs2[r] for r,s,w in terms])-mp.mpf(1)/(t+1))
   else:
    value=mp.fdot([w for r,s,w in terms],[c[r]*mp.conj(c[s]) for r,s,w in terms]);F.extend([sqrt2*mp.re(value),sqrt2*mp.im(value)])
  return F
 def jacobian(c):
  J=np.zeros(((t+1)**2,2*(n+1)));re=np.array([float(mp.re(x)) for x in c]);im=np.array([float(mp.im(x)) for x in c]);row=0
  for a,b,terms in pairs:
   if a==b:
    for r,s,w in terms:
     h=float(w);J[row,r]+=2*h*re[r];J[row,n+1+r]+=2*h*im[r]
    row+=1
   else:
    for r,s,w in terms:
     h=math.sqrt(2)*float(w)
     J[row,r]+=h*re[s];J[row,s]+=h*re[r];J[row,n+1+r]+=h*im[s];J[row,n+1+s]+=h*im[r]
     J[row+1,r]-=h*im[s];J[row+1,s]+=h*im[r];J[row+1,n+1+r]+=h*re[s];J[row+1,n+1+s]-=h*re[r]
    row+=2
  return J[:,active+[n+1+i for i in active]]
 return residual,jacobian

def refine(state,target=mp.mpf('1e-65'),maxiter=16):
 start=time.time();n=state['n'];t=state['reportedOrder'];c=[mp.mpc(*pair) for pair in state['originalCoefficients']];original=list(c);zeros=set(state['structuralZeroIndices']);active=[r for r in range(n+1) if r not in zeros]
 assert all(c[i]==0 for i in zeros);Ffun,Jfun=setup(n,t,active);F=Ffun(c);norm=lambda f:mp.sqrt(mp.fdot(f,f));before=norm(F);history=[mp.nstr(before,30)];ranks=[];lastrefresh=-99;Rmp=None;singular=None;family=None
 if t==2 and len(state['repeatedPoles'])==1:
  pole=state['repeatedPoles'][0];g=pole['multiplicity']
  if n==3*g+1:
   a,b=(g,n) if pole['pole']=='south' else (0,n-g)
   pa=(mp.mpf(b)-mp.mpf(n)/2)/(b-a);pb=1-pa
   c=[mp.mpc(0)]*(n+1)
   c[a]=mp.sqrt(pa)*original[a]/abs(original[a]);c[b]=mp.sqrt(pb)*original[b]/abs(original[b])
   F=Ffun(c);history.append(mp.nstr(norm(F),30))
   family={'name':'Moment-saturated second-order two-term family','support':[a,b],'probabilities':[mp.nstr(pa,90),mp.nstr(pb,90)],'phaseChoice':'Each nonzero coefficient retains its original phase; this is the closest vector in the forced two-term family in Euclidean coefficient distance.','additionalZeroIndices':[r for r in range(n+1) if r not in (a,b) and r not in zeros]}
 for iteration in range(maxiter):
  current=norm(F)
  if current<=target:break
  if Rmp is None or (iteration-lastrefresh>=4 and len(history)>1 and mp.mpf(history[-1])>=mp.mpf(history[-2])*mp.mpf('1e-6')):
   J=Jfun(c);U,S,Vh=svd(J,full_matrices=False,check_finite=True,lapack_driver='gesdd');threshold=max(J.shape)*np.finfo(float).eps*S[0]*10;rank=int(np.count_nonzero(S>threshold));R=(Vh[:rank].T/S[:rank])@U[:,:rank].T
   Rmp=[[mp.mpf(float(v)) for v in row] for row in R];ranks.append(rank);singular={'largest':float(S[0]),'smallestRetained':float(S[rank-1]),'cutoff':float(threshold),'fullResidualDimension':len(F),'realVariableCount':2*len(active)};lastrefresh=iteration
  delta=[-mp.fdot(row,F) for row in Rmp];m=len(active);candidate=list(c)
  for j,r in enumerate(active):candidate[r]+=mp.mpc(delta[j],delta[j+m])
  nextF=Ffun(candidate);nextnorm=norm(nextF)
  # The source is already close; refuse a large departure or an unstable step.
  if max(abs(a-b) for a,b in zip(candidate,original))>mp.mpf('1e-5'):
   break
  if nextnorm>current*mp.mpf('1.1'):
   Rmp=None
   if iteration>=maxiter-2:break
   continue
  c,F=candidate,nextF;history.append(mp.nstr(nextnorm,30))
 after=norm(F);assert all(c[i]==0 for i in zeros)
 serialized=strings(c);roundtrip=[mp.mpc(*pair) for pair in serialized];serializedResidual=norm(Ffun(roundtrip));normerror=abs(sum(abs(x)**2 for x in roundtrip)-1);delta=max(abs(a-b) for a,b in zip(c,original));passed=serializedResidual<=mp.mpf('1e-60')
 result={'id':state['id'],'n':n,'reportedOrder':t,'partition':state['partition'],'status':'target reached' if passed else 'target not reached','exactProof':False,'precisionDigits':100,'serializedDigits':85,'targetHSResidual':'1e-60','workingTargetHSResidual':mp.nstr(target),'beforeHSResidual':mp.nstr(before,40),'afterHSResidual':mp.nstr(after,40),'serializedHSResidual':mp.nstr(serializedResidual,40),'beforeSquaredResidual':mp.nstr(before**2,40),'afterSquaredResidual':mp.nstr(serializedResidual**2,40),'normalizedCoefficientNormError':mp.nstr(normerror,40),'maxComponentChange':mp.nstr(delta,40),'coefficientDistance':mp.nstr(mp.sqrt(sum(abs(a-b)**2 for a,b in zip(c,original))),40),'iterationResiduals':history,'acceptedSteps':len(history)-1,'jacobianRanks':ranks,'singularValueSummary':singular,'structuralZeroIndices':sorted(zeros),'repeatedPoles':state['repeatedPoles'],'refinedCoefficients':serialized,'original90DigitCoefficients':state['originalCoefficients'],'originalPoints':state['originalPoints'],'sourceCoordinateFilename':state['sourceCoordinateFilename'],'singularFamilyCorrection':family,'sourceCoefficientSelection':state['sourceCoefficientSelection'],'method':'Analytic real Jacobian and fixed approximate right inverse, with100-digit residuals and corrections; no higher-order search.','runtimeSeconds':time.time()-start}
 print(state['id'],'n',n,'t',t,'before',mp.nstr(before,5),'after',mp.nstr(serializedResidual,5),'steps',len(history)-1,'rank',ranks,'seconds',round(time.time()-start,2),flush=True);return result

def main():
 p=argparse.ArgumentParser();p.add_argument('--input',type=Path,default=HERE/'original-states.json');p.add_argument('--output',type=Path,default=HERE/'refined.json');p.add_argument('--only',default=None);p.add_argument('--resume',action='store_true');a=p.parse_args();states=json.loads(a.input.read_text())['records'];results=[];errors=[]
 if a.only:states=[x for x in states if x['id']==a.only]
 existing={}
 if a.resume and a.output.exists():existing={x['id']:x for x in json.loads(a.output.read_text())['records'] if x['status']=='target reached'}
 for state in states:
  try:
   if state['id'] in existing:
    result=existing[state['id']];result['sourceCoefficientSelection']=state['sourceCoefficientSelection'];result.setdefault('singularFamilyCorrection',None)
   else:result=refine(state)
   results.append(result)
  except Exception as e:errors.append({'id':state['id'],'error':repr(e)});print('ERROR',state['id'],repr(e),flush=True)
  summary={'requested':len(states),'completed':len(results),'targetReached':sum(x['status']=='target reached' for x in results),'failedTarget':sum(x['status']!='target reached' for x in results),'errors':len(errors)}
  a.output.write_text(json.dumps({'purpose':'Improve numerical accuracy at the original reported order. All results remain numerical. Original files are unchanged.','summary':summary,'records':results,'errors':errors},indent=2))
 print('DONE',summary,flush=True)
if __name__=='__main__':main()
