#!/usr/bin/env python3
"""Numerical audit of manifest-selected published Liège-library counterparts.
Never scans for additional states; numerical evidence is not an exact proof.
"""
import argparse,json,re,math,hashlib,time
from pathlib import Path
from collections import Counter
import mpmath as mp
from analyze_core import *
HERE=Path(__file__).resolve().parent
def localpath(value):
 p=Path(value)
 return p if p.is_absolute() else HERE/p
def portablepath(value):
 p=localpath(value)
 try:return str(p.relative_to(HERE))
 except ValueError:return str(p)
NUMBER=r'[+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?'
def clean(s):
 s=re.sub(r'`{1,2}[+-]?(?:\d+(?:\.\d*)?|\.\d+)','',s)
 return s.replace('*^','e').replace(' ','').strip()
def real_tokens(path):
 data=clean(path.read_text());tokens=data.split()
 # Preserve line boundaries, also accept whitespace/comma-separated exports.
 data=re.sub(r'`{1,2}[+-]?(?:\d+(?:\.\d*)?|\.\d+)','',path.read_text()).replace('*^','e')
 tokens=[x for x in re.split(r'[\s,;{}\[\]()]',data) if x]
 if not all(re.fullmatch(NUMBER,x) for x in tokens):raise ValueError(f'Unrecognized coordinate syntax: {path}')
 return [mp.mpf(x) for x in tokens],tokens

def read_complex(path):
 out=[]
 for original in path.read_text().splitlines():
  fields=re.sub(r'`{1,2}[+-]?(?:\d+(?:\.\d*)?|\.\d+)','',original).replace('*^','e').split()
  if len(fields)==2 and all(re.fullmatch(NUMBER,x) for x in fields):
   out.append(mp.mpc(mp.mpf(fields[0]),mp.mpf(fields[1])));continue
  text=clean(original)
  if not text:continue
  text=text.replace('*I','j').replace('I','j').replace('*','')
  if text.endswith('j'):
   body=text[:-1]
   if body in ('','+','-'):out.append(mp.mpc(0,-1 if body=='-' else 1));continue
   matches=re.findall(NUMBER,body)
   if ''.join(matches)!=body or len(matches) not in (1,2):raise ValueError(f'Unrecognized complex syntax: {original}')
   out.append(mp.mpc(0,mp.mpf(matches[0])) if len(matches)==1 else mp.mpc(mp.mpf(matches[0]),mp.mpf(matches[1])))
  elif re.fullmatch(NUMBER,text):out.append(mp.mpc(mp.mpf(text)))
  else:raise ValueError(f'Unrecognized complex syntax: {original}')
 return out

def compare_dicke(source,target):
 norm=sum(abs(x)**2 for x in source);normalized=[x/mp.sqrt(norm) for x in source];n=len(source)-1;best=None
 for reverse in (False,True):
  for conjugate in (False,True):
   for alternating in (False,True):
    cs=list(reversed(normalized)) if reverse else list(normalized)
    if conjugate:cs=[mp.conj(x) for x in cs]
    if alternating:cs=[(-1)**r*x for r,x in enumerate(cs)]
    overlap=sum(mp.conj(a)*b for a,b in zip(target,cs));infidelity=max(mp.mpf(0),1-abs(overlap)**2)
    if best is None or infidelity<best[0]:
     phase=mp.conj(overlap)/abs(overlap) if overlap else 1
     best=(infidelity,[x*phase for x in cs],{'reverseExcitationOrder':reverse,'complexConjugate':conjugate,'alternateSigns':alternating})
 return norm,best

def analyze(record):
 start=time.time();n=int(record['n']);t=int(record['t']);coord=localpath(record['coordinateLocalfile']);data,tokens=real_tokens(coord);assert len(data)==3*n,(coord,len(data),n);pts=[tuple(data[i:i+3]) for i in range(0,len(data),3)]
 c,unitpts,radialerror=reconstruct(pts);cluster=clustering(unitpts);coordinatecheck=rho_residual(c,t);chosen=c;coeffsource='high-precision physical-spinor reconstruction from coordinate counterpart';dickeobs=None
 dickepath=record.get('DickeLocalfile') or record.get('dickeLocalfile')
 if dickepath:
  raw=read_complex(localpath(dickepath));assert len(raw)==n+1,(dickepath,len(raw),n+1)
  norm,best=compare_dicke(raw,c);infidelity,mapped,convention=best;dickeres=rho_residual(mapped,t);match=infidelity<mp.mpf('1e-20')
  dickeobs={'path':portablepath(dickepath),'sourceNormError':mp.nstr(abs(norm-1),25),'normalizedACCheck':dickeres,'bestInfidelityToCoordinateState':mp.nstr(infidelity,25),'matchedByConventionAndGlobalPhase':bool(match),'conventionTransform':convention,'comparisonScope':'Global phase plus reversal, conjugation and alternating signs; arbitrary spatial rotation was not optimized. Failure to match does not prove the files describe inequivalent states.'}
  if match:chosen=mapped;coeffsource='normalized original Dicke counterpart, convention-matched to coordinate state'
 check=rho_residual(chosen,t) if chosen is not c else coordinatecheck;passed=mp.mpf(check['purityExcess'])<=mp.mpf('1e-20')
 textrows=[tuple(tokens[i:i+3]) for i in range(0,len(tokens),3)];exactgroups=Counter(textrows)
 visual=[]
 for group in cluster['groups']:
  avg=[sum(unitpts[i][k] for i in group)/len(group) for k in range(3)];norm=mp.sqrt(sum(x*x for x in avg));visual.append({'x':float(avg[0]/norm),'y':float(avg[1]/norm),'z':float(avg[2]/norm),'multiplicity':len(group)})
 ident=record.get('id') or f"liege-{record.get('sourcecategory','source')}-n{n}-t{t}-g{record.get('declaredMinimumMultiplicity',1)}"
 out={'id':ident,'n':n,'claimedOrder':t,'achievedOrder':t,'claimNumericallySupported':passed,'partition':cluster['partition'],'maxMultiplicity':cluster['maxMultiplicity'],'distinctPoints':cluster['distinctPoints'],'partitionStatus':'numerical clustering only','exactProof':False,'coefficients':[[float(mp.re(x)),float(mp.im(x))] for x in chosen],'highPrecisionCoefficients':coefficient_strings(chosen,35),'coefficientSource':coeffsource,'support':[r for r,x in enumerate(chosen) if abs(x)>mp.mpf('1e-30')],'supportThreshold':'1e-30 (numerical UI support only)','majoranaPoints':visual,'rawUnitPoints':[[mp.nstr(x,35) for x in p] for p in unitpts],'coordinates':{'path':portablepath(coord),'maximumRadiusErrorBeforeNormalization':mp.nstr(radialerror,25),'literalRepeatedRowPartition':sorted(exactgroups.values(),reverse=True),'ACCheckFromReconstruction':coordinatecheck},'dickeComparison':dickeobs,'ACCheck':check,'normalizedCoefficientNormError':mp.nstr(abs(sum(abs(x)**2 for x in chosen)-1),18),'criterion':{'metric':'squared Frobenius distance rho_t-I/(t+1)','acceptanceThreshold':'1e-20','interpretation':'numerical consistency only; does not prove exact anticoherence'},'clustering':cluster,'sourceRecord':record,'runtimeSeconds':time.time()-start}
 print(ident,'n',n,'t',t,'g',cluster['maxMultiplicity'],'ACsq',check['purityExcess'],'DickeMatch',None if dickeobs is None else dickeobs['matchedByConventionAndGlobalPhase'],flush=True)
 return out

def main():
 p=argparse.ArgumentParser();p.add_argument('--manifest',type=Path,default=HERE/'manifest.json');p.add_argument('--output',type=Path,default=HERE/'observations.json');args=p.parse_args();manifest=json.loads(args.manifest.read_text());records=manifest if isinstance(manifest,list) else manifest.get('records',manifest.get('states',[]));assert records,'No selected source records';results=[];errors=[];missing=[]
 for record in records:
  if not record.get('coordinateLocalfile'):
   missing.append({'id':record.get('id'),'n':record['n'],'claimedOrder':record['t'],'status':'missing coordinate counterpart','reason':'The published filename is listed, but no matching coordinate file was available locally and the host could not be reached. No state or multiplicity is inferred.','sourceRecord':record})
  else:
   try:results.append(analyze(record))
   except Exception as e:errors.append({'sourceRecord':record,'error':repr(e)});print('ERROR',record.get('id'),repr(e),flush=True)
  summary={'indexedRecords':len(records),'analyzedRecords':len(results),'missingRecords':len(missing),'errors':len(errors),'numericallyConsistentAtClaimedOrder':sum(x['claimNumericallySupported'] for x in results),'DickePairsCompared':sum(x['dickeComparison'] is not None for x in results),'DickePairsMatched':sum(x['dickeComparison'] is not None and x['dickeComparison']['matchedByConventionAndGlobalPhase'] for x in results),'unmatchedDickePairsUseCoordinateReconstruction':True}
  args.output.write_text(json.dumps({'method':'90-digit physical-spinor reconstruction; numerical AC and clustering checks. No exact-state proof is inferred.','sourceCaveat':'Files are local counterparts selected by filenames listed on the primary pages; online byte identity could not be checked while the host was unreachable.','summary':summary,'records':results,'missing':missing,'errors':errors},indent=2))
 print('DONE',len(results),'records',len(missing),'missing',len(errors),'errors',flush=True)
if __name__=='__main__':main()
