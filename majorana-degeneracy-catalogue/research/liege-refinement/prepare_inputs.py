#!/usr/bin/env python3
"""Freeze the previously selected numerical source states at100-digit precision."""
import sys,json,time
from pathlib import Path
HERE=Path(__file__).resolve().parent;SOURCE=HERE.parent/'liege-library'
if not SOURCE.exists():SOURCE=HERE.parent/'liege-ac-library'
sys.path.insert(0,str(SOURCE))
from analyze_records import real_tokens,read_complex,compare_dicke,localpath
from analyze_core import reconstruct,mp,coefficient_strings
mp.mp.dps=100
obs=json.loads((SOURCE/'observations.json').read_text());records=[]
for old in obs['records']:
 src=old['sourceRecord'];points,token=real_tokens(localpath(src['coordinateLocalfile']));points=[tuple(points[i:i+3]) for i in range(0,len(points),3)];c,pts,radialerror=reconstruct(points)
 if old['dickeComparison'] and old['dickeComparison']['matchedByConventionAndGlobalPhase']:
  _,best=compare_dicke(read_complex(localpath(src['DickeLocalfile'])),c);assert best[0]<mp.mpf('1e-20');c=best[1]
 poles=[];zeros=set()
 for group in old['clustering']['groups']:
  g=len(group)
  if g<2:continue
  x,y,z=pts[group[0]];assert x==0 and y==0 and abs(z)==1,(old['id'],g,x,y,z)
  indices=list(range(g)) if z<0 else list(range(old['n']-g+1,old['n']+1));zeros.update(indices);poles.append({'pole':'south' if z<0 else 'north','multiplicity':g,'indices':indices})
 assert all(abs(c[i])<mp.mpf('1e-60') for i in zeros)
 for i in zeros:c[i]=mp.mpc(0)
 records.append({'id':old['id'],'n':old['n'],'reportedOrder':old['claimedOrder'],'partition':old['partition'],'originalCoefficients':coefficient_strings(c,95),'originalPoints':[[mp.nstr(x,95) for x in p] for p in pts],'structuralZeroIndices':sorted(zeros),'repeatedPoles':poles,'sourceCategory':src['sourcecategory'],'sourceCoordinateFilename':Path(src['coordinateLocalfile']).name,'sourceCoefficientSelection':old['coefficientSource'],'originalReportedResidual':old['ACCheck']})
(HERE/'original-states.json').write_text(json.dumps({'precisionDigits':100,'source':'Frozen copies of the144 previously selected source records; originals are unchanged.','records':records},indent=2));print('Prepared',len(records),'records')
