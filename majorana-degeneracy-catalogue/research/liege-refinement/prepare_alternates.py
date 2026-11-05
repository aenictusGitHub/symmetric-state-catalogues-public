#!/usr/bin/env python3
"""Freeze unmatched original Dicke counterparts separately; no partition claim."""
import sys,json
from pathlib import Path
HERE=Path(__file__).resolve().parent;SOURCE=HERE.parent/'liege-library'
if not SOURCE.exists():SOURCE=HERE.parent/'liege-ac-library'
sys.path.insert(0,str(SOURCE))
from analyze_records import read_complex,localpath
from analyze_core import mp,coefficient_strings
mp.mp.dps=100;records=[]
for old in json.loads((SOURCE/'observations.json').read_text())['records']:
 if not old['dickeComparison'] or old['dickeComparison']['matchedByConventionAndGlobalPhase']:continue
 src=old['sourceRecord'];c=read_complex(localpath(src['DickeLocalfile']));norm=mp.sqrt(sum(abs(x)**2 for x in c));c=[x/norm for x in c];n=len(c)-1
 low=0
 while low<len(c) and c[low]==0:low+=1
 high=0
 while high<len(c) and c[-1-high]==0:high+=1
 poles=[];zeros=[]
 for name,g,indices in [('south',low,list(range(low))),('north',high,list(range(n-high+1,n+1)))]:
  if g>=2:poles.append({'pole':name,'multiplicity':g,'indices':indices});zeros+=indices
 records.append({'id':old['id']+'-original-dicke-alternate','sourceObservationId':old['id'],'n':n,'reportedOrder':old['claimedOrder'],'partition':None,'originalCoefficients':coefficient_strings(c,95),'originalPoints':[],'structuralZeroIndices':sorted(zeros),'repeatedPoles':poles,'sourceCategory':'unmatched original Dicke alternate','sourceCoordinateFilename':Path(src['DickeLocalfile']).name,'sourceCoefficientSelection':'The original Dicke file, normalized without changing its basis order or global phase. It is not assigned the coordinate counterpart’s partition or sphere plot.','originalReportedResidual':old['dickeComparison']['normalizedACCheck']})
assert len(records)==16
(HERE/'alternate-inputs.json').write_text(json.dumps({'precisionDigits':100,'note':'Separate original Dicke counterparts; no full-partition or coordinate-equivalence claim.','records':records},indent=2));print('Prepared',len(records),'alternates')
