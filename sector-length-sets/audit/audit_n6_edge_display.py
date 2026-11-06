"""Read-only diagnosis of high-zoom N6 edge-display artifacts.

Measures existing analytical paths; it does not alter any surface, path,
attainment flag, state, or manuscript source.
"""
from pathlib import Path
import json,math
atlas=Path(__file__).resolve().parents[1]
data=json.loads((atlas/'data/n6-smooth-geometry.json').read_text())
def world(p):
    return [2*(p[0]+1)/9-1,2*(p[2]+1)/23-1,2*(p[1]+1.5)/17.5-1]
rows=[]
for i,e in enumerate(data['edges']):
    points=[world(p) for p in e['points']]
    length=sum(math.dist(a,b) for a,b in zip(points,points[1:]))
    rows.append({'edgeIndex':i,'kind':e['kind'],'attained':e['attained'],
                 'worldArcLength':length,'wholeArcFitsOneLegacyDash':length<=.045})
width,height=460,520
span=max(1.72,1.8/(width/height))
legacy=[{'zoom':z,'dashPixels':.045*height*z/(2*span),
         'gapPixels':.035*height*z/(2*span)} for z in [1,3.5,10,80]]
# A minimal numerical depth example: smaller depth is nearer the camera.
faceDepth=5.;hiddenDepth=6.;laterRearDepth=5.8
hiddenPass=hiddenDepth>faceDepth
corruptedDepth=hiddenDepth if hiddenPass else faceDepth
report={
 'status':'DIAGNOSED',
 'pathCount':len(rows),
 'legacyWorldDash':.045,'legacyWorldGap':.035,
 'representativeViewport':{'width':width,'height':height,'halfSpan':span},
 'legacyDashScaling':legacy,
 'shortNontrivialGrayArcs':[r for r in rows if not r['attained'] and 1e-5<r['worldArcLength']<.045],
 'depthWriteCounterexample':{
    'opaqueFaceDepth':faceDepth,'hiddenLineDepth':hiddenDepth,
    'laterRearSolidLineDepth':laterRearDepth,
    'laterSolidPassesAfterHiddenDepthWrite':laterRearDepth<=corruptedDepth,
    'laterSolidPassesWithReadOnlyDepth':laterRearDepth<=faceDepth},
 'recommendations':[
    'Set depthWrite=false on both solid and dashed edge materials so opaque faces remain the authoritative depth buffer.',
    'Set dashed line dashSize/gapSize from current camera world units per CSS pixel on every render/resize/zoom. Six-pixel dashes and five-pixel gaps avoid high-zoom solid-looking hidden cap curves.',
    'The short null_upper/null_lower/flow paths are explicitly intended analytical envelope creases in compatibility_null_envelope_seams, compatibility_tiny_cap_seams, compatibility_research_internal_rims. Do not delete them merely for being short.',
    'Duplicate full body traces can be deduplicated for stroke consistency, but are distinct from the high-zoom solid rear-fan bug.'
 ],
 'knownDuplicateCoverage':[
    {'primary':[0,8],'duplicate':44,'description':'C5 locus re-emits spectral then orange exposed rim'},
    {'primary':37,'duplicate':47,'description':'octa-GHZ straight edge'},
    {'primary':38,'duplicate':50,'description':'W-product straight edge'},
    {'primary':39,'duplicate':49,'description':'GHZ-product straight edge'},
    {'primary':40,'duplicate':48,'description':'octa-W straight edge'},
    {'primary':4,'duplicate':52,'description':'paired-null attained rim, apart from endpoint tolerance'}
 ],
 'pathMeasurements':rows,
 'scientificDataChanges':0
}
(atlas/'audit/n6-edge-display-diagnosis.json').write_text(json.dumps(report,indent=2)+'\n')
print(json.dumps({k:v for k,v in report.items() if k!='pathMeasurements'},indent=2))
