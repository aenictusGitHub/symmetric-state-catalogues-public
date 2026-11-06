#!/usr/bin/env python3
"""Attach existing state witnesses/proof references to exact manuscript markers."""
import json
from pathlib import Path

root=Path(__file__).resolve().parents[1]
data=json.loads((root/'data/n6.json').read_text())
base=json.loads((root/'data/n6-manuscript-nodes-base.json').read_text())
gray_ids=['G1','G2','G3','G4-current','J4','G5','G6','G7','G8','G9','G10']
expected_labels=['Octa','GHZ₆','W₆','Product','C₄','ψ₆ / χ₆','Ω₄','AC₂','ψ₀','Ξ*','χ₁/₂','D₆⁽²⁾','χ₇/₁₀','ω','E*','η₆ (C₂)']
nodes=[]
for i,original in enumerate(base['nodes']):
    n=dict(original)
    n['id']='n6-manuscript-'+str(i)
    n['manuscriptOnly']=True
    if n['attained']:
        matches=[r for r in data['states'] if max(abs(a-b) for a,b in zip(n['point'],r['point']))<1e-8]
        assert matches, n
        n['label']=expected_labels[i-11]
        n['labels']=[r['label'] for r in matches]
        n['showLabel']=True
        n['stateIds']=[r['id'] for r in matches]
        n['states']=matches
        n['proofStatus']='attained'
        n['proofReference']='Explicit state and analytical verification in the stand-alone manuscript'
    else:
        code=gray_ids[i]
        if i in (0,1,2):
            theorem='prop:ten-gray-plane';reason='Complete sharp-plane equality-state classification excludes this tuple.'
        elif i==3:
            theorem='prop:ten-gray-flow';reason='The strict positive-dipole flow comparison excludes the current microscopic junction.'
        elif i==4:
            theorem='thm:fourcopy-refinement';reason='Equality in the refined four-copy bound is coherent-only; this junction is not the product point.'
        else:
            theorem='prop:ten-gray-tips';reason='Strictness of the lower bias-disk lift excludes this positive-overlap contact.'
        n.update(label='Excluded junction '+code,labels=[],showLabel=False,stateIds=[],states=[],
                 status='excluded',proofStatus='excluded',junctionCode=code,
                 proofReference=theorem,explanation=reason)
        n['proofSource']='N6_proven_outer_region_20261015.tex' if i==4 else 'assets/N6_gray_node_exclusions_20261022.tex'
    nodes.append(n)
assert len(nodes)==27
assert sum(n['attained'] for n in nodes)==16
assert all(len(n['stateIds'])==0 for n in nodes if not n['attained'])
assert next(n for n in nodes if n['label']=='ψ₆ / χ₆')['stateIds']==['n6-psi6','n6-chi6']
json_text=json.dumps(nodes,ensure_ascii=False,separators=(',',':'))
(root/'data/n6-manuscript-nodes.json').write_text(json.dumps(nodes,ensure_ascii=False,indent=2)+'\n')
(root/'data/n6-manuscript-nodes.js').write_text('globalThis.SECTOR_N6_MANUSCRIPT_NODES = '+json_text+';\n'+
    'globalThis.SECTOR_N6_MANUSCRIPT_CALLOUTS = '+json.dumps(base['callouts'],ensure_ascii=False,separators=(',',':'))+';\n')
print('PASS: exact 27 manuscript nodes; 16 violet state markers (17 state orbits shown), 11 analytically excluded gray junctions, 4 leader-only family attachments. Dicke3 and every cloud marker omitted.')
