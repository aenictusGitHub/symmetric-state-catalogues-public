#!/usr/bin/env node
/* Numerical point groups of the weighted representatives; no AC claims change. */
'use strict';
const fs=require('node:fs'),path=require('node:path'),crypto=require('node:crypto');
const root=path.resolve(__dirname,'..');
const detector=require(path.join(root,'majorana-symmetry.js'));
const rotations=require(path.join(root,'symmetry-groups.js'));
const points=require(path.join(root,'point-groups.js'));
const catalogue=JSON.parse(fs.readFileSync(path.join(root,'data/catalogue.json'),'utf8'));
const inputs=[...catalogue.states,...(catalogue.numericalStates||[])].filter(s=>s.majoranaPoints?.length);
const tolerance=1e-6;
const states=inputs.map(state=>{
  const detected=detector.detect(state.majoranaPoints,{tolerance});
  const rotationalGroup=rotations.classify(detected);
  const pointGroup=points.classify(state.majoranaPoints,detected,rotationalGroup);
  if(!rotationalGroup.symbol||!pointGroup.symbol)throw Error('Unclassified representative: '+state.id);
  return {id:state.id,n:state.n,partition:state.partition,achievedOrder:state.achievedOrder,
    pointGroup,rotationalGroup};
});
const hash=value=>crypto.createHash('sha256').update(value).digest('hex');
const result={schemaVersion:'1.0',method:'Numerical geometric classification of each weighted constellation, with permutation closure and improper-coset checks.',
  scope:'Groups belong to individual representatives, not to whole multiplicity partitions. Full point groups include improper symmetries; rotational subgroups contain only proper rotations. These checks do not certify exact algebraic symmetry or anticoherence.',
  coordinateTolerance:tolerance,constellationInputsSHA256:hash(JSON.stringify(inputs.map(s=>({id:s.id,majoranaPoints:s.majoranaPoints})))),
  implementationSHA256:Object.fromEntries(['majorana-symmetry.js','symmetry-groups.js','point-groups.js'].map(f=>[f,hash(fs.readFileSync(path.join(root,f)))])),
  stateCount:states.length,states};
fs.writeFileSync(path.join(root,'data/symmetry-groups.json'),JSON.stringify(result,null,2)+'\n');
console.log('Classified '+states.length+' weighted constellations.');
