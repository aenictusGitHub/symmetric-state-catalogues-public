/* Numeric cross-check of independently exact metadata against browser tensors. */
'use strict';
const fs=require('node:fs'),path=require('node:path'),assert=require('node:assert/strict');
const root=path.resolve(__dirname,'..'),O=require('../operator-math.js');
const catalogue=JSON.parse(fs.readFileSync(path.join(root,'data/catalogue.json'))),metadata=JSON.parse(fs.readFileSync(path.join(root,'data/operator-exact-weights.json')));
const records=new Map(catalogue.subspaces.map(r=>[r.id,r]));let checked=0,maxError=0;
for(const[id,state]of Object.entries(metadata.states)){
 const model=O.fromBasis(records.get(id));
 for(const[rank,w]of Object.entries(state.weights)){
  const error=Math.abs(model.multipoles[Number(rank)].weight-w.numericValue);maxError=Math.max(maxError,error);assert.ok(error<1e-8,`${id} rank${rank}: ${error}`);checked++;
 }
}
console.log(JSON.stringify({states:Object.keys(metadata.states).length,weights:checked,maxNumericalError:maxError,status:'PASS'}));
