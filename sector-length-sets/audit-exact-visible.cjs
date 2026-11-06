'use strict';
const fs=require('node:fs'),path=require('node:path'),vm=require('node:vm'),katex=require('./vendor/katex/katex.min.js');
const context={};vm.createContext(context);
for(const file of ['exact-small-n.js','exact-n6.js'])vm.runInContext(fs.readFileSync(path.join(__dirname,'data',file),'utf8'),context);
const exact={...context.EXACT_SMALL_EXPRESSIONS,...context.EXACT_N6_EXPRESSIONS};
const small=JSON.parse(fs.readFileSync(path.join(__dirname,'data/small-n.json'))),n6=JSON.parse(fs.readFileSync(path.join(__dirname,'data/n6.json'))),nodes=JSON.parse(fs.readFileSync(path.join(__dirname,'data/n6-manuscript-nodes.json')));
function sampleBoundary(d){const out=[];for(const curve of d.curves){const count=curve.points.length;for(let j=1;j<10;j++){const p=curve.points[Math.round(j*(count-1)/10)];if(!out.some(q=>Math.hypot(q.sectors[0]-p.sectors[0],q.sectors[1]-p.sectors[1])<.02)&&!d.knownStates.some(q=>Math.hypot(q.sectors[0]-p.sectors[0],q.sectors[1]-p.sectors[1])<.02))out.push({...p,curveId:curve.id,family:curve.label,label:p.label||curve.label});}}return out;}
const output={pass:true,visibleStates:0,byN:{},expressionsRendered:0,parameterBlocksRendered:0,sectorExpressionsRendered:0,missing:[],errors:[]};
function render(id,field,source){try{katex.renderToString(source,{throwOnError:true,strict:false});}catch(e){output.errors.push({id,field,message:e.message});}}
for(const n of [3,4,5,6]){
 const d=n===6?n6:small[n];const named=(d.states||d.knownStates).filter(s=>n===5?s.id!=='n5-dicke2':n!==6||nodes.some(p=>p.attained&&p.stateIds.includes(s.id)));const visible=[...named,...(n===5?sampleBoundary({...d,knownStates:named}):[])];output.byN[n]=visible.length;
 for(const state of visible){output.visibleStates++;const ex=exact[state.id];if(!ex){output.missing.push(state.id);continue;}render(state.id,'latex',ex.latex);output.expressionsRendered++;
  if(ex.parameterLatex){const source=ex.parameterLatex,parts=source.includes('begin{aligned}')?[source]:source.split(/,?\\qquad\s*|,\\quad\s*/);for(const p of parts){render(state.id,'parameterLatex displayed block',p);output.parameterBlocksRendered++;}render(state.id,'parameterLatex full',source);}
  if(ex.sectorLatex){render(state.id,'sectorLatex',ex.sectorLatex);output.sectorExpressionsRendered++;}
 }
}
output.pass=!output.missing.length&&!output.errors.length;
fs.writeFileSync(path.join(__dirname,'exact-visible-audit.json'),JSON.stringify(output,null,2)+'\n');console.log(JSON.stringify(output,null,2));if(!output.pass)process.exit(1);
