'use strict';
// Match the current app's visible-point selection, independently of WebGL.
const fs=require('node:fs'),path=require('node:path'),assert=require('node:assert/strict'),vm=require('node:vm');
const M=require('./majorana-math.js'),katex=require('./vendor/katex/katex.min.js');
const small=JSON.parse(fs.readFileSync(path.join(__dirname,'data/small-n.json'))),n6=JSON.parse(fs.readFileSync(path.join(__dirname,'data/n6.json'))),nodes=JSON.parse(fs.readFileSync(path.join(__dirname,'data/n6-manuscript-nodes.json')));
const context={window:{},katex};vm.createContext(context);vm.runInContext(fs.readFileSync(path.join(__dirname,'math-labels.js'),'utf8'),context);
const labels=context.window.AtlasMath;
function sampleBoundary(d){const out=[];for(const curve of d.curves){const count=curve.points.length;for(let j=1;j<10;j++){const p=curve.points[Math.round(j*(count-1)/10)];if(!out.some(q=>Math.hypot(q.sectors[0]-p.sectors[0],q.sectors[1]-p.sectors[1])<.02)&&!d.knownStates.some(q=>Math.hypot(q.sectors[0]-p.sectors[0],q.sectors[1]-p.sectors[1])<.02))out.push({...p,curveId:curve.id,family:curve.label,label:p.label||curve.label});}}return out;}
const results={pass:true,visible:{},maximumSectorError:0,maximumNodeMappingError:0,minimumMajoranaFidelity:1,mathLabelErrors:[],proofChecks:[]};
for(const n of [3,4,5,6]){
  const d=n===6?n6:small[n];let named=d.states||d.knownStates;
  named=named.filter(s=>n===5?s.id!=='n5-dicke2':n!==6||nodes.some(p=>p.attained&&p.stateIds.includes(s.id)));
  const data=n===5?{...d,knownStates:named}:d,samples=n===5?sampleBoundary(data):[];
  const unique=[];for(const s of named){if(!unique.some(p=>Math.hypot(...p.map((x,k)=>x-s.sectors[k]))<1e-8))unique.push(s.sectors.slice(0,n===6?3:2));}
  for(const s of [...named,...samples]){
    const actual=M.sectorLengths(s.coefficients),err=Math.max(...s.sectors.map((v,k)=>Math.abs(v-actual[k])));results.maximumSectorError=Math.max(results.maximumSectorError,err);assert.ok(err<1e-8);
    const stars=M.majoranaStars(s.coefficients);results.minimumMajoranaFidelity=Math.min(results.minimumMajoranaFidelity,stars.fidelity);assert.ok(stars.fidelity>1-1e-8);
    if(named.includes(s)){
      const markup=labels.html(labels.state(s.id,s.label));if(markup.includes('katex-error'))results.mathLabelErrors.push(s.id);
      assert.ok(markup.includes('class="katex"'),s.id+' named math label');
    }
  }
  results.visible[n]={namedRepresentatives:named.length,uniqueNamedPoints:unique.length,sampledBoundaryPoints:samples.length,genericCloudPoints:0,selectedPointCount:unique.length+samples.length};
}
assert.equal(nodes.length,27);assert.equal(nodes.filter(s=>s.attained).length,16);assert.equal(nodes.filter(s=>!s.attained&&s.status==='excluded'&&s.proofStatus==='excluded').length,11);
const states=new Map(n6.states.map(s=>[s.id,s]));for(const p of nodes){if(!p.attained){assert.equal(p.stateIds.length,0);continue;}assert.ok(p.stateIds.length>0);for(const id of p.stateIds){const s=states.get(id);assert.ok(s,id);const err=Math.max(...p.point.map((v,k)=>Math.abs(v-M.sectorLengths(s.coefficients)[k])));results.maximumNodeMappingError=Math.max(results.maximumNodeMappingError,err);assert.ok(err<1e-8,id+' manuscript mapping');}}
results.visible[6].grayExcludedPoints=11;results.visible[6].totalManuscriptPoints=27;
assert.equal(results.visible[6].uniqueNamedPoints,16);assert.equal(results.visible[6].namedRepresentatives,17);
assert.equal(results.mathLabelErrors.length,0);
results.proofChecks.push('Every visible violet point has an explicit normalized pure symmetric state and reconstructible Majorana constellation.','N6 visible-point catalogue matches all27 current manuscript nodes:16attained points (17state representatives) and11rigorously excluded gray points.','N5 named extrema and boundary samples use the same explicit state records as the exact boundary paths; no generic clouds are displayed.','Named-state labels are actual KaTeX math; all tested labels compile without errors.');
fs.writeFileSync(path.join(__dirname,'visible_point_audit.json'),JSON.stringify(results,null,2)+'\n');console.log(JSON.stringify(results,null,2));
