/* Exact N=3 interval witnesses and actual line-picking regressions. */
'use strict';
const assert=require('node:assert/strict'),fs=require('node:fs'),path=require('node:path'),vm=require('node:vm');
const root=path.resolve(__dirname,'../sector-length-sets');
const M=require(path.join(root,'majorana-math.js')),T=require(path.join(root,'vendor/three.min.js'));
const katex=require(path.join(root,'vendor/katex/katex.min.js'));
require(path.join(root,'data/small-n-constructors.js'));
const C=global.SectorConstructors;
const close=(a,b,tol=1e-10)=>assert.ok(Number.isFinite(a)&&Number.isFinite(b)&&Math.abs(a-b)<=tol,`${a} != ${b}`);
let randomSeed=39217;
const random=()=>{randomSeed=(1664525*randomSeed+1013904223)>>>0;return randomSeed/2**32;};
let witnesses=0,maximumSectorError=0,majoranaChecks=0;
function rawSectors(coefficients){
  const x=3*(2*M.reducedPurity(coefficients,1)-1);
  const y=3*(4*M.reducedPurity(coefficients,2)-1-2*x/3);
  const z=8*M.reducedPurity(coefficients,3)-1-x-y;
  return [x,y,z];
}
function checkState(x,withStars=false){
  const state=C.n3(x);assert.ok(state,'Every point of the closed interval has a witness');
  assert.equal(state.coefficients.length,4);
  close(state.coefficients.reduce((sum,[r,i])=>sum+r*r+i*i,0),1,1e-14);
  const expected=[x,3,4-x],physical=rawSectors(state.coefficients);
  physical.forEach((v,i)=>{maximumSectorError=Math.max(maximumSectorError,Math.abs(v-expected[i]));close(v,expected[i],2e-12);close(state.sectors[i],expected[i],1e-14);});
  if(withStars){const stars=M.majoranaStars(state.coefficients);assert.equal(stars.stars.length,3);assert.ok(stars.reliable);assert.ok(stars.fidelity>1-1e-10);majoranaChecks++;}
  witnesses++;return state;
}
for(let i=0;i<=300;i++)checkState(i/100);
for(let i=0;i<1000;i++)checkState(3*random());
for(const epsilon of [0,1e-15,1e-13,1e-10,1e-7,1e-4]){checkState(epsilon);checkState(3-epsilon);}
for(const x of [0,.1,1,1.5,2.5,3])checkState(x,true);
for(const x of [-1e-5,3.00001,NaN,Infinity,-Infinity])assert.equal(C.n3(x),null);

// Run the current production viewer, with an actual orthographic camera and
// ray/plane projection; no duplicate screen-coordinate implementation.
const context={window:{THREE:T,SectorConstructors:C},AtlasMath:{state:(_id,label)=>label||''}};
vm.runInNewContext(fs.readFileSync(path.join(root,'viewer.js'),'utf8'),context);
const View=context.window.SectorView,view=Object.create(View.prototype);
Object.assign(view,{n:3,dimension:1,data:{},ranges:[[0,3],[0,1]],width:800,height:600,
  raycaster:new T.Raycaster(),camera:new T.OrthographicCamera(-2,2,1.5,-1.5,.01,100),
  canvas:{getBoundingClientRect:()=>({left:5,top:7}),style:{}},
  namedWorld:[],namedProjections:[],projected:[],grayPoints:[]});
function camera(tx,ty,zoom){view.camera.zoom=zoom;view.camera.updateProjectionMatrix();view.camera.position.set(tx,ty,7);view.camera.lookAt(tx,ty,0);view.camera.updateMatrixWorld();}
let screenHits=0,offLineRejects=0;
for(const [tx,ty,zoom]of [[0,0,1],[.4,-.3,2],[-.7,.8,80]]){
  camera(tx,ty,zoom);
  for(const x of [0,1e-8,.1,1/3,1.25,2,3-1e-8,3]){
    const screen=view.project(view.map([x]));
    for(const offset of [0,-9,9]){
      const got=view.regionCoordinates(screen.x,screen.y+offset);assert.ok(got,`Line hit remains selectable at zoom ${zoom}`);assert.equal(got.length,1);close(got[0],x);screenHits++;
    }
    for(const offset of [-12,12]){assert.equal(view.regionCoordinates(screen.x,screen.y+offset),null,'The line hit target remains 10 screen pixels at every zoom');offLineRejects++;}
  }
  for(const x of [-.001,3.001]){const p=view.project(view.map([x]));assert.equal(view.regionCoordinates(p.x,p.y),null,'A click outside the closed interval is rejected');}
}
camera(0,0,1);
const screen=view.project(view.map([1.25]));let coordinates=null,selected=null,hovered='unchanged';
view.coordinateCallback=p=>coordinates=p;view.selectCallback=s=>selected=s;view.hoverCallback=s=>hovered=s;
const click={clientX:screen.x+5,clientY:screen.y+7};
view.pick(click,true);assert.ok(coordinates);assert.equal(coordinates.length,1);close(coordinates[0],1.25);
assert.equal(view.canvas.style.cursor,'pointer');
view.pick({clientX:screen.x+5,clientY:screen.y+19},true);assert.equal(selected,null);
const named={id:'known-n3-node'};view.namedWorld=[{state:named}];view.namedProjections=[screen];coordinates=null;
view.pick(click,true);assert.equal(selected,named);assert.equal(coordinates,null,'A named node wins over the interval underneath it');
view.pick(click,false);assert.equal(hovered,named);
view.namedWorld=[];view.namedProjections=[];
const sampled={id:'n3-boundary-sample'};view.projected=[{...screen,point:{state:sampled}}];selected=null;
view.pick(click,true);assert.equal(selected,sampled,'Existing sampled nodes keep priority');assert.equal(coordinates,null);
view.projected=[];

// Generic interval selection creates one transient marker, then removes it
// cleanly when a named state is selected. Persistent named markers survive.
const persistent={kind:'sector-violet',world:view.map([0])};
view.labels=[];view.markers=[persistent];view.selection={position:new T.Vector3()};
view.updateLabelVisibility=()=>{};view.render=()=>{};view.addLabel=()=>{};
view.setSelected({id:'n3-point:1',sectors:[1,3,3],interactivePoint:true});
view.setSelected({id:'n3-point:2',sectors:[2,3,2],interactivePoint:true});
assert.equal(view.markers.filter(m=>m.kind==='selected-point').length,1);
assert.ok(view.markers.includes(persistent));close(view.selectedWorld[0],view.map([2])[0]);
view.setSelected({id:'n3-ghz',sectors:[0,3,4]});assert.equal(view.markers.length,1);assert.equal(view.markers[0],persistent);

// Execute the application's real state/share-ID builders, rather than a test
// copy of their normalization or formula formatting.
const app=fs.readFileSync(path.join(root,'app.js'),'utf8');
const start=app.indexOf('  function threeQubitState('),end=app.indexOf('  function notify(');
assert.ok(start>=0&&end>start,'The N=3 application constructor is present');
const helpers={SectorConstructors:C};vm.runInNewContext(app.slice(start,end),helpers);
let restoredChecks=0;
for(const x of [0,1e-14,.1,1/3,1.25,2.9,3-1e-14,3,...Array.from({length:25},()=>3*random())]){
  const state=helpers.threeQubitState(x);assert.ok(state);assert.ok(state.id.startsWith('n3-point:'));assert.equal(state.interactivePoint,true);
  const restored=helpers.threeQubitStateFromId(state.id);assert.ok(restored);assert.equal(restored.id,state.id);
  assert.equal(JSON.stringify(restored.coefficients),JSON.stringify(state.coefficients));
  M.sectorLengths(state.coefficients).forEach((v,i)=>close(v,state.sectors[i],2e-12));
  for(const source of [state.mathLabel.slice(1,-1),state.analyticalExpression.latex,state.analyticalExpression.parameterLatex])if(source)katex.renderToString(source,{throwOnError:true});
  restoredChecks++;
}
for(const id of ['bad','n3-point:foo','n3-point:NaN','n3-point:Infinity','n3-point:-0.01','n3-point:3.01','n3-point:1:2','n4-point:1:3'])assert.equal(helpers.threeQubitStateFromId(id),null,'Invalid share IDs never produce a state');
for(const x of [-.01,3.01,NaN,Infinity])assert.equal(helpers.threeQubitState(x),null);
console.log(JSON.stringify({status:'PASS',witnesses,maximumSectorError,majoranaChecks,screenHits,offLineRejects,restoredChecks,nodePriority:'PASS',selectedMarkerLifecycle:'PASS',LaTeX:'PASS'}));
