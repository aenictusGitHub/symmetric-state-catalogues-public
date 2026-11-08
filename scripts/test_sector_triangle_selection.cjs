/* Exact N=4 inverse and screen-coordinate selection regression checks. */
const assert = require('node:assert/strict');
const fs = require('node:fs');
const path = require('node:path');
const vm = require('node:vm');
const root = path.resolve(__dirname, '../sector-length-sets');
const M = require(path.join(root, 'majorana-math.js'));
const T = require(path.join(root, 'vendor/three.min.js'));
const katex = require(path.join(root, 'vendor/katex/katex.min.js'));
require(path.join(root, 'data/small-n-constructors.js'));
const C = global.SectorConstructors;
const close = (a, b, tol = 1e-10) => assert.ok(Math.abs(a - b) < tol, `${a} != ${b}`);
function rawSectors(c) {
  const out = [1];
  for (let k = 1; k <= 4; k++) {
    let s = 2 ** k * M.reducedPurity(c, k);
    for (let j = 0; j < k; j++) s -= M.binomial(4 - j, k - j) / M.binomial(4, k) * out[j];
    out[k] = M.binomial(4, k) * s;
  }
  return out.slice(1);
}
let seed = 741, checked = 0, maximumError = 0;
function random() { seed = (1664525 * seed + 1013904223) >>> 0; return seed / 2 ** 32; }
function check(x, y, stars = false) {
  const state = C.n4(x, y); assert.ok(state);
  close(state.coefficients.reduce((s, [r, i]) => s + r*r + i*i, 0), 1, 1e-14);
  rawSectors(state.coefficients).forEach((s, i) => {
    maximumError = Math.max(maximumError, Math.abs(s - state.sectors[i]));
    close(s, state.sectors[i], 1e-10);
  });
  if (stars) {
    const result = M.majoranaStars(state.coefficients);
    assert.equal(result.stars.length, 4); assert.ok(result.reliable);
    assert.ok(result.fidelity > 1 - 1e-10);
  }
  checked++;
}
for (let i = 0; i <= 40; i++) for (let j = 0; j <= 40; j++) {
  const x = i/10, y = x + 2 + (4 - x)*j/40;
  check(x, y, j === 20 && i % 10 === 0);
}
for (let i = 0; i < 1000; i++) { const x = 4*random(); check(x, x+2+(4-x)*random()); }
for (const e of [0, 1e-15, 1e-13, 2.1e-12, 1e-10, 1e-7]) {
  for (const x of [0, 1e-14, .1, 1, 2.5, 3.9]) check(x, 6-e);
  for (const y of [2, 3, 5, 6]) check(e, y+e > 6 ? 6 : y+e);
}
for (const p of [[0,2],[0,6],[4,6],[1.25,4.5]]) check(...p,true);
for (const p of [[-1,4],[2,3],[1,7],[NaN,4],[1,Infinity]]) assert.equal(C.n4(...p),null);

const context = {window:{THREE:T}};
vm.runInNewContext(fs.readFileSync(path.join(root, 'viewer.js'),'utf8'),context);
const View = context.window.SectorView;
const view = Object.create(View.prototype);
Object.assign(view,{n:4,dimension:2,data:{},ranges:[[0,4],[2,6]],width:800,height:600,
  raycaster:new T.Raycaster(),camera:new T.OrthographicCamera(-2,2,1.5,-1.5,.01,100),
  canvas:{getBoundingClientRect:()=>({left:5,top:7}),style:{}},namedWorld:[],namedProjections:[],projected:[],grayPoints:[]});
for (const [tx,ty,zoom] of [[0,0,1],[.4,-.3,2],[-.7,.8,80]]) {
  view.camera.zoom=zoom;view.camera.updateProjectionMatrix();view.camera.position.set(tx,ty,7);
  view.camera.lookAt(tx,ty,0);view.camera.updateMatrixWorld();
  for (const xy of [[0,2],[0,6],[4,6],[1.25,4.5],[2.5,5]]) {
    const p=view.project(view.map(xy)), got=view.triangleCoordinates(p.x,p.y);
    assert.ok(got);close(got[0],xy[0]);close(got[1],xy[1]);
  }
  for (const xy of [[-1,4],[2,3],[2,7]]) {
    const p=view.project(view.map(xy));assert.equal(view.triangleCoordinates(p.x,p.y),null);
  }
}
view.camera.zoom=1;view.camera.position.set(0,0,7);view.camera.lookAt(0,0,0);view.camera.updateProjectionMatrix();view.camera.updateMatrixWorld();
const point=view.project(view.map([1.25,4.5]));let coordinates=null, selected=null;
view.coordinateCallback=p=>coordinates=p;view.selectCallback=s=>selected=s;view.hoverCallback=()=>{};
view.pick({clientX:point.x+5,clientY:point.y+7},true);assert.ok(coordinates);close(coordinates[0],1.25);
const named={id:'known'};view.namedWorld=[{state:named}];view.namedProjections=[point];coordinates=null;
view.pick({clientX:point.x+5,clientY:point.y+7},true);assert.equal(selected,named);assert.equal(coordinates,null);
for (const n of [3,5,6]) {view.n=n;assert.equal(view.triangleCoordinates(point.x,point.y),null);}
view.n=4;view.namedWorld=[];view.labels=[];view.markers=[];view.selection={position:new T.Vector3()};
view.updateLabelVisibility=()=>{};view.render=()=>{};view.addLabel=()=>{};
view.setSelected({sectors:[1,4],interactivePoint:true});view.setSelected({sectors:[2,5],interactivePoint:true});
assert.equal(view.markers.filter(m=>m.kind==='selected-point').length,1);
view.setSelected({sectors:[0,6]});assert.equal(view.markers.length,0);

const app=fs.readFileSync(path.join(root,'app.js'),'utf8');
const helpers={SectorConstructors:C};
vm.runInNewContext(app.slice(app.indexOf('  function triangleState('),app.indexOf('  function notify(')),helpers);
const state=helpers.triangleState(1.25,4.5), restored=helpers.triangleStateFromId(state.id);
assert.equal(restored.id,state.id);assert.equal(JSON.stringify(restored.coefficients),JSON.stringify(state.coefficients));
for (const source of [state.mathLabel.slice(1,-1),state.analyticalExpression.latex,state.analyticalExpression.parameterLatex]) katex.renderToString(source,{throwOnError:true});
for (const id of ['n4-point:foo:4','n4-point:1:Infinity','n4-point:3:2','bad']) assert.equal(helpers.triangleStateFromId(id),null);
console.log(`PASS: ${checked} N=4 witnesses; max sector error ${maximumError.toExponential(2)}; camera/pan/zoom, hit priority, marker lifecycle, share links and LaTeX.`);
