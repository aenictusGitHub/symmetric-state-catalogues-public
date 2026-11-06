/* Geometry/XML smoke checks for the standalone vector exporter. */
const fs=require('fs'),path=require('path'),vm=require('vm'),assert=require('assert');
const root=path.resolve(__dirname,'..'),T=require(path.join(root,'vendor/three.min.js'));
const scope={THREE:T,getComputedStyle:()=>({getPropertyValue:()=> '#faf9f5'}),XMLSerializer:class{},document:{documentElement:{}},console};
scope.window=scope;vm.createContext(scope);
vm.runInContext(fs.readFileSync(path.join(root,'data/export-fonts.js'),'utf8'),scope);
vm.runInContext(fs.readFileSync(path.join(root,'svg-export.js'),'utf8'),scope);
function view(){
 const scene=new T.Scene(),content=new T.Group(),camera=new T.OrthographicCamera(-2,2,1.6,-1.6,.01,100);
 camera.position.set(0,0,5);camera.lookAt(0,0,0);camera.updateMatrixWorld();scene.add(content);
 return{scene,content,camera,width:640,height:512,labels:[],svg:{innerHTML:''},markerSvg:{innerHTML:'<circle cx="320" cy="256" r="4.8" fill="#ad7549" stroke="#7f5233"/>'},render(){}};
}
function line(points,material){return new T.Line(new T.BufferGeometry().setFromPoints(points.map(p=>new T.Vector3(...p))),material||new T.LineBasicMaterial({color:0x75978e}));}
const sphere=view();sphere.sphere=new T.Group();sphere.hull=new T.Group();sphere.content.add(sphere.sphere,sphere.hull);
sphere.sphere.add(new T.Mesh(new T.SphereGeometry(1),new T.MeshBasicMaterial({color:0xe5ddc9,transparent:true,opacity:.18,depthWrite:false})));
const ring=Array.from({length:5},(_,i)=>[Math.cos(i*2*Math.PI/5),Math.sin(i*2*Math.PI/5),0]);
const fan=[];for(let i=1;i<4;i++)fan.push(...ring[0],...ring[i],...ring[i+1]);
const g=new T.BufferGeometry();g.setAttribute('position',new T.Float32BufferAttribute(fan,3));
sphere.hull.add(new T.Mesh(g,new T.MeshBasicMaterial({color:0x5f9588,transparent:true,opacity:.16,side:T.DoubleSide,depthWrite:false})));
for(let i=0;i<5;i++)sphere.hull.add(line([ring[i],ring[(i+1)%5]]));
const a=scope.AtlasSVG.serialize(sphere,{title:'Planar pentagon <&> test'});
assert.equal((a.match(/<polygon /g)||[]).length,1,'Coplanar hull must remain a single polygon fill');
assert.equal((a.match(/<path /g)||[]).length,5,'Only five polygon boundary edges');
assert(!a.includes('<image'),'Export must contain no raster image');
assert(!/url\((?!["']?(?:data:|#))/.test(a),'Font and image URLs must be embedded');
assert(a.includes('Planar pentagon &lt;&amp;&gt; test'),'XML title escaping');
fs.writeFileSync(path.join(__dirname,'svg-export-pentagon.svg'),a);

const sector=view();sector.surface=new T.Group();sector.edges=new T.Group();sector.axisGroup=new T.Group();sector.content.add(sector.surface,sector.edges,sector.axisGroup);
sector.camera.position.set(3,2,5);sector.camera.lookAt(0,0,0);sector.camera.updateMatrixWorld();
sector.surface.add(new T.Mesh(new T.BoxGeometry(1,1,1),new T.MeshBasicMaterial({color:0xb894d6})));
sector.edges.add(line([[-.4,0,-.5],[.4,0,-.5]]));
sector.edges.add(line([[-.4,.5,.5],[.4,.5,.5]]));
const hiddenMaterial=new T.LineDashedMaterial({color:0x75978e,depthFunc:T.GreaterDepth});
sector.edges.add(line([[-.4,0,-.5],[.4,0,-.5]],hiddenMaterial));
sector.axisGroup.add(line([[-1,0,-.6],[1,0,-.6]]));
const b=scope.AtlasSVG.serialize(sector,{title:'Opaque body with hidden edges'});
assert.equal((b.match(/<polygon /g)||[]).length,12);
assert(b.includes('stroke-dasharray="'),'Rear body edge must be dashed');
assert(!b.includes('<image'));
fs.writeFileSync(path.join(__dirname,'svg-export-sector.svg'),b);
const report={pass:true,planarHull:{polygons:1,trueEdges:5},sector:{triangles:12,hiddenEdgeDash:true},rasterImages:0,embeddedFonts:true,domLabelRendering:'Requires the root browser export check'};
fs.writeFileSync(path.join(__dirname,'svg-export-geometry.json'),JSON.stringify(report,null,2));
console.log(JSON.stringify(report,null,2));
