const fs=require('fs'),vm=require('vm'),path=require('path'),assert=require('assert/strict');
const base=path.resolve(__dirname,'../../..'),names=['symmetric-sector-catalogue','anticoherent-state-catalogue','ac-subspace-catalogue','majorana-degeneracy-catalogue'];
const out=process.argv[2];if(out)fs.mkdirSync(out,{recursive:true});
(async()=>{
 for(const name of names){
  const root=path.join(base,name);let blob=null,filename='',removed=false,revoked=false,timer=null,frame=null;
  const win={devicePixelRatio:2,addEventListener(){},removeEventListener(){}};
  const anchor={click(){filename=this.download;},remove(){removed=true;}};
  const context={window:win,Blob,URL:{createObjectURL(value){blob=value;return 'blob:local-test';},revokeObjectURL(url){assert.equal(url,'blob:local-test');revoked=true;}},document:{createElement:tag=>{assert.equal(tag,'a');return anchor;},body:{append(){}}},setTimeout:f=>{timer=f;},requestAnimationFrame:f=>(frame=f,1),cancelAnimationFrame(){frame=null;}};
  for(const file of ['majorana-hull.js','majorana-svg.js','majorana-symmetry.js'])vm.runInNewContext(fs.readFileSync(path.join(root,file),'utf8'),context);
  let source=fs.readFileSync(path.join(root,name==='majorana-degeneracy-catalogue'?'app.js':'majorana-viewer.js'),'utf8');
  if(name==='majorana-degeneracy-catalogue')source='(()=>{const symmetryColor=()=>"#34786d";'+source.slice(source.indexOf('  // Axes are unoriented.'),source.indexOf('  function renderProofs()'))+'window.test={attachSphere};})();';
  else source=source.replace('  window.MajoranaViewer=','  window.test={attachSphere};\n  window.MajoranaViewer=');
  vm.runInNewContext(source,context);
  const el=()=>({handlers:{},checked:false,value:'all',addEventListener(k,f){this.handlers[k]=f;},removeEventListener(k){delete this.handlers[k];}});
  const drawing=new Proxy({createRadialGradient:()=>({addColorStop(){}})},{get:(o,k)=>o[k]||(()=>{})});
  const canvas={...el(),width:1200,height:800,getContext:()=>drawing,getBoundingClientRect:()=>({width:600,height:400}),classList:{remove(){}},hasPointerCapture:()=>false};
  const controls={download:el(),exportLabel:'Example <state> & "test"',hull:el(),axes:el(),align:el(),reset:el(),symmetryToggle:el(),symmetryOrder:el()};
  const points=[[1,1,1],[1,-1,-1],[-1,1,-1],[-1,-1,1]].map(([x,y,z])=>({x:x/Math.sqrt(3),y:y/Math.sqrt(3),z:z/Math.sqrt(3),multiplicity:2}));
  controls.symmetry=win.MajoranaSymmetry.detect(points);controls.symmetryToggle.checked=true;
  let dispose=win.test.attachSphere(canvas,points,controls);
  const download=async()=>{controls.download.handlers.click();assert.equal(blob.type,'image/svg+xml;charset=utf-8');assert.ok(removed);assert.ok(filename.endsWith('.svg'));const s=await blob.text();assert.ok(!/NaN|Infinity|<image\b|<script\b|<foreignObject\b/.test(s));assert.ok(s.includes('viewBox="0 0 600 400"'));assert.equal(canvas.width,1200,'export does not resize the live canvas');timer();assert.ok(revoked);return s;};
  const diamondCentre=svg=>{const d=svg.match(/<path d="([^"]+)"[^>]*stroke="#a12b75"/);assert.ok(d,'barycentre diamond exported');const p=d[1].match(/^M([-\d.]+) ([-\d.]+)/);return [Number(p[1]),Number(p[2])+6];};
  const near=(a,b)=>assert.ok(Math.abs(a-b)<2e-6,`${a} != ${b}`);
  const initial=await download();near(diamondCentre(initial)[0],300);near(diamondCentre(initial)[1],200);assert.ok(initial.includes('stroke="#35424a"'),'sphere centre cross exported');assert.ok(initial.includes('Example &lt;state&gt; &amp; &quot;test&quot;'));assert.ok(initial.includes('rgb(58,126,119)'),'hull present');assert.ok(initial.includes('C₃'));assert.ok(initial.includes('>2</text>'));assert.ok(!initial.includes('>1</text>'));assert.ok(!initial.includes('>x</text>'));
  assert.ok(initial.includes('153.9'),'95% radius');
  const axisLines=[...initial.matchAll(/<path d="M300 200 L([-\d.]+) ([-\d.]+)"[^>]*stroke-width="1.8"/g)];
  assert.ok(axisLines.length>=6,'symmetry axes exported');
  for(const [,x,y]of axisLines)assert.ok(Math.abs(Math.hypot(Number(x)-300,Number(y)-200)-153.9*1.08)<2e-6,'both axis ends extend 8% past the sphere outline');
  assert.ok(!/rgba\(|feDropShadow|\sfr=/.test(initial),'no unsupported SVG colour/filter/radial-radius syntax');
  assert.ok(initial.includes('fill-opacity="0.16"')&&initial.includes('fill-opacity="0.07"'),'hull transparency explicit');
  assert.ok(initial.includes('stop-opacity="0.12"')&&initial.includes('stop-opacity="0"'),'ground shadow fades to transparent');
  assert.ok(initial.includes('flood-opacity="0.19"')&&initial.includes('<feGaussianBlur'),'SVG 1.1 star shadow');
  assert.ok(initial.includes('stroke-opacity="0.28"'),'rear hull edge alpha explicit');
  canvas.handlers.keydown({key:'ArrowRight',preventDefault(){}});const turned=await download();assert.notEqual(turned,initial,'export captures rotation even before the next animation frame');
  controls.axes.checked=true;const withAxes=await download();assert.ok(withAxes.includes('>x</text>'));
  controls.axes.checked=false;controls.hull.checked=false;controls.symmetryToggle.checked=false;const bare=await download();assert.ok(!bare.includes('rgb(58,126,119)')&&!bare.includes('C₃'),'hidden layers omitted');
  controls.hull.checked=true;controls.symmetryToggle.checked=true;controls.reset.handlers.click();assert.equal(await download(),initial,'reset export reproducible');
  if(out)fs.writeFileSync(path.join(out,name+'.svg'),initial);
  dispose();assert.ok(!controls.download.handlers.click,'download listener removed');
  // Continuous symmetry, multiple poles, and representative replacement.
  const poles=[{x:0,y:0,z:1,multiplicity:3},{x:0,y:0,z:-1,multiplicity:3}];controls.symmetry=win.MajoranaSymmetry.detect(poles);
  dispose=win.test.attachSphere(canvas,poles,controls);const polar=await download();assert.ok(polar.includes('C∞')&&polar.includes('>3</text>'));assert.ok(!polar.includes('rgb(58,126,119)'));dispose();
  // Count all stars, not just distinct positions: 3 north + 1 south => b_z=1/2.
  const weighted=[{x:0,y:0,z:1,multiplicity:3},{x:0,y:0,z:-1,multiplicity:1}],saved=JSON.stringify(weighted);
  controls.symmetry=win.MajoranaSymmetry.detect(weighted);dispose=win.test.attachSphere(canvas,weighted,controls);
  const weightedSVG=await download(),position=diamondCentre(weightedSVG);near(position[0],300);near(position[1],200-153.9*.5*Math.cos(.28));
  canvas.handlers.keydown({key:'ArrowRight',preventDefault(){}});const rotated=diamondCentre(await download());near(rotated[0],300+153.9*.5*Math.sin(.28)*Math.sin(.08));near(rotated[1],position[1]);
  assert.equal(JSON.stringify(weighted),saved,'marker computation preserves physical state');
  if(out)fs.writeFileSync(path.join(out,name+'-weighted.svg'),weightedSVG);dispose();
  console.log(name+': PASS — current orientation, visibility, multiplicities, XML escaping, vector-only output, download filename/MIME, cleanup and replacement.');
 }
})().catch(e=>{console.error(e);process.exitCode=1;});
