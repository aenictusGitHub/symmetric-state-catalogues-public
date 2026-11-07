/* Basis-independent projector Husimi and Hermitian-operator stellar viewer.
 * API: SubspaceViewer.open(catalogueSubspace), SubspaceViewer.close().
 * Numeric construction and signed symmetry verification: OperatorMath. */
(() => {
  'use strict';
  const esc=x=>String(x??'').replace(/[&<>"']/g,c=>({'&':'&amp;','<':'&lt;','>':'&gt;','"':'&quot;',"'":'&#39;'}[c]));
  const math=tex=>{try{return window.katex.renderToString(tex,{throwOnError:false,trust:false,strict:'ignore'});}catch{return esc(tex);}};
  const fmt=x=>!Number.isFinite(x)?'—':Math.abs(x)>0&&Math.abs(x)<1e-4?x.toExponential(3):Number(x.toPrecision(6)).toString();
  function exactWeight(subspace,rank,value){const state=window.OperatorExactWeights?.states?.[subspace.id],entry=state?.weights?.[rank];return state?.n===subspace.n&&state?.dimension===subspace.dimension&&state?.normalization==='rho=P_V/d'&&entry?.verified===true&&typeof entry.latex==='string'&&Number.isFinite(entry.numericValue)&&Math.abs(entry.numericValue-value)<=1e-8*Math.max(1,Math.abs(value))?entry:null;}
  const weightTex=(subspace,rank,value)=>exactWeight(subspace,rank,value)?.latex||fmt(value);
  const dot=(a,b)=>a.reduce((s,x,i)=>s+x*b[i],0),cross=(a,b)=>[a[1]*b[2]-a[2]*b[1],a[2]*b[0]-a[0]*b[2],a[0]*b[1]-a[1]*b[0]];
  const unit=a=>a.map(x=>x/Math.hypot(...a));
  const ident=()=>[[1,0,0],[0,1,0],[0,0,1]];
  const apply=(R,a)=>R.map(row=>dot(row,a));
  const mul=(A,B)=>A.map(row=>[0,1,2].map(j=>row.reduce((s,x,i)=>s+x*B[i][j],0)));
  function rotation(axis,angle){const [x,y,z]=unit(axis),c=Math.cos(angle),s=Math.sin(angle),v=1-c;return [[c+x*x*v,x*y*v-z*s,x*z*v+y*s],[y*x*v+z*s,c+y*y*v,y*z*v-x*s],[z*x*v-y*s,z*y*v+x*s,c+z*z*v]];}
  const axisColor=order=>({2:'#34786d',3:'#82609b',4:'#a45635',5:'#b38b23',continuous:'#356f94'}[order]||'#526b9b');
  const rankColor=rank=>['#a87038','#34786d','#82609b','#52799a','#a45668','#5c6850'][rank%6];
  const mixColor=(hex,target,t)=>'#'+[0,1,2].map(i=>{const a=parseInt(hex.slice(1+2*i,3+2*i),16),b=parseInt(target.slice(1+2*i,3+2*i),16);return Math.round(a+(b-a)*t).toString(16).padStart(2,'0');}).join('');
  const subscript=n=>String(n).replace(/\d/g,d=>'₀₁₂₃₄₅₆₇₈₉'[Number(d)]);
  let dialog=null,returnFocus=null,cleanup=()=>{},current=null,serial=0,palettePreference='viridis';
  const representationKey='ac-subspace-representation',representations=['husimi','wigner','operator','spheres'];
  let representationPreference='husimi';
  try{const saved=window.localStorage?.getItem(representationKey);if(representations.includes(saved))representationPreference=saved;}catch{}
  function rememberRepresentation(value){representationPreference=value;try{window.localStorage?.setItem(representationKey,value);}catch{}}
  // Keep only a few sampled meshes: catalogue records live for the whole page.
  const cache=new Map(),wignerMeshes=new WeakMap();
  const $=id=>dialog.querySelector('#sv-'+id);
  function ensure(){
    if(dialog)return;
    dialog=document.createElement('dialog');dialog.id='subspace-viewer';dialog.className='sv-dialog';dialog.setAttribute('aria-labelledby','sv-title');dialog.setAttribute('aria-describedby','sv-definition');
    dialog.innerHTML=`<div class="sv-heading"><div><p class="sv-eyebrow">Whole subspace</p><h2 id="sv-title">Husimi, Wigner & operator constellations</h2></div><button type="button" id="sv-close" aria-label="Close whole-subspace viewer">Close ×</button></div>
      <p id="sv-label" class="sv-label"></p><p id="sv-definition" class="sv-definition"></p>
      <div class="sv-tabs" role="group" aria-label="Whole-subspace representation"><button type="button" id="sv-husimi" aria-pressed="true">Husimi function</button><button type="button" id="sv-wigner" aria-pressed="false">Wigner function</button><button type="button" id="sv-operator" aria-pressed="false">Operator Majorana</button><button type="button" id="sv-spheres" aria-pressed="false" title="Show all clickable multipole spheres at their relative radii">Multipole spheres</button></div>
      <div id="sv-gallery-axis-legend" class="sv-gallery-axis-legend" aria-label="Symmetry axes shown on all multipole spheres" hidden></div>
      <div class="sv-workspace"><div class="sv-frame"><canvas id="sv-canvas" class="sv-canvas" tabindex="0" role="img" aria-label="Whole-subspace Husimi sphere"></canvas><div id="sv-math-labels" class="sv-math-labels"></div><span class="sv-frame-help">Drag to rotate · Scroll to zoom · Arrow keys rotate</span><div id="sv-status" class="sv-status" role="status">Computing the subspace projector…</div></div>
      <aside><div class="sv-controls"><div class="sv-turn" role="group" aria-label="Rotate whole-subspace view"><button type="button" id="sv-left" aria-label="Rotate left">←</button><button type="button" id="sv-up" aria-label="Rotate up">↑</button><button type="button" id="sv-down" aria-label="Rotate down">↓</button><button type="button" id="sv-right" aria-label="Rotate right">→</button><button type="button" id="sv-reset">Reset view</button></div>
      <label class="sv-mode-control" id="sv-shape-control">Display <select id="sv-shape"><option value="sphere">Colored sphere</option><option value="radial">Radial surface</option></select></label>
      <label class="sv-mode-control" id="sv-palette-control">Palette <select id="sv-palette"><option value="viridis">Violet &amp; teal</option><option value="plum">Plum &amp; gold</option><option value="ocean">Ocean &amp; coral</option></select></label>
      <label id="sv-contour-control"><input type="checkbox" id="sv-contours" checked> Contour lines</label>
      <label class="sv-mode-control" id="sv-separate-control" hidden><input type="checkbox" id="sv-separate">Separate all spheres</label>
      <label class="sv-mode-control" id="sv-rank-control" hidden>Multipole <select id="sv-rank" aria-label="Operator multipole rank"></select></label>
      <label class="sv-mode-control" id="sv-fit-control" hidden><input id="sv-fit" type="checkbox" checked> Fit selected sphere</label>
      <label><input id="sv-axes" type="checkbox" checked> ${math('C_n')} symmetry axes</label><label>Show <select id="sv-order"><option value="highest">Highest order</option><option value="all">All axes</option></select></label>
      <button id="sv-align" type="button">Align highest ${math('C_n')} with ${math('z')}</button>
      <label><input type="checkbox" id="sv-coordinate">Coordinate axes</label>
      </div><div id="sv-readout" class="sv-readout" aria-live="polite"></div><div id="sv-symmetry" class="sv-readout"></div><div class="sv-downloads"><button id="sv-svg" type="button">Export SVG ↓</button><button id="sv-json" type="button">Download operator ↓</button></div></aside></div>
      <details class="sv-method"><summary>What these pictures represent</summary><div id="sv-method"></div></details>
      <details class="sv-method"><summary>Multipoles & numerical verification</summary><div id="sv-verification"></div></details>`;
    document.body.append(dialog);$('close').addEventListener('click',close);dialog.addEventListener('close',()=>finish());
    dialog.addEventListener('click',e=>{if(e.target!==dialog)return;const b=dialog.getBoundingClientRect();if(e.clientX<b.left||e.clientX>b.right||e.clientY<b.top||e.clientY>b.bottom)close();});
  }
  function finish(){window.SubspaceSphereDetail?.close();serial++;cleanup();cleanup=()=>{};document.body.classList.remove('sv-open');if(returnFocus?.isConnected)returnFocus.focus({preventScroll:true});returnFocus=null;current=null;}
  function close(){window.SubspaceSphereDetail?.close();if(dialog?.open)dialog.close();else finish();}
  function message(text){$('status').textContent=text;$('status').hidden=!text;}
  function buildMesh(model){
    const rows=128,cols=256,vertices=[],faces=[];let qMin=Infinity,qMax=0;
    for(let i=0;i<=rows;i++){const t=Math.PI*i/rows,s=Math.sin(t),z=Math.cos(t);for(let j=0;j<=cols;j++){const p=2*Math.PI*j/cols,n=[s*Math.cos(p),s*Math.sin(p),z],q=window.OperatorMath.husimi(model,n);vertices.push({n,q});qMin=Math.min(qMin,q);qMax=Math.max(qMax,q);}}
    for(let i=0;i<rows;i++)for(let j=0;j<cols;j++){const a=i*(cols+1)+j,b=a+cols+1;faces.push([a,b,b+1,a+1]);}
    const mesh={vertices,faces,qMin,qMax,rows,cols};mesh.contours=window.HusimiRendering.contours(mesh);return mesh;
  }
  function buildWignerMesh(model){
    if(wignerMeshes.has(model))return wignerMeshes.get(model);
    const mesh=window.SpinWigner.prepare(model).grid(128,256);mesh.colorMax=Math.max(Math.abs(mesh.qMin),Math.abs(mesh.qMax),1e-15);
    const levels=window.HusimiRendering.defaultLevels(mesh);
    if(mesh.qMin<0&&mesh.qMax>0&&!levels.some(x=>Math.abs(x)<1e-14))levels.push(0);
    mesh.contours=window.HusimiRendering.contours(mesh,levels);wignerMeshes.set(model,mesh);return mesh;
  }
  function direction(p){return Array.isArray(p)?p.slice(0,3):[p.x,p.y,p.z];}
  function makeRenderer(model,mesh,symmetry,ui=$,onOpenRank=null){
    const canvas=ui('canvas'),ctx=canvas.getContext('2d'),texture=document.createElement('canvas'),textureContext=texture.getContext('2d');let R=mul(rotation([1,0,0],-.25),rotation([0,1,0],-.45)),alignment=ident(),zoom=1,mode='husimi',rank=1,record=null,drag=null,frame=0,hitCells=[];
    const maxWeight=Math.max(...model.multipoles.filter(m=>m.rank>0).map(m=>m.weight),1e-15);
    function redraw(){if(frame)return;frame=requestAnimationFrame(()=>{frame=0;render();});}
    const displayedAxes=()=>{if(!ui('axes').checked)return[];const axes=symmetry.axes||[],v=ui('order').value;if(v==='all')return axes;if(v==='highest'){const highest=Math.max(...axes.map(a=>a.order==null?Infinity:a.order));return axes.filter(a=>(a.order==null?Infinity:a.order)===highest);}return axes.filter(a=>String(a.order==null?'continuous':a.order)===v);};
    const colorBounds=()=>mode==='wigner'?[-mesh.colorMax,mesh.colorMax]:[mesh.qMin,mesh.qMax];
    const scalarPalette=()=>mode==='wigner'?'wigner':palettePreference;
    const scalarColor=q=>{const [lo,hi]=colorBounds();return window.HusimiRendering.color((q-lo)/Math.max(hi-lo,1e-14),scalarPalette());};
    function contourColor(level){
      const luminance=hex=>{const c=[1,3,5].map(i=>parseInt(hex.slice(i,i+2),16)/255).map(v=>v<=.04045?v/12.92:((v+.055)/1.055)**2.4);return .2126*c[0]+.7152*c[1]+.0722*c[2];};
      const background=luminance(scalarColor(level)),dark='#172a39',light='#cbd4d9';
      const contrast=ink=>{const l=luminance(ink);return (Math.max(background,l)+.05)/(Math.min(background,l)+.05);};
      return contrast(dark)>=contrast(light)?dark:light;
    }
    function heatmap(width,height,cx,cy,scale){
      const x0=Math.max(0,Math.floor(cx-scale)),y0=Math.max(0,Math.floor(cy-scale)),x1=Math.min(Math.ceil(width),Math.ceil(cx+scale)),y1=Math.min(Math.ceil(height),Math.ceil(cy+scale)),w=Math.max(1,x1-x0),h=Math.max(1,y1-y0),image=textureContext.createImageData(w,h),pixels=image.data,M=mul(R,alignment),bounds=colorBounds(),range=Math.max(bounds[1]-bounds[0],1e-14),stride=mesh.cols+1;
      for(let y=0;y<h;y++)for(let x=0;x<w;x++){
        const vx=(x+x0+.5-cx)/scale,vz=(cy-y-y0-.5)/scale,r2=vx*vx+vz*vz;if(r2>1)continue;const vy=Math.sqrt(Math.max(0,1-r2));
        const nx=M[0][0]*vx+M[1][0]*vy+M[2][0]*vz,ny=M[0][1]*vx+M[1][1]*vy+M[2][1]*vz,nz=M[0][2]*vx+M[1][2]*vy+M[2][2]*vz;
        const row=Math.acos(Math.max(-1,Math.min(1,nz)))*mesh.rows/Math.PI,col=((Math.atan2(ny,nx)/(2*Math.PI)+1)%1)*mesh.cols,i=Math.min(mesh.rows-1,Math.floor(row)),j=Math.floor(col),a=row-i,b=col-j,k=i*stride+j;
        const q=(1-a)*((1-b)*mesh.vertices[k].q+b*mesh.vertices[k+1].q)+a*((1-b)*mesh.vertices[k+stride].q+b*mesh.vertices[k+stride+1].q),t=Math.max(0,Math.min(1,(q-bounds[0])/range)),o=4*(y*w+x);
        const color=window.HusimiRendering.rgb(t,scalarPalette());
        pixels[o]=color[0];pixels[o+1]=color[1];pixels[o+2]=color[2];pixels[o+3]=255*Math.min(1,(1-Math.sqrt(r2))*scale);
      }
      texture.width=w;texture.height=h;textureContext.putImageData(image,0,0);return {type:'bitmap',canvas:texture,x:x0,y:y0};
    }
    function scene(width,height,vectorExport=false,cellView=null){
      const displayRecord=cellView?.record??record,displayRank=cellView?.rank??rank;
      if(!cellView&&mode==='operator'&&displayRank==='all'&&ui('separate').checked){
        const layout=window.SubspaceSphereLayout.layout(width,height,displayRecord?.records||[],maxWeight),items=[];hitCells=layout.cells;
        const shift=(p,x,y)=>p?{...p,x:p.x+x,y:p.y+y}:p;
        for(const cell of layout.cells){
          const rec=displayRecord.records.find(r=>r.rank===cell.rank);
          items.push({type:'clip-start',x:cell.x,y:cell.y,width:cell.width,height:cell.height});
          for(const o of scene(cell.width,cell.height,vectorExport,{record:rec,rank:cell.rank,scale:layout.maxDisplayRadius*zoom,cx:cell.cx-cell.x,cy:cell.cy-cell.y,separate:true})){
            items.push({...o,...(o.p?{p:shift(o.p,cell.x,cell.y)}:{}),...(o.a?{a:shift(o.a,cell.x,cell.y),b:shift(o.b,cell.x,cell.y)}:{}),...(o.points?{points:o.points.map(p=>shift(p,cell.x,cell.y))}:{}),...(o.cx!=null?{cx:o.cx+cell.x,cy:o.cy+cell.y}:{})});
          }
          const weightTex='w_{'+cell.rank+'}\\simeq '+fmt(cell.weight);
          items.push({type:'text',p:{x:cell.labelArea.x+cell.labelArea.width/2,y:cell.labelArea.baseline},content:weightTex,tex:weightTex,color:'#61566a',size:12,offset:[0,0],textAnchor:'middle'});
          items.push({type:'clip-end'});
        }
        return items;
      }
      if(!cellView)hitCells=[];
      const cx=cellView?.cx??width/2,cy=cellView?.cy??height/2,scale=cellView?.scale??Math.min(width*.39,height*.40)*zoom,objects=[],markerLayer=[],labelLayer=[];
      const projectReference=a=>{const v=apply(R,a);return {x:cx+v[0]*scale,y:cy-v[2]*scale,z:v[1]};};
      const project=a=>projectReference(apply(alignment,a));
      const line=(a,b,color,width=1,dash=[],opacity=1)=>objects.push({type:'line',a:project(a),b:project(b),color,width,dash,opacity});
      const circle=(p,r,fill,stroke,opacity=1)=>markerLayer.push({type:'circle',p:project(p),r,fill,stroke,opacity});
      const text=(p,content,color='#3f4241',size=12,offset=[6,-6])=>objects.push({type:'text',p:project(p),content,color,size,offset});
      if(mode==='husimi'||mode==='wigner'){
        const radial=ui('shape').value==='radial';let surfaceDepth=null;
        if(!radial&&!vectorExport)objects.push(heatmap(width,height,cx,cy,scale));
        else{const max=(mode==='wigner'?mesh.colorMax:mesh.qMax)||1,magnitude=q=>mode==='wigner'?Math.abs(q):q,verts=mesh.vertices.map(v=>({p:project(v.n.map(x=>x*(radial?magnitude(v.q)/max:1))),q:v.q}));
          if(radial&&ui('contours').checked)surfaceDepth=window.SubspaceRenderGeometry.surfaceDepth(verts,mesh.faces,width,height);
          const faces=mesh.faces.map(ids=>({ids,z:ids.reduce((s,i)=>s+verts[i].p.z,0)/ids.length})).sort((a,b)=>a.z-b.z);
          for(const f of faces){const vs=f.ids.map(i=>verts[i]),avg=vs.reduce((s,v)=>s+v.q,0)/4;objects.push({type:'poly',points:vs.map(v=>v.p),fill:scalarColor(avg),opacity:1});}
        }
        if(ui('contours').checked){
          const radiusFor=level=>radial?(mode==='wigner'?Math.abs(level)/mesh.colorMax:level/(mesh.qMax||1)):1,visibleContours=[];
          for(const segment of mesh.contours.segments){
            let a=segment.a.map(x=>x*radiusFor(segment.level)),b=segment.b.map(x=>x*radiusFor(segment.level));
            if(radial&&surfaceDepth){
              for(const [t0,t1]of window.SubspaceRenderGeometry.visibleIntervals(project(a),project(b),surfaceDepth))visibleContours.push([a.map((x,k)=>x+t0*(b[k]-x)),a.map((x,k)=>x+t1*(b[k]-x)),segment.level]);
            }else{
              const za=project(a).z,zb=project(b).z;
              if(za<0&&zb<0)continue;
              if(za<0||zb<0){const t=za/(za-zb),cut=a.map((x,k)=>x+t*(b[k]-x));if(za<0)a=cut;else b=cut;}
              visibleContours.push([a,b,segment.level]);
            }
          }
          // One solid stroke per contour; the level fixes its background color.
          const contourColors=new Map();
          for(const [a,b,level] of visibleContours){
            if(!contourColors.has(level))contourColors.set(level,contourColor(level));
            line(a,b,contourColors.get(level),mode==='wigner'&&Math.abs(level)<1e-14?1.35:1,[],1);
          }
        }
      }else{
        const records=displayRecord?.all?displayRecord.records:displayRecord?[displayRecord]:[],radiusOf=r=>!cellView?.separate&&ui('fit').checked&&!displayRecord?.all?1:r.weight/maxWeight,active=records.filter(r=>r.available).map(r=>({record:r,radius:radiusOf(r),color:displayRecord?.all||cellView?.separate?rankColor(r.rank):'#78518f'})).sort((a,b)=>b.radius-a.radius),outerRadius=active[0]?.radius||0;
        // The dotted circle is only a scale reference. Active multipoles alone
        // receive a shaded spherical shell and their true norm radius is kept.
        if(outerRadius<1-1e-8)objects.push({type:'ring',cx,cy,r:scale,color:'#b8b6ad',width:.75,dash:[1,5],opacity:.5});
        for(const shell of active){const edge=mixColor(shell.color,'#49434b',.25),soft=mixColor(shell.color,'#fffaf0',.68);
          objects.push({type:'shell',cx,cy,r:scale*shell.radius,opacity:displayRecord?.all ? .42 : 1,stops:[[0,'#ffffff',.05],[.38,'#ffffff',.025],[.72,soft,.07],[.93,shell.color,.16],[1,edge,.28]]});
          objects.push({type:'ring',cx,cy,r:scale*shell.radius,color:displayRecord?.all?shell.color:edge,width:displayRecord?.all ? .65 : .95,opacity:displayRecord?.all ? .55 : .7});
        }
        if(outerRadius>1e-12){
          const gridColor=mixColor(active[0].color,'#8b8790',.52),curves=[];
          for(const latitude of [-Math.PI/3,-Math.PI/6,0,Math.PI/6,Math.PI/3])curves.push(Array.from({length:181},(_,i)=>{const a=2*Math.PI*i/180;return [outerRadius*Math.cos(latitude)*Math.cos(a),outerRadius*Math.cos(latitude)*Math.sin(a),outerRadius*Math.sin(latitude)];}));
          for(let longitude=0;longitude<6;longitude++){const p=longitude*Math.PI/6;curves.push(Array.from({length:181},(_,i)=>{const a=2*Math.PI*i/180;return [outerRadius*Math.cos(a)*Math.cos(p),outerRadius*Math.cos(a)*Math.sin(p),outerRadius*Math.sin(a)];}));}
          // Continuous paths keep the dash spacing stable around each circle.
          for(const front of [false,true])for(const curve of curves){let run=[];const flush=()=>{if(run.length>1)objects.push({type:'curve',points:run,color:gridColor,width:front?.8:.65,dash:front?[]:[2.2,4],opacity:front?(displayRecord?.all ? .35 : .53):(displayRecord?.all ? .13 : .22)});run=[];};for(let i=1;i<curve.length;i++){const a=projectReference(curve[i-1]),b=projectReference(curve[i]);if(((a.z+b.z)>=0)===front){if(!run.length)run.push(a);run.push(b);}else flush();}flush();}
        }
        const points=[];
        for(const r of records.filter(x=>x.available)){let chosen=r.selectedPoints||[],other=r.oppositePoints||[];if(!chosen.length&&!other.length)chosen=r.points||[];const radius=radiusOf(r),color=displayRecord?.all||cellView?.separate?rankColor(r.rank):'#78518f';points.push(...[...chosen.map(p=>({p,chosen:true})),...other.map(p=>({p,chosen:false}))].map(x=>({...x,radius,color,depth:project(direction(x.p).map(v=>v*radius)).z})));}
        points.sort((a,b)=>a.depth-b.depth);
        for(const {p,chosen,depth,radius,color} of points){const n=direction(p).map(v=>v*radius),u=Math.max(0,Math.min(1,(depth/Math.max(radius,1e-15)+.12)/.24)),front=u*u*(3-2*u),base=chosen?color:'#b0733d',fill=mixColor(base,'#fff9ee',.62*(1-front)),stroke=mixColor(base,'#44372c',.3*front);circle(n,cellView?.separate?Math.min(4.2,Math.max(1,scale*radius*.11)):5.2,fill,stroke,1);if(p.multiplicity>1)markerLayer.push({type:'text',p:project(n),content:String(p.multiplicity),color:'#3d3445',size:11,offset:[7,-8]});}
      }
      for(const a of displayedAxes()){
        const color=axisColor(a.order==null?'continuous':a.order),guide=window.SubspaceRenderGeometry.axisGuide(a.axis||a.direction,project);
        for(const half of guide.halves){
          if(half.dashedOpacity>1e-4)line(half.a,half.b,color,1.15,[5,5],half.dashedOpacity);
          if(half.solidOpacity>1e-4)line(half.a,half.b,color,1.25,[],half.solidOpacity);
        }
        if(guide.label.opacity>1e-4&&!cellView?.separate){
          const p=project(guide.label.point),[dx,dy]=guide.label.radialOffset;
          p.x=Math.max(22,Math.min(width-22,p.x+dx));p.y=Math.max(15,Math.min(height-15,p.y+dy));
          labelLayer.push({type:'text',p,content:'C'+(a.order==null?'∞':subscript(a.order)),tex:'C_{'+(a.order==null?'\\infty':a.order)+'}',color,size:17,offset:[0,0],opacity:guide.label.opacity,textAnchor:guide.label.textAnchor});
        }
      }
      if(ui('coordinate').checked)for(const [v,label]of [[[1,0,0],'x'],[[0,1,0],'y'],[[0,0,1],'z']]){objects.push({type:'line',a:projectReference([0,0,0]),b:projectReference(v.map(x=>x*1.2)),color:'#8e8f89',width:.8,dash:[],opacity:.7});const p=projectReference(v.map(x=>x*1.24));p.x=Math.max(8,Math.min(width-(label==='z'?40:24),p.x));p.y=Math.max(25,Math.min(height-18,p.y));objects.push({type:'text',p,content:label,tex:label,color:'#62655f',size:14,offset:[label==='z'?24:6,-6]});}
      if(mode==='operator'&&(!displayRecord?.available&&!displayRecord?.all||displayRecord?.zero))objects.push({type:'message',cx:width/2,cy:height/2,maxWidth:width-30,text:displayRecord?.zero?'This multipole vanishes.':displayRecord?.error||'Choose a nonzero multipole rank.'});
      if(mode==='operator'){
        const values=displayRank==='all'?model.multipoles.filter(m=>m.rank>0&&!m.zero).map(m=>m.rank):[displayRank];
        const step=values.length>1?values[1]-values[0]:0,regular=values.length>6&&values.every((v,i)=>!i||v-values[i-1]===step);
        const members=regular?values.slice(0,2).join(',')+',\\ldots,'+values.at(-1):values.join(',');
        const tex=displayRank==='all'?'\\ell\\in\\{'+members+'\\}':'\\ell='+displayRank;
        labelLayer.push({type:'text',p:{x:16,y:28},content:tex,tex,color:'#493a56',size:Math.min(17,Math.max(10,(width-36)/Math.max(1,members.length)*1.5)),offset:[0,0],badge:true,detailRank:cellView?.separate?displayRank:null});
      }
      return objects.concat(markerLayer,labelLayer);
    }
    function render(){
      const separate=mode==='operator'&&rank==='all'&&ui('separate').checked,frameElement=canvas.parentElement;
      frameElement.classList.toggle('sv-separated',separate);
      const visibleWidth=canvas.getBoundingClientRect().width,columns=Math.max(1,Math.min(4,Math.round(visibleWidth/180)));
      const height=separate?Math.max(310,Math.ceil((record?.records?.length||1)/columns)*190)+'px':'';
      if(canvas.style.height!==height)canvas.style.height=height;
      const b=canvas.getBoundingClientRect();ui('math-labels').style.height=b.height+'px';
      frameElement.querySelector('.sv-frame-help').textContent=separate?'Drag to rotate · Scroll through multipoles · Ctrl+scroll to zoom':'Drag to rotate · Scroll to zoom · Arrow keys rotate';if(!b.width||!b.height)return;const dpr=Math.min(window.devicePixelRatio||1,2);canvas.width=Math.round(b.width*dpr);canvas.height=Math.round(b.height*dpr);ctx.setTransform(dpr,0,0,dpr,0,0);ctx.lineCap='round';ctx.clearRect(0,0,b.width,b.height);const objects=scene(b.width,b.height);paintCanvas(objects,ctx,b.width,b.height);paintMathLabels(objects);}
    const mathMarkup=new Map();
    function paintMathLabels(objects){
      const layer=ui('math-labels');layer.replaceChildren();const hosts=[{element:layer,x:0,y:0}];
      for(const o of objects){
        const host=hosts.at(-1);
        if(o.type==='clip-start'){
          const cell=document.createElement('div');cell.className='sv-math-cell';cell.style.cssText=`left:${o.x-host.x}px;top:${o.y-host.y}px;width:${o.width}px;height:${o.height}px`;host.element.append(cell);hosts.push({element:cell,x:o.x,y:o.y});continue;
        }
        if(o.type==='clip-end'){hosts.pop();continue;}
        if(!o.tex)continue;
        const key=o.tex;let html=mathMarkup.get(key);if(!html){html=math(key);mathMarkup.set(key,html);}
        const label=document.createElement(o.detailRank!=null&&typeof onOpenRank==='function'?'button':'span');label.className='sv-math-label'+(o.badge?' sv-rank-badge':'');label.innerHTML=html;
        if(o.detailRank!=null&&typeof onOpenRank==='function'){label.type='button';label.classList.add('sv-sphere-open');label.setAttribute('aria-label','Open enlarged multipole '+o.detailRank);label.addEventListener('click',()=>onOpenRank(o.detailRank));}else label.setAttribute('aria-hidden','true');
        label.style.cssText=`left:${o.p.x+o.offset[0]-host.x}px;top:${o.p.y+o.offset[1]-host.y}px;font-size:${o.size}px;color:${o.color};opacity:${o.opacity??1}`;
        // Match scene text baseline; badge padding remains outside the formula.
        label.style.transform=o.textAnchor==='middle'?'translate(-50%,-50%)':'translateY(-80%)';host.element.append(label);
      }
    }
    function rotate(dx,dy){R=mul(rotation([0,0,1],dx),mul(rotation([1,0,0],dy),R));redraw();}
    function reset(){zoom=1;if((symmetry.axes||[]).length){align();}else{R=mul(rotation([1,0,0],-.25),rotation([0,1,0],-.45));alignment=ident();redraw();}}
    function align(){const axes=symmetry.axes||[],a=axes.slice().sort((x,y)=>(y.order??Infinity)-(x.order??Infinity))[0];if(!a)return;const v=unit(a.axis||a.direction),target=[0,0,1],c=Math.max(-1,Math.min(1,dot(v,target)));alignment=c>1-1e-12?ident():c< -1+1e-12?rotation([1,0,0],Math.PI):rotation(cross(v,target),Math.acos(c));R=rotation([1,0,0],-.20);redraw();}
    const onDown=e=>{drag={x:e.clientX,y:e.clientY,startX:e.clientX,startY:e.clientY,moved:false};canvas.setPointerCapture?.(e.pointerId);canvas.classList.add('dragging');};
    const onMove=e=>{if(!drag)return;drag.moved||=Math.hypot(e.clientX-drag.startX,e.clientY-drag.startY)>4;if(drag.moved)rotate((e.clientX-drag.x)*.008,(e.clientY-drag.y)*.008);drag.x=e.clientX;drag.y=e.clientY;};
    const onUp=e=>{const clicked=drag&&!drag.moved&&e?.type!=='pointercancel';drag=null;canvas.classList.remove('dragging');if(clicked&&typeof onOpenRank==='function'&&mode==='operator'&&rank==='all'&&ui('separate').checked){const b=canvas.getBoundingClientRect(),x=e.clientX-b.left,y=e.clientY-b.top,cell=hitCells.find(c=>x>=c.x&&x<=c.x+c.width&&y>=c.y&&y<=c.y+c.height);if(cell)onOpenRank(cell.rank);}};
    const onWheel=e=>{if(mode==='operator'&&rank==='all'&&ui('separate').checked&&!e.ctrlKey&&!e.metaKey)return;e.preventDefault();zoom=Math.max(.35,Math.min(5,zoom*Math.exp(-e.deltaY*.001)));redraw();};
    const onKey=e=>{const actions={ArrowLeft:()=>rotate(-.14,0),ArrowRight:()=>rotate(.14,0),ArrowUp:()=>rotate(0,-.14),ArrowDown:()=>rotate(0,.14),Home:reset,'+':()=>{zoom=Math.min(5,zoom*1.15);redraw();},'-':()=>{zoom=Math.max(.35,zoom/1.15);redraw();}};if(actions[e.key]){e.preventDefault();actions[e.key]();}};
    canvas.addEventListener('pointerdown',onDown);canvas.addEventListener('pointermove',onMove);canvas.addEventListener('pointerup',onUp);canvas.addEventListener('pointercancel',onUp);canvas.addEventListener('wheel',onWheel,{passive:false});canvas.addEventListener('keydown',onKey);
    const handlers={left:()=>rotate(-.14,0),right:()=>rotate(.14,0),up:()=>rotate(0,-.14),down:()=>rotate(0,.14),reset,align,shape:redraw,fit:redraw,palette:()=>{palettePreference=ui('palette').value;redraw();},contours:redraw,axes:redraw,order:redraw,coordinate:redraw};
    for(const [id,fn]of Object.entries(handlers))ui(id).addEventListener(['shape','fit','palette','contours','axes','order','coordinate'].includes(id)?'change':'click',fn);
    const observer=new ResizeObserver(redraw);observer.observe(canvas);
    return {align,set(next,r,c,nextMesh){mode=next;rank=r;record=c;if(nextMesh)mesh=nextMesh;canvas.setAttribute('aria-label',mode==='husimi'?'Whole-subspace Husimi sphere':mode==='wigner'?'Whole-subspace Wigner sphere':(rank==='all'?(ui('separate').checked?'All nonzero ranks, separate spheres':'All nonzero ranks'): 'Rank '+rank)+' operator Majorana constellation');redraw();},exportSvg(){const b=canvas.getBoundingClientRect(),description=mode==='operator'?(record?.all?'Multipole sphere radii are proportional to their norms.':ui('fit').checked?'Selected sphere fitted to unit display radius for inspecting star directions; true multipole norm '+(record?.weight??'')+'.':'Selected sphere radius is proportional to its multipole norm.'):mode==='wigner'?'Unit-integral signed spin Wigner function with a zero-centered blue-white-orange scale; '+(ui('contours').checked?'signed contours.':'no contours.'):'Husimi function with '+palettePreference+' palette and '+(ui('contours').checked?'contours.':'no contours.');return sceneSvg(scene(b.width,b.height,true),b.width,b.height).replace('</title>','</title><desc>'+esc(description)+'</desc>');},dispose(){observer.disconnect();cancelAnimationFrame(frame);for(const [id,fn]of Object.entries(handlers))ui(id).removeEventListener(['shape','fit','palette','contours','axes','order','coordinate'].includes(id)?'change':'click',fn);canvas.removeEventListener('pointerdown',onDown);canvas.removeEventListener('pointermove',onMove);canvas.removeEventListener('pointerup',onUp);canvas.removeEventListener('pointercancel',onUp);canvas.removeEventListener('wheel',onWheel);canvas.removeEventListener('keydown',onKey);},redraw,getDisplayedAxes:displayedAxes};
  }
  function paintCanvas(objects,ctx,width,height){
    for(const o of objects){
      ctx.globalAlpha=o.opacity??1;ctx.setLineDash(o.dash||[]);
      if(o.type==='clip-start'){ctx.save();ctx.beginPath();ctx.rect(o.x,o.y,o.width,o.height);ctx.clip();}
      else if(o.type==='clip-end')ctx.restore();
      else if(o.type==='bitmap')ctx.drawImage(o.canvas,o.x,o.y);
      else if(o.type==='poly'){ctx.beginPath();o.points.forEach((p,i)=>i?ctx.lineTo(p.x,p.y):ctx.moveTo(p.x,p.y));ctx.closePath();ctx.fillStyle=o.fill;ctx.fill();ctx.strokeStyle=o.fill;ctx.lineWidth=.7;ctx.stroke();}
      else if(o.type==='shell'){
        const gradient=ctx.createRadialGradient(o.cx-o.r*.36,o.cy-o.r*.42,0,o.cx,o.cy,o.r);
        for(const [offset,color,alpha]of o.stops){const rgb=[1,3,5].map(i=>parseInt(color.slice(i,i+2),16));gradient.addColorStop(offset,`rgba(${rgb.join(',')},${alpha})`);}
        ctx.beginPath();ctx.arc(o.cx,o.cy,o.r,0,2*Math.PI);ctx.fillStyle=gradient;ctx.fill();
      }
      else if(o.type==='curve'){ctx.beginPath();o.points.forEach((p,i)=>i?ctx.lineTo(p.x,p.y):ctx.moveTo(p.x,p.y));ctx.strokeStyle=o.color;ctx.lineWidth=o.width;ctx.stroke();}
      else if(o.type==='line'){ctx.beginPath();ctx.moveTo(o.a.x,o.a.y);ctx.lineTo(o.b.x,o.b.y);ctx.strokeStyle=o.color;ctx.lineWidth=o.width;ctx.stroke();}
      else if(o.type==='circle'){ctx.beginPath();ctx.arc(o.p.x,o.p.y,o.r,0,2*Math.PI);ctx.fillStyle=o.fill;ctx.fill();ctx.strokeStyle=o.stroke;ctx.lineWidth=1;ctx.stroke();}
      else if(o.type==='ring'){ctx.beginPath();ctx.arc(o.cx,o.cy,o.r,0,2*Math.PI);ctx.strokeStyle=o.color;ctx.lineWidth=o.width;ctx.stroke();}
      else if(o.type==='text'&&!o.tex){ctx.fillStyle=o.color;ctx.textAlign=o.textAnchor==='middle'?'center':'left';ctx.textBaseline=o.textAnchor==='middle'?'middle':'alphabetic';ctx.font=`${o.size}px "Latin Modern Roman", Georgia, serif`;ctx.fillText(o.content,o.p.x+o.offset[0],o.p.y+o.offset[1]);}
      else if(o.type==='message'){ctx.fillStyle='#6a5d70';ctx.font='13px "Latin Modern Roman", Georgia, serif';ctx.textAlign='center';const words=o.text.split(' ');let line='',row=0;for(const word of words){if(ctx.measureText(line+' '+word).width>(o.maxWidth??width-65)&&line){ctx.fillText(line,o.cx??width/2,(o.cy??height/2)+row*20);line=word;row++;}else line+=(line?' ':'')+word;}ctx.fillText(line,o.cx??width/2,(o.cy??height/2)+row*20);ctx.textAlign='left';}
    }
    ctx.globalAlpha=1;ctx.textAlign='left';ctx.textBaseline='alphabetic';ctx.setLineDash([]);
  }
  function sceneSvg(objects,width,height){
    const n=x=>Number(x.toFixed(3)),gradientIds=new Map(),prefix='sv-'+Date.now().toString(36)+'-'+(sceneSvg.sequence=(sceneSvg.sequence||0)+1)+'-shell-';let nextId=0;
    const clipIds=new Map();let clipIndex=0;
    const clips=objects.filter(o=>o.type==='clip-start').map(o=>{const id=prefix+'clip-'+clipIndex++;clipIds.set(o,id);return `<clipPath id="${id}"><rect x="${n(o.x)}" y="${n(o.y)}" width="${n(o.width)}" height="${n(o.height)}"/></clipPath>`;}).join('');
    const defs=clips+objects.filter(o=>o.type==='shell').map(o=>{const id=prefix+nextId++;gradientIds.set(o,id);return `<radialGradient id="${id}" gradientUnits="userSpaceOnUse" cx="${n(o.cx)}" cy="${n(o.cy)}" r="${n(o.r)}" fx="${n(o.cx-o.r*.36)}" fy="${n(o.cy-o.r*.42)}">${o.stops.map(([offset,color,alpha])=>`<stop offset="${offset}" stop-color="${color}" stop-opacity="${alpha}"/>`).join('')}</radialGradient>`;}).join('');
    return `<svg xmlns="http://www.w3.org/2000/svg" width="${n(width)}" height="${n(height)}" viewBox="0 0 ${n(width)} ${n(height)}" stroke-linecap="round"><title>Whole-subspace visualization</title><defs>${defs}</defs><rect width="100%" height="100%" fill="#faf8f1"/>`+objects.map(o=>{
      const style=`opacity="${o.opacity??1}"`;
      if(o.type==='clip-start')return `<g clip-path="url(#${clipIds.get(o)})">`;
      if(o.type==='clip-end')return '</g>';
      if(o.type==='poly')return `<polygon points="${o.points.map(p=>n(p.x)+','+n(p.y)).join(' ')}" fill="${o.fill}" stroke="${o.fill}" stroke-width=".7" ${style}/>`;
      if(o.type==='shell')return `<circle cx="${n(o.cx)}" cy="${n(o.cy)}" r="${n(o.r)}" fill="url(#${gradientIds.get(o)})" ${style}/>`;
      if(o.type==='curve')return `<polyline points="${o.points.map(p=>n(p.x)+','+n(p.y)).join(' ')}" fill="none" stroke="${o.color}" stroke-width="${o.width}" stroke-dasharray="${(o.dash||[]).join(' ')}" ${style}/>`;
      if(o.type==='line')return `<path d="M${n(o.a.x)} ${n(o.a.y)}L${n(o.b.x)} ${n(o.b.y)}" stroke="${o.color}" stroke-width="${o.width}" stroke-dasharray="${(o.dash||[]).join(' ')}" ${style}/>`;
      if(o.type==='circle')return `<circle cx="${n(o.p.x)}" cy="${n(o.p.y)}" r="${o.r}" fill="${o.fill}" stroke="${o.stroke}" ${style}/>`;
      if(o.type==='ring')return `<circle cx="${n(o.cx)}" cy="${n(o.cy)}" r="${n(o.r)}" fill="none" stroke="${o.color}" stroke-width="${o.width}" stroke-dasharray="${(o.dash||[]).join(' ')}" ${style}/>`;
      if(o.type==='text'&&o.tex)return `<g ${style}>${window.SubspaceSvgExport.mathLabel(o.tex,o.p.x+o.offset[0],o.p.y+o.offset[1]+(o.textAnchor==='middle'?o.size*.35:0),o.size,o.color,o.textAnchor||'start')}</g>`;
      if(o.type==='text')return `<text x="${n(o.p.x+o.offset[0])}" y="${n(o.p.y+o.offset[1])}" font-family="Latin Modern Roman,Georgia,serif" font-size="${o.size}" fill="${o.color}" text-anchor="${o.textAnchor||'start'}"${o.textAnchor==='middle'?' dominant-baseline="middle"':''} ${style}>${esc(o.content)}</text>`;
      return '';
    }).join('')+'</svg>';
  }
    function download(name,text,type){const url=URL.createObjectURL(new Blob([text],{type})),a=document.createElement('a');a.href=url;a.download=name;a.click();setTimeout(()=>URL.revokeObjectURL(url),1000);}
  function method(model){return `<p>The <strong>Wigner function</strong> uses the all-positive Stratonovich–Weyl spin kernel, scaled to ${math('\\int W_V\\,d\\Omega=1')}. In the tensor convention used here, ${math('W_V(\\mathbf n)=\\sqrt{(N+1)/(4\\pi)}\\sum_{\\ell=0}^{N}\\sum_{m=-\\ell}^{\\ell}\\rho_{\\ell m}Y_{\\ell m}(\\mathbf n)')}. The zero-centered color scale preserves negative values. In radial magnitude mode the radius is ${math('|W_V|/\\max|W_V|')}, with its sign shown by color. See <a href="https://doi.org/10.3390/app15095155" target="_blank" rel="noopener">Sánchez-Soto et al., <em>Phase Space Insights</em> (2025)</a> for the spin kernel.</p><p>${math('P_V=\\sum_{a=1}^{d}|u_a\\rangle\\langle u_a|')} and ${math('\\rho_V=P_V/d')} describe the whole subspace independently of its orthonormal basis. The plotted function is ${math('Q_V(\\mathbf n)=\\langle\\mathbf n^{\\otimes N}|\\rho_V|\\mathbf n^{\\otimes N}\\rangle')}. Its spherical mean is ${math('1/(N+1)')}; the projector overlap is ${math('dQ_V')}. Each selectable palette uses this subspace’s sampled minimum and maximum; changing colors preserves Q values and contour levels. Contours mark eight equally spaced interior Q values and rotate with the surface. Pale strokes on dark regions and dark strokes on bright regions improve contrast; the line color has no physical meaning. Radial mode uses radius ${math('Q_V/\\max Q_V')}.</p><p>The operator view follows <a href="https://doi.org/10.1103/PhysRevA.101.022332" target="_blank" rel="noopener">Serrano-Ensástiga and Braun, <em>Majorana representation for mixed states</em></a>. A nonzero rank-${math('\\ell')} multipole gives ${math('2\\ell')} antipodally paired stars with norm ${math('w_\\ell')}. With “Fit selected sphere” off, the shell radius is ${math('w_\\ell/\\max_{k\\ge1}w_k')}. Fitting a selected sphere uses a unit display radius for inspection; the true multipole norm remains listed. The all-rank view always compares relative norms. The “Multipole spheres” button places all multipoles side by side on a common radius scale. The selected representation is remembered when reopening the viewer or changing subspace. Dotted circles mark the maximum norm; scrolling reveals every rank. Short analytical weight expressions are displayed when verified from the exact source basis; other weights remain numerical. <a href="data/operator-exact-weights.json" download>Download exact weight formulas</a>. Zero multipoles have no constellation. Violet/rank-colored stars choose a signed subconstellation, and amber stars are its antipodal partners. Both are filled: dark shading indicates the front and light shading the rear. Swapping an even number of selected stars with their antipodal partners preserves the signed block; an odd number reverses its sign. Signed tensor coefficients are also retained in the operator download.</p><p>Symmetry axes are tested against the complete signed operator, rather than one basis vector or an unsigned constellation. Detection uses numerical tolerances, and is not an exact proof. Low multipoles vanishing in ${math('\\rho_V')} are a consequence of AC; the whole-subspace condition also involves off-diagonal basis transitions.</p>`;}
  async function open(subspace){
    window.SubspaceSphereDetail?.close();ensure();if(dialog.open){cleanup();cleanup=()=>{};}returnFocus=document.activeElement;current=subspace;const ticket=++serial;
    $('label').textContent=subspace.label||subspace.id||'Anticoherent subspace';$('definition').innerHTML=math('\\rho_V=P_V/d')+' · Basis-independent whole-subspace view';$('readout').innerHTML='';$('math-labels').replaceChildren();$('symmetry').innerHTML='';$('symmetry').hidden=false;$('gallery-axis-legend').hidden=true;$('gallery-axis-legend').innerHTML='';$('method').innerHTML='';$('verification').innerHTML='';$('husimi').setAttribute('aria-pressed','true');$('wigner').setAttribute('aria-pressed','false');$('operator').setAttribute('aria-pressed','false');$('spheres').setAttribute('aria-pressed','false');$('shape-control').hidden=false;$('rank-control').hidden=true;$('shape').value='sphere';$('palette').value=palettePreference;$('contours').checked=true;$('fit').checked=true;$('separate').checked=false;$('separate-control').hidden=true;$('axes').checked=true;$('coordinate').checked=false;
    if(!dialog.open)dialog.showModal();document.body.classList.add('sv-open');message('Computing the subspace projector…');
    await new Promise(resolve=>requestAnimationFrame(()=>setTimeout(resolve,0)));if(ticket!==serial)return;
    let model,mesh,symmetry;
    try{let saved=cache.get(subspace);if(!saved){model=window.OperatorMath.fromBasis(subspace);mesh=buildMesh(model);symmetry=window.OperatorMath.symmetry(model);saved={model,mesh,symmetry};}cache.delete(subspace);cache.set(subspace,saved);while(cache.size>4)cache.delete(cache.keys().next().value);({model,mesh,symmetry}=saved);}catch(e){message('The whole-subspace view is unavailable: '+e.message);return;}
    if(ticket!==serial)return;
    $('definition').innerHTML=math(`N=${model.n},\\quad d=${model.dimension},\\quad\\rho_V=P_V/${model.dimension}`);
    const ranks=model.multipoles.filter(m=>m.rank>0),first=ranks.find(m=>!m.zero)?.rank||1;let rank=first,mode='husimi',activeRecord=null,wignerMesh=null;
    $('rank').innerHTML='<option value="all">All nonzero ranks</option>'+ranks.map(m=>`<option value="${m.rank}">${m.rank}${m.zero?' — vanishes':''}</option>`).join('');$('rank').value=String(rank);
    const axes=symmetry.axes||[],orders=[...new Set(axes.map(a=>a.order==null?'continuous':String(a.order)))].sort((a,b)=>a==='continuous'?-1:b==='continuous'?1:Number(b)-Number(a));
    $('order').innerHTML='<option value="highest">Highest order</option><option value="all">All axes</option>'+orders.map(o=>`<option value="${esc(o)}">${o==='continuous'?'Continuous':o+'-fold'}</option>`).join('');$('order').value='highest';$('axes').disabled=!axes.length;$('order').disabled=!axes.length;$('align').disabled=!axes.length;
    const counts=orders.map(o=>{const count=axes.filter(a=>(a.order==null?'continuous':String(a.order))===o).length;return `<span class="sv-axis-key"><i style="background:${axisColor(o)}"></i>${math('C_{'+(o==='continuous'?'\\infty':o)+'}')} : ${count} ${count===1?'axis':'axes'}</span>`;}).join('');
    $('symmetry').innerHTML=symmetry.isotropic?'<p>Rotationally invariant operator: every spatial axis.</p>':counts||'<p>No nontrivial rotational symmetry detected.</p>';
    if(symmetry.perpendicularHalfTurns)$('symmetry').innerHTML+='<p class="sv-note">Every axis perpendicular to '+math('C_\\infty')+' is also a '+math('C_2')+' axis.</p>';
    $('symmetry').innerHTML+=`<p class="sv-note">Whole-projector rotations · numerical tolerance ${esc(fmt(symmetry.tolerance||1e-6))}</p>`;
    $('method').innerHTML=method(model);
    const low=model.multipoles.filter(m=>m.rank>0&&m.rank<=(subspace.t||0)).reduce((s,m)=>s+m.weight*m.weight,0);
    $('verification').innerHTML=`<p>${math('\\operatorname{Tr}\\rho_V^2=1/d')} (a normalized rank-${model.dimension} projector). Multipole ranks run from ${math('0')} to ${math(String(model.n))}; the scalar rank ${math('0')} is fixed by the trace.</p><p>Reported AC order ${esc(subspace.t??'—')}; low-rank squared weight ${esc(fmt(low))}. Coefficient diagnostics: ${esc(JSON.stringify(model.diagnostics||{}))}.</p><div id="sv-bars" class="sv-multipole-bar" aria-label="Multipole norms"></div><p class="sv-bar-label">Multipole norms, rank 1 → ${model.n}; select a bar to inspect its constellation.</p><div id="sv-tensor" class="sv-tensor-values"></div>`;
    const max=Math.max(...ranks.map(m=>m.weight),1e-15);$('bars').innerHTML=ranks.map(m=>`<button type="button" data-rank="${m.rank}" style="height:${Math.max(2,58*m.weight/max)}px" title="Rank ${m.rank}: ${esc(fmt(m.weight))}" aria-label="Rank ${m.rank}, norm ${esc(fmt(m.weight))}" aria-pressed="false"${m.zero?' disabled':''}></button>`).join('');
    const openRank=selected=>{const m=ranks.find(x=>x.rank===selected);if(!m||m.zero)return;const rec=window.OperatorMath.constellation(model,selected);window.SubspaceSphereDetail.open({model,mesh,symmetry,record:rec,rank:selected,weightTex:weightTex(subspace,selected,m.weight),exactWeight:!!exactWeight(subspace,selected,m.weight),weightNumeric:m.weight,subspaceId:subspace.id});};
    const renderer=makeRenderer(model,mesh,symmetry,$,openRank);
    if(axes.length)renderer.align();
    function updateGalleryAxisLegend(){
      const legend=$('gallery-axis-legend'),gallery=mode==='operator'&&rank==='all'&&$('separate').checked;
      const shown=gallery?renderer.getDisplayedAxes():[],groups=new Map();
      for(const axis of shown){const order=axis.order==null?'continuous':String(axis.order);groups.set(order,(groups.get(order)||0)+1);}
      legend.hidden=!groups.size;
      legend.innerHTML=groups.size?'<span class="sv-gallery-axis-label">Axes on all spheres:</span>'+[...groups].map(([order,count])=>`<span class="sv-gallery-axis-item"><i style="background:${axisColor(order)}" aria-hidden="true"></i>${math('C_{'+(order==='continuous'?'\\infty':order)+'}')}<span class="sv-gallery-axis-count">(${count} ${count===1?'axis':'axes'})</span></span>`).join(''):'';
      $('symmetry').hidden=!legend.hidden;
    }
    function info(){
      updateGalleryAxisLegend();
      if(mode==='husimi')$('readout').innerHTML=`<p>${math('Q_V=\\langle\\mathbf n^{\\otimes N}|\\rho_V|\\mathbf n^{\\otimes N}\\rangle')}</p><div class="sv-color-key"><span>${fmt(mesh.qMin)}</span><i class="sv-color-ramp" style="background:${window.HusimiRendering.gradient(palettePreference)}" aria-hidden="true"></i><span>${fmt(mesh.qMax)}</span></div><p class="sv-note">Sampled angular range · ${$('shape').value==='radial'?'Radius proportional to Q':'Color represents Q'}${$('contours').checked&&mesh.contours.levels.length?' · '+mesh.contours.levels.length+' equally spaced Q contours':''}</p>${$('contours').checked&&mesh.contours.levels.length?'<p class="sv-note">Contours join equal values of '+math('Q_V')+'. Light and dark strokes are chosen only for contrast.</p>':''}`;
      else if(mode==='wigner'){
        const bound=wignerMesh.colorMax;$('readout').innerHTML=`<p>${math('W_V(\\mathbf n),\\quad\\int W_V\\,d\\Omega=1')}</p><div class="sv-color-key"><span>${fmt(-bound)}</span><i class="sv-color-ramp" style="background:${window.HusimiRendering.gradient('wigner')}" aria-hidden="true"></i><span>${fmt(bound)}</span></div><p class="sv-note">Blue: negative · White: zero · Orange: positive.</p><p class="sv-note">Sampled range: ${fmt(wignerMesh.qMin)} to ${fmt(wignerMesh.qMax)}. ${$('shape').value==='radial'?'Radius proportional to '+math('|W_V|')+'; color retains the sign.':'Color represents the signed Wigner function.'}</p>`;
      }
      else if(rank==='all'){$('readout').innerHTML=($('separate').checked?'<p class="sv-note">Separate spheres share one radius scale. Dotted circles show the maximum norm; drag rotates them together. Click a sphere or its rank label to enlarge it independently. Scroll to see every rank.</p>':'')+(activeRecord?.records||[]).map(r=>`<span class="sv-axis-key"><i style="background:${rankColor(r.rank)}"></i>${math('\\ell='+r.rank+',\\ w_\\ell='+weightTex(subspace,r.rank,r.weight))}${r.available?'':' — roots unavailable'}</span>`).join('')+'<p class="sv-note">Shell radius follows multipole norm. Violet/rank colors mark the signed subconstellation; amber marks its antipodal partners. Dark stars are in front; light stars are at the rear. Signed coefficients remain in the operator download. Maximum signed-block reconstruction residual '+fmt(Math.max(0,...(activeRecord?.records||[]).filter(r=>r.available).map(r=>r.reconstructionResidual||0)))+'.</p>';}
      else{const m=ranks.find(x=>x.rank===rank);$('readout').innerHTML=`<p>${math('\\ell='+rank+',\\quad w_\\ell='+weightTex(subspace,rank,m?.weight||0))}</p>`+(exactWeight(subspace,rank,m?.weight||0)?`<p class="sv-note">${math('w_\\ell\\simeq '+fmt(m?.weight||0))}</p>`:'')+(m?.zero?'<p>This multipole vanishes; it has no stars.</p>':activeRecord?.available?`<p>${2*rank} stars, counted with multiplicity.</p><p class="sv-note">${$('fit').checked?'Sphere fitted for inspection; the multipole norm above is unchanged.':'Shaded sphere: active multipole.'+((m.weight/max)<1-1e-8?' Outer dotted circle: maximum multipole norm.':'')}</p><span class="sv-sign-key"><i></i>Selected subconstellation</span><span class="sv-sign-key"><i class="opposite"></i>Antipodal partners</span><p class="sv-note">Dark stars: front · Light stars: rear.</p><p class="sv-note">${activeRecord.signReliable?'The marked subconstellation retains the signed block.':'Signed coefficients are retained below and in the download.'} Root residual ${fmt(activeRecord.rootResidual)}; signed-block reconstruction residual ${fmt(activeRecord.reconstructionResidual)}.</p>`:`<p>${esc(activeRecord?.error||'Constellation unavailable at this rank.')}</p><p class="sv-note">The signed coefficients remain available below and in the download.</p>`);}
      const inspected=rank==='all'?ranks.filter(m=>!m.zero):ranks.filter(m=>m.rank===rank);$('tensor').textContent=inspected.map(m=>`Rank ${m.rank}, m=${m.rank}…${-m.rank}: `+(m.coefficients||[]).map(v=>Array.isArray(v)?fmt(v[0])+(v[1]<0?' − ':' + ')+fmt(Math.abs(v[1]))+'i':String(v)).join('; ')).join('\n\n');
      $('bars').querySelectorAll('button').forEach(b=>b.setAttribute('aria-pressed',String(mode==='operator'&&(rank==='all'&&!b.disabled||Number(b.dataset.rank)===rank))));
    }
    function change(next){if(next==='wigner'&&!wignerMesh)wignerMesh=buildWignerMesh(model);if(rank!=='all')$('separate').checked=false;const enteringOperator=next==='operator'&&mode!=='operator',enteringWigner=next==='wigner'&&mode!=='wigner';mode=next;const gallery=mode==='operator'&&rank==='all'&&$('separate').checked;rememberRepresentation(gallery?'spheres':mode);$('spheres').setAttribute('aria-pressed',String(gallery));$('husimi').setAttribute('aria-pressed',String(mode==='husimi'));$('operator').setAttribute('aria-pressed',String(mode==='operator'&&!gallery));$('wigner').setAttribute('aria-pressed',String(mode==='wigner'));$('shape-control').hidden=mode==='operator';$('shape').querySelector('option[value="radial"]').textContent=mode==='wigner'?'Radial magnitude':'Radial surface';$('rank-control').hidden=mode!=='operator';$('separate-control').hidden=mode!=='operator';$('contour-control').hidden=mode==='operator';$('palette-control').hidden=mode!=='husimi';$('fit-control').hidden=mode!=='operator'||rank==='all';if(mode==='operator'){if(rank==='all'){activeRecord={all:true,records:ranks.filter(m=>!m.zero).map(m=>{try{return window.OperatorMath.constellation(model,m.rank);}catch(e){return {rank:m.rank,weight:m.weight,available:false,error:e.message};}})};}else{const m=ranks.find(x=>x.rank===rank);try{activeRecord=m?.zero?{zero:true}:window.OperatorMath.constellation(model,rank);}catch(e){activeRecord={available:false,error:e.message};}}}renderer.set(mode,rank,activeRecord,mode==='wigner'?wignerMesh:mesh);if(enteringOperator||enteringWigner)renderer.align();info();}
    const onH=()=>change('husimi'),onW=()=>change('wigner'),onO=()=>{rank='all';$('rank').value='all';$('separate').checked=false;change('operator');},onSpheres=()=>{rank='all';$('rank').value='all';$('separate').checked=true;change('operator');},onSeparate=()=>{if($('separate').checked){rank='all';$('rank').value='all';}change('operator');},onRank=()=>{rank=$('rank').value==='all'?'all':Number($('rank').value);change('operator');},onAxisLegend=()=>updateGalleryAxisLegend(),onShape=()=>info(),onFit=()=>info(),onPalette=()=>info(),onContours=()=>info(),onBars=e=>{const b=e.target.closest('[data-rank]');if(b){rank=Number(b.dataset.rank);$('rank').value=String(rank);change('operator');}},onSvg=()=>{
      const included=mode!=='operator'?[]:rank!=='all'?ranks.filter(m=>m.rank===rank):ranks.filter(m=>!m.zero);
      const formulaRows=included.map(m=>({tex:'\\ell='+m.rank+',\\quad w_{'+m.rank+'}='+weightTex(subspace,m.rank,m.weight),color:mode==='operator'&&rank==='all'?rankColor(m.rank):null,note:exactWeight(subspace,m.rank,m.weight)?'≈ '+fmt(m.weight):'Numerical value'}));
      if(mode==='operator'&&!formulaRows.length)formulaRows.push({tex:'w_\\ell=0\\quad(1\\le\\ell\\le '+model.n+')'});
      const symmetryRows=orders.map(o=>{const count=axes.filter(a=>(a.order==null?'continuous':String(a.order))===o).length;return {tex:'C_{'+(o==='continuous'?'\\infty':o)+'}:\\;'+count+'\\,\\text{'+(count===1?'axis':'axes')+'}',color:axisColor(o)};});
      const notes=[];
      if(!axes.length)notes.push(symmetry.isotropic?'Rotationally invariant: every spatial axis.':'No nontrivial rotational symmetry detected.');
      if(symmetry.perpendicularHalfTurns)notes.push('Every perpendicular axis is also a C₂ axis.');
      notes.push('Whole-projector symmetries; numerical tolerance '+fmt(symmetry.tolerance||1e-6)+'.');
      if(mode==='husimi'&&$('contours').checked&&mesh.contours.levels.length)notes.push('Contours join equal Husimi values. Light and dark strokes are chosen only for contrast.');
      if(mode==='operator')notes.push(rank==='all'?($('separate').checked?'Separate spheres share one common multipole-norm radius scale.':'Sphere radii follow multipole norms.'):$('fit').checked?'Selected sphere fitted for inspection; its weight is unchanged.':'Sphere radius follows multipole norm.');
      const colorScale=mode==='wigner'?{title:'Wigner color scale',mid:'0',stops:window.HusimiRendering.palettes.wigner.map(p=>({t:p.t,color:window.HusimiRendering.color(p.t,'wigner')})),min:fmt(-wignerMesh.colorMax),max:fmt(wignerMesh.colorMax),label:'W_V',note:'Blue: negative; white: zero; orange: positive. The scale is centered at zero.'}:mode==='husimi'?{stops:window.HusimiRendering.palettes[palettePreference].map(p=>({t:p.t,color:window.HusimiRendering.color(p.t,palettePreference)})),min:fmt(mesh.qMin),max:fmt(mesh.qMax),label:'Q_V',note:'Sampled angular range'+($('contours').checked?' · '+mesh.contours.levels.length+' contours':'')}:null;
      const svg=window.SubspaceSvgExport.compose(renderer.exportSvg(),{title:mode==='operator'?'Operator Majorana constellations':mode==='wigner'?'Whole-subspace Wigner function':'Whole-subspace Husimi function',subtitle:subspace.label||subspace.id,contextTex:'N='+model.n+',\\quad d='+model.dimension+',\\quad\\rho_V=P_V/d'+(mode==='wigner'?',\\quad\\int W_V\\,d\\Omega=1':''),formulaRows,colorScale,symmetryRows,notes});
      download((subspace.id||'subspace')+'-'+mode+'.svg',svg,'image/svg+xml');
    },onJson=()=>download((subspace.id||'subspace')+'-operator.json',JSON.stringify({subspace:subspace.id,n:model.n,dimension:model.dimension,normalization:'rho=P_V/d',phaseSpace:{representation:mode,wignerConvention:'All-positive Stratonovich–Weyl spin kernel, integral W dOmega = 1',wignerReference:'https://doi.org/10.3390/app15095155'},rho:model.rho,multipoles:model.multipoles,exactMultipoleWeights:window.OperatorExactWeights?.states?.[subspace.id]?.weights||{},constellations:mode==='operator'?(activeRecord?.all?activeRecord.records:[activeRecord]):[],symmetry,diagnostics:model.diagnostics,reference:'E. Serrano-Ensástiga and D. Braun, Phys. Rev. A 101, 022332 (2020), doi:10.1103/PhysRevA.101.022332'},null,2),'application/json');
    $('husimi').addEventListener('click',onH);$('operator').addEventListener('click',onO);$('wigner').addEventListener('click',onW);$('spheres').addEventListener('click',onSpheres);$('rank').addEventListener('change',onRank);$('separate').addEventListener('change',onSeparate);$('axes').addEventListener('change',onAxisLegend);$('order').addEventListener('change',onAxisLegend);$('shape').addEventListener('change',onShape);$('palette').addEventListener('change',onPalette);$('fit').addEventListener('change',onFit);$('contours').addEventListener('change',onContours);$('bars').addEventListener('click',onBars);$('svg').addEventListener('click',onSvg);$('json').addEventListener('click',onJson);
    cleanup=()=>{renderer.dispose();$('husimi').removeEventListener('click',onH);$('operator').removeEventListener('click',onO);$('wigner').removeEventListener('click',onW);$('spheres').removeEventListener('click',onSpheres);$('rank').removeEventListener('change',onRank);$('separate').removeEventListener('change',onSeparate);$('axes').removeEventListener('change',onAxisLegend);$('order').removeEventListener('change',onAxisLegend);$('shape').removeEventListener('change',onShape);$('palette').removeEventListener('change',onPalette);$('fit').removeEventListener('change',onFit);$('contours').removeEventListener('change',onContours);$('bars').removeEventListener('click',onBars);$('svg').removeEventListener('click',onSvg);$('json').removeEventListener('click',onJson);};
    message('');if(representationPreference==='spheres'||representationPreference==='operator'){rank='all';$('rank').value='all';$('separate').checked=representationPreference==='spheres';change('operator');}else change(representationPreference);renderer.redraw();
  }
  window.SubspaceSphereRenderer=makeRenderer;
  window.SubspaceViewer={open,close};
})();
