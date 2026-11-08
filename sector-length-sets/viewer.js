/* Original, build-free viewers. Three.js is vendored; no network requests. */
(function () {
  'use strict';
  const T = window.THREE;
  const colors = {violet:0x7946b0, pale:0xc3a5df, gray:0xb5b6ba, edge:0x67696b, grid:0xdeddd6};
  const v3 = p => new T.Vector3(...p);
  const svgNS = 'http://www.w3.org/2000/svg';
  function line(points, color, dashed=false, hidden=false) {
    const g=new T.BufferGeometry().setFromPoints(points.map(v3));
    const m=dashed ? new T.LineDashedMaterial({color,dashSize:.045,gapSize:.035,transparent:true,opacity:.42,depthWrite:false,depthFunc:hidden?T.GreaterDepth:T.LessEqualDepth}) : new T.LineBasicMaterial({color,transparent:true,opacity:.9,depthWrite:false});
    const out=new T.Line(g,m);out.computeLineDistances();out.renderOrder=3;return out;
  }
  function disposeGroup(group) {group.traverse(o=>{o.geometry?.dispose();if(o.material){(Array.isArray(o.material)?o.material:[o.material]).forEach(m=>m.dispose());}});group.clear();}
  class BaseView {
    constructor(canvas, overlay, options={}) {
      this.canvas=canvas;this.overlay=overlay;this.container=canvas.parentElement;this.labels=[];this.markers=[];
      this.renderer=new T.WebGLRenderer({canvas,antialias:true,alpha:true,preserveDrawingBuffer:true});this.renderer.setPixelRatio(Math.min(devicePixelRatio||1,2));this.renderer.setClearColor(0xfaf9f5,0);
      this.scene=new T.Scene();this.camera=new T.OrthographicCamera(-2,2,2,-2,.01,100);this.target=new T.Vector3();this.yaw=.8;this.pitch=.45;this.zoom=1;this.rotatable=true;
      this.scene.add(new T.AmbientLight(0xffffff,2.1));const sun=new T.DirectionalLight(0xffffff,1.4);sun.position.set(-3,6,4);this.scene.add(sun);
      this.content=new T.Group();this.scene.add(this.content);this.svg=document.createElementNS(svgNS,'svg');this.overlay.append(this.svg);this.markerSvg=document.createElementNS(svgNS,'svg');this.markerSvg.classList.add('screen-markers');this.markerSvg.setAttribute('aria-hidden','true');this.overlay.append(this.markerSvg);
      this._resize=new ResizeObserver(()=>this.resize());this._resize.observe(this.container);
      this.drag=null;this.moved=false;
      canvas.addEventListener('pointerdown',e=>this.beginPointerGesture(e));
      canvas.addEventListener('pointermove',e=>this.updatePointerGesture(e));
      canvas.addEventListener('pointerup',e=>this.endPointerGesture(e));
      canvas.addEventListener('contextmenu',e=>{if(this.panAllowed)e.preventDefault();});
      canvas.addEventListener('pointercancel',e=>this.cancelPointerGesture(e));
      canvas.addEventListener('lostpointercapture',e=>this.cancelPointerGesture(e));
      canvas.addEventListener('pointerleave',()=>{if(!this.drag)this.onLeave?.();});
      canvas.addEventListener('wheel',e=>{
        e.preventDefault();
        if(this.n===6&&this.drag?.point)return;
        // Keep the pointed-at part of a sector set fixed while magnifying it.
        const rect=canvas.getBoundingClientRect(),ndc=new T.Vector3(2*(e.clientX-rect.left)/rect.width-1,1-2*(e.clientY-rect.top)/rect.height,0);
        const anchor=this.panAllowed?ndc.clone().unproject(this.camera):null;
        this.zoom=Math.max(.55,Math.min(this.maxZoom||3.5,this.zoom*Math.exp(-e.deltaY*.001)));
        this.camera.zoom=this.zoom;this.camera.updateProjectionMatrix();
        if(anchor)this.target.add(anchor.sub(ndc.unproject(this.camera)));
        this.resize();
      },{passive:false});
      this.resize();
    }
    beginPointerGesture(event){
      if(this.drag||![0,1,2].includes(event.button)||(!this.panAllowed&&event.button!==0))return false;
      const point=event.button===0&&!event.shiftKey?this.beginPointDrag?.(event):null;
      const pan=!point&&this.panAllowed&&(!this.rotatable||event.shiftKey||event.button!==0);
      this.drag={x:event.clientX,y:event.clientY,pointerId:event.pointerId,yaw:this.yaw,pitch:this.pitch,
        pan,point,target:this.target.clone(),right:new T.Vector3().setFromMatrixColumn(this.camera.matrixWorld,0),up:new T.Vector3().setFromMatrixColumn(this.camera.matrixWorld,1)};
      this.moved=false;if(point)this.canvas.style.cursor='grabbing';
      try{this.canvas.setPointerCapture?.(event.pointerId);}catch{}
      return true;
    }
    updatePointerGesture(event){
      const drag=this.drag;
      if(!drag){this.onHover?.(event);return false;}
      if(drag.pointerId!==undefined&&event.pointerId!==undefined&&drag.pointerId!==event.pointerId)return false;
      const dx=event.clientX-drag.x,dy=event.clientY-drag.y;
      if(Math.hypot(dx,dy)>4)this.moved=true;
      if(!this.moved)return true;
      if(drag.point){this.updatePointDrag(event,drag);return true;}
      if(drag.pan){
        const sx=(this.camera.right-this.camera.left)/(this.camera.zoom*this.width),sy=(this.camera.top-this.camera.bottom)/(this.camera.zoom*this.height);
        this.target.copy(drag.target).addScaledVector(drag.right,-dx*sx).addScaledVector(drag.up,dy*sy);this.updateCamera();this.render();
      }else if(this.rotatable){
        this.yaw=drag.yaw-dx*.008;this.pitch=Math.max(-1.3,Math.min(1.3,drag.pitch+dy*.008));this.updateCamera();this.render();
      }
      return true;
    }
    releaseGestureCapture(pointerId){
      try{if(!this.canvas.hasPointerCapture||this.canvas.hasPointerCapture(pointerId))this.canvas.releasePointerCapture?.(pointerId);}catch{}
    }
    endPointerGesture(event){
      const drag=this.drag;if(!drag)return false;
      if(drag.pointerId!==undefined&&event.pointerId!==undefined&&drag.pointerId!==event.pointerId)return false;
      if(Math.hypot(event.clientX-drag.x,event.clientY-drag.y)>4)this.moved=true;
      const moved=this.moved;
      if(drag.point&&moved)this.finishPointDrag(event,drag,false);
      this.drag=null;this.moved=false;this.releaseGestureCapture(drag.pointerId);
      if(drag.point)this.canvas.style.cursor='grab';
      if(!moved&&event.button===0)this.onClick?.(event);
      return true;
    }
    cancelPointerGesture(event){
      const drag=this.drag;if(!drag)return false;
      if(event?.pointerId!==undefined&&drag.pointerId!==undefined&&event.pointerId!==drag.pointerId)return false;
      if(drag.point)this.finishPointDrag(null,drag,true);
      this.drag=null;this.moved=false;this.releaseGestureCapture(drag.pointerId);
      if(drag.point)this.canvas.style.cursor='grab';
      return true;
    }
    resize(){const r=this.container.getBoundingClientRect();this.width=Math.max(1,r.width);this.height=Math.max(1,r.height);this.renderer.setSize(this.width,this.height,false);const aspect=this.width/this.height,span=Math.max(this.span||1.55,(this.minHalfWidth||0)/aspect);this.camera.left=-span*aspect;this.camera.right=span*aspect;this.camera.top=span;this.camera.bottom=-span;this.camera.zoom=this.zoom;this.camera.updateProjectionMatrix();this.svg.setAttribute('viewBox',`0 0 ${this.width} ${this.height}`);this.markerSvg.setAttribute('viewBox',`0 0 ${this.width} ${this.height}`);this.updateCamera();this.render();}
    updateCamera(){const r=7;this.camera.position.set(r*Math.cos(this.pitch)*Math.sin(this.yaw),r*Math.sin(this.pitch),r*Math.cos(this.pitch)*Math.cos(this.yaw)).add(this.target);this.camera.lookAt(this.target);this.camera.updateMatrixWorld();}
    project(p){const v=v3(p).project(this.camera);return {x:(v.x+1)*this.width/2,y:(1-v.y)*this.height/2,z:v.z};}
    clearLabels(){this.labels.forEach(l=>l.el.remove());this.labels=[];this.markers=[];this.svg.replaceChildren();this.markerSvg.replaceChildren();}
    addLabel(text,position,type='plot-label',extra={}){const el=document.createElement('span');el.className=type;el.innerHTML=window.AtlasMath?AtlasMath.html(AtlasMath.axis(text)):text;this.overlay.append(el);const label={el,text,position,...extra};this.labels.push(label);return label;}
    layoutLabels(){this.svg.replaceChildren();const rects=[];for(const l of this.labels){if(l.el.hidden)continue;const p=this.project(l.position);if(p.z>1||p.x<0||p.y<0||p.x>this.width||p.y>this.height){l.el.style.display='none';continue;}l.el.style.display='';let x=p.x+(l.offset?.[0]||0),y=p.y+(l.offset?.[1]||0);const w=l.el.offsetWidth||40,h=l.el.offsetHeight||16;if(l.avoid){const candidates=l.fixedOffset?[l.fixedOffset]:[[17,-15],[-20,-16],[20,18],[-25,20],[0,-35],[35,0],[-45,0],[0,38],[25,-50],[-30,50],[55,-20],[-55,25]];let best=Infinity;for(const[dx,dy]of candidates){let xx=Math.max(w/2+4,Math.min(this.width-w/2-4,p.x+dx)),yy=Math.max(46+h/2,Math.min(this.height-h/2-20,p.y+dy));let penalty=0;for(const r of rects)if(Math.abs(xx-r.x)<(w+r.w)/2+5&&Math.abs(yy-r.y)<(h+r.h)/2+3)penalty+=1000;for(const n of this.namedProjections||[])if(Math.hypot(xx-n.x,yy-n.y)<Math.min(w/2,18)+8)penalty+=100;penalty+=Math.hypot(dx,dy);if(penalty<best){best=penalty;x=xx;y=yy;}}rects.push({x,y,w,h});const dx=x-p.x,dy=y-p.y,dist=Math.hypot(dx,dy);if(dist>8){const endRatio=Math.max(0,1-((Math.abs(dx)>Math.abs(dy)?w/2:h/2)+2)/dist),node=5.5/dist;const lead=document.createElementNS(svgNS,'line');lead.setAttribute('x1',p.x+dx*node);lead.setAttribute('y1',p.y+dy*node);lead.setAttribute('x2',p.x+dx*endRatio);lead.setAttribute('y2',p.y+dy*endRatio);lead.setAttribute('stroke','#6d716c');lead.setAttribute('stroke-width','.7');this.svg.append(lead);}}else rects.push({x,y,w,h});l.el.style.left=x+'px';l.el.style.top=y+'px';}}
    render(){if(!this.renderer)return;
      // Dash spacing stays readable at any zoom; faces alone determine depth.
      const unit=(this.camera.top-this.camera.bottom)/(this.camera.zoom*this.height);
      this.scene.traverse(o=>{if(o.isLine){for(const m of (Array.isArray(o.material)?o.material:[o.material]))if(m?.isLineDashedMaterial){m.dashSize=6*unit;m.gapSize=5*unit;m.depthWrite=false;}}});
      this.renderer.render(this.scene,this.camera);this.layoutLabels();this.renderMarkers();}
    renderMarkers(){
      this.markerSvg.replaceChildren();
      const projected=this.markers.map(m=>({...m,screen:this.project(m.world)}))
        .filter(m=>m.screen.z>=-1&&m.screen.z<=1&&m.screen.x>=-8&&m.screen.x<=this.width+8&&m.screen.y>=-8&&m.screen.y<=this.height+8)
        .sort((a,b)=>Number(b.rear)-Number(a.rear)||b.screen.z-a.screen.z);
      for(const m of projected){
        const {x,y}=m.screen;
        if(this.selectedWorld&&Math.hypot(...m.world.map((v,i)=>v-this.selectedWorld[i]))<1e-8){
          const halo=document.createElementNS(svgNS,'circle');
          for(const [key,value] of Object.entries({cx:x,cy:y,r:m.radius+2.2,fill:'none',stroke:'#315d51','stroke-width':1.1,class:'marker-halo'}))halo.setAttribute(key,value);
          this.markerSvg.append(halo);
        }
        const circle=document.createElementNS(svgNS,'circle');
        for(const [key,value] of Object.entries({cx:x,cy:y,r:m.radius,fill:m.fill,stroke:m.stroke,'stroke-width':.95,'vector-effect':'non-scaling-stroke',class:'point-marker','data-kind':m.kind,'data-rear':String(!!m.rear)}))circle.setAttribute(key,value);
        this.markerSvg.append(circle);
      }
    }
    exportImage(name){const a=document.createElement('a');a.download=name;a.href=this.canvas.toDataURL('image/png');a.click();}
  }
  class SectorView extends BaseView {
    constructor(canvas,overlay){super(canvas,overlay);this.panAllowed=true;this.maxZoom=80;this.surface=new T.Group();this.samples=new T.Group();this.edges=new T.Group();this.nodes=new T.Group();this.content.add(this.surface,this.samples,this.edges,this.nodes);this.raycaster=new T.Raycaster();this.hoverCallback=()=>{};this.selectCallback=()=>{};this.onHover=e=>this.pick(e,false);this.onClick=e=>this.pick(e,true);this.onLeave=()=>this.hoverCallback(null);}
    map(p){if(this.n===6)return [2*(p[0]+1)/9-1,2*(p[2]+1)/23-1,2*(p[1]+1.5)/17.5-1];const xr=this.ranges[0],yr=this.ranges[1];return [((p[0]-xr[0])/(xr[1]-xr[0])-.5)*2.65,this.dimension===1?0:((p[1]-yr[0])/(yr[1]-yr[0])-.5)*2,0];}
    setData(d){this.cancelPointerGesture();this.data=d;this.n=d.n||d.N;this.dimension=d.dimension;this.ranges=this.n===3?[[0,3],[0,1]]:this.n===4?[[0,4],[2,6]]:[[0,5],[3,10]];this.clearLabels();this.occlusionIndex=null;this.focusedStateId=null;[this.surface,this.samples,this.edges,this.nodes].forEach(disposeGroup);this.witnesses=d.boundarySamples||[];this.named=d.states||d.knownStates;this.grayPoints=[];this.allPoints=this.witnesses.map(s=>({state:s,world:this.map(s.sectors)}));this.rotatable=this.n===6;this.zoom=1;this.span=this.n===6?1.72:1.5;this.minHalfWidth=this.n===6?1.8:1.67;this.yaw=this.n===6?-1.04:0;this.pitch=this.n===6?.22:0;this.target.set(0,0,0);
      if(this.n===6){for(const patch of d.outer.patches)this.addPatch(patch);for(const edge of d.outer.edges){const p=edge.points.map(x=>this.map(x));const c=edge.attained?colors.violet:colors.edge;this.edges.add(line(p,c),line(p,c,true,true));}for(const path of d.attainedPaths||[]){const p=path.points.map(x=>this.map(x));const front=line(p,colors.violet),rear=line(p,colors.violet,true,true);front.userData.kind=rear.userData.kind='attained-state-locus';front.userData.pathId=rear.userData.pathId=path.id;this.edges.add(front,rear);}for(const node of (d.plotNodes||d.outer.nodes).filter(p=>!p.attained)){const gray={state:{...node,outerOnly:true,sectors:node.point,label:node.label||'Outer-region junction'},world:this.map(node.point)};this.grayPoints.push(gray);gray.marker={world:gray.world,radius:4.6,fill:'#989b9f',stroke:'#686c70',kind:'sector-gray'};this.markers.push(gray.marker);}}
      else{if(this.dimension===2){const shape=new T.Shape(d.polygon.map(p=>{const q=this.map(p);return new T.Vector2(q[0],q[1]);}));const geom=new T.ShapeGeometry(shape);const mesh=new T.Mesh(geom,new T.MeshBasicMaterial({color:0xdac7eb,side:T.DoubleSide}));mesh.position.z=-.02;this.surface.add(mesh);}for(const c of d.curves||[])this.edges.add(line(c.points.map(p=>this.map(p.sectors)),colors.violet));}
      // Manuscript nodes and the requested explicit boundary samples are plotted.
      for(const p of this.allPoints){p.marker={world:p.world,radius:3.4,fill:'#8854ba',stroke:'#64368e',kind:'boundary-sample'};this.markers.push(p.marker);}
      this.namedWorld=[];for(const s of this.named){const p=this.map(s.sectors),same=this.namedWorld.find(q=>Math.hypot(...q.world.map((v,i)=>v-p[i]))<1e-8);if(same){same.aliases.push(s);same.caption+=', '+AtlasMath.state(s.id,s.label);}else this.namedWorld.push({state:s,world:p,aliases:[s],caption:AtlasMath.state(s.id,s.label)});}for(const p of this.namedWorld){p.marker={world:p.world,radius:4.6,fill:'#7946b0',stroke:'#4d2676',kind:'sector-violet'};this.markers.push(p.marker);}
      this.selection=new T.Mesh(new T.RingGeometry(.019,.029,32),new T.MeshBasicMaterial({color:0x244e43,side:T.DoubleSide,depthTest:false}));this.selection.renderOrder=15;this.nodes.add(this.selection);this.selection.visible=false;this.selectedWorld=null;
      if(this.n===6&&window.StaticTriangleOccluder)this.occlusionIndex=new StaticTriangleOccluder(this.surface.children.map(m=>m.geometry.getAttribute('position').array));
      this.drawAxes();for(const p of this.namedWorld){p.label=this.addLabel(this.shortLabel(p.caption),p.world,'plot-label',{avoid:true,named:true,stateIds:p.aliases.map(s=>s.id),keyState:this.n<=5||p.aliases.some(s=>['n6-product','n6-ghz','n6-w','n6-octa'].includes(s.id)),fixedOffset:this.n===3?({'n3-ghz':[-18,-28],'n3-w':[18,43],'n3-product':[-18,-28]}[p.state.id]):null});}
      if(this.n===6)for(const callout of window.SECTOR_N6_MANUSCRIPT_CALLOUTS||[])this.addLabel('$'+callout.manuscriptLabel+'$',this.map(callout.point),'plot-label',{avoid:true,named:true,family:true});this.resize();this.applyOptions();}
    shortLabel(s){return s.replace(' (C₅ junction)','').replace(' (sharp C₂ endpoint)','').replace(' (pyramid)','').replace('Product','|0⟩⊗⁶').replace('Octa','|octa⟩');}
    addPatch(patch){
      const mapped=patch.positions.map(p=>this.map(p));
      let normals=patch.worldNormals;
      if(!normals){const indexed=new T.BufferGeometry();indexed.setAttribute('position',new T.Float32BufferAttribute(mapped.flat(),3));indexed.setIndex(patch.triangles.flat());indexed.computeVertexNormals();const a=indexed.getAttribute('normal');normals=mapped.map((_,i)=>[a.getX(i),a.getY(i),a.getZ(i)]);indexed.dispose();}
      const pos=[],col=[],normal=[];
      for(let i=0;i<patch.triangles.length;i++){
        const c=new T.Color(patch.attained[i]?0xb894d6:0xb9babd);
        for(const j of patch.triangles[i]){pos.push(...mapped[j]);col.push(c.r,c.g,c.b);normal.push(...normals[j]);}
      }
      const g=new T.BufferGeometry();g.setAttribute('position',new T.Float32BufferAttribute(pos,3));g.setAttribute('color',new T.Float32BufferAttribute(col,3));g.setAttribute('normal',new T.Float32BufferAttribute(normal,3));
      const mesh=new T.Mesh(g,new T.MeshLambertMaterial({vertexColors:true,side:T.DoubleSide,polygonOffset:true,polygonOffsetFactor:1,polygonOffsetUnits:1}));mesh.userData.kind=patch.kind;mesh.userData.patch=patch;this.surface.add(mesh);
    }
    drawAxes(){if(this.axisGroup){this.content.remove(this.axisGroup);disposeGroup(this.axisGroup);}this.axisGroup=new T.Group();this.content.add(this.axisGroup);const add=(p,c=colors.grid)=>this.axisGroup.add(line(p,c));if(this.n===6){for(const a of [0,1,2])for(const b of [-1,1])for(const c of [-1,1]){const p=[0,0,0],q=[0,0,0];p[a]=-1.04;q[a]=1.04;p[(a+1)%3]=q[(a+1)%3]=b*1.04;p[(a+2)%3]=q[(a+2)%3]=c*1.04;add([p,q]);}const axes=[{a:0,values:[0,2,4,6,8],limits:[-1,8],name:'S₁',point:v=>[v,-1.04,1.04]},{a:2,values:[0,4,8,12,16],limits:[-1.5,16],name:'S₂',point:v=>[-1.04,-1.04,v]},{a:1,values:[0,8,16,22],limits:[-1,22],name:'S₃',point:v=>[-1.04,v,-1.04]}];for(const ax of axes){const [lo,hi]=ax.limits;for(const val of ax.values){const p=ax.point(2*(val-lo)/(hi-lo)-1);this.addLabel(String(val),p,'plot-label axis-tick',{offset:ax.a===1?[-14,0]:[0,14]});}const p=ax.point(1.25);this.addLabel(ax.name,p,'plot-label axis-label',{offset:[10,0]});}}else{const xlo=this.ranges[0][0],xhi=this.ranges[0][1],ylo=this.ranges[1][0],yhi=this.ranges[1][1];if(this.dimension===2){for(let x=xlo;x<=xhi;x++){const p=this.map([x,ylo]),q=this.map([x,yhi]);p[2]=q[2]=-.03;add([p,q]);this.addLabel(String(x),p,'plot-label axis-tick',{offset:[0,17]});}const step=this.n===5?2:1;for(let y=ylo;y<=yhi;y+=step){const p=this.map([xlo,y]),q=this.map([xhi,y]);p[2]=q[2]=-.03;add([p,q]);this.addLabel(String(y),p,'plot-label axis-tick',{offset:[-18,0]});}this.addLabel('S₂',[-1.48,1.15,0],'plot-label axis-label');}else{add([[-1.325,0,0],[1.325,0,0]],colors.edge);for(let x=xlo;x<=xhi;x++){const p=this.map([x]);this.addLabel(String(x),p,'plot-label axis-tick',{offset:[0,25]});}}this.addLabel('S₁',[1.52,this.dimension===2?-1:0,0],'plot-label axis-label');}}
    updateLabelVisibility(){
      const mode=document.getElementById('show-labels').value;
      this.labels.forEach(l=>{if(l.named)l.el.hidden=mode==='none'||(mode!=='all'&&!l.keyState&&!l.stateIds?.includes(this.focusedStateId));});
    }
    applyOptions(){if(!this.surface)return;this.surface.visible=document.getElementById('show-outer').checked;this.samples.visible=this.n===5;this.updateLabelVisibility();this.render();}
    setPreset(name){const presets={main:[-1.04,.22],alternative:[-1.95,.50],rear:[2.10,.377]};[this.yaw,this.pitch]=presets[name]||presets.main;this.zoom=1;this.target.set(0,0,0);this.resize();}
    resetView(){if(this.n===6)this.setPreset('main');else{this.target.set(0,0,0);this.zoom=1;this.resize();}}
    setSelected(state){
      if(!this.selection)return;
      this.selection.position.copy(v3(this.map(state.sectors)));this.selection.visible=false;this.selectedWorld=this.map(state.sectors);this.focusedStateId=state.id;
      this.markers=this.markers.filter(m=>m.kind!=='selected-point');
      if(state.interactivePoint)this.markers.push({world:this.selectedWorld,radius:4.6,fill:'#7946b0',stroke:'#4d2676',kind:'selected-point'});
      for(const l of this.labels.filter(l=>l.focusOnly))l.el.remove();this.labels=this.labels.filter(l=>!l.focusOnly);
      if(state.id&&!this.namedWorld.some(p=>p.aliases.some(s=>s.id===state.id)))this.addLabel(AtlasMath.state(state.id,state.mathLabel||state.label||state.family),this.selectedWorld,'plot-label',{avoid:true,named:true,focusOnly:true,stateIds:[state.id]});
      this.updateLabelVisibility();this.render();
    }
    render(){if(this.selection)this.selection.quaternion.copy(this.camera.quaternion);if(this.namedWorld){this.namedProjections=this.namedWorld.map(p=>this.project(p.world)).concat((this.grayPoints||[]).map(p=>this.project(p.world)));for(const p of this.namedWorld){const rear=this.n===6&&this.occluded(p.world);p.marker.rear=rear;p.marker.fill=rear?'#b49acb':'#7946b0';p.marker.stroke=rear?'#8c72a4':'#4d2676';p.label&&(p.label.el.style.opacity=rear?'.65':'1');}for(const p of this.allPoints){const rear=this.n===6&&this.occluded(p.world);p.marker.rear=rear;p.marker.fill=rear?'#b49acb':'#8854ba';p.marker.stroke=rear?'#8c72a4':'#64368e';}for(const p of this.grayPoints){const rear=this.occluded(p.world);p.marker.rear=rear;p.marker.fill=rear?'#bec1c3':'#989b9f';p.marker.stroke=rear?'#a0a4a8':'#686c70';}for(const marker of this.markers.filter(m=>m.kind==='selected-point')){const rear=this.n===6&&this.occluded(marker.world);marker.rear=rear;marker.fill=rear?'#b49acb':'#7946b0';marker.stroke=rear?'#8c72a4':'#4d2676';this.labels.filter(l=>l.focusOnly).forEach(l=>l.el.style.opacity=rear?'.65':'1');}this.projected=this.allPoints.map(p=>({...this.project(p.world),point:p}));}super.render();}
    occluded(world){if(!this.surface.visible||this.n!==6)return false;const q=v3(world).project(this.camera);this.raycaster.setFromCamera(new T.Vector2(q.x,q.y),this.camera);const distance=v3(world).sub(this.raycaster.ray.origin).dot(this.raycaster.ray.direction);if(this.occlusionIndex)return this.occlusionIndex.occludes(this.raycaster.ray,distance-.025);const hits=this.raycaster.intersectObjects(this.surface.children,false);return hits.length>0&&hits[0].distance<distance-.025;}

    selectedPointAt(x,y){
      if(![3,4,5,6].includes(this.n)||!this.selectedWorld)return false;
      const callback=this.n===6?this.surfaceDragCallback:this.pointDragCallback;
      if(typeof callback!=='function')return false;
      const p=this.project(this.selectedWorld);
      if(p.z < -1||p.z>1||Math.hypot(p.x-x,p.y-y)>13)return false;
      if(this.n!==6)return true;
      // A rear marker is still selectable, but cannot drag through the front
      // body. At a visible selected point its own first face must be violet.
      if(!this.surface?.visible||this.occluded(this.selectedWorld))return false;
      const hit=this.attainedSurfaceAt(p.x,p.y);if(!hit)return false;
      const separation=v3(this.map(hit.point)).distanceTo(v3(this.selectedWorld));
      // This allowance is only for analytical-graph/triangle chord error;
      // unlike the 13px grab target it must not grow when zooming out.
      return separation<=.002;
    }
    beginPointDrag(event){
      const rect=this.canvas.getBoundingClientRect(),x=event.clientX-rect.left,y=event.clientY-rect.top;
      if(!this.selectedPointAt(x,y))return null;
      const p=this.project(this.selectedWorld);
      if(this.n===6){
        const hit=this.attainedSurfaceAt(p.x,p.y);if(!hit)return null;
        const [wx,wy,wz]=this.selectedWorld,coordinates=[4.5*(wx+1)-1,8.75*(wz+1)-1.5,11.5*(wy+1)-1];
        return {n:6,offset:[p.x-x,p.y-y],initial:coordinates,lastCoordinates:coordinates.slice(),
          initialMeta:hit,lastMeta:hit,acceptedScreen:[p.x,p.y],acceptedWorld:this.selectedWorld.slice(),
          started:false,latest:null,frame:null};
      }
      const coordinates=this.regionCoordinates(p.x,p.y,{clampToSet:true});
      return coordinates?{n:this.n,offset:[p.x-x,p.y-y],initial:coordinates,lastCoordinates:coordinates.slice(),started:false,latest:null,frame:null}:null;
    }
    emitPointDrag(point,coordinates,phase,meta){
      const callback=point.n===6?this.surfaceDragCallback:this.pointDragCallback;
      return callback?.call(this,coordinates,{...meta,phase});
    }
    updatePointDrag(event,drag){
      const point=drag.point;
      if(!point.started){point.started=true;this.emitPointDrag(point,point.initial.slice(),'start',point.initialMeta);}
      if(this.drag!==drag)return;
      point.latest={clientX:event.clientX,clientY:event.clientY};
      if(point.frame===null)point.frame=requestAnimationFrame(()=>{point.frame=null;if(this.drag===drag)this.flushPointDrag(drag);});
    }
    traceSurfaceDrag(point,x,y){
      // Preserve the latest valid connected screen path. A fast pointer move
      // must not skip a gray stripe to land on an unrelated violet patch.
      const start=point.acceptedScreen,end=[Math.max(0,Math.min(this.width,x)),Math.max(0,Math.min(this.height,y))];
      const dx=end[0]-start[0],dy=end[1]-start[1],distance=Math.hypot(dx,dy);
      if(distance<1e-8)return null;
      const steps=Math.min(1024,Math.ceil(distance/2)),pixel=(this.camera.top-this.camera.bottom)/(this.camera.zoom*this.height);
      let previous=point.acceptedWorld,previousScreen=start,last=null;
      for(let i=1;i<=steps;i++){
        const fraction=Math.min(1,2*i/distance),screen=[start[0]+fraction*dx,start[1]+fraction*dy];
        const hit=this.attainedSurfaceAt(screen[0],screen[1]);if(!hit)break;
        const world=this.map(hit.point),screenStep=Math.hypot(screen[0]-previousScreen[0],screen[1]-previousScreen[1]);
        // Adjacent facets can meet at a crease. A large depth discontinuity
        // instead exposes a different front sheet and ends this drag path.
        if(v3(world).distanceTo(v3(previous))>Math.max(.025,12*pixel*Math.max(1,screenStep)))break;
        last={hit,screen,world};previous=world;previousScreen=screen;
      }
      return last;
    }
    flushPointDrag(drag){
      const point=drag.point,event=point.latest;point.latest=null;
      if(!event||point.n!==this.n)return;
      const rect=this.canvas.getBoundingClientRect(),x=event.clientX-rect.left+point.offset[0],y=event.clientY-rect.top+point.offset[1];
      if(point.n===6){
        const result=this.traceSurfaceDrag(point,x,y);if(!result)return;
        const coordinates=result.hit.point.slice();
        // The app owns the physical inverse. A failed state reconstruction
        // never advances this accepted position, even on a violet mesh cell.
        if(this.emitPointDrag(point,coordinates,'move',result.hit)===true){
          point.lastCoordinates=coordinates;point.lastMeta=result.hit;
          point.acceptedScreen=result.screen;point.acceptedWorld=result.world;
        }
        return;
      }
      const coordinates=this.regionCoordinates(x,y,{clampToSet:true});
      if(coordinates){point.lastCoordinates=coordinates.slice();this.emitPointDrag(point,coordinates,'move');}
    }
    finishPointDrag(event,drag,cancel){
      const point=drag.point;
      if(point.frame!==null){cancelAnimationFrame(point.frame);point.frame=null;}
      if(cancel){point.latest=null;if(point.started)this.emitPointDrag(point,point.lastCoordinates.slice(),'cancel',point.lastMeta);return;}
      if(!point.started){point.started=true;this.emitPointDrag(point,point.initial.slice(),'start',point.initialMeta);}
      point.latest={clientX:event.clientX,clientY:event.clientY};this.flushPointDrag(drag);
      this.emitPointDrag(point,point.lastCoordinates.slice(),'end',point.lastMeta);
    }
    regionCoordinates(x,y,{clampToSet=false}={}){
      if(this.n!==3&&this.n!==4&&this.n!==5)return null;
      this.raycaster.setFromCamera(new T.Vector2(2*x/this.width-1,1-2*y/this.height),this.camera);
      const p=this.raycaster.ray.intersectPlane(new T.Plane(new T.Vector3(0,0,1),0),new T.Vector3());
      if(!p)return null;
      const [xmin,xmax]=this.ranges[0],tol=1e-10;
      const sx=(p.x/2.65+.5)*(xmax-xmin)+xmin;
      if(this.n===3){
        // A pixel-sized hit band stays easy to select at every zoom level.
        const linePoint=this.project([p.x,0,0]);
        if(!clampToSet&&(Math.abs(linePoint.y-y)>10||sx<xmin-tol||sx>xmax+tol))return null;
        return [Math.max(xmin,Math.min(xmax,sx))];
      }
      const [ymin,ymax]=this.ranges[1],sy=(p.y/2+.5)*(ymax-ymin)+ymin;
      if(!clampToSet&&(sx < xmin-tol || sx > xmax+tol || sy > ymax+tol))return null;
      const a=Math.max(xmin,Math.min(xmax,sx));
      const bottom=this.n===4?a+2:window.SectorConstructors?.n5Lower(a)?.sectors[1];
      if(bottom===undefined||(!clampToSet&&sy < bottom-tol))return null;
      return [a,Math.max(bottom,Math.min(ymax,sy))];
    }
    triangleCoordinates(x,y){return this.n===4?this.regionCoordinates(x,y):null;}
    attainedSurfaceAt(x,y){
      if(this.n!==6||!this.surface.visible)return null;
      this.raycaster.setFromCamera(new T.Vector2(2*x/this.width-1,1-2*y/this.height),this.camera);
      let hit,mesh;
      if(this.occlusionIndex?.closestHit){hit=this.occlusionIndex.closestHit(this.raycaster.ray);mesh=hit&&this.surface.children[hit.source];}
      else{hit=this.raycaster.intersectObjects(this.surface.children,false)[0];mesh=hit?.object;}
      if(!hit||!mesh)return null;
      const patch=mesh.userData.patch,face=hit.faceIndex;
      // The first hit occludes all rear faces, including when it is gray.
      if(!patch?.attained[face])return null;
      const attr=mesh.geometry.getAttribute('position'),offset=face*3;
      const a=new T.Vector3().fromBufferAttribute(attr,offset),b=new T.Vector3().fromBufferAttribute(attr,offset+1),c=new T.Vector3().fromBufferAttribute(attr,offset+2);
      const weights=T.Triangle.getBarycoord(hit.point,a,b,c,new T.Vector3());
      if(!weights)return null;
      const indices=patch.triangles[face],w=[weights.x,weights.y,weights.z];
      const point=[0,1,2].map(k=>indices.reduce((sum,j,i)=>sum+w[i]*patch.positions[j][k],0));
      const anchorPoint=[0,1,2].map(k=>indices.reduce((sum,j)=>sum+patch.positions[j][k],0)/3);
      return {point,anchorPoint,surface:patch.surface,kind:patch.kind,faceIndex:face};
    }
    pick(event,select){if(!this.data)return;const rect=this.canvas.getBoundingClientRect(),x=event.clientX-rect.left,y=event.clientY-rect.top;let best=null,dist=13*13;for(let i=0;i<this.namedWorld.length;i++){const p=this.namedProjections[i],d=(p.x-x)**2+(p.y-y)**2;if(d<dist){dist=d;best=this.namedWorld[i].state;}}if(this.projected.length){for(const p of this.projected){const d=(p.x-x)**2+(p.y-y)**2;if(d<dist){dist=d;best=p.point.state;}}}if(!best){for(const p of this.grayPoints||[]){const q=this.project(p.world),d=(q.x-x)**2+(q.y-y)**2;if(d<dist){dist=d;best=p.state;}}}
      const coordinates=best?null:this.regionCoordinates(x,y);
      const surface=best||coordinates?null:this.attainedSurfaceAt(x,y);
      this.canvas.style.cursor=this.selectedPointAt(x,y)?'grab':best||coordinates||surface?'pointer':this.rotatable?'grab':'default';
      if(select&&surface)this.surfaceCallback?.(surface.point,surface);
      else if(select&&coordinates)this.coordinateCallback?.(coordinates);
      else (select?this.selectCallback:this.hoverCallback)(best,{x,y});}
  }
  // Keep polygon facet boundaries separate from the triangulation used for fill.
  function hullTopology(points){
    const triangles=[],edgeKeys=new Set(),faceKeys=new Set(),eps=1e-7;
    if(points.length===2)return {triangles,edges:[[0,1]]};
    for(let i=0;i<points.length;i++)for(let j=i+1;j<points.length;j++)for(let k=j+1;k<points.length;k++){
      const a=v3(points[i]),normal=v3(points[j]).sub(a).cross(v3(points[k]).sub(a));
      if(normal.length()<eps)continue;normal.normalize();
      const distances=points.map(p=>normal.dot(v3(p).sub(a)));
      if(distances.some(d=>d>eps)&&distances.some(d=>d< -eps))continue;
      const ids=distances.map((d,id)=>Math.abs(d)<eps?id:-1).filter(id=>id>=0),key=ids.join(',');
      if(faceKeys.has(key))continue;faceKeys.add(key);
      const center=new T.Vector3();ids.forEach(id=>center.add(v3(points[id])));center.divideScalar(ids.length);
      if(normal.dot(center)<0)normal.negate();
      const u=v3(points[ids[0]]).sub(center).normalize(),v=normal.clone().cross(u);
      ids.sort((x,y)=>{const xx=v3(points[x]).sub(center),yy=v3(points[y]).sub(center);return Math.atan2(xx.dot(v),xx.dot(u))-Math.atan2(yy.dot(v),yy.dot(u));});
      for(let h=0;h<ids.length;h++)edgeKeys.add([ids[h],ids[(h+1)%ids.length]].sort((a,b)=>a-b).join(','));
      for(let h=1;h<ids.length-1;h++)triangles.push([ids[0],ids[h],ids[h+1]]);
    }
    return {triangles,edges:[...edgeKeys].map(key=>key.split(',').map(Number))};
  }
  class MajoranaView extends BaseView {
    constructor(canvas,overlay){super(canvas,overlay);this.displayRotation=new T.Quaternion();this.span=1.3;this.yaw=.55;this.pitch=.25;this.sphere=new T.Group();this.hull=new T.Group();this.stars=new T.Group();this.symmetryAxesGroup=new T.Group();this.content.add(this.sphere,this.hull,this.stars,this.symmetryAxesGroup);const sphere=new T.Mesh(new T.SphereGeometry(1,48,32),new T.MeshBasicMaterial({color:0xe5ddc9,transparent:true,opacity:.18,depthWrite:false}));this.sphere.add(sphere);for(const lat of [-Math.PI/3,-Math.PI/6,0,Math.PI/6,Math.PI/3]){const points=[];for(let t=0;t<=128;t++){const a=t*Math.PI/64;points.push([Math.cos(lat)*Math.cos(a),Math.sin(lat),Math.cos(lat)*Math.sin(a)]);}this.sphere.add(line(points,0xd5d0c2));}for(let phi=0;phi<Math.PI;phi+=Math.PI/4){const points=[];for(let t=0;t<=128;t++){const a=t*Math.PI/64;points.push([Math.cos(a)*Math.cos(phi),Math.sin(a),Math.cos(a)*Math.sin(phi)]);}this.sphere.add(line(points,0xd5d0c2));}this.resize();}
    setState(coefficients,{preserveOrientation=false}={}){const data=SectorMath.majoranaStars(coefficients);this.majorana=data;this.clearLabels();if(!preserveOrientation)this.displayRotation.identity();this.aligned=false;this.starLabelMarkers=[];this.hull.quaternion.copy(this.displayRotation);this.symmetryAxesGroup.quaternion.copy(this.displayRotation);document.getElementById('alignment-note').hidden=true;disposeGroup(this.stars);disposeGroup(this.hull);this.groups=data.groups;const points=data.groups.map(p=>[p.x,p.z,-p.y]);this.starMarkers=points.map((world,i)=>{const marker={world,baseWorld:world.slice(),radius:4.8,fill:'#ad7549',stroke:'#7f5233',kind:'majorana'};this.markers.push(marker);if(data.groups[i].multiplicity>1){const label=this.addLabel(String(data.groups[i].multiplicity),world,'plot-label',{offset:[12,-10]});this.starLabelMarkers.push({label,baseWorld:world.slice()});}return marker;});if(points.length>=2){const topology=hullTopology(points),array=[];topology.triangles.forEach(f=>f.forEach(i=>array.push(...points[i])));if(array.length){const g=new T.BufferGeometry();g.setAttribute('position',new T.Float32BufferAttribute(array,3));g.computeVertexNormals();this.hull.add(new T.Mesh(g,new T.MeshBasicMaterial({color:0x5f9588,transparent:true,opacity:.16,side:T.DoubleSide,depthWrite:false})));}for(const[a,b]of topology.edges)this.hull.add(line([points[a],points[b]],0x75978e));}

      this.setupSymmetry(data.groups);this.hull.visible=document.getElementById('show-hull').checked;this.resize();return data;}
    setupSymmetry(groups){
      disposeGroup(this.symmetryAxesGroup);this.symmetry=MajoranaSymmetry.analyze(groups);this.axisSegments=[];this.symmetryAxisLabels=[];
      const colors={2:'#789080',3:'#8d77aa',4:'#b17d5c',5:'#657d9c',6:'#a86d4f',infinity:'#9c7555'};
      for(const axis of this.symmetry.axes){
        const d=[axis.direction[0],axis.direction[2],-axis.direction[1]],color=colors[axis.order]||'#7c8490';
        for(const sign of [-1,1]){
          const segment=[[0,0,0],d.map(v=>v*sign*1.12)],solid=line(segment,color),dashed=line(segment,color,true,false),occludedSolid=line(segment,color),occludedDashed=line(segment,color,true,false);
          for(const [object,base,occluded]of [[solid,.75,false],[dashed,.45,false],[occludedSolid,.75,true],[occludedDashed,.45,true]]){object.material.opacity=base*(occluded?.22:1);object.renderOrder=5;object.userData.symmetryAxis=true;object.userData.axisOccluded=occluded;this.symmetryAxesGroup.add(object);}
          this.axisSegments.push({direction:d,order:axis.order,sign,solid,dashed,occludedSolid,occludedDashed});
        }
        const tex=axis.order==='infinity'?'C_{\\infty}':'C_{'+axis.order+'}';const label=this.addLabel('$'+tex+'$',d.map(v=>v*1.18),'plot-label axis-symmetry-label',{offset:[7,-5]});label.el.style.color=color;this.symmetryAxisLabels.push({direction:d,order:axis.order,label,tex});
      }
      const select=document.getElementById('symmetry-order'),previous=select.value;select.replaceChildren();
      for(const [value,text]of [['highest','Highest order'],['all','All orders']]){const o=document.createElement('option');o.value=value;o.textContent=text;select.append(o);}
      const orders=[...new Set(this.symmetry.axes.map(a=>a.order))].sort((a,b)=>(b==='infinity'?Infinity:b)-(a==='infinity'?Infinity:a));this.highestAxisOrder=orders[0];
      for(const order of orders){const o=document.createElement('option');o.value=String(order);o.textContent=order==='infinity'?'Infinite-fold axis':order+'-fold axes';select.append(o);}select.value=[...select.options].some(o=>o.value===previous)?previous:'highest';select.disabled=!orders.length;document.getElementById('align-symmetry').disabled=!orders.length;
      const summary=document.getElementById('symmetry-summary');const counts=orders.map(order=>{const count=this.symmetry.axes.filter(a=>a.order===order).length;return '$C_{'+(order==='infinity'?'\\infty':order)+'}$: '+count+' '+(count===1?'axis':'axes');});let text=counts.join(' · ')||'$C_1$: no nontrivial rotational axis';if(this.symmetry.continuous?.transverseC2)text+=' · Continuous transverse $C_2$ axes';summary.innerHTML=AtlasMath.html(text);
    }
    updateSymmetry(){if(!this.axisSegments)return;const enabled=document.getElementById('show-symmetry').checked,filter=document.getElementById('symmetry-order').value;const visible=order=>enabled&&(filter==='all'||String(order)===filter||(filter==='highest'&&order===this.highestAxisOrder));this.symmetryAxesGroup.visible=enabled;
      const cameraDirection=this.camera.position.clone().sub(this.target).normalize(),extent=1.12;
      const section=(object,part,lo,hi,on)=>{
        object.visible=on&&hi-lo>1e-10;if(!object.visible)return;
        const position=object.geometry.getAttribute('position'),distance=object.geometry.getAttribute('lineDistance');
        for(const [i,r]of [[0,lo],[1,hi]]){position.setXYZ(i,...part.direction.map(v=>v*part.sign*r));distance.setX(i,r);}position.needsUpdate=distance.needsUpdate=true;object.frustumCulled=false;
      };
      for(const part of this.axisSegments){
        const depth=v3(part.direction).applyQuaternion(this.displayRotation).dot(cameraDirection)*part.sign,front=depth>=0,on=visible(part.order);
        const boundary=Math.min(extent,front?1:1/Math.sqrt(Math.max(1e-30,1-depth*depth)));
        section(part.occludedSolid,part,0,boundary,on&&front);section(part.occludedDashed,part,0,boundary,on&&!front);
        section(part.solid,part,boundary,extent,on&&front);section(part.dashed,part,boundary,extent,on&&!front);
      }
      for(const a of this.symmetryAxisLabels){const direction=v3(a.direction).applyQuaternion(this.displayRotation),sign=direction.dot(this.camera.position)>=0?1:-1;a.label.position=direction.multiplyScalar(sign*1.18).toArray();a.label.el.hidden=!visible(a.order);const aligned=this.aligned&&Math.abs(v3(a.direction).applyQuaternion(this.displayRotation).y)>1-1e-8;a.label.el.innerHTML=AtlasMath.html('$'+a.tex+(aligned?'\\parallel z':'')+'$');}
    }
    alignHighest(){const axis=this.symmetry?.axes.find(a=>a.order===this.highestAxisOrder);if(!axis)return;const direction=new T.Vector3(axis.direction[0],axis.direction[2],-axis.direction[1]);this.displayRotation.setFromUnitVectors(direction.normalize(),new T.Vector3(0,1,0));this.hull.quaternion.copy(this.displayRotation);this.symmetryAxesGroup.quaternion.copy(this.displayRotation);this.aligned=true;document.getElementById('alignment-note').hidden=false;this.yaw=.55;this.pitch=.25;this.zoom=1;this.resize();}
    render(){this.updateSymmetry();if(this.starMarkers){this.starMarkers.forEach(marker=>{const position=v3(marker.baseWorld).applyQuaternion(this.displayRotation),front=position.dot(this.camera.position)>0;marker.world=position.toArray();marker.rear=!front;marker.fill=front?'#ad7549':'#ddc6ad';marker.stroke=front?'#7f5233':'#ad8f72';});for(const a of this.starLabelMarkers||[])a.label.position=v3(a.baseWorld).applyQuaternion(this.displayRotation).toArray();}super.render();}

    reset(){this.displayRotation.identity();this.hull.quaternion.identity();this.symmetryAxesGroup.quaternion.identity();this.aligned=false;document.getElementById('alignment-note').hidden=true;this.yaw=.55;this.pitch=.25;this.zoom=1;this.resize();}
  }
  window.SectorView=SectorView;window.MajoranaView=MajoranaView;
})();
