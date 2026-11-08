/* Shared offline Majorana viewer. Coordinates are supplied by the catalogue.
 * API: MajoranaViewer.open({n, label, points:[{x,y,z,multiplicity}], evidence})
 * Symmetry modules are optional; KaTeX is optional. All data remain local. */
(() => {
  'use strict';
  const escape=value=>String(value??'').replace(/[&<>"']/g,c=>({'&':'&amp;','<':'&lt;','>':'&gt;','"':'&quot;',"'":'&#39;'}[c]));
  const math=tex=>{try{return window.katex.renderToString(tex,{throwOnError:false,trust:false,strict:'ignore'});}catch{return escape(tex);}};
  const symmetryColor=order=>({2:'#34786d',3:'#82609b',4:'#a45635',5:'#b38b23',continuous:'#356f94'}[order]||'#526b9b');
  const unknownGroup={symbol:null,name:'Not determined',order:null};
  let dialog=null,returnFocus=null,dispose=()=>{};
  const symmetryCache=new WeakMap();
  const $=id=>dialog.querySelector('#mv-'+id);
  function ensureDialog(){
    if(dialog)return;
    dialog=document.createElement('dialog');
    dialog.id='majorana-viewer';dialog.className='mv-dialog';
    dialog.setAttribute('aria-labelledby','mv-title');dialog.setAttribute('aria-describedby','mv-help');
    dialog.innerHTML=`<div class="mv-heading"><div><p class="mv-eyebrow">Majorana representation</p><h2 id="mv-title">State constellation</h2></div><button id="mv-close" type="button" class="mv-close" aria-label="Close Majorana viewer">Close ×</button></div>
      <p id="mv-label" class="mv-label"></p><p id="mv-summary" class="mv-summary"></p><p id="mv-evidence" class="mv-note"></p>
      <div class="mv-workspace"><div class="mv-sphere-frame"><canvas id="mv-sphere" class="mv-sphere" tabindex="0" aria-label="Interactive Majorana sphere" aria-describedby="mv-help"></canvas></div>
        <div class="mv-sidebar"><div class="mv-controls">
          <div class="mv-turn" role="group" aria-label="Rotate constellation"><button type="button" id="mv-left" aria-label="Rotate left">←</button><button type="button" id="mv-up" aria-label="Rotate up">↑</button><button type="button" id="mv-down" aria-label="Rotate down">↓</button><button type="button" id="mv-right" aria-label="Rotate right">→</button><button type="button" id="mv-reset">Reset view</button></div>
          <button type="button" id="mv-align" aria-label="Align highest-order symmetry axis to z">Align to z</button>
          <button type="button" id="mv-download-svg">Download SVG</button>
          <label><input type="checkbox" id="mv-symmetry">Symmetry axes</label><label for="mv-symmetry-order">Show <select id="mv-symmetry-order"><option value="all">All orders</option></select></label>
          <label><input type="checkbox" id="mv-hull">Convex hull</label>
          <label><input type="checkbox" id="mv-axes">Coordinate axes</label>
        </div><div id="mv-groups" class="mv-groups"></div><div id="mv-symmetry-summary" class="mv-symmetry-summary" aria-live="polite"></div>
        <div class="mv-legend"><span><i class="mv-dot" aria-hidden="true"></i>Front hemisphere</span><span><i class="mv-dot mv-back" aria-hidden="true"></i>Back hemisphere</span><span><b class="mv-centre-key" aria-hidden="true">+</b>Sphere centre</span><span><b class="mv-barycentre-key" aria-hidden="true">◇</b><span id="mv-barycentre-label">Barycentre</span></span><span id="mv-multiplicity-note" hidden style="display:none">Numbers beside points = multiplicities · Unlabelled = 1</span></div>
      </div></div>
      <details class="mv-method"><summary>How to use the viewer</summary><p id="mv-help" class="mv-help">Drag the sphere or use the rotation buttons. With the sphere focused, arrow keys rotate and Home resets. Escape closes. Alignment chooses a displayed coordinate frame with the highest-order symmetry axis along z. Reset restores the default aligned view. Saved coordinates and state coefficients are unchanged.</p></details>
      <details class="mv-method"><summary>About the representation and symmetry</summary><p id="mv-method-copy" class="mv-note"></p><p id="mv-coordinate-note" class="mv-note" hidden></p></details>`;
    document.body.append(dialog);
    $('close').addEventListener('click',close);
    dialog.addEventListener('close',()=>{if(!dialog.open)cleanup();});
    dialog.addEventListener('click',event=>{if(event.target!==dialog)return;const box=dialog.getBoundingClientRect();if(event.clientX<box.left||event.clientX>box.right||event.clientY<box.top||event.clientY>box.bottom)close();});
  }
  function cleanup(){
    dispose();dispose=()=>{};document.body.classList.remove('mv-open');
    if(returnFocus?.isConnected&&!returnFocus.disabled)returnFocus.focus({preventScroll:true});
    returnFocus=null;
  }
  function close(){if(dialog?.open)dialog.close();cleanup();}
  function groupCard(group,title,unit){
    const size=group.continuous?'Continuous':Number.isInteger(group.order)?`${group.order} ${unit}`:'Not determined';
    return `<div class="mv-group"><span class="mv-group-kind">${title}</span><strong>${group.symbol?math(group.symbol):'—'}</strong><span class="mv-group-name">${escape(group.name)}</span><span class="mv-group-size">${escape(size)}</span></div>`;
  }
  function symmetryFor(points){
    if(Array.isArray(points)&&symmetryCache.has(points))return symmetryCache.get(points);
    let result={status:'invalid',axes:[]},rotational=unknownGroup,point=unknownGroup;
    try{result=window.MajoranaSymmetry?.detect(points,{tolerance:1e-6})||result;}catch{}
    try{rotational=window.MajoranaGroups?.classify(result)||rotational;point=window.MajoranaPointGroups?.classify(points,result,rotational)||point;}catch{}
    const groups={result,rotational,point};
    if(Array.isArray(points))symmetryCache.set(points,groups);
    return groups;
  }
  function groupLabel(state){
    if(!state?.points?.length)return '';
    return symmetryFor(state.points).point.label||'Not determined';
  }
  function groupSummary(state,compact=false){
    if(!state?.points?.length)return '';
    const {point,rotational}=symmetryFor(state.points);
    const item=(group,label)=>`<span class="representative-group" title="${escape(group.name)} · Numerically identified"><span>${label}</span><strong>${group.symbol?math(group.symbol):'Not determined'}</strong></span>`;
    return `<div class="representative-symmetry${compact?' compact':''}" data-state-id="${escape(state.id||'')}" aria-label="Symmetry of this representative">${item(point,'Point group')}${item(rotational,'Rotations')}${compact?'':'<span class="representative-symmetry-note">Numerical geometry · This representative</span>'}</div>`;
  }
  function prepareSymmetry(points){
    const {result,rotational,point}=symmetryFor(points);
    $('groups').innerHTML=groupCard(point,'Point group','symmetries')+groupCard(rotational,'Rotational subgroup','rotations')+'<p class="mv-group-evidence">Numerical geometry · This representative only</p>';
    const counts=new Map();for(const axis of result.axes){const key=axis.order==null?'continuous':String(axis.order);counts.set(key,(counts.get(key)||0)+1);}
    const keys=[...counts.keys()].sort((a,b)=>a==='continuous'?-1:b==='continuous'?1:Number(b)-Number(a));
    if(result.continuousFamilies?.length&&!keys.includes('2'))keys.push('2');
    $('symmetry-order').innerHTML='<option value="all">All orders</option>'+keys.map(key=>`<option value="${escape(key)}">${key==='continuous'?'Continuous axis':escape(key)+'-fold axes'}</option>`).join('');
    $('symmetry-order').value=result.axes.length>8?keys[0]:'all';
    const available=Boolean(result.axes.length||result.continuousFamilies?.length);
    $('symmetry').checked=available;$('symmetry').disabled=!available;$('symmetry-order').disabled=!available;
    const family=result.continuousFamilies?.length?`<span class="mv-symmetry-family">Every axis in the highlighted perpendicular plane is ${math('C_2')}.</span>`:'';
    $('symmetry-summary').innerHTML=available?keys.filter(key=>counts.has(key)).map(key=>`<span class="mv-symmetry-key"><i style="background:${symmetryColor(key)}" aria-hidden="true"></i>${math('C_{'+(key==='continuous'?'\\infty':key)+'}')}<span>${counts.get(key)} ${counts.get(key)===1?'axis':'axes'}</span></span>`).join('')+family:result.status==='invalid'?'Symmetry classification unavailable.':'No nontrivial rotational symmetry detected.';
    $('method-copy').innerHTML=`The ${math('n')} Majorana points describe the selected pure symmetric state. In the normalized Dicke basis, <span class="mv-formula">${math('P(X,Y)=\\sum_{r=0}^{n}(-1)^r\\sqrt{\\binom{n}{r}}c_rX^{n-r}Y^r')}.</span> A root ${math('z=Y/X')} corresponds to the spinor ratio ${math('b/a=1/z')}: zero is the south pole and infinity the north pole. Displayed multiplicities follow the coordinate record; numerical grouping alone does not certify an exact degeneracy. Coordinates and symmetry labels are numerical renderings; they do not establish an anticoherence or optimality claim. ${math('C_m')} marks a rotation axis of order ${math('m')}; opposite ends denote one axis. Symmetries preserve the whole constellation, including multiplicities. The full point group includes any reflections and inversion; the rotational subgroup contains proper rotations only. These operations are checked geometrically with unit-sphere distance tolerance ${math('10^{-6}')}.`;
    return result;
  }
  function open(state){
    if(!state||!Number.isInteger(state.n)||state.n<1||!Array.isArray(state.points)||!state.points.length)return false;
    const points=state.points.map(p=>({x:Number(p.x),y:Number(p.y),z:Number(p.z),multiplicity:p.multiplicity??1}));
    if(points.some(p=>![p.x,p.y,p.z].every(Number.isFinite)||!Number.isInteger(p.multiplicity)||p.multiplicity<1||Math.abs(Math.hypot(p.x,p.y,p.z)-1)>1e-5)||points.reduce((sum,p)=>sum+p.multiplicity,0)!==state.n)return false;
    ensureDialog();dispose();dispose=()=>{};
    if(!dialog.open){returnFocus=document.activeElement;dialog.showModal();}
    document.body.classList.add('mv-open');
    $('title').textContent=`${state.n} qubits · Majorana sphere`;
    $('label').textContent=state.label||'Selected representative';
    $('summary').textContent=`${state.n} points · ${points.length} distinct ${points.length===1?'position':'positions'}`;
    const repeated=points.some(point=>point.multiplicity>1);
    $('multiplicity-note').hidden=!repeated;
    $('multiplicity-note').setAttribute('style',repeated?'':'display:none');
    $('barycentre-label').textContent=repeated?'Barycentre (with multiplicities)':'Barycentre';
    let caption=state.evidence||'Numerical coordinates of the selected representative.';
    const coordinateNotes=[];
    for(const note of ['Coordinates are a numerical rendering of the saved state.','Coordinates use the saved decimal coefficients; coincident finite roots are not inferred by clustering.']){
      if(caption.includes(note)){coordinateNotes.push(note);caption=caption.replaceAll(note,'');}
    }
    $('evidence').textContent=caption.replace(/\s+/g,' ').trim();
    $('evidence').hidden=!$('evidence').textContent;
    $('coordinate-note').textContent=coordinateNotes.join(' ');
    $('coordinate-note').hidden=!coordinateNotes.length;
    $('sphere').setAttribute('aria-label',`Majorana sphere: ${state.n} qubits, ${points.length} distinct positions. ${state.label||'Selected representative'}`);
    $('axes').checked=false;$('hull').checked=true;
    const symmetry=prepareSymmetry(state.points);
    dispose=attachSphere($('sphere'),points,{download:$('download-svg'),exportLabel:state.label||state.id||'Majorana constellation',reset:$('reset'),align:$('align'),hull:$('hull'),axes:$('axes'),symmetry,symmetryToggle:$('symmetry'),symmetryOrder:$('symmetry-order'),turns:[[$('left'),-.2,0],[$('right'),.2,0],[$('up'),0,-.2],[$('down'),0,.2]]});
    dialog.scrollTop=0;$('sphere').focus({preventScroll:true});return true;
  }
  // Axes are unoriented. Prefer the northern end, then y/x at the equator.
  // Equal orders retain the detector's stable ordering; a continuous axis wins.
  function principalAxis(symmetry){
    let chosen=null,best=1;
    for(const a of symmetry?.axes||[]){
      const rank=a.order==null?Infinity:Number(a.order),v=a.direction;
      if(!Array.isArray(v)||v.length!==3||!v.every(Number.isFinite)||Math.hypot(...v)<1e-12||rank<=best)continue;
      const norm=Math.hypot(...v),direction=v.map(x=>x/norm);
      const sign=[direction[2],direction[1],direction[0]].find(x=>Math.abs(x)>1e-12)<0?-1:1;
      chosen={order:a.order,direction:direction.map(x=>sign*x)};best=rank;
    }
    return chosen;
  }
  function axisToZ(direction){
    const norm=Math.hypot(...direction),[x,y,z]=direction.map(v=>v/norm);
    // Shortest proper rotation from direction to +z; handle the antipode.
    if(z<-1+1e-12)return [0,1,0,0];
    const q=[1+z,y,-x,0],length=Math.hypot(...q);
    return q.map(v=>v/length);
  }
  function attachSphere(canvas,points,controls={}) {
    const canvasContext=canvas.getContext('2d');if(!canvasContext)return ()=>{};
    const dense=points.length>30;
    const weight=points.reduce((sum,p)=>sum+(p.multiplicity??1),0);
    const barycentre=Object.fromEntries(['x','y','z'].map(axis=>[axis,points.reduce((sum,p)=>sum+(p.multiplicity??1)*p[axis],0)/weight]));
    const hull=window.MajoranaHull?.build(points);
    if(controls.hull){controls.hull.disabled=!hull?.faces.length;controls.hull.checked=Boolean(hull?.faces.length);controls.hull.title=hull?.faces.length?'Convex hull of the distinct Majorana positions':'At least three non-collinear positions are needed';}
    const multiply=(a,b)=>[a[0]*b[0]-a[1]*b[1]-a[2]*b[2]-a[3]*b[3],a[0]*b[1]+a[1]*b[0]+a[2]*b[3]-a[3]*b[2],a[0]*b[2]-a[1]*b[3]+a[2]*b[0]+a[3]*b[1],a[0]*b[3]+a[1]*b[2]-a[2]*b[1]+a[3]*b[0]];
    const axis=(x,y,z,angle)=>[Math.cos(angle/2),x*Math.sin(angle/2),y*Math.sin(angle/2),z*Math.sin(angle/2)];
    const principal=principalAxis(controls.symmetry);
    const initial=principal?axis(1,0,0,-Math.PI/2+.28):multiply(axis(1,0,0,-.28),axis(0,1,0,.42));
    const initialFrame=principal?axisToZ(principal.direction):[1,0,0,0];
    // Choose one end in the aligned state frame, independent of the camera.
    const axisDirections=new Map();
    for(const a of controls.symmetry?.axes||[]){
      const v=multiply(multiply(initialFrame,[0,...a.direction]),[initialFrame[0],-initialFrame[1],-initialFrame[2],-initialFrame[3]]);
      const sign=[v[3],v[2],v[1]].find(x=>Math.abs(x)>1e-12)<0?-1:1;
      axisDirections.set(a.direction,a.direction.map(x=>sign*x));
    }
    let orientation=[...initial],coordinateFrame=[...initialFrame],drag=null,frame=0,disposed=false;
    const listeners=[];
    const listen=(element,event,handler)=>{if(element){element.addEventListener(event,handler);listeners.push(()=>element.removeEventListener(event,handler));}};
    const rotate=(point,inFrame=true)=>{
      const q=inFrame?multiply(orientation,coordinateFrame):orientation;
      const v=multiply(multiply(q,[0,point.x,point.y,point.z]),[q[0],-q[1],-q[2],-q[3]]);
      return {x:v[1],y:v[2],z:v[3],multiplicity:point.multiplicity};
    };
    function turn(dx,dy){orientation=multiply(multiply(axis(0,1,0,dx),axis(1,0,0,dy)),orientation);const norm=Math.hypot(...orientation);orientation=orientation.map(v=>v/norm);schedule();}
    function schedule(){if(!frame&&!disposed)frame=requestAnimationFrame(()=>{frame=0;draw();});}
    function draw(context=canvasContext,exportBox=null){
      const box=exportBox||canvas.getBoundingClientRect();if(!box.width||!box.height)return;
      const dpr=exportBox?1:Math.min(window.devicePixelRatio||1,2);if(!exportBox){canvas.width=Math.round(box.width*dpr);canvas.height=Math.round(box.height*dpr);}context.setTransform(dpr,0,0,dpr,0,0);
      const width=box.width,height=box.height,cx=width/2,cy=height/2,radius=Math.min(width*.365,height*.405)*.95;
      context.clearRect(0,0,width,height);
      context.save();context.translate(cx,cy+radius*.95);context.scale(1,.12);const shadow=context.createRadialGradient(0,0,1,0,0,radius);shadow.addColorStop(0,'rgba(83,61,40,.12)');shadow.addColorStop(1,'rgba(83,61,40,0)');context.fillStyle=shadow;context.beginPath();context.arc(0,0,radius,0,Math.PI*2);context.fill();context.restore();
      const fill=context.createRadialGradient(cx-radius*.32,cy-radius*.38,radius*.06,cx,cy,radius);fill.addColorStop(0,'#fffefa');fill.addColorStop(.75,'#f3eee4');fill.addColorStop(1,'#e8dfd1');context.fillStyle=fill;context.beginPath();context.arc(cx,cy,radius,0,Math.PI*2);context.fill();
      const project=p=>({x:cx+radius*p.x,y:cy-radius*p.y,z:p.z});
      const circles=[];
      for(const latitude of [-Math.PI/3,-Math.PI/6,0,Math.PI/6,Math.PI/3])circles.push(Array.from({length:181},(_,i)=>{const a=i*Math.PI/90;return rotate({x:Math.cos(latitude)*Math.cos(a),y:Math.sin(latitude),z:Math.cos(latitude)*Math.sin(a)},false);}));
      for(const longitude of [0,Math.PI/4,Math.PI/2,3*Math.PI/4])circles.push(Array.from({length:181},(_,i)=>{const a=i*Math.PI/90;return rotate({x:Math.cos(a)*Math.cos(longitude),y:Math.sin(a),z:Math.cos(a)*Math.sin(longitude)},false);}));
      for(const front of [false,true]){
        context.strokeStyle=front?'rgba(133,112,88,.28)':'rgba(133,112,88,.16)';context.lineWidth=.8;context.setLineDash(front?[]:[2.5,4]);
        for(const circle of circles){context.beginPath();for(let i=1;i<circle.length;i++){const a=circle[i-1],b=circle[i];if((.5*(a.z+b.z)>=0)!==front)continue;const pa=project(a),pb=project(b);context.moveTo(pa.x,pa.y);context.lineTo(pb.x,pb.y);}context.stroke();}
      }
      context.setLineDash([]);context.strokeStyle='#c8b9a6';context.lineWidth=1;context.beginPath();context.arc(cx,cy,radius,0,Math.PI*2);context.stroke();
      if(controls.hull?.checked&&hull?.faces.length){
        const vertices=hull.vertices.map(point=>rotate(point));
        const faces=hull.faces.map(ids=>({ids,depth:ids.reduce((s,i)=>s+vertices[i].z,0)/ids.length})).sort((a,b)=>a.depth-b.depth);
        context.save();context.setLineDash([]);
        for(const face of faces){
          context.fillStyle=face.depth>=0?'rgba(58,126,119,.16)':'rgba(58,126,119,.07)';
          context.beginPath();face.ids.forEach((i,j)=>{const p=project(vertices[i]);if(j)context.lineTo(p.x,p.y);else context.moveTo(p.x,p.y);});context.closePath();context.fill();
        }
        context.lineWidth=1.1;
        for(const [a,b] of hull.edges){
          context.strokeStyle=(vertices[a].z+vertices[b].z)/2>=0?'rgba(45,111,105,.7)':'rgba(45,111,105,.28)';
          const p=project(vertices[a]),q=project(vertices[b]);context.beginPath();context.moveTo(p.x,p.y);context.lineTo(q.x,q.y);context.stroke();
        }
        context.restore();
      }
      const axisLabels=[];
      function axisLine(vector,label,color,unoriented=true,order=null){
        // Keep the saved end of an axis: changing its sign with the viewing
        // depth makes the label jump to the opposite end during rotation.
        const direction=unoriented?(axisDirections.get(vector)||vector):vector;
        const endpoint=rotate({x:direction[0],y:direction[1],z:direction[2]},unoriented);
        const p=project(endpoint),dx=p.x-cx,dy=p.y-cy,length=Math.hypot(dx,dy);
        // Extend in 3D, then project. An axis pointing towards the camera
        // must foreshorten with the constellation, not span a full diameter.
        const scale=unoriented?1.08:1.05,extent=length*scale;
        context.save();context.strokeStyle=color;context.lineWidth=unoriented?1.8:1;context.globalAlpha=.8;
        for(const sign of [-1,1]){
          const depth=sign*endpoint.z,rear=depth<0;
          context.setLineDash(rear?[4,4]:[]);
          // A rear ray stays behind the sphere until it clears the silhouette.
          // Coordinate guides keep their original uninterrupted appearance.
          const boundary=unoriented?(rear?1/Math.sqrt(Math.max(1e-30,1-depth*depth)):1):0;
          const cut=Math.min(scale,boundary),spans=unoriented?[[0,cut,.22],[cut,scale,1]]:[[0,scale,1]];
          for(const [lo,hi,fade]of spans){
            if(hi-lo<1e-12)continue;
            context.globalAlpha=.8*fade;context.lineDashOffset=-length*lo;
            context.beginPath();context.moveTo(cx+sign*dx*lo,cy+sign*dy*lo);context.lineTo(cx+sign*dx*hi,cy+sign*dy*hi);context.stroke();
          }
        }
        const endOn=unoriented?Math.max(0,1-length/10):0;
        if(endOn>0){context.globalAlpha=.8*.22*endOn;context.setLineDash([]);context.lineDashOffset=0;context.beginPath();context.arc(cx,cy,4,0,2*Math.PI);context.stroke();}
        context.restore();
        // A continuous screen-space offset stays readable near an end-on
        // axis without dividing by its almost-zero projected length.
        const padding=unoriented?12:28,offset=padding/(extent+padding);
        axisLabels.push({label,color,order,x:Math.max(24,Math.min(width-24,cx+dx*scale*(1+offset))),y:Math.max(13,Math.min(height-13,cy+dy*scale*(1+offset)-padding*offset))});
      }
      if(controls.axes?.checked)for(const [name,vector,color] of [['x',[1,0,0],'#a15435'],['y',[0,1,0],'#41776e'],['z',[0,0,1],'#76648a']])axisLine(vector,name,color,false);
      if(controls.symmetryToggle?controls.symmetryToggle.checked:Boolean(controls.symmetry?.axes?.length)){
        const order=controls.symmetryOrder?.value||'all';
        for(const a of controls.symmetry.axes){const key=a.order==null?'continuous':String(a.order);if(order!=='all'&&order!==key)continue;const label=key==='continuous'?'C∞':'C'+key.replace(/\d/g,n=>'₀₁₂₃₄₅₆₇₈₉'[Number(n)]);axisLine(a.direction,label,symmetryColor(key),true,key);}
        if(controls.symmetry.continuousFamilies?.length&&(order==='all'||order==='2')){
          const normal=controls.symmetry.axes.find(a=>a.order==null)?.direction;
          if(normal){
            const helper=Math.abs(normal[0])<.8?[1,0,0]:[0,1,0];
            const u=[normal[1]*helper[2]-normal[2]*helper[1],normal[2]*helper[0]-normal[0]*helper[2],normal[0]*helper[1]-normal[1]*helper[0]],length=Math.hypot(...u);for(let i=0;i<3;i++)u[i]/=length;
            const v=[normal[1]*u[2]-normal[2]*u[1],normal[2]*u[0]-normal[0]*u[2],normal[0]*u[1]-normal[1]*u[0]];
            context.save();context.strokeStyle=symmetryColor('2');context.lineWidth=2.2;
            for(let i=0;i<180;i++){const point=angle=>rotate({x:u[0]*Math.cos(angle)+v[0]*Math.sin(angle),y:u[1]*Math.cos(angle)+v[1]*Math.sin(angle),z:u[2]*Math.cos(angle)+v[2]*Math.sin(angle)}),a=point(i*Math.PI/90),b=point((i+1)*Math.PI/90),pa=project(a),pb=project(b);context.globalAlpha=a.z<0?.22:1;context.setLineDash(a.z<0?[3,4]:[]);context.beginPath();context.moveTo(pa.x,pa.y);context.lineTo(pb.x,pb.y);context.stroke();}context.restore();
          }
        }
      }
      const ordered=points.map(point=>rotate(point)).sort((a,b)=>a.z-b.z);
      for(const point of ordered){
        const p=project(point),front=point.z>=0,labelled=Number(point.multiplicity)>1,nodeRadius=4.2+(point.z+1)*.6;
        context.save();context.globalAlpha=front?1:.85;context.lineWidth=front?1.3:1.4;context.setLineDash([]);
        context.fillStyle=front?'#a45a34':'#ddb392';context.strokeStyle=front?'#814321':'#8a4d2c';
        if(front){context.shadowColor='rgba(80,43,23,.19)';context.shadowBlur=3;context.shadowOffsetY=1;}
        context.beginPath();context.arc(p.x,p.y,nodeRadius,0,Math.PI*2);context.fill();context.stroke();context.shadowBlur=0;context.shadowOffsetY=0;context.setLineDash([]);
        if(labelled){const text=String(point.multiplicity),x=p.x+nodeRadius+5,y=p.y-nodeRadius-3;context.textAlign='left';context.textBaseline='middle';context.font='600 12px system-ui,-apple-system,sans-serif';context.lineJoin='round';context.lineWidth=2.5;context.strokeStyle='rgba(255,253,247,.85)';context.strokeText(text,x,y);context.fillStyle=front?'#703b20':'#744b33';context.fillText(text,x,y);}context.restore();
      }
      for(const label of axisLabels){if(window.MajoranaLatex&&(label.order!=null||/^[xyz]$/.test(label.label))){window.MajoranaLatex.draw(context,label.order??label.label,label.x,label.y,16,label.color);continue;}context.save();context.font='600 13px system-ui,-apple-system,sans-serif';context.textAlign='center';context.textBaseline='middle';context.lineJoin='round';context.lineWidth=3;context.strokeStyle='rgba(255,253,247,.82)';context.strokeText(label.label,label.x,label.y);context.fillStyle=label.color;context.fillText(label.label,label.x,label.y);context.restore();}
      // Draw both markers on top of the translucent geometry. The cross
      // remains visible inside the diamond when their projections coincide.
      context.save();context.globalAlpha=1;context.shadowBlur=0;context.setLineDash([]);context.lineJoin='round';
      context.beginPath();context.moveTo(cx-4,cy);context.lineTo(cx+4,cy);context.moveTo(cx,cy-4);context.lineTo(cx,cy+4);
      context.strokeStyle='#fffdf7';context.lineWidth=4;context.stroke();context.strokeStyle='#35424a';context.lineWidth=1.7;context.stroke();
      const bary=project(rotate(barycentre));
      context.beginPath();context.moveTo(bary.x,bary.y-6);context.lineTo(bary.x+6,bary.y);context.lineTo(bary.x,bary.y+6);context.lineTo(bary.x-6,bary.y);context.closePath();
      context.strokeStyle='#fffdf7';context.lineWidth=4;context.stroke();context.strokeStyle='#a12b75';context.lineWidth=1.8;context.stroke();context.restore();
    }
    const down=event=>{if(event.button!==0&&event.pointerType==='mouse')return;canvas.focus({preventScroll:true});drag={id:event.pointerId,x:event.clientX,y:event.clientY};canvas.setPointerCapture(event.pointerId);canvas.classList.add('dragging');};
    const move=event=>{if(!drag||drag.id!==event.pointerId)return;const dx=event.clientX-drag.x,dy=event.clientY-drag.y;drag.x=event.clientX;drag.y=event.clientY;turn(dx*.008,dy*.008);};
    const up=event=>{if(drag&&drag.id===event.pointerId){drag=null;canvas.classList.remove('dragging');if(canvas.hasPointerCapture(event.pointerId))canvas.releasePointerCapture(event.pointerId);}};
    const reset=()=>{orientation=[...initial];coordinateFrame=[...initialFrame];schedule();};
    if(controls.align){
      controls.align.disabled=!principal;
      controls.align.title=principal?`Align ${principal.order==null?'the continuous':`a ${principal.order}-fold`} symmetry axis with z in the displayed frame`:'No rotational symmetry axis detected';
    }
    listen(controls.align,'click',()=>{
      if(!principal)return;
      coordinateFrame=axisToZ(principal.direction);
      // z is almost vertical, slightly tilted towards the viewer for depth.
      orientation=axis(1,0,0,-Math.PI/2+.28);
      if(controls.symmetryToggle)controls.symmetryToggle.checked=true;
      if(controls.symmetryOrder)controls.symmetryOrder.value=principal.order==null?'continuous':String(principal.order);
      schedule();
    });
    if(controls.download)controls.download.disabled=!window.MajoranaSVG;
    listen(controls.download,'click',()=>{
      const box=canvas.getBoundingClientRect();if(!box.width||!box.height)return;
      const svg=new window.MajoranaSVG.Context();draw(svg,box);
      const n=points.reduce((sum,p)=>sum+p.multiplicity,0);
      const label=controls.exportLabel||'Majorana constellation';
      const filename=`majorana-n${n}-${label.normalize('NFKD').replace(/[^a-zA-Z0-9]+/g,'-').replace(/^-|-$/g,'').slice(0,80)||'state'}.svg`;
      window.MajoranaSVG.download(svg.toSVG(box.width,box.height,`${n} qubits · ${label}`,'Majorana constellation in the current viewing orientation. Visible hull and axes follow the viewer settings. Numbers denote multiplicities greater than one. The dark cross marks the sphere centre; the magenta diamond marks the barycentre, averaging all Majorana vectors with their multiplicities.'),filename);
    });
    const key=event=>{const step=event.shiftKey ? 0.25 : 0.08;const directions={ArrowLeft:[-step,0],ArrowRight:[step,0],ArrowUp:[0,-step],ArrowDown:[0,step]};if(event.key==='Home'){event.preventDefault();reset();}else if(directions[event.key]){event.preventDefault();turn(...directions[event.key]);}};
    listen(canvas,'pointerdown',down);listen(canvas,'pointermove',move);listen(canvas,'pointerup',up);listen(canvas,'pointercancel',up);listen(canvas,'lostpointercapture',up);listen(canvas,'keydown',key);listen(controls.reset,'click',reset);
    listen(controls.hull,'change',schedule);listen(controls.axes,'change',schedule);
    listen(controls.symmetryToggle,'change',schedule);listen(controls.symmetryOrder,'change',schedule);
    for(const [button,dx,dy] of controls.turns||[])listen(button,'click',()=>turn(dx,dy));
    const observer=typeof ResizeObserver!=='undefined'?new ResizeObserver(schedule):null;observer?.observe(canvas);window.addEventListener('resize',schedule);schedule();
    return ()=>{disposed=true;observer?.disconnect();window.removeEventListener('resize',schedule);if(frame)cancelAnimationFrame(frame);listeners.forEach(remove=>remove());if(drag&&canvas.hasPointerCapture(drag.id))canvas.releasePointerCapture(drag.id);drag=null;canvas.classList.remove('dragging');};
  }

  window.MajoranaViewer=Object.freeze({open,close,symmetryFor,groupSummary,groupLabel});
})();
