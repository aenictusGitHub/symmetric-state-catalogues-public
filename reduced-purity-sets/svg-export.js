/* Stand-alone vector export of the currently displayed orthographic view. */
(function(global){
  'use strict';
  const T=global.THREE;
  const NS='http://www.w3.org/2000/svg', XHTML='http://www.w3.org/1999/xhtml';
  const escape=value=>String(value).replace(/[&<>"']/g,c=>({'&':'&amp;','<':'&lt;','>':'&gt;','"':'&quot;',"'":'&#39;'}[c]));
  const number=value=>Number(value.toFixed(3)).toString();
  const finite=p=>Number.isFinite(p.x)&&Number.isFinite(p.y)&&Number.isFinite(p.z);
  const inside=(object,ancestor)=>{for(let p=object;p;p=p.parent)if(p===ancestor)return true;return false;};
  const project=(position,camera,width,height)=>{const p=position.clone().project(camera);return {x:(p.x+1)*width/2,y:(1-p.y)*height/2,z:p.z};};
  const colorHex=color=>'#'+color.getHexString();
  const pathPoints=points=>points.map((p,i)=>(i?'L':'M')+number(p.x)+' '+number(p.y)).join('');
  function triangleColor(mesh,a,b,c){
    const material=Array.isArray(mesh.material)?mesh.material[0]:mesh.material;
    const color=material.color?.clone()||new T.Color(0xffffff);
    const attribute=mesh.geometry.getAttribute('color');
    if(material.vertexColors&&attribute){
      const mix=new T.Color(
        (attribute.getX(a)+attribute.getX(b)+attribute.getX(c))/3,
        (attribute.getY(a)+attribute.getY(b)+attribute.getY(c))/3,
        (attribute.getZ(a)+attribute.getZ(b)+attribute.getZ(c))/3);
      color.multiply(mix);
    }
    return colorHex(color);
  }
  function collect(view){
    const triangles=[],lines=[];let globe=null;
    const width=view.width,height=view.height,camera=view.camera;
    view.scene.updateMatrixWorld(true);camera.updateMatrixWorld(true);
    view.content.traverseVisible(object=>{
      const geometry=object.geometry,material=Array.isArray(object.material)?object.material[0]:object.material;
      if(!geometry||!material||material.visible===false)return;
      const positions=geometry.getAttribute('position');if(!positions)return;
      if(object.isMesh){
        if(geometry.type==='RingGeometry')return;
        if(view.sphere&&inside(object,view.sphere)&&geometry.type==='SphereGeometry'){
          const center=object.getWorldPosition(new T.Vector3());
          const scale=object.getWorldScale(new T.Vector3());
          const right=new T.Vector3(1,0,0).applyQuaternion(camera.quaternion)
            .multiplyScalar((geometry.parameters?.radius||1)*Math.max(scale.x,scale.y,scale.z));
          const p=project(center,camera,width,height),q=project(center.clone().add(right),camera,width,height);
          globe={x:p.x,y:p.y,r:Math.hypot(q.x-p.x,q.y-p.y),fill:colorHex(material.color),opacity:material.opacity};
          return;
        }
        const index=geometry.getIndex(),count=index?index.count:positions.count;
        for(let i=0;i+2<count;i+=3){
          const ids=[0,1,2].map(k=>index?index.getX(i+k):i+k);
          const world=ids.map(k=>new T.Vector3().fromBufferAttribute(positions,k).applyMatrix4(object.matrixWorld));
          const points=world.map(p=>project(p,camera,width,height));
          if(!points.every(finite)||points.every(p=>p.z>1)||points.every(p=>p.z< -1))continue;
          const area=(points[1].x-points[0].x)*(points[2].y-points[0].y)-(points[1].y-points[0].y)*(points[2].x-points[0].x);
          if(Math.abs(area)<1e-8)continue;
          triangles.push({points,world,depth:points.reduce((s,p)=>s+p.z,0)/3,
            fill:triangleColor(object,...ids),opacity:material.opacity??1,
            opaque:!material.transparent&&(material.opacity??1)>.98,
            hull:!!view.hull&&inside(object,view.hull)});
        }
      }else if(object.isLine){
        // The displayed body has a second, depth-GREATER pass for hidden dashes.
        // Reconstruct its visible/hidden runs once instead of exporting duplicates.
        if(material.depthFunc===T.GreaterDepth)return;
        const index=geometry.getIndex(),count=index?index.count:positions.count;
        const points=Array.from({length:count},(_,i)=>{
          const j=index?index.getX(i):i;
          return project(new T.Vector3().fromBufferAttribute(positions,j).applyMatrix4(object.matrixWorld),camera,width,height);
        });
        if(points.length<2||!points.every(finite))return;
        const bodyEdge=!!view.edges&&inside(object,view.edges);
        const grid=!!view.axisGroup&&inside(object,view.axisGroup)||!!view.sphere&&inside(object,view.sphere);
        const pixelsPerWorld=height*camera.zoom/(camera.top-camera.bottom),axis=!!object.userData?.symmetryAxis,distances=geometry.getAttribute('lineDistance');
        const distanceStart=distances?.getX(0)||0,distanceEnd=distances?.getX(count-1)||0;
        const axisScale=axis&&distanceEnd>distanceStart?Math.hypot(points.at(-1).x-points[0].x,points.at(-1).y-points[0].y)/(distanceEnd-distanceStart):pixelsPerWorld;
        lines.push({points,axis,dashOffset:axis?-distanceStart*axisScale:0,color:colorHex(material.color||new T.Color(0x777777)),opacity:material.opacity??1,
          width:material.linewidth||1,dashed:material.isLineDashedMaterial||material.type==='LineDashedMaterial',
          dash:bodyEdge?[6,5]:[(material.dashSize||6/pixelsPerWorld)*axisScale/(material.scale||1),(material.gapSize||5/pixelsPerWorld)*axisScale/(material.scale||1)],
          bodyEdge,grid,segments:!!object.isLineSegments,loop:!!object.isLineLoop,
          priority:grid?0:bodyEdge?2:1});
      }
    });
    return {triangles,lines,globe};
  }
  function hullPolygons(triangles){
    const groups=new Map(),others=[];
    for(const triangle of triangles){
      if(!triangle.hull){others.push(triangle);continue;}
      const [a,b,c]=triangle.world,normal=b.clone().sub(a).cross(c.clone().sub(a)).normalize();
      const component=[normal.x,normal.y,normal.z].find(v=>Math.abs(v)>1e-7)||1;
      if(component<0)normal.negate();
      const key=[normal.x,normal.y,normal.z,normal.dot(a)].map(v=>Math.round(v*1e6)).join(',')+triangle.fill;
      if(!groups.has(key))groups.set(key,{...triangle,normal,vertices:new Map()});
      const group=groups.get(key);
      triangle.world.forEach((p,i)=>group.vertices.set([p.x,p.y,p.z].map(v=>Math.round(v*1e7)).join(','),{world:p,point:triangle.points[i]}));
    }
    for(const group of groups.values()){
      const vertices=[...group.vertices.values()],center=new T.Vector3();
      vertices.forEach(v=>center.add(v.world));center.divideScalar(vertices.length);
      const u=vertices[0].world.clone().sub(center).normalize(),v=group.normal.clone().cross(u);
      vertices.sort((a,b)=>{
        const p=a.world.clone().sub(center),q=b.world.clone().sub(center);
        return Math.atan2(p.dot(v),p.dot(u))-Math.atan2(q.dot(v),q.dot(u));
      });
      others.push({...group,points:vertices.map(v=>v.point),depth:vertices.reduce((s,v)=>s+v.point.z,0)/vertices.length});
    }
    return others.sort((a,b)=>b.depth-a.depth);
  }
  function depthLookup(triangles,width,height){
    // Numerical visibility calculation only; no raster data is put in the SVG.
    const factor=Math.min(2,1024/Math.max(width,height)),w=Math.max(1,Math.ceil(width*factor)),h=Math.max(1,Math.ceil(height*factor));
    const depth=new Float32Array(w*h);depth.fill(Infinity);
    for(const triangle of triangles){
      if(!triangle.opaque)continue;
      const [a,b,c]=triangle.points.map(p=>({x:p.x*factor,y:p.y*factor,z:p.z}));
      const denominator=(b.y-c.y)*(a.x-c.x)+(c.x-b.x)*(a.y-c.y);
      if(Math.abs(denominator)<1e-10)continue;
      const xmin=Math.max(0,Math.floor(Math.min(a.x,b.x,c.x))),xmax=Math.min(w-1,Math.ceil(Math.max(a.x,b.x,c.x)));
      const ymin=Math.max(0,Math.floor(Math.min(a.y,b.y,c.y))),ymax=Math.min(h-1,Math.ceil(Math.max(a.y,b.y,c.y)));
      for(let y=ymin;y<=ymax;y++)for(let x=xmin;x<=xmax;x++){
        const px=x+.5,py=y+.5;
        const u=((b.y-c.y)*(px-c.x)+(c.x-b.x)*(py-c.y))/denominator;
        const v=((c.y-a.y)*(px-c.x)+(a.x-c.x)*(py-c.y))/denominator,wgt=1-u-v;
        if(u< -1e-6||v< -1e-6||wgt< -1e-6)continue;
        const z=u*a.z+v*b.z+wgt*c.z,index=y*w+x;
        if(z<depth[index])depth[index]=z;
      }
    }
    return point=>{
      const x=Math.floor(point.x*factor),y=Math.floor(point.y*factor);
      return x>=0&&x<w&&y>=0&&y<h&&point.z>depth[y*w+x]+.0006;
    };
  }
  function lineSVG(line,isHidden){
    const parts=[],points=line.loop?[...line.points,line.points[0]]:line.points;
    let run=null;
    for(let i=0;i+1<points.length;i+=line.segments?2:1){
      const a=points[i],b=points[i+1],steps=Math.max(1,Math.ceil(Math.hypot(b.x-a.x,b.y-a.y)/2));
      if(line.segments)run=null;
      for(let j=0;j<steps;j++){
        const interpolate=t=>({x:a.x+(b.x-a.x)*t,y:a.y+(b.y-a.y)*t,z:a.z+(b.z-a.z)*t});
        const hidden=isHidden(interpolate((j+.5)/steps)),start=interpolate(j/steps),end=interpolate((j+1)/steps);
        if(hidden&&!line.bodyEdge){run=null;continue;}
        if(!run||run.hidden!==hidden){run={hidden,points:[start,end]};parts.push(run);}else run.points.push(end);
      }
    }
    return parts.map(part=>{
      const dashed=part.hidden||line.dashed;
      return `<path d="${pathPoints(part.points)}" fill="none" stroke="${line.color}" stroke-width="${number(line.width)}" stroke-linecap="round" stroke-linejoin="round" opacity="${number(part.hidden?Math.min(line.opacity,.38):line.opacity)}"${dashed?' stroke-dasharray="'+line.dash.map(number).join(' ')+'"'+(line.axis?' stroke-dashoffset="'+number(line.dashOffset)+'"':''):''}/>`;
    }).join('');
  }
  function symmetryLabelSVG(order,tex,x,y,size,color,opacity){
    const math=global.MajoranaLatex;
    const axis=math.layout(order,size);
    if(!tex.includes('\\parallel')){
      return '<g opacity="'+escape(opacity)+'">'+math.toSVG(order,x,y,size,color)+'</g>';
    }
    const z=math.layout('z',size),gap=size*.55;
    const left=x-(axis.width+gap+z.width)/2, mid=left+axis.width+gap/2;
    const baseline=y+axis.baselineOffset;
    const marks=[mid-size*.06,mid+size*.06].map(px=>
      'M'+number(px)+' '+number(baseline-size*.7)+'v'+number(size*.87)).join('');
    return '<g opacity="'+escape(opacity)+'" aria-label="'+escape(tex)+'">'+
      math.toSVG(order,left+axis.width/2,y,size,color)+
      '<path d="'+marks+'" fill="none" stroke="'+escape(color)+'" stroke-width="'+number(size*.05)+'"/>'+
      math.toSVG('z',left+axis.width+gap+z.width/2,baseline-z.baselineOffset,size,color)+'</g>';
  }
  function labelSVG(view){
    const serializer=new XMLSerializer(),output=[];
    const axes=new Map((view.symmetryAxisLabels||[]).map(axis=>[axis.label,axis]));
    for(const label of view.labels||[]){
      const element=label.el,style=getComputedStyle(element);
      if(element.hidden||style.display==='none'||style.visibility==='hidden')continue;
      const left=parseFloat(element.style.left),top=parseFloat(element.style.top);
      const width=element.offsetWidth||element.getBoundingClientRect().width,height=element.offsetHeight||element.getBoundingClientRect().height;
      if(![left,top,width,height].every(Number.isFinite)||width<=0||height<=0)continue;
      const axis=axes.get(label);
      if(axis&&global.MajoranaLatex){
        const tex=element.querySelector('annotation[encoding="application/x-tex"]')?.textContent||axis.tex;
        output.push(symmetryLabelSVG(axis.order,tex,left,top,parseFloat(style.fontSize),style.color,style.opacity));
        continue;
      }
      const wrapper=document.createElementNS(XHTML,'div'),clone=element.cloneNode(true);
      wrapper.setAttribute('style','margin:0;padding:3px;white-space:nowrap;overflow:visible;');
      for(const property of ['font-family','font-size','font-weight','font-style','line-height','color','opacity','text-shadow'])clone.style.setProperty(property,style.getPropertyValue(property));
      clone.style.position='static';clone.style.left='auto';clone.style.top='auto';clone.style.transform='none';clone.style.whiteSpace='nowrap';clone.style.display='inline-block';clone.removeAttribute('id');
      wrapper.append(clone);
      output.push(`<foreignObject x="${number(left-width/2-3)}" y="${number(top-height/2-3)}" width="${number(width+6)}" height="${number(height+6)}">${serializer.serializeToString(wrapper)}</foreignObject>`);
    }
    return output.join('');
  }
  function serialize(view,{title='Sector-length geometry'}={}){
    if(!T||!view?.scene||!view.camera)throw new Error('A rendered atlas view is required.');
    if(typeof global.ATLAS_EXPORT_FONT_CSS!=='string')throw new Error('The embedded SVG export fonts have not loaded.');
    view.render?.();
    const width=Math.max(1,view.width),height=Math.max(1,view.height);
    const {triangles,lines,globe}=collect(view),hidden=depthLookup(triangles,width,height);
    const background=getComputedStyle(document.documentElement).getPropertyValue('--paper').trim()||'#faf9f5';
    const fonts=global.ATLAS_EXPORT_FONT_CSS.replace(/\]\]>/g,']]]]><![CDATA[>');
    const body=hullPolygons(triangles).map(face=>`<polygon points="${face.points.map(p=>number(p.x)+','+number(p.y)).join(' ')}" fill="${face.fill}" fill-opacity="${number(face.opacity)}" stroke="none"/>`).join('');
    const sphere=globe?`<circle cx="${number(globe.x)}" cy="${number(globe.y)}" r="${number(globe.r)}" fill="${globe.fill}" opacity="${number(globe.opacity??.18)}"/>`:'';
    const geometryLines=lines.sort((a,b)=>a.priority-b.priority).map(line=>lineSVG(line,hidden)).join('');
    return `<?xml version="1.0" encoding="UTF-8"?>\n<svg xmlns="${NS}" xmlns:xlink="http://www.w3.org/1999/xlink" width="${number(width)}" height="${number(height)}" viewBox="0 0 ${number(width)} ${number(height)}" role="img" aria-labelledby="atlas-export-title"><title id="atlas-export-title">${escape(title)}</title><desc>Vector projection of the current view, including geometric edges, point markers, labels and embedded LaTeX fonts.</desc><defs><clipPath id="atlas-export-clip"><rect width="${number(width)}" height="${number(height)}"/></clipPath><style><![CDATA[${fonts}]]></style></defs><rect width="100%" height="100%" fill="${escape(background)}"/><g clip-path="url(#atlas-export-clip)">${sphere}<g data-layer="faces">${body}</g><g data-layer="edges">${geometryLines}</g><g data-layer="leaders">${view.svg?.innerHTML||''}</g><g data-layer="markers">${view.markerSvg?.innerHTML||''}</g><g data-layer="labels">${labelSVG(view)}</g></g></svg>`;
  }
  function download(view,{filename='sector-length-view.svg',title}={}){
    const svg=serialize(view,{title}),blob=new Blob([svg],{type:'image/svg+xml;charset=utf-8'}),url=URL.createObjectURL(blob),anchor=document.createElement('a');
    anchor.href=url;anchor.download=filename.endsWith('.svg')?filename:filename+'.svg';anchor.hidden=true;document.body.append(anchor);anchor.click();anchor.remove();setTimeout(()=>URL.revokeObjectURL(url),60000);return svg;
  }
  global.AtlasSVG={serialize,download};
})(window);
