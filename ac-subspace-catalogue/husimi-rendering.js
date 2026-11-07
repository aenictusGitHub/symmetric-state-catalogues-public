/* Shared Husimi palette and physical-sphere isocontours.
 * Pure geometry: contours approximate the sampled field, never alter Q values.
 * CommonJS export is supplied for reproducible tests; no browser dependencies. */
(function(root){
  'use strict';
  const palettes={
    wigner:[{t:0,rgb:[43,88,149]},{t:.25,rgb:[126,174,205]},{t:.5,rgb:[251,250,244]},{t:.75,rgb:[226,156,112]},{t:1,rgb:[177,65,28]}],
    plum:[{t:0,rgb:[34,13,71]},{t:.2,rgb:[87,24,107]},{t:.4,rgb:[143,39,109]},{t:.6,rgb:[203,70,84]},{t:.8,rgb:[241,124,63]},{t:1,rgb:[253,218,142]}],
    ocean:[{t:0,rgb:[18,35,85]},{t:.22,rgb:[36,99,176]},{t:.44,rgb:[56,183,196]},{t:.66,rgb:[243,232,162]},{t:.84,rgb:[241,148,82]},{t:1,rgb:[190,54,61]}],
    viridis:[{t:0,rgb:[68,1,84]},{t:.25,rgb:[59,82,139]},{t:.5,rgb:[33,145,140]},{t:.75,rgb:[94,201,98]},{t:1,rgb:[253,231,37]}]
  };
  const paletteStops=palettes.viridis;
  const stopsFor=name=>palettes[name]||paletteStops;
  const clamp=x=>Math.max(0,Math.min(1,Number.isFinite(x)?x:0));
  const hex=a=>'#'+a.map(x=>Math.round(x).toString(16).padStart(2,'0')).join('');
  function rgb(t,palette='viridis'){
    const stops=stopsFor(palette),u=clamp(t);let k=0;
    while(k<stops.length-2&&u>stops[k+1].t)k++;
    const a=stops[k],b=stops[k+1],f=(u-a.t)/(b.t-a.t);
    return a.rgb.map((v,i)=>v+(b.rgb[i]-v)*f);
  }
  const color=(t,palette='viridis')=>hex(rgb(t,palette));
  const gradient=(palette='viridis')=>`linear-gradient(90deg,${stopsFor(palette).map(s=>hex(s.rgb)+' '+100*s.t+'%').join(',')})`;
  const norm=v=>Math.hypot(...v);
  const normalize=v=>{const n=norm(v);return n>1e-14?v.map(x=>x/n):null;};
  function defaultLevels(mesh,count=8){
    const lo=Number(mesh.qMin),hi=Number(mesh.qMax),range=hi-lo;
    if(!Number.isFinite(range)||range<=1e-12*Math.max(1,Math.abs(lo),Math.abs(hi)))return[];
    const n=Math.max(0,Math.min(32,Math.floor(count)));
    return Array.from({length:n},(_,i)=>lo+range*(i+1)/(n+1));
  }
  function contours(mesh,requestedLevels){
    const vertices=mesh?.vertices||[],faces=mesh?.faces||[];
    const levels=(Array.isArray(requestedLevels)?requestedLevels:defaultLevels(mesh)).filter(Number.isFinite).sort((a,b)=>a-b);
    const regular=mesh?.values&&Number.isInteger(mesh.rows)&&Number.isInteger(mesh.cols)&&mesh.rows>0&&mesh.cols>0&&mesh.values.length===(mesh.rows+1)*(mesh.cols+1);
    if((!regular&&(!vertices.length||!faces.length))||!levels.length)return{segments:[],levels:[]};
    const range=mesh.qMax-mesh.qMin;
    if(!(range>1e-12*Math.max(1,Math.abs(mesh.qMin),Math.abs(mesh.qMax))))return{segments:[],levels:[]};
    const segments=[],seen=new Set(),epsilon=1e-10;
    const key=p=>p.map(v=>Math.round(v/epsilon)).join(',');
    function triangle(a,b,c,level,t){
      const v=[a,b,c],p=[];
      for(let i=0;i<3;i++){
        const x=v[i],y=v[(i+1)%3];
        if((x.q>level)===(y.q>level))continue;
        const alpha=clamp((level-x.q)/(y.q-x.q));
        const point=normalize(x.n.map((z,k)=>z+(y.n[k]-z)*alpha));
        if(point&&!p.some(q=>norm(point.map((z,k)=>z-q[k]))<epsilon))p.push(point);
      }
      if(p.length!==2||norm(p[0].map((z,k)=>z-p[1][k]))<epsilon)return;
      const keys=p.map(key).sort(),id=level+':'+keys.join('|');
      if(seen.has(id))return;seen.add(id);
      segments.push({a:p[0],b:p[1],level,t});
    }
    // Four triangles around the bilinear center avoid choosing one persistent
    // diagonal through every cell. They also resolve four-edge saddle cases
    // into a deterministic contour topology without connecting across a face.
    function cell(corners){
      if(corners.length<3||corners.some(v=>!v?.n||!Number.isFinite(v.q)))return;
      const low=Math.min(...corners.map(v=>v.q)),high=Math.max(...corners.map(v=>v.q));
      if(low===high)return;
      if(corners.length===3){for(const level of levels)if(level>=low&&level<=high)triangle(...corners,level,(level-mesh.qMin)/range);return;}
      const centerN=normalize([0,1,2].map(k=>corners.reduce((s,v)=>s+v.n[k],0)/corners.length));
      if(!centerN)return;
      const center={n:centerN,q:corners.reduce((s,v)=>s+v.q,0)/corners.length};
      for(const level of levels){
        if(level<low||level>high)continue;
        const t=(level-mesh.qMin)/range;
        for(let k=0;k<corners.length;k++)triangle(corners[k],corners[(k+1)%corners.length],center,level,t);
      }
    }
    if(regular){
      // Keep only two rows of directions: a fine contour field does not need
      // hundreds of thousands of retained surface vertices and face objects.
      const {rows,cols,values}=mesh,stride=cols+1;
      const azimuth=Array.from({length:cols},(_,j)=>[Math.cos(2*Math.PI*j/cols),Math.sin(2*Math.PI*j/cols)]);
      const row=i=>{const theta=Math.PI*i/rows,z=Math.cos(theta),r=i===0||i===rows?0:Math.sin(theta);return Array.from({length:cols+1},(_,j)=>{const [x,y]=azimuth[j===cols?0:j];return {n:[r*x,r*y,z],q:values[i*stride+j]};});};
      let upper=row(0);
      for(let i=0;i<rows;i++){const lower=row(i+1);for(let j=0;j<cols;j++)cell([upper[j],upper[j+1],lower[j+1],lower[j]]);upper=lower;}
    }else for(const ids of faces)cell(ids.map(i=>vertices[i]));
    return{segments,levels};
  }
  // Join only degree-two vertices. Saddles and other junctions remain separate;
  // simplification selects existing points and has a bounded chord deviation.
  function contourPaths(result,options={}){
    const epsilon=options.epsilon??1e-9,tolerance=options.tolerance??0;
    const byLevel=new Map(),paths=[];
    const key=p=>p.map(v=>Math.round(v/epsilon)).join(',');
    for(const s of result?.segments||[]){
      let graph=byLevel.get(s.level);if(!graph){graph={nodes:new Map(),edges:[]};byLevel.set(s.level,graph);}
      const ids=[s.a,s.b].map(p=>{const id=key(p);if(!graph.nodes.has(id))graph.nodes.set(id,{p,edges:[]});return id;});
      if(ids[0]===ids[1])continue;
      const edge={ids,used:false,t:s.t};const index=graph.edges.push(edge)-1;
      for(const id of ids)graph.nodes.get(id).edges.push(index);
    }
    const simplify=points=>{
      if(!(tolerance>0)||points.length<3)return points;
      const keep=new Uint8Array(points.length);keep[0]=keep[points.length-1]=1;const stack=[[0,points.length-1]],threshold=tolerance*tolerance;
      while(stack.length){const [lo,hi]=stack.pop(),a=points[lo],b=points[hi],ab=b.map((v,k)=>v-a[k]),length=ab.reduce((sum,v)=>sum+v*v,0);let maximum=threshold,index=-1;
        for(let i=lo+1;i<hi;i++){const ap=points[i].map((v,k)=>v-a[k]),u=length?clamp(ap.reduce((sum,v,k)=>sum+v*ab[k],0)/length):0;const distance=ap.reduce((sum,v,k)=>sum+(v-u*ab[k])**2,0);if(distance>maximum){maximum=distance;index=i;}}
        if(index>=0){keep[index]=1;stack.push([lo,index],[index,hi]);}
      }
      return points.filter((_,i)=>keep[i]);
    };
    for(const [level,graph]of byLevel){
      function walk(start,index){const points=[graph.nodes.get(start).p],first=graph.edges[index];let current=start;
        while(index!=null){const edge=graph.edges[index];if(edge.used)break;edge.used=true;current=edge.ids[0]===current?edge.ids[1]:edge.ids[0];const node=graph.nodes.get(current);points.push(node.p);if(current===start||node.edges.length!==2)break;index=node.edges.find(i=>!graph.edges[i].used);}
        if(points.length>1)paths.push({points:simplify(points),level,t:first.t,closed:current===start});
      }
      for(const [id,node]of graph.nodes)if(node.edges.length!==2)for(const index of node.edges)if(!graph.edges[index].used)walk(id,index);
      for(let i=0;i<graph.edges.length;i++)if(!graph.edges[i].used)walk(graph.edges[i].ids[0],i);
    }
    return paths;
  }
  const api={color,rgb,gradient,paletteStops,palettes,defaultLevels,contours,contourPaths,version:'1.1.0'};
  root.HusimiRendering=api;
  if(typeof module==='object'&&module.exports)module.exports=api;
})(typeof window==='object'?window:globalThis);
