/* Shared Husimi palette and physical-sphere isocontours.
 * Pure geometry: contours approximate the sampled field, never alter Q values.
 * CommonJS export is supplied for reproducible tests; no browser dependencies. */
(function(root){
  'use strict';
  const palettes={
    plum:[{t:0,rgb:[34,13,71]},{t:.2,rgb:[87,24,107]},{t:.4,rgb:[143,39,109]},{t:.6,rgb:[203,70,84]},{t:.8,rgb:[241,124,63]},{t:1,rgb:[253,218,142]}],
    ocean:[{t:0,rgb:[18,35,85]},{t:.22,rgb:[36,99,176]},{t:.44,rgb:[56,183,196]},{t:.66,rgb:[243,232,162]},{t:.84,rgb:[241,148,82]},{t:1,rgb:[190,54,61]}],
    viridis:[{t:0,rgb:[68,1,84]},{t:.25,rgb:[59,82,139]},{t:.5,rgb:[33,145,140]},{t:.75,rgb:[94,201,98]},{t:1,rgb:[253,231,37]}]
  };
  const paletteStops=palettes.plum;
  const stopsFor=name=>palettes[name]||paletteStops;
  const clamp=x=>Math.max(0,Math.min(1,Number.isFinite(x)?x:0));
  const hex=a=>'#'+a.map(x=>Math.round(x).toString(16).padStart(2,'0')).join('');
  function rgb(t,palette='plum'){
    const stops=stopsFor(palette),u=clamp(t);let k=0;
    while(k<stops.length-2&&u>stops[k+1].t)k++;
    const a=stops[k],b=stops[k+1],f=(u-a.t)/(b.t-a.t);
    return a.rgb.map((v,i)=>v+(b.rgb[i]-v)*f);
  }
  const color=(t,palette='plum')=>hex(rgb(t,palette));
  const gradient=(palette='plum')=>`linear-gradient(90deg,${stopsFor(palette).map(s=>hex(s.rgb)+' '+100*s.t+'%').join(',')})`;
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
    if(!vertices.length||!faces.length||!levels.length)return{segments:[],levels:[]};
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
    for(const ids of faces){
      const corners=ids.map(i=>vertices[i]);
      if(corners.length<3||corners.some(v=>!v?.n||!Number.isFinite(v.q)))continue;
      const low=Math.min(...corners.map(v=>v.q)),high=Math.max(...corners.map(v=>v.q));
      if(low===high)continue;
      if(corners.length===3){for(const level of levels)if(level>=low&&level<=high)triangle(...corners,level,(level-mesh.qMin)/range);continue;}
      const centerN=normalize([0,1,2].map(k=>corners.reduce((s,v)=>s+v.n[k],0)/corners.length));
      if(!centerN)continue;
      const center={n:centerN,q:corners.reduce((s,v)=>s+v.q,0)/corners.length};
      for(const level of levels){
        if(level<low||level>high)continue;
        const t=(level-mesh.qMin)/range;
        for(let k=0;k<corners.length;k++)triangle(corners[k],corners[(k+1)%corners.length],center,level,t);
      }
    }
    return{segments,levels};
  }
  const api={color,rgb,gradient,paletteStops,palettes,defaultLevels,contours,version:'1.1.0'};
  root.HusimiRendering=api;
  if(typeof module==='object'&&module.exports)module.exports=api;
})(typeof window==='object'?window:globalThis);
