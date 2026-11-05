/* Convex hull of distinct displayed points. Numerical geometry only.
 * Incremental 3D hull; coplanar triangles are merged into polygonal faces.
 * Repeated positions do not change the hull. No external dependencies. */
(() => {
  'use strict';
  const EPS=1e-12;
  const sub=(a,b)=>a.map((x,i)=>x-b[i]);
  const dot=(a,b)=>a.reduce((s,x,i)=>s+x*b[i],0);
  const cross=(a,b)=>[a[1]*b[2]-a[2]*b[1],a[2]*b[0]-a[0]*b[2],a[0]*b[1]-a[1]*b[0]];
  const norm=a=>Math.hypot(...a);
  function polygon(ids,vertices,normal){
    const helper=Math.abs(normal[0])<.8?[1,0,0]:[0,1,0];
    let u=cross(normal,helper);u=u.map(x=>x/norm(u));
    const v=cross(normal,u),points=ids.map(id=>({id,x:dot(vertices[id],u),y:dot(vertices[id],v)})).sort((a,b)=>a.x-b.x||a.y-b.y);
    const turn=(a,b,c)=>(b.x-a.x)*(c.y-a.y)-(b.y-a.y)*(c.x-a.x);
    const half=list=>{const out=[];for(const p of list){while(out.length>1&&turn(out.at(-2),out.at(-1),p)<=EPS)out.pop();out.push(p);}return out;};
    const lower=half(points),upper=half([...points].reverse());lower.pop();upper.pop();return [...lower,...upper].map(p=>p.id);
  }
  function build(points){
    const vertices=[];
    for(const p of points){const v=[p.x,p.y,p.z];if(v.every(Number.isFinite)&&!vertices.some(q=>norm(sub(v,q))<EPS))vertices.push(v);}
    const finish=(dimension,faces)=>{
      const edges=new Map();for(const face of faces)for(let i=0;i<face.length;i++){const a=face[i],b=face[(i+1)%face.length],key=[Math.min(a,b),Math.max(a,b)].join(',');edges.set(key,[a,b]);}
      return {dimension,vertices:vertices.map(([x,y,z])=>({x,y,z})),faces,edges:[...edges.values()]};
    };
    if(vertices.length<2)return finish(0,[]);
    const farthest=score=>vertices.reduce((best,v,i)=>score(v)>score(vertices[best])?i:best,0);
    const a=0,b=farthest(v=>norm(sub(v,vertices[a]))),line=sub(vertices[b],vertices[a]);
    const c=farthest(v=>norm(cross(line,sub(v,vertices[a]))));
    if(norm(cross(line,sub(vertices[c],vertices[a])))<EPS)return finish(1,[]);
    let normal=cross(line,sub(vertices[c],vertices[a]));normal=normal.map(x=>x/norm(normal));
    const d=farthest(v=>Math.abs(dot(normal,sub(v,vertices[a]))));
    if(Math.abs(dot(normal,sub(vertices[d],vertices[a])))<EPS)return finish(2,[polygon(vertices.map((_,i)=>i),vertices,normal)]);
    const inside=[0,1,2].map(i=>(vertices[a][i]+vertices[b][i]+vertices[c][i]+vertices[d][i])/4);
    const face=(i,j,k)=>{
      let ids=[i,j,k],n=cross(sub(vertices[j],vertices[i]),sub(vertices[k],vertices[i])),length=norm(n);
      if(length<EPS)return null;
      if(dot(n,sub(inside,vertices[i]))>0){ids=[i,k,j];n=n.map(x=>-x);}
      n=n.map(x=>x/length);return {ids,normal:n,offset:dot(n,vertices[i])};
    };
    let faces=[[a,b,c],[a,d,b],[b,d,c],[c,d,a]].map(ids=>face(...ids));
    for(let i=0;i<vertices.length;i++){
      if([a,b,c,d].includes(i))continue;
      const visible=faces.filter(f=>dot(f.normal,vertices[i])-f.offset>EPS);
      if(!visible.length)continue;
      const boundary=new Map();
      for(const f of visible)for(let e=0;e<3;e++){const u=f.ids[e],v=f.ids[(e+1)%3],key=[Math.min(u,v),Math.max(u,v)].join(',');if(boundary.has(key))boundary.delete(key);else boundary.set(key,[u,v]);}
      const removed=new Set(visible);faces=faces.filter(f=>!removed.has(f));
      for(const [u,v] of boundary.values()){const f=face(u,v,i);if(f)faces.push(f);}
    }
    const groups=[];
    for(const f of faces){let g=groups.find(g=>norm(sub(g.normal,f.normal))<1e-10&&Math.abs(g.offset-f.offset)<EPS);if(!g){g={normal:f.normal,offset:f.offset,ids:new Set()};groups.push(g);}f.ids.forEach(i=>g.ids.add(i));}
    return finish(3,groups.map(g=>polygon([...g.ids],vertices,g.normal)));
  }
  window.MajoranaHull=Object.freeze({build});
})();
