/* Common-scale layout for separated operator multipole spheres.
 * Every radius is maxDisplayRadius * weight / maxWeight. In particular, tiny
 * multipoles are never enlarged independently. Coordinates are canvas pixels.
 * This helper does not rotate, project, or change any multipole coefficients.
 */
(function(root){
  'use strict';
  function layout(width,height,ranks,maxWeight){
    width=Number(width);height=Number(height);
    if(!Number.isFinite(width)||!Number.isFinite(height)||width<0||height<0)throw new RangeError('Finite nonnegative viewport dimensions are required.');
    if(!Array.isArray(ranks))throw new TypeError('The multipole ranks must be an array.');
    const entries=ranks.map((entry,index)=>{
      const rank=Number(entry.rank),weight=Number(entry.weight);
      if(!Number.isInteger(rank)||rank<1||!Number.isFinite(weight)||weight<0)throw new RangeError('Each multipole needs a positive integer rank and a finite nonnegative weight.');
      return {rank,weight,index};
    }).sort((a,b)=>a.rank-b.rank||a.index-b.index);
    const largest=Math.max(0,...entries.map(entry=>entry.weight));
    const supplied=Number(maxWeight);
    // A stale caller maximum cannot make a sphere overflow its cell. Normally
    // the caller supplies the same global maximum used by the nested view.
    const normalization=Math.max(largest,Number.isFinite(supplied)&&supplied>0?supplied:0);
    const padding=Math.min(16,Math.min(width,height)*.035);
    const gap=Math.min(12,Math.min(width,height)*.025);
    const empty={cells:[],columns:0,rows:0,maxDisplayRadius:0,maxWeight:normalization,padding,gap,labelHeight:0,topLabelHeight:0};
    if(!entries.length||!width||!height)return empty;
    const innerWidth=Math.max(0,width-2*padding),innerHeight=Math.max(0,height-2*padding);
    let best=null;
    // Test every rectangular grid and use the one that gives the largest
    // common radius. Prefer fewer empty cells if two layouts are equivalent.
    for(let columns=1;columns<=entries.length;columns++){
      const rows=Math.ceil(entries.length/columns);
      const cellWidth=(innerWidth-(columns-1)*gap)/columns;
      const cellHeight=(innerHeight-(rows-1)*gap)/rows;
      if(cellWidth<=0||cellHeight<=0)continue;
      const topLabelHeight=Math.min(22,cellHeight*.22);
      const labelHeight=Math.min(28,Math.max(16,cellHeight*.23),cellHeight*.28);
      const contentHeight=cellHeight-topLabelHeight-labelHeight;
      const inset=Math.min(8,cellWidth*.06,contentHeight*.06);
      const radius=Math.max(0,Math.min(cellWidth-2*inset,contentHeight-2*inset)/2);
      const unused=columns*rows-entries.length;
      if(!best||radius>best.radius+1e-9||Math.abs(radius-best.radius)<=1e-9&&unused<best.unused){
        best={columns,rows,cellWidth,cellHeight,labelHeight,topLabelHeight,inset,radius,unused};
      }
    }
    if(!best)return empty;
    const cells=entries.map((entry,index)=>{
      const row=Math.floor(index/best.columns),column=index%best.columns;
      const rowCount=Math.min(best.columns,entries.length-row*best.columns);
      const rowWidth=rowCount*best.cellWidth+(rowCount-1)*gap;
      const x=(width-rowWidth)/2+column*(best.cellWidth+gap);
      const y=padding+row*(best.cellHeight+gap);
      const labelY=y+best.cellHeight-best.labelHeight;
      const relativeRadius=normalization>0?entry.weight/normalization:0;
      return {
        rank:entry.rank,weight:entry.weight,row,column,x,y,
        width:best.cellWidth,height:best.cellHeight,
        cx:x+best.cellWidth/2,cy:y+best.topLabelHeight+(best.cellHeight-best.topLabelHeight-best.labelHeight)/2,
        radius:best.radius*relativeRadius,relativeRadius,
        maxDisplayRadius:best.radius,
        topLabelArea:{x:x+best.inset,y,width:best.cellWidth-2*best.inset,height:best.topLabelHeight,baseline:y+Math.min(18,best.topLabelHeight*.8)},
        labelArea:{x:x+best.inset,y:labelY,width:best.cellWidth-2*best.inset,height:best.labelHeight,baseline:labelY+Math.min(18,best.labelHeight*.8)}
      };
    });
    return {cells,columns:best.columns,rows:best.rows,maxDisplayRadius:best.radius,maxWeight:normalization,padding,gap,labelHeight:best.labelHeight,topLabelHeight:best.topLabelHeight};
  }
  const api={layout};
  if(typeof module==='object'&&module.exports)module.exports=api;
  root.SubspaceSphereLayout=api;
})(typeof window==='object'?window:globalThis);
