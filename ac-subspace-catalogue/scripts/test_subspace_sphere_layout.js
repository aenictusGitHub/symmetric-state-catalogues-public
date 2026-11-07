/* Reproducible geometry checks; no browser or graphics library is required. */
'use strict';
const assert=require('node:assert/strict');
const {layout}=require('../subspace-sphere-layout.js');
const near=(a,b,eps=1e-9)=>assert.ok(Math.abs(a-b)<=eps*Math.max(1,Math.abs(a),Math.abs(b)),`${a} != ${b}`);
const overlap=(a,b)=>a.x<b.x+b.width-1e-9&&b.x<a.x+a.width-1e-9&&a.y<b.y+b.height-1e-9&&b.y<a.y+a.height-1e-9;
let cases=0;
for(const [width,height]of [[320,310],[486,437],[900,600],[390,844],[1000,310],[70,60]]){
  for(let count=1;count<=42;count++){
    const ranks=Array.from({length:count},(_,i)=>({rank:count-i,weight:(1+(i*17)%43)/43}));
    const result=layout(width,height,ranks,1),cells=result.cells;
    assert.equal(cells.length,count);
    assert.deepEqual(cells.map(c=>c.rank),Array.from({length:count},(_,i)=>i+1));
    for(const c of cells){
      assert.ok(c.x>=-1e-9&&c.y>=-1e-9&&c.x+c.width<=width+1e-9&&c.y+c.height<=height+1e-9,'Cell outside viewport');
      assert.ok(c.cx-c.radius>=c.x-1e-9&&c.cx+c.radius<=c.x+c.width+1e-9,'Sphere outside cell horizontally');
      assert.ok(c.cy-c.radius>=c.topLabelArea.y+c.topLabelArea.height-1e-9&&c.cy+c.radius<=c.labelArea.y+1e-9,'Sphere overlaps label or cell border');
      assert.ok(c.cx-result.maxDisplayRadius>=c.x-1e-9&&c.cx+result.maxDisplayRadius<=c.x+c.width+1e-9,'Maximum-norm reference outside cell');
      assert.ok(c.cy-result.maxDisplayRadius>=c.topLabelArea.y+c.topLabelArea.height-1e-9&&c.cy+result.maxDisplayRadius<=c.labelArea.y+1e-9,'Maximum-norm reference overlaps labels');
      assert.ok(c.labelArea.y+c.labelArea.height<=c.y+c.height+1e-9,'Label outside cell');
      near(c.radius,result.maxDisplayRadius*c.weight);
    }
    for(let i=0;i<cells.length;i++)for(let j=i+1;j<cells.length;j++)assert.equal(overlap(cells[i],cells[j]),false,'Cells overlap');
    const lastRow=cells.filter(c=>c.row===result.rows-1);
    near(lastRow[0].x,width-(lastRow.at(-1).x+lastRow.at(-1).width));
    cases++;
  }
}
const tiny=layout(600,400,[{rank:2,weight:1},{rank:4,weight:1e-12}],1);
near(tiny.cells[1].radius/tiny.cells[0].radius,1e-12,1e-20);
assert.ok(tiny.cells[1].radius<1e-9,'Tiny weights must not be independently enlarged');
const common=layout(600,400,[{rank:2,weight:.1},{rank:4,weight:.2}],1);
near(common.cells[1].radius/common.cells[0].radius,2);
near(common.cells[1].relativeRadius,.2);
assert.deepEqual(layout(300,300,[],1).cells,[]);
assert.deepEqual(layout(0,300,[{rank:2,weight:1}],1).cells,[]);
const zeros=layout(300,300,[{rank:2,weight:0}],0);assert.equal(zeros.cells[0].radius,0);
assert.throws(()=>layout(NaN,100,[],1),RangeError);
assert.throws(()=>layout(100,100,[{rank:2,weight:-1}],1),RangeError);
assert.throws(()=>layout(100,100,[{rank:1.2,weight:1}],1),RangeError);
console.log(`PASS: ${cases} responsive layouts; rank order, common-scale radii, no overlaps, labels, last-row centering and tiny-weight preservation.`);
