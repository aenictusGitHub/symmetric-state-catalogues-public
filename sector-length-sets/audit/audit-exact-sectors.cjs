'use strict';
const fs=require('node:fs'),path=require('node:path');
const root=path.resolve(__dirname,'..');
const katex=require(path.join(root,'vendor/katex/katex.min.js'));
const exact=JSON.parse(fs.readFileSync(path.join(root,'data/exact-sector-values.json'),'utf8'));
let values=0,definitions=0;
for(const [id,r]of Object.entries(exact)){
  if(!Array.isArray(r.values)||!r.values.length)throw Error('Missing tuple: '+id);
  for(const value of r.values){katex.renderToString(value,{throwOnError:true,strict:false});values++;}
  if(r.definitionsLatex){katex.renderToString(r.definitionsLatex,{throwOnError:true,strict:false});definitions++;}
}
const result={status:'PASS',entries:Object.keys(exact).length,valuesRendered:values,definitionsRendered:definitions};
fs.writeFileSync(path.join(__dirname,'exact-sector-typesetting-audit.json'),JSON.stringify(result,null,2)+'\n');
process.stdout.write(JSON.stringify(result,null,2)+'\n');
