'use strict';
const fs=require('node:fs'),path=require('node:path');
const root=path.resolve(__dirname,'..');
const M=require(path.join(root,'majorana-math.js'));
const katex=require(path.join(root,'vendor/katex/katex.min.js'));
const data=JSON.parse(fs.readFileSync(path.join(root,'data/n6-boundary-samples.json'),'utf8'));
let expressions=0,minFidelity=1,maxSectorError=0;
const entries=[];
for(const state of data.samples){
  const math=[state.analyticalExpression.latex,state.analyticalExpression.parameterLatex,...state.exactSectorValues.values,state.exactSectorValues.definitionsLatex,state.mathLabel.slice(1,-1)].filter(Boolean);
  for(const source of math){
    if(/\\\\(?:frac|sqrt|begin|sum|mathrm)/.test(source))throw Error('Double-escaped command: '+state.id);
    katex.renderToString(source,{throwOnError:true,strict:false});expressions++;
  }
  const rendered=katex.renderToString(state.analyticalExpression.latex,{throwOnError:true});
  if(!rendered.includes('class="mfrac"'))throw Error('Missing real typeset fraction: '+state.id);
  const sectors=M.sectorLengths(state.coefficients);
  const error=Math.max(...sectors.map((v,k)=>Math.abs(v-state.sectors[k])));
  const stars=M.majoranaStars(state.coefficients);
  if(error>1e-9||stars.fidelity<1-1e-10||stars.stars.length!==6)throw Error('State or constellation mismatch: '+state.id);
  minFidelity=Math.min(minFidelity,stars.fidelity);maxSectorError=Math.max(maxSectorError,error);
  entries.push({id:state.id,sectorError:error,majoranaFidelity:stars.fidelity,starCount:stars.stars.length});
}
const result={status:'PASS',states:entries.length,expressionsRendered:expressions,minimumMajoranaFidelity:minFidelity,maximumBrowserSectorError:maxSectorError,checks:entries};
fs.writeFileSync(path.join(__dirname,'n6-boundary-sample-ui-math.json'),JSON.stringify(result,null,2)+'\n');
process.stdout.write(JSON.stringify({...result,checks:undefined},null,2)+'\n');
