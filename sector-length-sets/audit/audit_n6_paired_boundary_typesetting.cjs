'use strict';
const fs=require('node:fs'),path=require('node:path');
const root=path.resolve(__dirname,'..'),katex=require(path.join(root,'vendor/katex/katex.min.js'));
const M=require(path.join(root,'majorana-math.js'));
const payload=JSON.parse(fs.readFileSync(path.join(root,'data/n6-paired-boundary-samples.json'),'utf8'));
let count=0,minFidelity=1,maxBrowserSectorError=0;
for(const state of payload.samples){
 const sources=[state.analyticalExpression.latex,state.analyticalExpression.parameterLatex,...state.exactSectorValues.values,state.mathLabel.slice(1,-1)];
 for(const source of sources){
  if(/\\\\(?:frac|sqrt|begin|sum|mathrm)/.test(source))throw Error('Double-escaped command: '+state.id);
  katex.renderToString(source,{throwOnError:true,strict:false});count++;
 }
 const formula=katex.renderToString(state.analyticalExpression.latex,{throwOnError:true});
 if(!formula.includes('class="mfrac"'))throw Error('Fractions not actually typeset: '+state.id);
 const stars=M.majoranaStars(state.coefficients),S=M.sectorLengths(state.coefficients);
 if(stars.stars.length!==6||stars.fidelity<1-1e-9)throw Error('Constellation failure: '+state.id);
 minFidelity=Math.min(minFidelity,stars.fidelity);
 maxBrowserSectorError=Math.max(maxBrowserSectorError,...S.map((x,k)=>Math.abs(x-state.sectors[k])));
}
const report={status:'PASS',states:payload.samples.length,expressionsRendered:count,minimumMajoranaFidelity:minFidelity,maximumBrowserSectorError:maxBrowserSectorError};
fs.writeFileSync(path.join(__dirname,'n6-paired-boundary-typesetting-audit.json'),JSON.stringify(report,null,2)+'\n');
process.stdout.write(JSON.stringify(report,null,2)+'\n');
