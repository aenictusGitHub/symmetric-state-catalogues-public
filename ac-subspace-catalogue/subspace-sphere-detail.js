/* Enlarged operator-Majorana view. Each dialog has its own camera/renderer.
 * The tensor coefficients and true multipole weight are never rescaled.
 * API: SubspaceSphereDetail.open({model,mesh,symmetry,record,rank,weightTex,
 *     exactWeight,weightNumeric,subspaceId,onClose}); close().
 */
(function(global){
  'use strict';
  const esc=x=>String(x??'').replace(/[&<>"']/g,c=>({'&':'&amp;','<':'&lt;','>':'&gt;','"':'&quot;',"'":'&#39;'}[c]));
  const math=tex=>{try{return global.katex.renderToString(String(tex),{throwOnError:false,trust:false,strict:'ignore'});}catch{return esc(tex);}};
  const number=x=>Number.isFinite(x)?Number(x.toPrecision(8)).toString():'—';
  const numericTex=x=>{const value=number(x),parts=value.split('e');return parts.length===2?parts[0]+'\\times10^{'+Number(parts[1])+'}':value;};
  const axisColor=order=>({2:'#34786d',3:'#82609b',4:'#a45635',5:'#b38b23',continuous:'#356f94'}[order]||'#526b9b');
  let dialog=null,renderer=null,returnFocus=null,settings=null;
  const ui=id=>dialog.querySelector('#ssd-'+id);
  function ensure(){
    if(dialog)return;
    dialog=document.createElement('dialog');dialog.id='subspace-sphere-detail';dialog.className='ssd-dialog';dialog.setAttribute('aria-labelledby','ssd-title');dialog.setAttribute('aria-describedby','ssd-description');
    dialog.innerHTML=`<div class="ssd-heading"><div><p class="ssd-eyebrow">Individual multipole</p><h2 id="ssd-title">Operator Majorana</h2></div><button type="button" id="ssd-close" aria-label="Close enlarged sphere">Close ×</button></div>
      <p id="ssd-description" class="ssd-description">Drag to rotate this sphere independently · Scroll to zoom · Arrow keys rotate</p>
      <div id="ssd-weight" class="ssd-equation"></div><p id="ssd-approx" class="ssd-caption"></p>
      <div class="ssd-workspace"><div class="sv-frame ssd-frame"><canvas id="ssd-canvas" class="sv-canvas" role="img" tabindex="0" aria-label="Enlarged operator Majorana sphere"></canvas><div id="ssd-math-labels" class="sv-math-labels" aria-hidden="true"></div><span class="sv-frame-help">Drag to rotate · Scroll to zoom · Arrow keys rotate</span></div>
      <aside class="ssd-controls"><div class="ssd-turn" role="group" aria-label="Rotate enlarged sphere"><button id="ssd-left" type="button" aria-label="Rotate enlarged sphere left">←</button><button id="ssd-up" type="button" aria-label="Rotate enlarged sphere up">↑</button><button id="ssd-down" type="button" aria-label="Rotate enlarged sphere down">↓</button><button id="ssd-right" type="button" aria-label="Rotate enlarged sphere right">→</button><button id="ssd-reset" type="button">Reset view</button></div>
      <label><input type="checkbox" id="ssd-fit" checked>Fit sphere to window</label><p class="ssd-note">Fitting enlarges the display; the multipole weight above is unchanged.</p>
      <label><input type="checkbox" id="ssd-axes" checked>${math('C_n')} symmetry axes</label><label>Show <select id="ssd-order" aria-label="Enlarged sphere symmetry axes"></select></label>
      <button id="ssd-align" type="button">Align highest ${math('C_n')} with ${math('z')}</button>
      <label><input type="checkbox" id="ssd-coordinate">Coordinate axes</label>
      <div class="ssd-stars-key"><span><i></i>Selected subconstellation</span><span><i class="opposite"></i>Antipodal partners</span><p>Dark stars: front · Light stars: rear.</p></div>
      <div id="ssd-symmetry" class="ssd-symmetry"></div><button type="button" id="ssd-svg">Export SVG ↓</button><p id="ssd-error" class="ssd-error" role="alert" hidden></p></aside></div>
      <div hidden><select id="ssd-shape"><option value="sphere">Sphere</option></select><select id="ssd-palette"><option value="viridis">Violet & teal</option></select><input id="ssd-contours" type="checkbox"><input id="ssd-separate" type="checkbox"></div>`;
    document.body.append(dialog);
    ui('close').addEventListener('click',close);
    dialog.addEventListener('cancel',event=>event.stopPropagation());
    dialog.addEventListener('close',()=>{if(!dialog.open&&(renderer||settings))finish();});
    dialog.addEventListener('click',event=>{
      if(event.target!==dialog)return;
      const b=dialog.getBoundingClientRect();if(event.clientX<b.left||event.clientX>b.right||event.clientY<b.top||event.clientY>b.bottom)close();
    });
    ui('svg').addEventListener('click',exportSvg);
  }
  function finish(){
    renderer?.dispose();renderer=null;
    const focus=returnFocus,callback=settings?.onClose;returnFocus=null;settings=null;
    document.body.classList.remove('ssd-open');
    const parent=focus?.closest?.('dialog');
    if(focus?.isConnected&&(!parent||parent.open))focus.focus({preventScroll:true});
    callback?.();
  }
  function close(){if(dialog?.open)dialog.close();if(renderer||settings)finish();}
  function open(options){
    if(!options?.model||!Number.isInteger(Number(options.rank))||Number(options.rank)<1)throw new Error('A nonzero multipole rank is required.');
    if(typeof global.SubspaceSphereRenderer!=='function')throw new Error('The sphere renderer has not loaded.');
    ensure();if(dialog.open)close();
    settings={...options,rank:Number(options.rank)};returnFocus=document.activeElement;
    const {model,rank,symmetry={}}=settings,axes=symmetry.axes||[];
    const weight=Number.isFinite(settings.weightNumeric)?settings.weightNumeric:Number(settings.record?.weight||0);
    settings.weightNumeric=weight;
    const value=settings.exactWeight&&typeof settings.weightTex==='string'?settings.weightTex:numericTex(weight);
    settings.displayWeightTex=value;
    ui('title').innerHTML='Operator Majorana · '+math('\\ell='+rank);
    ui('weight').innerHTML=math('w_{'+rank+'}='+value);
    ui('approx').innerHTML=settings.exactWeight?math('w_{'+rank+'}\\simeq '+numericTex(weight)):'Numerically evaluated multipole weight.';
    ui('fit').checked=true;ui('axes').checked=true;ui('coordinate').checked=false;ui('separate').checked=false;ui('error').hidden=true;ui('math-labels').replaceChildren();
    const orders=[...new Set(axes.map(a=>a.order==null?'continuous':String(a.order)))].sort((a,b)=>a==='continuous'?-1:b==='continuous'?1:Number(b)-Number(a));
    ui('order').innerHTML='<option value="highest">Highest order</option><option value="all">All axes</option>'+orders.map(o=>`<option value="${esc(o)}">${o==='continuous'?'Continuous':o+'-fold'}</option>`).join('');ui('order').value='highest';
    for(const id of ['axes','order','align'])ui(id).disabled=!axes.length;
    ui('symmetry').innerHTML=orders.map(o=>{const count=axes.filter(a=>(a.order==null?'continuous':String(a.order))===o).length;return `<span class="ssd-axis-key"><i style="background:${axisColor(o)}"></i>${math('C_{'+(o==='continuous'?'\\infty':o)+'}')} : ${count} ${count===1?'axis':'axes'}</span>`;}).join('')||(symmetry.isotropic?'<p>Every spatial axis is a symmetry axis.</p>':'<p>No nontrivial rotational symmetry detected.</p>');
    if(symmetry.perpendicularHalfTurns)ui('symmetry').innerHTML+='<p class="ssd-note">Every perpendicular axis is also a '+math('C_2')+' axis.</p>';
    ui('symmetry').innerHTML+='<p class="ssd-note">Axes describe the whole projector.</p>';
    dialog.showModal();document.body.classList.add('ssd-open');
    try{
      renderer=global.SubspaceSphereRenderer(model,settings.mesh,symmetry,ui);
      renderer.set('operator',rank,settings.record);if(axes.length)renderer.align();renderer.redraw();
    }catch(error){close();throw error;}
  }
  function exportSvg(){
    if(!renderer||!settings)return;
    try{
      const {rank,model,symmetry={},displayWeightTex,weightNumeric,exactWeight}=settings;
      const counts=new Map();for(const axis of symmetry.axes||[]){const order=axis.order==null?'continuous':String(axis.order);counts.set(order,(counts.get(order)||0)+1);}
      const symmetryRows=[...counts].sort((a,b)=>a[0]==='continuous'?-1:b[0]==='continuous'?1:Number(b[0])-Number(a[0])).map(([o,count])=>({tex:'C_{'+(o==='continuous'?'\\infty':o)+'}:\\;'+count+'\\,\\text{'+(count===1?'axis':'axes')+'}',color:axisColor(o)}));
      const svg=global.SubspaceSvgExport.compose(renderer.exportSvg(),{title:'Operator Majorana · multipole '+rank,contextTex:'N='+model.n+',\\quad d='+model.dimension+',\\quad\\rho_V=P_V/d',formulaRows:[{tex:'\\ell='+rank+',\\quad w_{'+rank+'}='+displayWeightTex,note:exactWeight?'≈ '+number(weightNumeric):'Numerical value'}],symmetryRows,notes:[ui('fit').checked?'Sphere fitted for inspection; its multipole weight is unchanged.':'Sphere radius follows the multipole norm.','Whole-projector symmetry axes.']});
      const url=URL.createObjectURL(new Blob([svg],{type:'image/svg+xml'})),link=document.createElement('a');link.href=url;link.download=(settings.subspaceId||'operator')+'-multipole-'+rank+'.svg';link.click();setTimeout(()=>URL.revokeObjectURL(url),1000);
    }catch(error){ui('error').textContent='The SVG export could not be created. '+error.message;ui('error').hidden=false;}
  }
  global.SubspaceSphereDetail={open,close};
})(window);
