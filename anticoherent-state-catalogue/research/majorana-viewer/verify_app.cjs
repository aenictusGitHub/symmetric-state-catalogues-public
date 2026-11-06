const fs=require('fs'),vm=require('vm'),path=require('path');
const catalogueRoot=path.resolve(__dirname,'../..'),base=path.dirname(catalogueRoot);let assertions=0;
const assert=(v,m)=>{assertions++;if(!v)throw Error(m);};
function session(name,withCoordinates=true){
 const root=path.join(base,name),isAC=Object.hasOwn(JSON.parse(fs.readFileSync(path.join(root,'data/catalogue.json'),'utf8')),'states'),elements=new Map();let doc,opened=null,lastBlob=null,closes=0;
 const classes=()=>({add(){},remove(){},contains(){return false;}});
 class El{
  constructor(id){this.id=id;this.handlers={};this.listeners=new Map();this.value='';this.classList=classes();this.dataset={};this.attributes={};this.owned=[];this._html='';this.isConnected=true;}
  set innerHTML(v){for(const id of this.owned){elements.get(id).isConnected=false;elements.delete(id);}this.owned=[];this._html=v;for(const m of v.matchAll(/\bid="([^"]+)"/g)){const e=new El(m[1]);elements.set(e.id,e);this.owned.push(e.id);}}
  get innerHTML(){return this._html;}querySelectorAll(){return [];}querySelector(){return null;}
  addEventListener(k,f){this.handlers[k]=f;}removeEventListener(){}setAttribute(k,v){this.attributes[k]=v;}focus(){doc.activeElement=this;}click(){this.handlers.click?.({target:this});}remove(){this.isConnected=false;}scrollIntoView(){}
 }
 for(const m of fs.readFileSync(path.join(root,'index.html'),'utf8').matchAll(/\bid="([^"]+)"/g))elements.set(m[1],new El(m[1]));
 doc={activeElement:null,getElementById:id=>elements.get(id)||null,querySelectorAll:()=>[],querySelector:()=>null,addEventListener(){},createTextNode:v=>v,createElement:tag=>tag==='template'?{content:{querySelectorAll:()=>[]},innerHTML:''}:new El(tag),body:{append(){}}};
 class FakeURL extends URL {static createObjectURL(blob){lastBlob=blob;return 'blob:test';}static revokeObjectURL(){}}
 const hashListeners=[],loc={hash:'',href:'file:///offline/'+name+'/index.html'},win={katex:require(path.join(root,'vendor/katex/katex.min.js')),addEventListener(type,listener){if(type==='hashchange')hashListeners.push(listener);},MajoranaViewer:{open(r){opened=r;return true;},close(){opened=null;closes++;}}};
 const ctx={window:win,document:doc,location:loc,history:{replaceState(a,b,h){loc.hash=h;}},URLSearchParams,URL:FakeURL,Blob,console,setTimeout:()=>0,clearTimeout(){},navigator:{},matchMedia:()=>({matches:false})};
 vm.runInNewContext(fs.readFileSync(path.join(root,'data/catalogue.js'),'utf8'),ctx);
 if(withCoordinates)vm.runInNewContext(fs.readFileSync(path.join(root,'data/majorana.js'),'utf8'),ctx);
 const modalInterceptors=win.MajoranaViewer;
 for(const file of ['majorana-symmetry.js','symmetry-groups.js','point-groups.js','majorana-viewer.js'])vm.runInNewContext(fs.readFileSync(path.join(root,file),'utf8'),ctx,{filename:file});
 win.MajoranaViewer={...win.MajoranaViewer,...modalInterceptors};
 const expose=isAC?'window.__audit={select,getSelected:()=>selected,renderDetail};':'window.__audit={select(n,k,m){selectedN=n;selectedK=k;mode=m;renderTable();renderDetail();updateHash();},getSelected:selected,renderDetail};';
 const source=fs.readFileSync(path.join(root,'app.js'),'utf8').replace(/\}\)\(\);\s*$/,expose+'})();');vm.runInNewContext(source,ctx);
 return {elements,api:win.__audit,win,changeHash(hash){loc.hash=hash;for(const listener of hashListeners)listener();},get opened(){return opened;},get blob(){return lastBlob;},get closes(){return closes;},data:isAC?win.AC_CATALOGUE:win.SECTOR_CATALOGUE};
}
(async()=>{
 const summaries=[];
 for(const name of [path.basename(catalogueRoot)]){
  const env=session(name),ac=Object.hasOwn(env.data,'states'),records=env.win.MAJORANA_DATA.records;let selectedCount=0,emptyCount=0,downloads=0,switches=0;
  if(ac){
   const comparison=env.data.sphericalDesignComparison;
   assert(comparison.rows.length===100,'spherical-design comparison covers all100 point counts');
   assert(JSON.stringify(comparison.mismatchNumbers)==='[29,87,94]','only three recorded-order differences');
   const matrix=env.elements.get('state-table').innerHTML;
   const badges=[...matrix.matchAll(/<button[^>]*data-n="(\d+)" data-t="(\d+)"[^>]*>([\s\S]*?)<\/button>/g)].filter(m=>m[3].includes('design-gap-badge'));
   assert(JSON.stringify(badges.map(m=>`${m[1]}:${m[2]}`))==='["29:5","87:11","94:12"]','comparison badges occur only on differing row maxima');
   for(const badge of badges){
    const row=comparison.rows.find(row=>row.n===Number(badge[1]));
    assert(badge[3].includes(`Design order ${row.designOrder}`),'badge names the design reference order in plain language');
    assert(badge[0].includes(`for ${row.n} points`),'badge description identifies the point count');
    assert(badge[0].includes('not the design order of this state'),'badge distinguishes the selected constellation');
    assert(badge[0].includes('is conjectured'),'badge preserves maximum uncertainty');
   }
   assert(env.elements.get('design-comparison').hidden===false,'comparison panel available');
   const index=fs.readFileSync(path.join(catalogueRoot,'index.html'),'utf8');
   assert(!/<details\b[^>]*id="design-comparison"[^>]*\bopen\b/.test(index),'comparison starts collapsed');
   assert(env.elements.get('design-comparison-count').textContent==='3 differences in reported orders','plain comparison summary');
   const body=env.elements.get('design-comparison-body').innerHTML;
   assert(body.includes('true maxima are not proved to differ'),'comparison does not claim unequal true maxima');
   assert(body.includes('https://neilsloane.com/sphdesigns/index.html'),'comparison links primary benchmark');
   assert(body.includes('Numerical record; rigorous bounds 7–16')&&body.includes('Numerical record; rigorous bounds 9–17'),'numerical AC record bounds remain distinct');
   for(const row of comparison.rows.filter(row=>!row.ordersMatch)){
    assert(body.includes(`href="#n=${row.n}&amp;t=${row.acOrder}"`),'comparison links row maximum');
    env.changeHash(`#n=${row.n}&t=${row.acOrder}`);
    assert(env.api.getSelected().n===row.n && env.api.getSelected().t===row.acOrder,'comparison hash selects correct AC cell');
    const detail=env.elements.get('result-detail').innerHTML;
    assert(detail.includes(`With ${row.n} points on the sphere`)&&detail.includes(`order ${row.designOrder}</strong>`),'selected result explains the design badge');
    assert(detail.includes('numerical construction')&&detail.includes('independent benchmark'),'selected result scopes the numerical design record');
   }
   env.changeHash('#design-comparison');
   assert(env.elements.get('design-comparison').open===true,'detail comparison link opens its explanation');
  }
  for(const entry of env.data.entries){
   const candidateIds=ac?entry.stateIds||[]:[entry.id],cases=candidateIds.length?candidateIds:[null];
   for(const id of cases){
    const previous=env.closes;ac?env.api.select(entry.n,entry.t,id):env.api.select(entry.n,entry.k,entry.mode);
    assert(env.closes>previous,'rerender closes viewer');
    const button=env.elements.get('view-majorana');const allowed=!ac||['exists','numerical'].includes(entry.status)&&id;
    if(allowed){assert(button,'viewer trigger for '+id);button.click();assert(env.opened===records[id],'trigger passes exact representative '+id);assert(env.opened.n===entry.n,'qubit count consistent');selectedCount++;}
    else{assert(!button,'no phantom trigger '+entry.id);assert(env.opened===null,'no stale viewer');emptyCount++;}
    const detail=env.elements.get('result-detail').innerHTML;
    if(allowed){
     if(ac){
      const st=env.data.states.find(s=>s.id===id);
      const maximum=st.maximumACOrder;
      if(maximum){
       assert(detail.includes(`${maximum.status==='proven'?'Proved':'Conjectured'} maximum AC order: <strong>${maximum.order}</strong>`),'maximum order status shown '+id);
      }else assert(!detail.includes('<p class="maximum-order-tag'),'no maximum badge on lower/superseded state '+id);
      if(st.status==='numerical'){
       assert(detail.includes('Numerical') && detail.includes('residual'),'numerical evidence visible '+id);
       assert(!detail.includes('Exact reduced-state verification'),'numerical witness never labelled exact '+id);
       assert(!detail.includes('Why the AC order is exactly'),'no exact next-order claim '+id);
       assert(!detail.includes('Exact amplitude'),'numerical coefficients not called exact '+id);
      }
     }
     assert(detail.includes(env.win.MajoranaViewer.groupSummary(records[id])),'matching full/proper group summary beside representative '+id);
     assert(detail.includes('Point group')&&detail.includes('Rotations'),'both group roles beside representative '+id);
     assert(detail.includes('This representative'),'group scoped to representative '+id);
     if(ac){
      const options=[...detail.matchAll(/<option value="([^"]+)"[^>]*>([^<]*)<\/option>/g)];
      for(const candidate of candidateIds){const option=options.find(m=>m[1]===candidate);assert(option&&option[2].includes(env.win.MajoranaViewer.groupLabel(records[candidate])),'AC option includes own group '+candidate);}
     }
    }else assert(!detail.includes('Point group'),'no group assigned without representative '+entry.id);
    env.elements.get(ac?'download-record':'download-result').click();const payload=JSON.parse(await env.blob.text());
    assert(JSON.stringify(payload.majoranaRepresentation)===JSON.stringify(allowed?records[id]:null),'download coordinates '+entry.id);downloads++;
    if(ac&&candidateIds.length>1){const alternate=candidateIds.find(x=>x!==id),selector=env.elements.get('witness-select');selector.handlers.change({target:{value:alternate}});env.elements.get('view-majorana').click();assert(env.opened===records[alternate],'representative switch '+alternate);assert(env.api.getSelected().stateId===alternate,'selected state updates');assert(env.elements.get('result-detail').innerHTML.includes(env.win.MajoranaViewer.groupSummary(records[alternate])),'group summary switches with representative '+alternate);switches++;}
   }
  }
  let constructionCardsChecked=0;
  if(ac)for(const match of env.elements.get('construction-list').innerHTML.matchAll(/<article class="construction-card">([\s\S]*?)<\/article>/g)){const id=match[1].match(/data-state-link="([^"]+)"/)?.[1];assert(id&&records[id],'construction links a coordinate record');assert(match[1].includes(env.win.MajoranaViewer.groupSummary(records[id])),'construction card shows matching symmetry '+id);constructionCardsChecked++;}
  const offline=session(name,false);assert(!offline.elements.get('view-majorana'),'missing supplemental data no broken trigger');
  summaries.push({catalogue:name,entries:env.data.entries.length,representativeSelections:selectedCount,emptyCellsChecked:emptyCount,coordinateDownloads:downloads,alternativeSwitches:switches,representativeGroupSummariesChecked:selectedCount,ACGroupOptionChecks:ac,constructionCardsChecked,worksWithoutSupplementalData:true});
 }
 const report={passed:true,assertions,catalogues:summaries,checks:['every catalogue entry selects without runtime errors','viewer receives exact selected state coordinate record','alternative AC representative selection updates viewer','exclusion and unresolved cells have no phantom stars','downloaded record includes matching Majorana representation','selection change disposes prior viewer','file:// paths supported and optional coordinates missing handled','full and proper groups shown beside every representative','AC options show each alternative group and cards update on selection','spherical-design comparison badges only on three differing row maxima','comparison starts collapsed and links correct cells','reported orders remain distinct from true maxima and selected constellation design strength'],limitation:'Minimal DOM harness; browser visual QA separate.'};
 fs.writeFileSync(path.join(__dirname,'app-audit.json'),JSON.stringify(report,null,2)+'\n');console.log(JSON.stringify(report,null,2));
})().catch(e=>{console.error(e);process.exitCode=1;});
