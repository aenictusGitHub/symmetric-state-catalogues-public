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
 const loc={hash:'',href:'file:///offline/'+name+'/index.html'},win={katex:require(path.join(root,'vendor/katex/katex.min.js')),addEventListener(){},MajoranaViewer:{open(r){opened=r;return true;},close(){opened=null;closes++;}}};
 const ctx={window:win,document:doc,location:loc,history:{replaceState(a,b,h){loc.hash=h;}},URLSearchParams,URL:FakeURL,Blob,console,setTimeout:()=>0,clearTimeout(){},navigator:{},matchMedia:()=>({matches:false})};
 vm.runInNewContext(fs.readFileSync(path.join(root,'data/catalogue.js'),'utf8'),ctx);
 if(withCoordinates)vm.runInNewContext(fs.readFileSync(path.join(root,'data/majorana.js'),'utf8'),ctx);
 const modalInterceptors=win.MajoranaViewer;
 for(const file of ['majorana-symmetry.js','symmetry-groups.js','point-groups.js','majorana-viewer.js'])vm.runInNewContext(fs.readFileSync(path.join(root,file),'utf8'),ctx,{filename:file});
 win.MajoranaViewer={...win.MajoranaViewer,...modalInterceptors};
 const expose=isAC?'window.__audit={select,getSelected:()=>selected,renderDetail};':'window.__audit={select(n,k,m){selectedN=n;selectedK=k;mode=m;renderTable();renderDetail();updateHash();},getSelected:selected,renderDetail};';
 const source=fs.readFileSync(path.join(root,'app.js'),'utf8').replace(/\}\)\(\);\s*$/,expose+'})();');vm.runInNewContext(source,ctx);
 return {elements,api:win.__audit,win,get opened(){return opened;},get blob(){return lastBlob;},get closes(){return closes;},data:isAC?win.AC_CATALOGUE:win.SECTOR_CATALOGUE};
}
(async()=>{
 const summaries=[];
 for(const name of [path.basename(catalogueRoot)]){
  const env=session(name),ac=Object.hasOwn(env.data,'states'),records=env.win.MAJORANA_DATA.records;let selectedCount=0,emptyCount=0,downloads=0,switches=0;
  for(const entry of env.data.entries){
   const candidateIds=ac?entry.stateIds||[]:[entry.id],cases=candidateIds.length?candidateIds:[null];
   for(const id of cases){
    const previous=env.closes;ac?env.api.select(entry.n,entry.t,id):env.api.select(entry.n,entry.k,entry.mode);
    assert(env.closes>previous,'rerender closes viewer');
    const button=env.elements.get('view-majorana');const allowed=!ac||entry.status==='exists'&&id;
    if(allowed){assert(button,'viewer trigger for '+id);button.click();assert(env.opened===records[id],'trigger passes exact representative '+id);assert(env.opened.n===entry.n,'qubit count consistent');selectedCount++;}
    else{assert(!button,'no phantom trigger '+entry.id);assert(env.opened===null,'no stale viewer');emptyCount++;}
    const detail=env.elements.get('result-detail').innerHTML;
    if(allowed){
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
 const report={passed:true,assertions,catalogues:summaries,checks:['every catalogue entry selects without runtime errors','viewer receives exact selected state coordinate record','alternative AC representative selection updates viewer','exclusion and unresolved cells have no phantom stars','downloaded record includes matching Majorana representation','selection change disposes prior viewer','file:// paths supported and optional coordinates missing handled','full and proper groups shown beside every representative','AC options show each alternative group and cards update on selection'],limitation:'Minimal DOM harness; browser visual QA separate.'};
 fs.writeFileSync(path.join(__dirname,'app-audit.json'),JSON.stringify(report,null,2)+'\n');console.log(JSON.stringify(report,null,2));
})().catch(e=>{console.error(e);process.exitCode=1;});
