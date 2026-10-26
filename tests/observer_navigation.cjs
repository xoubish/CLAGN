// Execute the generated page's real JavaScript with a minimal DOM. Catch eager
// rendering of the complete spectral catalogue, which previously froze All.
const fs=require('fs'),vm=require('vm'),assert=require('assert'),path=require('path');
const html=fs.readFileSync(path.join(__dirname,'../docs/index.html'),'utf8');
const payload=html.match(/<script id="candidate-data" type="application\/json">([\s\S]*?)<\/script>/)[1];
const ids=new Set([...html.matchAll(/\bid="([^"]+)"/g)].map(m=>m[1]));
const elements=new Map();
function element(id){
  if(!ids.has(id)&&!elements.get('cards')?.innerHTML.includes(`id="${id}"`))return null;
  if(!elements.has(id)){
    const classes=new Set();
    elements.set(id,{textContent:id==='candidate-data'?payload:'',innerHTML:'',value:'',dataset:{},
      classList:{add(c){classes.add(c);},remove(c){classes.delete(c);},toggle(c,on){on?classes.add(c):classes.delete(c);},contains(c){return classes.has(c);}},
      insertAdjacentHTML(position,s){this.innerHTML+=s;},setAttribute(){},addEventListener(){},scrollIntoView(){}});
  }
  return elements.get(id);
}
const document={getElementById:element,querySelectorAll:()=>[],querySelector:()=>({addEventListener(){}}),
  addEventListener(){},documentElement:{dataset:{theme:'dark'}}};
const location={hash:''};
const context=vm.createContext({document,window:{addEventListener(){}},location,
  history:{pushState(){},replaceState(){}},console,URL,Blob,setTimeout});
const script=[...html.matchAll(/<script(?:\s[^>]*)?>([\s\S]*?)<\/script>/g)]
  .map(m=>m[1]).filter(s=>!s.trimStart().startsWith('{')).join('\n');
vm.runInContext(script,context);
const evaluate=s=>vm.runInContext(s,context);
assert.equal(evaluate('state.view'),'about');
assert.equal(element('cards').innerHTML,'');
assert.equal((element('table').innerHTML.match(/data-catalog-target=/g)||[]).length,evaluate('D.all_targets.length'));
const nav=html.match(/<nav[\s\S]*?<\/nav>/)[0];
assert.deepEqual([...nav.matchAll(/data-view="([^"]+)"/g)].map(m=>m[1]),['about','sep23','observed','oct26','oct27']);
for(const [tab,count] of [['sep23',14],['observed',18],['oct26',40],['oct27',40],['about',0]]){
  context.setView(tab);
  assert.equal((element('cards').innerHTML.match(/<article /g)||[]).length,count,tab);
}
context.revealTarget('P9694');
assert.equal((element('cards').innerHTML.match(/<article /g)||[]).length,1);
assert(element('cards').innerHTML.includes('id="card-P9694"'));
context.revealTarget('P8544');
assert.equal((element('cards').innerHTML.match(/<article /g)||[]).length,1);
assert(element('cards').innerHTML.includes('id="card-P8544"'));
evaluate("state.q='p9694';refresh()");
assert.equal((element('table').innerHTML.match(/data-catalog-target=/g)||[]).length,1);
assert.equal((element('bigmap').innerHTML.match(/data-name=/g)||[]).length,
  evaluate('D.targets.filter(t=>Number.isFinite(t.ux)&&Number.isFinite(t.uy)).length'));
location.hash='#all';assert.equal(context.tabFromHash(),'about');
location.hash='#pool';assert.equal(context.tabFromHash(),'about');
console.log('PASS: default About, five tabs, complete catalogue and manifold, one target card at a time, search and legacy links.');

for(const [night,n] of [['sep23',14],['oct26',40],['oct27',40]]){
  assert.equal((context.renderVisibility(night).match(/data-visibility-target=/g)||[]).length,n);
  assert(context.visibilityReadout(night,0).includes('0 of '+n));
}
evaluate("state.q='jp9694';refresh()");
assert.equal((element('table').innerHTML.match(/data-catalog-target=/g)||[]).length,1);
assert(element('table').innerHTML.includes('>JP9694</a>'));
assert(element('cards').innerHTML.includes('id="card-P9694"') || element('cards').innerHTML==='');
console.log('PASS: three visibility charts and both original and prefixed name search.');

if(evaluate("(D.october.manual_visits||[]).length")){
  context.setView('oct26');
  assert(element('october').innerHTML.includes('data-manual-visit="Z16 Cyg A/B/C"'));
  assert(element('october').innerHTML.includes('not included in the executable NGPS CSV'));
  assert(element('october').innerHTML.includes('20:03–20:17'));
  assert(evaluate("D.visibility.oct26.standards.some(r=>r.name==='Z16 Cyg A/B/C'&&r.kind==='manual visit')"));
  assert(evaluate("D.files['oct26_ngps.csv'].includes('PAUSE after this visit')"));
  console.log('PASS: reserved manual acquisition, explicit exposure warning, pause note and blocked visibility interval.');
}

// Optional local archive additions must render with unknown redshift/spectra.
if(evaluate("D.targets.some(t=>t.name==='H5MUSES-327')")){
  context.revealTarget('H5MUSES-327');
  assert(element('cards').innerHTML.includes('id="card-H5MUSES-327"'));
  assert(element('cards').innerHTML.includes('Archive review'));
  assert(element('cards').innerHTML.includes('stellar'));
  assert(!element('cards').innerHTML.includes('NaN'));
  assert(element('bigmap').innerHTML.includes('data-name="H5MUSES-327"'));
  evaluate("state.q='5muses';refresh()");
  assert.equal((element('table').innerHTML.match(/data-catalog-target=/g)||[]).length,1);
  console.log('PASS: requested target, unknown redshift, archive review, manifold marker and search.');
}

if(evaluate("D.targets.some(t=>t.name==='HUGC-03097')")){
  for(const name of ['HNGC-1377','HIC-1953','HUGC-03097','HIRAS-F23290-1252','HNGC-0814']){
    context.revealTarget(name);
    assert(element('cards').innerHTML.includes(`id="card-${name}"`));
    assert(element('cards').innerHTML.includes('Legacy Surveys'));
    assert(element('cards').innerHTML.includes('No saved W1 projection'));
    assert(!element('cards').innerHTML.includes('NaN'));
  }
  context.revealTarget('HUGC-03097');
  assert(element('cards').innerHTML.includes('No Moon-safe window'));
  assert(!evaluate("D.october.rows.some(r=>r.name==='HUGC-03097')"));
  console.log('PASS: five proposal cards, unavailable projections, actual image surveys and Moon-blocked target.');
}

if(evaluate("D.targets.some(t=>t.name==='HNGC-1377'&&t.ztf_archive_matches)")){
  for(const name of ['HNGC-1377','HIC-1953','HUGC-03097','HIRAS-F23290-1252','HNGC-0814']){
    const t=JSON.parse(evaluate(`JSON.stringify(D.targets.find(t=>t.name==='${name}'))`));
    assert(t.ztf.g.length>0&&t.ztf.r.length>0);
    assert(Object.values(t.ztf_archive_matches).every(m=>m.separation_arcsec<=2));
    assert(t.archive_review.note.includes('intrinsic variability has not been established'));
  }
  context.revealTarget('HIRAS-F23290-1252');
  assert(element('cards').innerHTML.includes('was rejected'));
  assert.equal(evaluate("D.targets.find(t=>t.name==='HIRAS-F23290-1252').spec"),null);
  context.revealTarget('HNGC-1377');
  assert(element('cards').innerHTML.includes('Flux / own median (relative)'));
  assert(element('cards').innerHTML.includes('unverified observing dates'));
  assert(!element('cards').innerHTML.includes('NaN'));
  evaluate("spectralSelection.set('HNGC-1377',new Set([1]))");
  const calibrated=evaluate("specPanel(D.targets.find(t=>t.name==='HNGC-1377'))");
  assert(calibrated.includes('10⁻¹⁷ erg'));
  assert(!calibrated.includes('Flux / own median'));
  const sdss=JSON.parse(evaluate("JSON.stringify(D.targets.find(t=>t.name==='HNGC-0814'))"));
  assert.equal(sdss.n_spec,1);
  assert.equal(sdss.spec.epochs.length,1);
  assert.equal(sdss.spec.epochs[0].date,'2024-12-08');
  assert(Math.abs(sdss.spec.epochs[0].meta.z-sdss.z)<.001);
  console.log('PASS: ZTF matches, rejected DESI neighbor, honest spectral dates, relative/physical units and SDSS-V association.');
}
