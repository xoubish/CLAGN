const TAB_KEYS=['all','sep23','observed','oct26','oct27','about'];
const TAB_LABELS={all:'All',sep23:'Sep23-targets',observed:'Observed',oct26:'Oct26',oct27:'Oct27',about:'About'};
function tabFromHash(){const key=location.hash.slice(1);return ({pool:'all',sequence:'sep23',october:'oct26','oct26-sequence':'oct26','oct27-sequence':'oct27'})[key]||(TAB_KEYS.includes(key)?key:'all');}
function setView(view){
  if(!TAB_KEYS.includes(view))view='all';
  state.view=view;state.night=['sep23','oct26','oct27'].includes(view)?view:'all';state.show='all';
  state.q='';$('search').value='';refresh();history.pushState(null,'','#'+view);
  $('tab-content').scrollIntoView({behavior:'smooth',block:'start'});
}
function revealTarget(name){
  if(!D.all_targets.some(t=>t.name===name))return;
  if(!document.getElementById('card-'+name)){
    state.view='all';state.night='all';state.q='';$('search').value='';refresh();
  }
  document.getElementById('card-'+name)?.scrollIntoView({behavior:'smooth',block:'start'});
}
function catalogRows(){
  return D.all_targets.filter(t=>!state.q||`${t.name} ${t.jwst_target} ${t.jwst_id}`.toLowerCase().includes(state.q))
    .sort((a,b)=>state.sort==='brightness'?(a.r??99)-(b.r??99)||a.ra-b.ra:a.ra-b.ra);
}
function extraCard(t){
  const r=(D.october?.rows||[]).find(r=>r.name===t.name);
  return `<article class="card" id="card-${esc(t.name)}"><div class="card-title"><div><h2>${esc(t.name)} <span class="small muted">${esc(t.jwst_target)}</span></h2><div class="mono small">${esc(t.ra_hms)} · ${esc(t.dec_dms)}</div></div><span class="badge seq">${esc(t.jwst_id||'October target')}</span></div><div class="detail-note"><p>z = ${fmt(t.z,4)} · cached r = ${fmt(t.r,2)}</p>${r?queueTelescope(r):''}<p class="plot-note">This target was added to the October queue. Detailed light curves, imaging and archival spectra are not included in the cached explorer for this object.</p></div></article>`;
}
function queueTelescope(r){
  const backups=(D.october?.backups||[]).filter(b=>b.night===r.night&&b.replaces_name===r.name);
  return `<div class="tele"><strong>${r.night==='oct26'?'October 26':'October 27'} · visit ${r.order}${r.jwst_id?' · JWST '+esc(r.jwst_id):''}</strong><p>${r.start_pdt.slice(11,16)}–${r.end_pdt.slice(11,16)} PDT · 2 × 300 s · 12 min including overhead. Latest feasible start ${r.latest_start_pdt.slice(11,16)} PDT.</p><p>Max X ${fmt(+r.airmass_max_actual,2)} · Moon ≥${fmt(+r.moon_min,1)}°.</p><p>${esc(r.priority_reason)}. ${esc(r.flags)}</p>${backups.length?`<details><summary>Alternatives for this time slot</summary><p>${backups.map(b=>`${esc(b.name)} (r ${fmt(+b.r_adopted,2)})`).join(' · ')}</p><p class="plot-note">Skip any already observed. A replacement does not cover a skipped JWST visit.</p></details>`:''}</div>`;
}
function targetCard(name){
  const t=D.targets.find(t=>t.name===name);
  return t?card(t):extraCard(D.all_targets.find(t=>t.name===name));
}
function renderAll(){
  const list=catalogRows(),names=new Set(list.map(t=>t.name)),detailed=D.targets.filter(t=>names.has(t.name));
  $('pool-night-label').textContent=`${list.length} of ${D.all_targets.length} candidates`;
  $('count').textContent=`${list.length} shown`;
  $('bigmap').innerHTML=bigMap(detailed);
  $('night-summary').innerHTML=`<p>${D.all_targets.length} unique candidates, including all 24 JWST proposal targets.</p><p class="small muted">${D.targets.length} have detailed cached cards. Additional October targets are listed with their coordinates and queue notes. The manifold shows only targets with cached projections.</p><p class="small muted">Sep23-targets is the planned list; Observed contains the 18 targets actually observed. Oct26 and Oct27 show the proposed 40-target queues.</p>`;
  $('table').innerHTML=`<table><thead><tr>${['TARGET','RA (H:M:S)','DEC (D:M:S)','CACHED r','z','JWST','SEP 23 PLAN','OBSERVED','OCTOBER QUEUE'].map(h=>`<th>${h}</th>`).join('')}</tr></thead><tbody>${list.map(t=>`<tr data-catalog-target="${esc(t.name)}"><td><a href="#card-${esc(t.name)}">${esc(t.name)}</a></td><td class="mono">${esc(t.ra_hms)}</td><td class="mono">${esc(t.dec_dms)}</td><td>${fmt(t.r,2)}</td><td>${fmt(t.z,4)}</td><td>${esc(t.jwst_id||'—')}</td><td>${t.sep23?'Yes':'—'}</td><td>${t.observed?'Yes':'—'}</td><td>${t.scheduled_night?TAB_LABELS[t.scheduled_night]:'—'}</td></tr>`).join('')}</tbody></table>`;
  $('cards').innerHTML=list.map(t=>targetCard(t.name)).join('');
}
function refresh(){
  const view=state.view;
  ['pool','sequence','observed','october','about','cards'].forEach(id=>$(id).classList.add('hidden'));
  $('tab-content').setAttribute('aria-labelledby','tab-'+view);
  document.querySelectorAll('[data-view]').forEach(b=>{const active=b.dataset.view===view;b.setAttribute('aria-selected',active);b.setAttribute('aria-pressed',active);b.tabIndex=active?0:-1;});
  $('tab-summary').classList.toggle('hidden',view==='about');
  const counts={all:D.all_targets.length,sep23:Object.keys(SEQ).length,observed:D.observed.length,oct26:40,oct27:40};
  const descriptions={all:'Complete candidate catalogue, including the JWST additions. The CSV contains the full catalogue with default settings; use a night tab for an observing sequence.',sep23:'The 14 science targets planned for September 23. Standards and slot alternatives are shown below; the CSV contains the science targets only. Actual observations are in Observed.',observed:'The 18 science targets actually observed on September 23, in exposure order. The CSV is an observation log with actual exposure totals and measured S/N.',oct26:'40 science targets in observing order, including 11 JWST targets. The CSV is the science NGPS list; standards and alternatives are noted below.',oct27:'40 science targets in observing order, including 12 JWST targets. The CSV is the science NGPS list; standards and alternatives are noted below.'};
  if(view!=='about')$('tab-summary').innerHTML=`<div class="tab-heading"><div><h2>${TAB_LABELS[view]} <span class="count">${counts[view]} targets</span></h2><p class="muted">${descriptions[view]}</p></div><button class="tab-download" data-download="${D.tab_downloads[view]}">DOWNLOAD ${TAB_LABELS[view].toUpperCase()} CSV</button></div>`;
  if(view==='about'){$('about').classList.remove('hidden');return;}
  $('cards').classList.remove('hidden');
  if(view==='observed'){$('observed').classList.remove('hidden');$('cards').innerHTML=observedTargets().map(observedCard).join('');return;}
  if(view==='sep23'){$('sequence').classList.remove('hidden');$('cards').innerHTML=Object.entries(SEQ).sort((a,b)=>a[1].rank-b[1].rank).map(([name])=>targetCard(name)).join('');return;}
  if(view==='oct26'||view==='oct27'){$('october').classList.remove('hidden');fillOctober(view);$('cards').innerHTML=D.october.rows.filter(r=>r.night===view).map(r=>targetCard(r.name)).join('');return;}
  $('pool').classList.remove('hidden');renderAll();
}
