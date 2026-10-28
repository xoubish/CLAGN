const TAB_KEYS=['about','sep23','observed','oct26','observed_oct26','oct27','observed_oct27'];
const TAB_LABELS={sep23:'Sep23-targets',observed:'Sep23-observed',oct26:'Oct26',observed_oct26:'Oct26-observed',oct27:'Oct27',observed_oct27:'Oct27-observed',about:'About'};
function tabFromHash(){const key=location.hash.slice(1);return ({all:'about',pool:'about',sequence:'sep23',october:'oct26','oct26-sequence':'oct26','oct27-sequence':'oct27'})[key]||(TAB_KEYS.includes(key)?key:'about');}
function setView(view){
  if(!TAB_KEYS.includes(view))view='about';
  state.view=view;state.night=['sep23','oct26','oct27'].includes(view)?view:'all';state.show='all';
  state.q='';$('search').value='';refresh();history.pushState(null,'','#'+view);
  $('tab-content').scrollIntoView({behavior:'smooth',block:'start'});
}
function revealTarget(name){
  if(!D.all_targets.some(t=>t.name===name))return;
  if(state.view==='about'||!document.getElementById('card-'+name)){
    state.selected=name;
    state.view='about';state.night='all';state.q='';$('search').value='';refresh();
  }
  document.getElementById('card-'+name)?.scrollIntoView({behavior:'smooth',block:'start'});
}
function catalogRows(){
  return D.all_targets.filter(t=>!state.q||`${t.name} ${t.display_name} ${t.target_label}`.toLowerCase().includes(state.q))
    .sort((a,b)=>state.sort==='brightness'?(a.r??99)-(b.r??99)||a.ra-b.ra:a.ra-b.ra);
}
function extraCard(t){
  const r=(D.october?.rows||[]).find(r=>r.name===t.name);
  return `<article class="card" id="card-${esc(t.name)}"><div class="card-title"><div><h2>${esc(displayName(t.name))} <span class="small muted">${esc(t.target_label)}</span></h2><div class="mono small">${esc(t.ra_hms)} · ${esc(t.dec_dms)}</div></div><span class="badge seq">October target</span></div><div class="detail-note"><p>z = ${fmt(t.z,4)} · cached r = ${fmt(t.r,2)}</p>${r?queueTelescope(r):''}<p class="plot-note">This target was added to the October queue. Detailed light curves, imaging and archival spectra are not included in the cached explorer for this object.</p></div></article>`;
}
function queueTelescope(r){
  const backups=(D.october?.backups||[]).filter(b=>b.night===r.night&&b.replaces_name===r.name);
  return `<div class="tele"><strong>${r.night==='oct26'?'October 26':'October 27'} · visit ${r.order}</strong><p>${r.start_pdt.slice(11,16)}–${r.end_pdt.slice(11,16)} PDT · ${r.exposures} × ${r.seconds_each} s · ${r.visit_minutes} min including overhead. Latest feasible start ${r.latest_start_pdt.slice(11,16)} PDT.</p><p>Max X ${fmt(+r.airmass_max_actual,2)} · Moon ≥${fmt(+r.moon_min,1)}°.</p><p>${esc(r.priority_reason)}. ${esc(r.science_question||'')} ${esc(r.flags)}</p>${backups.length?`<details><summary>Alternatives for this time slot</summary><p>${backups.map(b=>`${esc(b.name)} (r ${fmt(+b.r_adopted,2)})`).join(' · ')}</p><p class="plot-note">Skip any already observed. Review priorities before replacing a target; preserve requested H targets and the P8544 / HIRAS-F23290-1252 exposure depth.</p></details>`:''}</div>`;
}
function targetCard(name){
  const field=(D.field_cards||[]).find(t=>t.name===name);
  if(field)return stellarCard(field);
  const t=D.targets.find(t=>t.name===name);
  return t?card(t):extraCard(D.all_targets.find(t=>t.name===name));
}
function renderCatalogue(){
  const list=catalogRows(),names=new Set(list.map(t=>t.name));
  const detailed=D.targets.filter(t=>Number.isFinite(t.ux)&&Number.isFinite(t.uy));
  $('pool-night-label').textContent=`${list.length} of ${D.all_targets.length} candidates`;
  $('count').textContent=`${list.length} shown`;
  $('bigmap').innerHTML=bigMap(detailed);
  $('night-summary').innerHTML=`<p>${D.all_targets.length} unique candidates.</p><p class="small muted">The manifold shows all ${detailed.length} targets with saved projections. The remaining ${D.all_targets.length-detailed.length} targets are in the catalogue below. Click a point or target name to open its details.</p><p class="small muted">Sep23-targets is the planned list; Sep23-observed, Oct26-observed and Oct27-observed contain the targets actually observed, with their new spectra. Oct26 and Oct27 show the observing queues in scheduled order.</p>`;
  $('table').innerHTML=`<table><thead><tr>${['TARGET','RA (H:M:S)','DEC (D:M:S)','CACHED r','z','SEP 23 PLAN','OBSERVED','OCTOBER QUEUE'].map(h=>`<th>${h}</th>`).join('')}</tr></thead><tbody>${list.map(t=>`<tr data-catalog-target="${esc(t.name)}"><td><a href="#card-${esc(t.name)}">${esc(displayName(t.name))}</a></td><td class="mono">${esc(t.ra_hms)}</td><td class="mono">${esc(t.dec_dms)}</td><td>${fmt(t.r,2)}</td><td>${fmt(t.z,4)}</td><td>${t.sep23?'Yes':'—'}</td><td>${t.observed?'Yes':'—'}</td><td>${t.scheduled_night?TAB_LABELS[t.scheduled_night]:'—'}</td></tr>`).join('')}</tbody></table>`;
  // Rendering all detailed spectra and SVG plots at once freezes the browser.
  // Keep the entire catalogue, but create only the explicitly selected card.
  $('cards').innerHTML=state.selected&&names.has(state.selected)?targetCard(state.selected):'';
}
function refresh(){
  const view=state.view;
  ['pool','sequence','observed','october','about','cards'].forEach(id=>$(id).classList.add('hidden'));
  $('tab-content').setAttribute('aria-labelledby','tab-'+view);
  document.querySelectorAll('[data-view]').forEach(b=>{const active=b.dataset.view===view;b.setAttribute('aria-selected',active);b.setAttribute('aria-pressed',active);b.tabIndex=active?0:-1;});
  $('tab-summary').classList.toggle('hidden',view==='about');
  const counts={all:D.all_targets.length,sep23:Object.keys(SEQ).length,observed:D.observed.length,observed_oct26:(D.observed_oct26||[]).length,observed_oct27:(D.observed_oct27||[]).length,oct26:40,oct27:40};
  const descriptions={sep23:'The 14 science targets planned for September 23. Standards and slot alternatives are shown below; the CSV contains the science targets only. Actual observations are in Observed.',observed:'The 18 science targets actually observed on September 23, in exposure order. The CSV is an observation log with actual exposure totals and measured S/N.',observed_oct26:'The science targets observed on October 26, plus Z16 Cyg, in exposure order. The CSV is an observation log with actual exposure totals and measured S/N.',observed_oct27:'The science targets observed on October 27, in exposure order. The CSV is an observation log with actual exposure totals and measured S/N.',oct26:'',oct27:''};
  for(const night of ['oct26','oct27']){
    const rows=[...(D.october?.rows||[]),...(D.october?.manual_visits||[])].filter(r=>r.night===night);
    counts[night]=rows.length;
    descriptions[night]=`${rows.length} science targets in observing order. Reviewed ${D.october?.review_date||''}. The CSV includes the opening standard, all science targets, and the closing standard in order.`;
  }
  if(view!=='about')$('tab-summary').innerHTML=`<div class="tab-heading"><div><h2>${TAB_LABELS[view]} <span class="count">${counts[view]} targets</span></h2><p class="muted">${descriptions[view]}</p></div><button class="tab-download" data-download="${D.tab_downloads[view]}">DOWNLOAD ${TAB_LABELS[view].toUpperCase()} CSV</button></div>`;
  if(view==='about'){
    ['pool','about','cards'].forEach(id=>$(id).classList.remove('hidden'));
    renderCatalogue();return;
  }
  $('cards').classList.remove('hidden');
  if(view==='observed'||view==='observed_oct26'||view==='observed_oct27'){$('observed').classList.remove('hidden');fillObserved(view);$('cards').innerHTML=observedCards(view);return;}
  if(view==='sep23'){$('sequence').classList.remove('hidden');$('cards').innerHTML=Object.entries(SEQ).sort((a,b)=>a[1].rank-b[1].rank).map(([name])=>targetCard(name)).join('');return;}
  if(view==='oct26'||view==='oct27'){$('october').classList.remove('hidden');fillOctober(view);const fields=D.field_cards||[],manual=(D.october.manual_visits||[]).map(r=>({...r,name:r.upload_name||r.name})).filter(r=>!fields.some(f=>f.night===r.night&&f.name===r.name));$('cards').innerHTML=[...D.october.rows,...fields,...manual].filter(r=>r.night===view).sort((a,b)=>a.start_pdt.localeCompare(b.start_pdt)).map(r=>targetCard(r.name)).join('');return;}
}

function stellarCard(t){
  return `<article class="card" id="card-${esc(t.name)}">
<div class="card-title"><div><h2>${esc(t.label)}</h2><div class="mono small">${esc(t.ra_hms)} · ${esc(t.dec_dms)}</div></div><div><div class="time mono">${t.start_pdt.slice(11,16)}–${t.end_pdt.slice(11,16)} PDT</div><span class="badge seq">October 26</span></div></div>
<div class="stellar-grid"><section><h3>At the telescope</h3><dl class="facts"><dt>Exposure</dt><dd>${t.exposures} × ${t.seconds_each} s · maximum ${t.max_exposure_seconds} s each</dd><dt>Slit / binning</dt><dd>${t.slit_width_arcsec}″ per slice · ${t.binspat} × ${t.binspect}</dd><dt>Fixed sky PA</dt><dd>${fmt(t.pa_deg,2)}° east of north</dd><dt>Max airmass</dt><dd>${fmt(t.airmass_max_actual,2)}</dd><dt>Moon distance</dt><dd>${fmt(t.moon_min,1)}°</dd><dt>Upload name</dt><dd>${esc(t.name)}</dd></dl><h3 style="margin-top:20px">Components</h3><table><thead><tr><th>STAR</th><th>GAIA G</th><th>SLICE</th></tr></thead><tbody>${t.components.map(s=>`<tr><td>Z16 Cyg ${s.name}</td><td>${fmt(s.Gmag,2)}</td><td>${s.name==='C'?'Opposite outer':'One outer'}</td></tr>`).join('')}</tbody></table><p class="plot-note" style="margin-top:10px">One pointing for A/B/C; C is the priority. Same slit width and binning as the other targets.</p></section>
<aside class="imaging"><h3>${esc(t.image_source)} · ${t.field_arcsec}″ field</h3><a href="october/16Cyg_field.svg" target="_blank" rel="noopener" aria-label="Enlarge 16 Cyg image and slit overlay">${t.field_svg}</a><p class="plot-note">Gold outlines: three 50″ × 1.5″ slices at PA ${fmt(t.pa_deg,2)}°. A/B/C markers predict October 26, 2026 positions; proper motion shifts them from the archival image.</p><p class="plot-note">${esc(t.image_status)} Confirm C and centering in the acquisition view.</p></aside></div>${stellarArchive(t.archive_review)}</article>`;
}

function stellarArchive(review){
  if(!review)return '';
  return `<section class="spectra stellar-archive"><h3>Archival optical spectra · checked ${esc(review.checked_date)}</h3><div class="table-scroll"><table><thead><tr><th>COMPONENT</th><th>CONFIRMED SPECTRUM</th><th>DOWNLOAD</th></tr></thead><tbody>${review.components.map(r=>`<tr><td>Z16 Cyg ${r.component}</td><td>CFHT / ESPaDOnS · ${r.date} · ${Math.round(r.wavelength_min)}–${Math.round(r.wavelength_max)} Å</td><td><a href="${esc(r.fits_url)}" target="_blank" rel="noopener">FITS</a> · <a href="${esc(r.gaia_url)}" download>Gaia BP/RP</a></td></tr>`).join('')}<tr><td>Z16 Cyg C</td><td colspan="2">No securely identified optical spectrum found in this search.</td></tr></tbody></table></div><p class="plot-note" style="margin-top:10px">CFHT curves are continuum-normalized spectra of A and B; Gaia curves are low-resolution mean spectra. No spectrum of A is assigned to close companion C.</p><details><summary>View A and B spectra</summary><img src="${esc(review.plot)}" alt="Confirmed optical spectra of 16 Cyg A and B: CFHT overview and H-alpha region, plus Gaia BP/RP mean spectra" loading="lazy" style="display:block;width:100%;height:auto;margin-top:12px"></details><details><summary>Archive search coverage</summary><p class="plot-note">${esc(review.search_summary)}</p><p class="plot-note">${esc(review.limitations)}</p><p class="small">${review.sources.map(s=>`<a href="${esc(s.url)}" target="_blank" rel="noopener">${esc(s.label)}</a>`).join(' · ')}</p></details></section>`;
}
