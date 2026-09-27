/* October queues use the new 12-minute geometry; the older pool stays a diagnostic. */
function fillOctober(){
  const O=D.october;
  if(!O){$('october').classList.add('hidden');return;}
  const download=(name,label)=>`<button data-download="${name}">${label}</button>`;
  const summary=`<h2>October 26–27 · 40 targets per night</h2><p>80 distinct targets · <b>2 × 300 s · 1.5″ slit · 2 × 3 binning · 2 min overhead</b> · Moon ≥40° throughout every visit. Science takes 8 hours per night; standards reserve another 20 minutes.</p><p class="plot-note">23 of 24 JWST targets are scheduled below, including four September repeats. P11530 is Moon-blocked and retains September coverage. Cached magnitudes generally end in 2025; faint JWST targets are flagged. The dawn sequences are tight: protect JWST visits if running late. NGPS files specify order and settings; use the scheduled times and windows to preserve gaps.</p><div class="actions">${download('all_80_targets.csv','ALL 80 TARGETS CSV')}${download('jwst_coverage.csv','JWST COVERAGE CSV')}${download('october_visibility_windows.csv','VISIBILITY WINDOWS CSV')}${download('free_time.csv','FREE TIME CSV')}${download('slot_backups.csv','SLOT BACKUPS CSV')}</div>`;
  const nights=['oct26','oct27'].map(n=>{
    const rows=O.rows.filter(r=>r.night===n),meta=O.nights[n],jwst=rows.filter(r=>r.jwst_id).length;
    const visits=O.standards.filter(r=>r.night===n).map(r=>`${esc(r.name)} ${r.start_pdt.slice(11,16)}–${r.end_pdt.slice(11,16)} (${r.exposures} × ${r.seconds_each} s initial)`).join('; ');
    return `<section id="${n}-sequence" class="october-night"><h3>${n==='oct26'?'October 26':'October 27'} · 40 targets · ${jwst} JWST</h3><p class="small muted">PDT ${meta.dusk_pdt.slice(11,16)}–${meta.dawn_pdt.slice(11,16)} · Moon ${Math.round(meta.moon_illumination_percent)}% · ${Math.round(meta.remaining_minutes)} min spare after standards. Standards: ${visits}.</p><div class="actions">${download(n+'_ngps.csv','SCIENCE NGPS CSV')}${download(n+'_sequence.csv','TIMING + TARGET DETAILS CSV')}${download(n+'_standards_ngps.csv','STANDARDS NGPS CSV')}${download(n+'_backups_ngps.csv','BACKUPS NGPS CSV')}${download(n+'_coordinates.csv','COORDINATES CSV · NO HEADER')}</div><div class="table-scroll"><table><thead><tr>${['#','VISIT (PDT)','UTC START','TARGET','JWST','RA (H:M:S)','DEC (D:M:S)','CACHED r','MAX X','MIN MOON','PRIORITY / NOTES'].map(h=>`<th>${h}</th>`).join('')}</tr></thead><tbody>${rows.map(r=>{
      const inPool=D.targets.some(t=>t.name===r.name);
      const name=inPool?`<a href="#card-${esc(r.name)}">${esc(r.name)}</a>`:esc(r.name);
      return `<tr data-october-target="${esc(r.name)}" class="${r.jwst_id?'seq':''}"><td>${r.order}</td><td class="mono">${r.start_pdt.slice(11,16)}–${r.end_pdt.slice(11,16)}</td><td>${r.start_utc.slice(11,16)}</td><td>${name}</td><td>${esc(r.jwst_id||'—')}</td><td class="mono">${esc(r.ra_hms)}</td><td class="mono">${esc(r.dec_dms)}</td><td>${fmt(+r.r_adopted,2)}</td><td>${fmt(+r.airmass_max_actual,2)}</td><td>${fmt(+r.moon_min,1)}°</td><td class="october-notes"><b>${esc(r.priority_reason)}</b>${r.flags?'<br>'+esc(r.flags):''}</td></tr>`;
    }).join('')}</tbody></table></div></section>`;
  }).join('');
  $('october').innerHTML=summary+nights;
}
