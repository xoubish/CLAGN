/* Every available epoch is selectable. Labels are dates, never survey names. */
const spectralSelection = new Map();
function spectralColor(i) { return `hsl(${(210+i*137.508)%360} 58% ${document.documentElement.dataset.theme==='dark'?65:40}%)`; }
function selectedSpectra(t) {
  const all=(t.spec?.epochs||[]).map((_,i)=>i);
  return spectralSelection.has(t.name)?all.filter(i=>spectralSelection.get(t.name).has(i)):all;
}
function specPanel(t) {
  const epochs=t.spec?.epochs||[], selected=selectedSpectra(t);
  const dates=new Set(epochs.map(e=>e.epoch_day).filter(d=>d!=null));
  const missing=t.epochs.filter(e=>!e.file_available).length;
  let summary=`<p class="small spectrum-count"><b>${selected.length} of ${epochs.length} available spectra shown</b> · ${dates.size} dates with files · ${t.n_spec} dates identified.</p>`;
  if(missing)summary+=`<p class="small muted">${missing} identified ${missing===1?'date has':'dates have'} no spectrum file loaded yet.</p>`;
  if(!epochs.length)return summary+'<div class="empty">No spectrum files loaded for this source yet.</div>';
  let controls=`<div class="spectrum-actions"><button data-spec-action="all">Show all</button><button data-spec-action="ends">First &amp; last</button><button data-spec-action="none">Clear</button></div><div class="spectrum-choices">`;
  controls+=epochs.map((e,i)=>`<label title="${esc([e.quality_note,e.coadd?'Combined spectrum; date may represent an observing interval.':''].filter(Boolean).join(' '))}"><input type="checkbox" data-spectrum-index="${i}" ${selected.includes(i)?'checked':''}><i style="background:${spectralColor(i)}"></i><span>${esc(e.label)}</span>${e.quality_flag?'<span class="quality-flag">check quality</span>':''}</label>`).join('')+'</div>';
  if(!selected.length)return summary+controls+'<div class="empty">Select a date to plot its spectrum.</div>';
  const traces=selected.map(i=>({i,e:epochs[i]}));
  const values=traces.flatMap(({e})=>e.flux.filter(Number.isFinite)).sort((a,b)=>a-b);
  const waves=traces.flatMap(({e})=>[e.wave[0],e.wave[e.wave.length-1]]);
  if(!values.length)return summary+controls+'<div class="empty">The selected spectra have no finite flux values.</div>';
  let lo=Math.min(0,values[Math.floor(values.length*.01)]),hi=values[Math.min(values.length-1,Math.floor(values.length*.995))];
  if(hi<=lo)hi=lo+1;
  const pad=(hi-lo)*.05;lo-=pad;hi+=pad;
  const W=900,H=295,L=53,R=15,T=40,B=37;
  const x0=Math.min(...waves),x1=Math.max(...waves),X=lin(x0,x1,L,W-R),Y=lin(lo,hi,H-B,T);
  const clipId=`spectral-clip-${t.name}`;
  let svg=`<svg viewBox="0 0 ${W} ${H}" width="100%" role="img" aria-label="${selected.length} selected spectra of ${esc(t.name)}"><defs><clipPath id="${esc(clipId)}"><rect x="${L}" y="${T}" width="${W-L-R}" height="${H-B-T}"/></clipPath></defs>`;
  niceTicks(lo,hi,4).forEach(v=>{svg+=`<line class="grid" x1="${L}" x2="${W-R}" y1="${Y(v)}" y2="${Y(v)}"/><text x="${L-6}" y="${Y(v)+4}" text-anchor="end">${Math.abs(v)<10?v.toFixed(1):v.toFixed(0)}</text>`});
  if(t.z!=null)[['Mg II',2798,0],['Hβ',4861.33,1],['[O III]',5006.84,0],['Hα',6562.8,1]].forEach(([name,rest,row])=>{
    const obs=rest*(1+t.z);if(obs<x0||obs>x1)return;
    const x=X(obs);svg+=`<line x1="${x}" x2="${x}" y1="${T}" y2="${H-B}" stroke="var(--hair2)"/><text x="${x}" y="${row?T-5:T-20}" text-anchor="middle">${name}</text>`;
  });
  svg+=`<g clip-path="url(#${esc(clipId)})">`;
  for(const {i,e} of traces){let path='',pen=false;e.flux.forEach((f,k)=>{if(!Number.isFinite(f)||!Number.isFinite(e.wave[k])){pen=false;return}path+=`${pen?'L':'M'}${X(e.wave[k]).toFixed(1)},${Y(f).toFixed(1)}`;pen=true});svg+=`<path data-epoch-trace="${i}" d="${path}" fill="none" stroke="${spectralColor(i)}" stroke-width="1.15" opacity="${selected.length>8?.65:.9}"><title>${esc(e.label)}</title></path>`}
  svg+='</g>';
  niceTicks(x0,x1,7).forEach(v=>{svg+=`<text x="${X(v)}" y="${H-B+17}" text-anchor="middle">${v}</text>`});
  svg+=`<text x="${W-R}" y="${H-3}" text-anchor="end">Observed wavelength (Å)</text><text x="${L+8}" y="${T+15}">Fλ · 10⁻¹⁷ erg s⁻¹ cm⁻² Å⁻¹</text></svg>`;
  return summary+controls+svg+'<p class="plot-note">Dates are UTC observing dates; exact exposure times are not available for every combined spectrum. Date labels do not imply that spectra have been matched in aperture or flux calibration.</p>';
}
function redrawSpectra(name=state.selected) {
  const t=D.targets.find(t=>t.name===name), panel=document.querySelector(`[data-spectrum-target="${CSS.escape(name)}"] .spectrum-panel`)||$('spectrum-panel');
  if(t&&panel)panel.innerHTML=specPanel(t);
}
document.addEventListener('change',e=>{
  const box=e.target.closest('[data-spectrum-index]');if(!box)return;
  const name=box.closest('[data-spectrum-target]')?.dataset.spectrumTarget||state.selected;
  const t=D.targets.find(t=>t.name===name);if(!t)return;
  const chosen=new Set(selectedSpectra(t)),index=Number(box.dataset.spectrumIndex);
  if(box.checked)chosen.add(index);else chosen.delete(index);
  spectralSelection.set(t.name,chosen);redrawSpectra(t.name);
});
document.addEventListener('click',e=>{
  const button=e.target.closest('[data-spec-action]');if(!button)return;
  const name=button.closest('[data-spectrum-target]')?.dataset.spectrumTarget||state.selected;
  const t=D.targets.find(t=>t.name===name);if(!t)return;
  const n=t.spec?.epochs.length||0,action=button.dataset.specAction;
  spectralSelection.set(t.name,new Set(action==='all'?Array.from({length:n},(_,i)=>i):action==='ends'&&n?[0,n-1]:[]));
  redrawSpectra(t.name);
});
