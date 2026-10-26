/* Queue-ordered, full-visit geometry; times and ranges are precomputed locally. */
function displayName(name){return D.display_names?.[name]||name;}
const VIS_COLORS=['#677584','#f3cc5c','#5cace2','#aa8bd4'];
const VIS_LABELS=['Blocked','Best','Usable','Higher X'];
function visibilityStatus(V,r,index){
  const sample=r.samples[index];
  if(!sample)return {tier:0,text:'Too late to finish'};
  const [x,moon,tier]=sample;
  const start=V.sample_start+index*V.step;
  let reason=VIS_LABELS[tier];
  if(!tier){
    if(V.standards.some(s=>!(r.manual&&s.name===r.name)&&start<s.end&&start+r.duration*60000>s.start))reason='Reserved block overlap';
    else if(moon<40.02)reason='Moon too close';
    else reason=x===null?'Below horizon / too low':'Above X limit';
  }
  return {tier,text:`${reason} · X ${x===null?'—':x.toFixed(2)} · ☾ ${moon.toFixed(1)}°`};
}
function visibilitySVG(night,index){
  const V=D.visibility[night],W=1280,L=225,R=260,T=34,B=30,H=T+B+V.rows.length*25;
  const X=lin(V.start,V.end,L,W-R),time=V.sample_start+index*V.step;
  let s=`<svg viewBox="0 0 ${W} ${H}" width="100%" data-visibility-chart="${night}" role="img" aria-label="Full-visit start windows in observing order; all times PDT"><text x="8" y="18">ORDER · TARGET</text><text x="${W-R+14}" y="18">IF STARTED AT ${pdt(time)} PDT</text>`;
  V.standards.forEach(b=>s+=`<rect x="${X(b.start)}" y="${T-8}" width="${X(b.end)-X(b.start)}" height="${H-T-B+12}" fill="#80909b" opacity=".16"><title>${esc(b.name)} ${b.kind||'standard'}: ${pdt(b.start)}–${pdt(b.end)} PDT</title></rect>`);
  for(let t=Math.ceil(V.start/3600000)*3600000;t<=V.end;t+=3600000){
    s+=`<line x1="${X(t)}" x2="${X(t)}" y1="${T-8}" y2="${H-B}" stroke="var(--hair2)"/><text x="${X(t)}" y="18" text-anchor="middle">${pdt(t)}</text><text x="${X(t)}" y="${H-7}" text-anchor="middle">${pdt(t)}</text>`;
  }
  V.rows.forEach((r,i)=>{
    const y=T+i*25,status=visibilityStatus(V,r,index),label=displayName(r.name);
    const targetLink=r.manual?`href="#${night}" data-manual-target="${esc(r.name)}"`:`href="#card-${esc(r.name)}"`;
    s+=`<g data-visibility-target="${esc(r.name)}"><line x1="0" x2="${W}" y1="${y+21}" y2="${y+21}" stroke="var(--hair)"/><a ${targetLink}><text x="8" y="${y+13}"><title>${esc(label)}${r.manual?' · manual acquisition · ≤10 s per exposure':''}${label!==r.name?' · alias '+esc(r.name):''} · ${r.duration} min · X ≤${r.airmass_limit}</title>${String(r.order).padStart(2,'0')} · ${esc(label.length>25?label.slice(0,24)+'…':label)}</text></a>`;
    r.windows.forEach(w=>s+=`<rect data-visibility-tier="${w.tier}" x="${X(w.start)}" y="${y+3}" width="${Math.max(1,X(w.end)-X(w.start))}" height="12" rx="1" fill="${VIS_COLORS[w.tier]}"><title>${esc(label)} · ${VIS_LABELS[w.tier]} starts ${pdt(w.start)}–${pdt(w.end)} PDT; last finish ${pdt(w.end+r.duration*60000)}. ${r.duration} min including overhead; Moon ≥40° throughout.</title></rect>`);
    s+=`<rect x="${X(r.planned_start)}" y="${y}" width="${X(r.planned_end)-X(r.planned_start)}" height="18" fill="none" stroke="var(--ink)" stroke-width="1.5"><title>Planned ${esc(label)} visit ${pdt(r.planned_start)}–${pdt(r.planned_end)} PDT</title></rect><text data-visibility-status="${status.tier}" x="${W-R+14}" y="${y+13}" style="fill:${status.tier?VIS_COLORS[status.tier]:'var(--ink3)'}">${esc(status.text)}</text></g>`;
  });
  return s+`<line x1="${X(time)}" x2="${X(time)}" y1="${T-10}" y2="${H-B}" stroke="var(--ink)" stroke-dasharray="3 3" pointer-events="none"/></svg>`;
}
function visibilityReadout(night,index){
  const V=D.visibility[night],n=V.rows.filter(r=>r.samples[index]?.[2]>0).length;
  return `${pdt(V.sample_start+index*V.step)} PDT · ${n} of ${V.rows.length} targets fit a complete visit`;
}
function renderVisibility(night){
  const V=D.visibility?.[night];if(!V)return '';
  const index=Math.max(0,Math.round((V.rows[0].planned_start-V.sample_start)/V.step));
  return `<section class="visibility"><h3>Visibility windows · observing order</h3><p class="plot-note">Colored bars show when a <b>complete visit can start</b>, including overhead, with Moon separation ≥40° throughout. Yellow: best, X ≤1.5. Blue: usable, X ≤1.8. Purple: X ≤2.0, only where the target’s limit permits. Blank: blocked. Outlined boxes show planned visits; shaded columns reserve standards and manual visits.</p><div class="visibility-controls"><label for="visibility-time-${night}">Try a start time (PDT)</label><input id="visibility-time-${night}" data-visibility-time="${night}" type="range" min="0" max="${V.count-1}" step="1" value="${index}" aria-valuetext="${pdt(V.sample_start+index*V.step)} PDT"><output id="visibility-readout-${night}" for="visibility-time-${night}">${visibilityReadout(night,index)}</output></div><div class="table-scroll visibility-plot" id="visibility-plot-${night}">${visibilitySVG(night,index)}</div><p class="plot-note">Drag the time slider or click the chart’s time axis to compare swaps. Right-hand values are maximum airmass and minimum Moon separation over the whole visit. Other science visits are not reserved here: check queue conflicts and avoid repeats when swapping. Palomar geometry on a 30-second grid, rounded inward at limits; confirm current conditions and sequencer limits at the telescope.</p></section>`;
}
function updateVisibility(night,index){
  const slider=$('visibility-time-'+night);
  index=Math.max(0,Math.min(D.visibility[night].count-1,Math.round(index)));
  slider.value=index;slider.setAttribute('aria-valuetext',pdt(D.visibility[night].sample_start+index*D.visibility[night].step)+' PDT');
  $('visibility-readout-'+night).textContent=visibilityReadout(night,index);
  $('visibility-plot-'+night).innerHTML=visibilitySVG(night,index);
}
document.addEventListener('input',e=>{if(e.target.dataset.visibilityTime)updateVisibility(e.target.dataset.visibilityTime,+e.target.value);});
document.addEventListener('click',e=>{
  const svg=e.target.closest('[data-visibility-chart]');if(!svg||e.target.closest('a'))return;
  const point=new DOMPoint(e.clientX,e.clientY).matrixTransform(svg.getScreenCTM().inverse());
  if(point.x<225||point.x>1020)return;
  const night=svg.dataset.visibilityChart,V=D.visibility[night];
  updateVisibility(night,((point.x-225)/795*(V.end-V.start)+V.start-V.sample_start)/V.step);
});
