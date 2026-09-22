// SDSS images are north up / east left. SVG coordinates below are arcseconds.
// A positive sky PA therefore requires a NEGATIVE SVG rotation.
const fieldState={target:null,night:null,size:'30',sample:0,showSlit:true};
function fieldNight(t){return state.night==='all'?nightsOf(t)[0]:state.night;}
function fieldSample(t,k){const p=t.slit_preview?.[k];return p?.samples[p.default_index];}
function fieldCaption(t,k,sample=fieldSample(t,k)){
  if(!sample)return 'No observable PA preview available';
  const p=t.slit_preview[k],label=sample[3]||(p.scheduled?'other time':'unscheduled preview');
  return `${D.nights[k].label} · ${pdt(mjdMs(sample[0]))} PDT · ${label} · PA ${sample[1]==null?'undefined at zenith':fmt(sample[1],1)+'° E of N'} · X ${fmt(sample[2])}`;
}
function fieldSvg(t,size,k,sample=fieldSample(t,k),showSlit=true){
  const wide=size==='wide'||!t.cut40,src=wide?t.cut?.image:t.cut40;
  const sourceSize=wide?Number(t.cut?.field_arcsec):40;
  if(!src||!Number.isFinite(sourceSize)||sourceSize<=0)return '<div class="empty">No scaled cutout available.</div>';
  const fov=wide?sourceSize:Number(size),h=fov/2,pa=sample?.[1];
  const caption=`${t.name}, ${fov} arcsecond SDSS field. North up, east left. ${fieldCaption(t,k,sample)}`;
  let svg=`<svg class="cutout field-svg" viewBox="${-h} ${-h} ${fov} ${fov}" role="img" aria-label="${esc(caption)}" style="overflow:hidden"><image href="${src}" x="${-sourceSize/2}" y="${-sourceSize/2}" width="${sourceSize}" height="${sourceSize}"/>`;
  if(showSlit&&Number.isFinite(pa))svg+=`<rect class="slit-footprint" x="-0.75" y="-30" width="1.5" height="60" transform="rotate(${-pa})" fill="none" stroke="#ffd480" stroke-width="1" vector-effect="non-scaling-stroke"><title>1.5″ × 60″ slit, centered on catalogue coordinates; PA ${fmt(pa,1)}° east of north</title></rect>`;
  // Compass and scale have their own pixel-like frame so they stay legible.
  const bar=fov>40?30:5,barPixels=bar/fov*240;
  svg+=`<svg x="${-h}" y="${-h}" width="${fov}" height="${fov}" viewBox="0 0 240 240"><g fill="none" stroke="#fff" stroke-width="1.2"><path d="M36 40V16m0 24H12M36 16l-3 5m3-5 3 5M12 40l5-3m-5 3 5 3"/><path d="M14 222h${barPixels}m-${barPixels} -3v6m${barPixels} -6v6"/></g><g style="paint-order:stroke;stroke:#000;stroke-width:2px;stroke-linejoin:round"><text x="32" y="12" style="fill:white;font-size:10px">N</text><text x="4" y="44" style="fill:white;font-size:10px">E</text><text x="14" y="216" style="fill:white;font-size:11px">${bar}″</text></g></svg></svg>`;
  return svg;
}
function fieldCard(t){
  const k=fieldNight(t),size=t.cut40?'30':'wide',fov=t.cut40?30:t.cut?.field_arcsec;
  return `<h3>SDSS · ${fov??'—'}″ × ${fov??'—'}″</h3><div><button class="field-button" data-field="${esc(t.name)}" data-size="${size}" aria-label="Enlarge field and inspect slit for ${esc(t.name)}">${fieldSvg(t,size,k)}</button><p class="plot-note"><span style="color:#ffd480">1.5″ slit</span> · North up · east left<br>${esc(fieldCaption(t,k))}</p><button class="small" data-field="${esc(t.name)}" data-size="${size}">Inspect slit / change time</button>${t.cut?.image?` <button class="small" data-field="${esc(t.name)}" data-size="wide">Wider field (${t.cut.field_arcsec}″)</button>`:''}</div>`;
}
function renderField(){
  const t=fieldState.target,k=fieldState.night,p=t.slit_preview?.[k],sample=p?.samples[fieldState.sample];
  const fov=fieldState.size==='wide'||!t.cut40?t.cut?.field_arcsec:fieldState.size;
  $('field-title').textContent=`${t.name} · ${fov}″ × ${fov}″`;
  $('field-image').innerHTML=fieldSvg(t,fieldState.size,k,sample,fieldState.showSlit);
  $('field-caption').textContent=fieldCaption(t,k,sample);
  $('field-size').value=fieldState.size;
  $('field-slit').checked=fieldState.showSlit;
  $('field-note').textContent=`The outline shows a 1.5″ × 60″ slit centered on the catalogue coordinates${Number(fov)<60?'; its ends lie outside this view':''}. PA is a planning estimate at the displayed time, modulo 180°. Refresh the NGPS prediction for the actual exposure. Seeing spreads light across the slit edges; use the acquisition image to check centering and companions.`;
}
function fieldTimes(){
  const t=fieldState.target,p=t.slit_preview?.[fieldState.night];
  fieldState.sample=p?.default_index??0;
  $('field-time').innerHTML=p?p.samples.map((s,i)=>`<option value="${i}">${pdt(mjdMs(s[0]))}${s[3]?' · '+esc(s[3]):''} · X ${fmt(s[2])}</option>`).join(''):'<option>No observable preview</option>';
  $('field-time').value=String(fieldState.sample);$('field-time').disabled=!p;
}
function openField(t,size){
  fieldState.target=t;fieldState.night=fieldNight(t);fieldState.size=t.cut40?size:'wide';
  $('field-night').innerHTML=nightsOf(t).map(k=>`<option value="${k}">${esc(D.nights[k].label)}</option>`).join('');
  $('field-night').value=fieldState.night;
  $('field-size').innerHTML=(t.cut40?'<option value="30">30″ × 30″</option><option value="40">40″ × 40″</option>':'')+(t.cut?.image?`<option value="wide">${t.cut.field_arcsec}″ · wide</option>`:'');
  fieldTimes();renderField();$('field-dialog').showModal();
}
