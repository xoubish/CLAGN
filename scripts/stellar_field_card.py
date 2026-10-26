"""Attach the scheduled 16 Cyg field with a WCS-projected slit overlay."""
import base64
import json
import warnings
from pathlib import Path

import astropy.units as u
import numpy as np
from astropy.coordinates import SkyCoord
from astropy.io import fits
from astropy.wcs import WCS, FITSFixedWarning

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / 'data/additional_targets/field_2954592_505204'


def attach_stellar_field_cards(payload):
    visits = payload.get('october', {}).get('manual_visits', [])
    payload['field_cards'] = []
    if not visits or not (DATA/'ps1_card.fits').exists():return
    stars = json.loads((DATA/'stars_2026.json').read_text())
    with warnings.catch_warnings():
        warnings.simplefilter('ignore', FITSFixedWarning)
        wcs = WCS(fits.getheader(DATA/'ps1_card.fits'))
    height, width = wcs.array_shape
    # FITS pixel centers are zero-based; JPEG rows run in the opposite direction.
    def pixel(coord):
        x,y = wcs.world_to_pixel(coord)
        return float(x+.5), float(height-y-.5)
    def points(coords):return ' '.join(f'{x:.3f},{y:.3f}' for x,y in map(pixel,coords))
    image = 'data:image/jpeg;base64,' + base64.b64encode((DATA/'ps1_card.jpg').read_bytes()).decode()
    for visit in visits:
        if visit.get('upload_name') != 'Z16CygABC':continue
        center = SkyCoord(visit['ra']*u.deg,visit['dec']*u.deg)
        angle = np.deg2rad(visit['pa_deg'])
        along = np.array([np.sin(angle),np.cos(angle)])
        across = np.array([np.cos(angle),-np.sin(angle)])
        polygons = []
        for k in [-1,0,1]:
            offsets = [a*along+(k*1.5+b)*across for a,b in [(-25,-.75),(25,-.75),(25,.75),(-25,.75)]]
            coords = [center.spherical_offsets_by(e*u.arcsec,n*u.arcsec) for e,n in offsets]
            polygons.append(f'<polygon class="stellar-slice" points="{points(coords)}" fill="#ffd480" fill-opacity=".08" stroke="#ffd480" stroke-width="1.1"/>')
        markers = []
        for name in 'ABC':
            star = stars[name]
            x,y = pixel(SkyCoord(star['ra']*u.deg,star['dec']*u.deg))
            dx,dy = {'A':(-29,-30),'B':(-30,30),'C':(42,23)}[name]
            markers.append(f'<g class="stellar-component" data-component="{name}"><circle cx="{x:.3f}" cy="{y:.3f}" r="3" fill="none" stroke="white" stroke-width="1"/><path d="M{x+dx} {y+dy}L{x} {y}" stroke="white" stroke-width=".8"/><text x="{x+dx}" y="{y+dy-4}" fill="white" font-size="14" font-family="sans-serif" font-weight="bold" style="paint-order:stroke;stroke:#000;stroke-width:3px">{name}</text></g>')
        cx,cy=pixel(center)
        # Draw the compass with the same WCS as the image, preserving its slight rotation.
        north=pixel(center.spherical_offsets_by(0*u.arcsec,6*u.arcsec))
        east=pixel(center.spherical_offsets_by(6*u.arcsec,0*u.arcsec))
        compass=''
        for label,(px,py) in [('N',north),('E',east)]:
            ex,ey=40+px-cx,45+py-cy
            compass+=f'<path d="M40 45L{ex:.2f} {ey:.2f}" stroke="white" stroke-width="1.2"/><text x="{ex-4:.2f}" y="{ey-5:.2f}" fill="white" font-size="11">{label}</text>'
        svg=f'''<svg xmlns="http://www.w3.org/2000/svg" class="cutout stellar-field" viewBox="0 0 {width} {height}" role="img" aria-label="16 Cyg A B and C: Pan-STARRS image with three 1.5 arcsecond by 50 arcsecond NGPS slices at PA 133.04 degrees. Markers predict October 26 positions; C is hidden by saturated A in the archive image.">
<title>16 Cyg A/B/C — Pan-STARRS with NGPS slit overlay</title>
<image href="{image}" width="{width}" height="{height}"/>
<g class="stellar-slit">{''.join(polygons)}</g>
{''.join(markers)}
<path d="M{cx-4} {cy}h8M{cx} {cy-4}v8" stroke="#ffd480" stroke-width="1"/>
<g font-family="sans-serif" style="paint-order:stroke;stroke:#000;stroke-width:2px">{compass}<path d="M16 340h40" stroke="white" stroke-width="2"/><text x="16" y="332" fill="white" font-size="11">10″</text></g>
</svg>'''
        dest=ROOT/'docs/october/16Cyg_field.svg';dest.write_text(svg)
        payload['field_cards'].append(dict(visit,name=visit['upload_name'],label=visit['name'],
            field_svg=svg,image_source='Pan-STARRS gri',field_arcsec=90,
            image_status='SDSS: outside imaging footprint. A/B are saturated; C is not resolved in this image.',
            components=[dict(stars[n],name=n) for n in 'ABC'],
            archive_review=json.loads((DATA/'spectra/archive_review.json').read_text())
                if (DATA/'spectra/archive_review.json').exists() else None))
