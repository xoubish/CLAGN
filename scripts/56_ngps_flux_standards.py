"""Build explicitly central-slice sensitivity functions for the same-night standards."""
import importlib
import sys
from concurrent.futures import ThreadPoolExecutor
import numpy as np
env = importlib.import_module('54_reduce_ngps_sep23')
from pypeit import specobjs

def channel(ch):
    sci=env.NIGHT/'central_exposures'/ch
    assert sci.exists(), 'Prepare corrected central spectra first.'
    out=env.BASE/'audit/standards'/ch
    out.mkdir(parents=True,exist_ok=True)
    for p in sorted(sci.glob('spec1d*fits')):
        if not any(s in p.name for s in ['P330E','BD284211']):
            continue
        objs=specobjs.SpecObjs.from_fitsfile(str(p))
        slit={'u':550,'g':240,'r':243,'i':236}[ch]
        ix=[i for i,o in enumerate(objs) if o.SLITID==slit]
        assert len(ix)==1, (p,ix)
        central=specobjs.SpecObjs(specobjs=np.array([objs[ix[0]]]),header=objs.header)
        cp=out/p.name
        central.write_to_fits(objs.header,str(cp),overwrite=True)
        sens=out/('sens_'+p.name)
        if not sens.exists():
            env.run([env.BIN/'pypeit_sensfunc',cp,'--algorithm','UVIS','--extr','OPT',
                     '-o',sens],f'sens_{ch}_{p.name}.log',cwd=out)
        print('Sensitivity ready:',sens.name,flush=True)

if __name__=='__main__':
    with ThreadPoolExecutor(max_workers=4) as pool:
        list(pool.map(channel,sys.argv[1:] or list('ugri')))
