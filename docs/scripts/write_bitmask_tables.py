#!/usr/bin/env python3

"""
Dynamically build the rst documentation of the bitmasks.
"""

import time
from importlib import resources

#-----------------------------------------------------------------------------

def write_bitmask(bitmask_class, opath, class_link=True):
    ofile = opath / f'{bitmask_class.__name__.lower()}.rst'
    lines = bitmask_class().to_rst_table(header=False, class_link=class_link)
    with open(ofile, 'w') as f:
        f.write('\n'.join(lines))

if __name__ == '__main__':
    t = time.perf_counter()

    root = resources.files('mangadap').parent
    path = root / 'docs' / 'tables'
    if not path.is_dir():
        path.mkdir(parents=True)

    # Tables to write:
    #
    from mangadap.util.drpbitmask import DRPQuality3DBitMask
    write_bitmask(DRPQuality3DBitMask, path)
    from mangadap.util.drpbitmask import DRPFitsBitMask
    write_bitmask(DRPFitsBitMask, path)

    from mangadap.dapfits import DAPQualityBitMask
    write_bitmask(DAPQualityBitMask, path)
    from mangadap.dapfits import DAPMapsBitMask
    write_bitmask(DAPMapsBitMask, path)
    from mangadap.dapfits import DAPCubeBitMask
    write_bitmask(DAPCubeBitMask, path)

    from mangadap.proc.spatiallybinnedspectra import SpatiallyBinnedSpectraBitMask
    write_bitmask(SpatiallyBinnedSpectraBitMask, path, class_link=False)
    from mangadap.proc.templatelibrary import TemplateLibraryBitMask
    write_bitmask(TemplateLibraryBitMask, path, class_link=False)
    from mangadap.proc.stellarcontinuummodel import StellarContinuumModelBitMask
    write_bitmask(StellarContinuumModelBitMask, path, class_link=False)
    from mangadap.proc.emissionlinemoments import EmissionLineMomentsBitMask
    write_bitmask(EmissionLineMomentsBitMask, path, class_link=False)
    from mangadap.proc.emissionlinemodel import EmissionLineModelBitMask
    write_bitmask(EmissionLineModelBitMask, path, class_link=False)
    from mangadap.proc.spectralindices import SpectralIndicesBitMask
    write_bitmask(SpectralIndicesBitMask, path, class_link=False)

    print('Elapsed time: {0} seconds'.format(time.perf_counter() - t))



