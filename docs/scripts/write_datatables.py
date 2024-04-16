#!/usr/bin/env python3

"""
Dynamically build the rst documentation of the bitmasks.
"""

import time
from importlib import resources

#-----------------------------------------------------------------------------

def write_datatable(datatable_class, opath, class_link=True):
    ofile = opath / f'{datatable_class.__name__.lower()}.rst'
    lines = datatable_class().to_rst_table(header=False, class_link=class_link)
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
    from mangadap.proc.reductionassessments import ReductionAssessmentDataTable
    write_datatable(ReductionAssessmentDataTable, path, class_link=False)
    from mangadap.proc.spatiallybinnedspectra import SpatiallyBinnedSpectraDataTable
    write_datatable(SpatiallyBinnedSpectraDataTable, path, class_link=False)
    from mangadap.proc.spectralfitting import StellarKinematicsFitDataTable
    write_datatable(StellarKinematicsFitDataTable, path, class_link=False)
    from mangadap.par.emissionmomentsdb import EmissionMomentsDefinitionTable
    write_datatable(EmissionMomentsDefinitionTable, path, class_link=False)
    from mangadap.proc.emissionlinemoments import EmissionLineMomentsDataTable
    write_datatable(EmissionLineMomentsDataTable, path, class_link=False)
    from mangadap.par.emissionlinedb import EmissionLineDefinitionTable
    write_datatable(EmissionLineDefinitionTable, path, class_link=False)
    from mangadap.proc.spectralfitting import EmissionLineFitDataTable
    write_datatable(EmissionLineFitDataTable, path, class_link=False)
    from mangadap.proc.spectralindices import SpectralIndicesDefinitionTable
    write_datatable(SpectralIndicesDefinitionTable, path, class_link=False)
    from mangadap.proc.spectralindices import SpectralIndicesDataTable
    write_datatable(SpectralIndicesDataTable, path, class_link=False)

    print('Elapsed time: {0} seconds'.format(time.perf_counter() - t))



