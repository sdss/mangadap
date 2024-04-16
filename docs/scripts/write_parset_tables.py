#!/usr/bin/env python3

"""
Dynamically build the rst documentation of the bitmasks.
"""

import time
from importlib import resources

#-----------------------------------------------------------------------------

def write_parset(parset_class, opath, class_link=True):
    ofile = opath / f'{parset_class.__name__.lower()}.rst'
    lines = parset_class().to_rst_table(header=False, class_link=class_link, nested=False)
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
#    from mangadap.config.analysisplan import AnalysisPlan
#    write_parset(AnalysisPlan, path, class_link=False)

    from mangadap.par.artifactdb import ArtifactPar
    write_parset(ArtifactPar, path, class_link=False)
    from mangadap.par.emissionlinedb import EmissionLinePar
    write_parset(EmissionLinePar, path, class_link=False)

    from mangadap.proc.reductionassessments import ReductionAssessmentDef
    write_parset(ReductionAssessmentDef, path, class_link=False)

    from mangadap.proc.spatiallybinnedspectra import SpatiallyBinnedSpectraDef
    write_parset(SpatiallyBinnedSpectraDef, path, class_link=False)
    from mangadap.proc.spatialbinning import RadialBinningPar, VoronoiBinningPar, SquareBinningPar
    write_parset(RadialBinningPar, path, class_link=False)
    write_parset(VoronoiBinningPar, path, class_link=False)
    write_parset(SquareBinningPar, path, class_link=False)
    from mangadap.proc.spectralstack import SpectralStackPar
    write_parset(SpectralStackPar, path, class_link=False)

    from mangadap.proc.templatelibrary import TemplateLibraryDef
    write_parset(TemplateLibraryDef, path, class_link=False)
    from mangadap.proc.stellarcontinuummodel import StellarContinuumModelDef
    write_parset(StellarContinuumModelDef, path, class_link=False)
    from mangadap.proc.ppxffit import PPXFFitPar
    write_parset(PPXFFitPar, path, class_link=False)

    from mangadap.proc.emissionlinemoments import EmissionLineMomentsDef
    write_parset(EmissionLineMomentsDef, path, class_link=False)
    from mangadap.proc.bandpassfilter import BandPassFilterPar
    write_parset(BandPassFilterPar, path, class_link=False)

    from mangadap.proc.emissionlinemodel import EmissionLineModelDef
    write_parset(EmissionLineModelDef, path, class_link=False)
#    from mangadap.proc.elric import ElricPar
#    write_parset(ElricPar, path, class_link=False)
    from mangadap.proc.sasuke import SasukePar
    write_parset(SasukePar, path, class_link=False)

    from mangadap.proc.spectralindices import SpectralIndicesDef
    write_parset(SpectralIndicesDef, path, class_link=False)

    print('Elapsed time: {0} seconds'.format(time.perf_counter() - t))



