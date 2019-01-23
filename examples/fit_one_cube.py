#!/usr/bin/env python3

import time
import os
import time
import astropy.constants
from astropy.io import fits
from mangadap.survey.manga_dap import manga_dap
from mangadap.par.obsinput import ObsInputPar
from mangadap.par.analysisplan import AnalysisPlan, AnalysisPlanSet

#-----------------------------------------------------------------------------
def get_obsinput(plt, ifu, drpall_file=None):
    """
    Grab the input parameters the DAP requires for each observation to
    fit a cube.  If the drpall file is None, use the default path.
    """
    hdu = fits.open(os.path.join(os.environ['MANGA_SPECTRO_REDUX'], os.environ['MANGADRP_VER'],
                                 'drpall-{0}.fits'.format(os.environ['MANGADRP_VER']))) \
                if drpall_file is None else fits.open(drpall_file)
    indx = hdu[1].data['PLATEIFU'] == '{0}-{1}'.format(plt, ifu)

    return ObsInputPar(plt, ifu, mode='CUBE',
                       vel=astropy.constants.c.to('km/s').value*hdu[1].data['NSA_Z'][indx][0],
                       vdisp=100.,
                       ell=1-hdu[1].data['NSA_ELPETRO_BA'][indx][0],
                       pa=hdu[1].data['NSA_ELPETRO_PHI'][indx][0],
                       reff=hdu[1].data['NSA_ELPETRO_TH50_R'][indx][0])


def fit_one_cube(plt, ifu, drpall_file=None, directory_path=None, analysis_path=None):
    # Grab the required input parameters
    obs = get_obsinput(plt, ifu, drpall_file=drpall_file)

    # Define how you want to analyze the data
    plan = AnalysisPlanSet([ AnalysisPlan(drpqa_key='SNRG',
                                          bin_key='HYB10',
                                          continuum_key='GAU-MILESHC',
                                          elmom_key='EMOMM',
                                          elfit_key='EFITMB',
                                          spindex_key='INDXEN') ])

#    plan = AnalysisPlanSet([ AnalysisPlan(drpqa_key='SNRG',
#                                          bin_key='VOR10',
#                                          continuum_key='GAU-MILESHC',
#                                          elmom_key='EMOMM',
#                                          elfit_key='EFITM',
#                                          spindex_key='INDXEN') ])

    # Run it!
    return manga_dap(obs, plan, verbose=2, directory_path=directory_path,
                     analysis_path=analysis_path)


#-----------------------------------------------------------------------------

if __name__ == '__main__':
    t = time.clock()

#    fit_one_cube(7815, 1902, drpall_file='./data/drpall-v2_5_3.fits', directory_path='./data',
#                 analysis_path='./output')

    fit_one_cube(7815, 3702, analysis_path='./output')

    print('Elapsed time: {0} seconds'.format(time.clock() - t))

