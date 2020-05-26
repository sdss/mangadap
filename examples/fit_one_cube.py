#!/usr/bin/env python3

import os
import numpy

import astropy.constants
from astropy.io import fits

from mangadap.datacube import MaNGADataCube
from mangadap.survey.manga_dap import manga_dap
from mangadap.par.analysisplan import AnalysisPlan, AnalysisPlanSet

#-----------------------------------------------------------------------------
def get_config(plt, ifu, config_file, drpall_file=None):
    if drpall_file is None:
        drpall_file = os.path.join(os.environ['MANGA_SPECTRO_REDUX'], os.environ['MANGADRP_VER'],
                                   'drpall-{0}.fits'.format(os.environ['MANGADRP_VER']))

    # Use the DRPall file
    with fits.open(drpall_file) as hdu:
        indx = numpy.where(hdu['MANGA'].data['PLATEIFU'] == '{0}-{1}'.format(plt, ifu))[0]
        if len(indx) != 1:
            raise ValueError('{0}-{1} either does not exist or has more than one match!'.format(
                             plt, ifu))

        MaNGADataCube.write_config(config_file, plt, ifu, z=hdu[1].data['z'][indx[0]],
                                   ell=1-hdu[1].data['nsa_elpetro_ba'][indx[0]],
                                   pa=hdu[1].data['nsa_elpetro_phi'][indx[0]],
                                   reff=hdu[1].data['nsa_elpetro_th50_r'][indx[0]],
                                   overwrite=True)


def fit_one_cube(plt, ifu, drpall_file=None, directory_path=None, analysis_path=None):
    # Grab the required input parameters
    config_file = '{0}-{1}.cfg'.format(plt, ifu)
    get_config(plt, ifu, config_file, drpall_file=drpall_file)

    # Read the datacube
    cube = MaNGADataCube.from_config(config_file, directory_path=directory_path)

    # Define how you want to analyze the data
    plan = AnalysisPlanSet([ AnalysisPlan(drpqa_key='SNRG',
                                          bin_key='VOR10', #'HYB10',
                                          continuum_key='MILESHCMPL10',
                                          elmom_key='EMOMMPL10',
                                          elfit_key='EFITMPL10', #'EFITMPL9DB',
                                          spindex_key='INDXEN') ])

    # Run it!
    return manga_dap(cube, plan, verbose=2, directory_path=directory_path,
                     analysis_path=analysis_path)
#-----------------------------------------------------------------------------


if __name__ == '__main__':
    drpver = 'v3_0_1'
    directory_path = os.path.join(os.environ['MANGADAP_DIR'], 'mangadap', 'data', 'remote')
    drpall_file = os.path.join(directory_path, 'drpall-{0}.fits'.format(drpver))
    fit_one_cube(7815, 3702, drpall_file=drpall_file, directory_path=directory_path,
                 analysis_path='./output')

