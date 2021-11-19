#!/usr/bin/env python3

import os
import numpy

import astropy.constants
from astropy.io import fits

from mangadap.datacube import MaNGADataCube
from mangadap.datacube import MUSEDataCube
from mangadap.survey.manga_dap import manga_dap
from mangadap.par.analysisplan import AnalysisPlan, AnalysisPlanSet
from mangadap.proc.reductionassessments import available_reduction_assessments
from mangadap.scripts.spotcheck_dap_maps import spotcheck_images
from mangadap.scripts.ppxffit_qa import ppxffit_qa_plot
from mangadap.scripts.fit_residuals_muse import fit_residuals_muse
from mangadap.scripts.manga_dap_inspector import manga_dap_inspector
from IPython import embed

#-----------------------------------------------------------------------------


def test_read_muse(cubefil, sresfil, directory_path=None, analysis_path=None,
                      redshift=None, vdisp=None, ellipticity=None, pa=None, reff=None):

    ifile = directory_path + cubefil

    # Read the datacube
    cube = MUSEDataCube(ifile, z=redshift, vdisp=vdisp, ell=ellipticity, pa=pa,
                        reff=reff, sres_ifile=sresfil, sres_fill=True, covar_ext=None)


    #assert cube.file_name == MaNGADataCube.build_file_name(cube.plate, cube.ifudesign,
    #                log=cube.log), 'Name mismatch'
    assert cube.log, 'Should read the log-binned version by default.'
    assert cube.wcs is not None, 'WCS should be defined.'
    assert cube.shape[:2] == cube.spatial_shape, 'Spatial shape should be first two axes.'
    assert cube.nspec == numpy.prod(cube.spatial_shape), 'Definition of number of spectra changed.'
    assert cube.sres is not None, 'Spectral resolution data was not constructed.'
    # Is it ok to ignore this? -- yes, but probably we're using pixelized LSF
    assert cube.sres_ext == 'LSFPRE', 'Should default to LSFPRE extension.'
    assert abs(cube.pixelscale - cube._get_pixelscale()) < 1e-6, 'Bad match in pixel scale.'
    # NOTE: This is worse than it should be because of how the WCS in MaNGA is defined.

    # MUSE is failing this!!!!
    assert numpy.all(numpy.absolute(cube.wave - cube._get_wavelength_vector(cube.nwave)) < 2e-4), \
            'Bad calculation of wavelength vector.'
    assert cube.covar is None, 'Covariance should not have been read'



def test_wcs_muse(cubefil, sresfil, directory_path=None, analysis_path=None,
                      redshift=None, objra=None, objdec=None, vdisp=None, ellipticity=None, pa=None, reff=None):
    ifile = directory_path + cubefil

    # Read the datacube
    cube = MUSEDataCube(ifile, objra=objra, objdec=objdec, z=redshift, vdisp=vdisp, ell=ellipticity, pa=pa,
                        reff=reff, sres_ifile=sresfil, sres_fill=True, covar_ext=None)

    x, y = cube.mean_sky_coordinates(offset=None)
    assert x[0,0] > x[-1,0], 'RA should increase from large to small indices'
    assert y[0,0] < y[0,-1], 'DEC should increase from small to small indices'

    print(numpy.unravel_index(numpy.argmin( numpy.square(x - cube.prihdr['OBJRA'])
                                            + numpy.square(y - cube.prihdr['OBJDEC'])), x.shape))
    assert numpy.unravel_index(numpy.argmin( numpy.square(x - cube.prihdr['OBJRA'])
                                            + numpy.square(y - cube.prihdr['OBJDEC'])), x.shape) \
                == (14,14), 'Object should be at cube center.'
    x, y = cube.mean_sky_coordinates(center_coo=(x[0,0], y[0,0]))
    assert numpy.isclose(x[0,0], 0.0) and numpy.isclose(y[0,0], 0.0), 'Offset incorrect'
    x, y = cube.mean_sky_coordinates()

    ## The below breaks with my by-eye estimated object centroid
    print(abs(x[14,14]) < 1e-2 and abs(y[14,14]))
    embed()
    assert abs(x[14,14]) < 1e-2 and abs(y[14,14]) < 1e-2, 'Offset incorrect'




def test_copyto_muse(cubefil, sresfil, directory_path=None, analysis_path=None,
                      redshift=None, objra=None, objdec=None, vdisp=None, ellipticity=None, pa=None, reff=None):
    ifile = directory_path + cubefil

    # Read the datacube
    cube = MUSEDataCube(ifile, objra=objra, objdec=objdec, z=redshift, vdisp=vdisp, ell=ellipticity, pa=pa,
                        reff=reff, sres_ifile=sresfil, sres_fill=True, covar_ext=None)

    flux = cube.copy_to_array()
    assert not isinstance(flux, numpy.ma.MaskedArray), 'Should output normal array'
    assert flux.shape[0] == cube.nspec, 'Should be flattened into a 2D array.'
    assert flux.shape[1] == cube.nwave, 'Should be flattened into a 2D array.'

    # Apply a wavelength mask
    waverange = [5000, 7000]
    flux = cube.copy_to_array(waverange=waverange)
    indx = (cube.wave > waverange[0]) & (cube.wave < waverange[1])
    assert flux.shape[1] == numpy.sum(indx), 'Wavelength range masking failed'

    # Find the spaxels with non-zero signal
    methods = available_reduction_assessments()
    i = numpy.where([m['key'] == 'SNRG' for m in methods])[0]
    assert len(i) == 1, 'Could not find correct reduction assessment definition.'

    ## Issue here with datacube.py, in copy_to_masked_array
    sig, var, snr = cube.flux_stats(response_func=methods[i[0]]['response_func'])
    indx = ((sig > 0) & numpy.invert(numpy.ma.getmaskarray(sig))).data.ravel()
    ngood = numpy.sum(indx)

    # Select the spaxels with non-zero signal
    flux = cube.copy_to_array(waverange=waverange, select_bins=indx)
    assert flux.shape[0] == ngood, 'Bin selection failed'

    # Get the masked array
    flux = cube.copy_to_masked_array()
    assert isinstance(flux, numpy.ma.MaskedArray), 'Should output a masked array'
    assert flux.shape[0] == cube.nspec, 'Should be flattened into a 2D array.'
    assert flux.shape[1] == cube.nwave, 'Should be flattened into a 2D array.'

    # Select the spaxels with non-zero signal
    flux = cube.copy_to_masked_array(select_bins=indx)
    assert flux.shape[0] == ngood, 'Bin selection failed'

    # Try to get the inverse variance
    i = cube.nspec//2 + cube.spatial_shape[1]//2
    ivar = cube.copy_to_masked_array(attr='ivar')
    assert ivar.shape == (cube.nspec, cube.nwave), 'Bad ivar shape'
    assert numpy.array_equal(cube.ivar[numpy.unravel_index(i, cube.spatial_shape)],
                             ivar[i].data), 'Did not pull ivar data.'

    # Try to get the spectral resolution
    sres = cube.copy_to_masked_array(attr='sres')
    assert sres.shape == (cube.nspec, cube.nwave), 'Bad sres shape'
    assert numpy.array_equal(cube.sres[numpy.unravel_index(i, cube.spatial_shape)],
                             sres[i].data), 'Did not pull sres data.'



def test_stats_muse(cubefil, sresfil, directory_path=None, analysis_path=None,
                      redshift=None, objra=None, objdec=None, vdisp=None, ellipticity=None, pa=None, reff=None):
    ifile = directory_path + cubefil

    # Read the datacube
    cube = MUSEDataCube(ifile, objra=objra, objdec=objdec, z=redshift, vdisp=vdisp, ell=ellipticity, pa=pa,
                        reff=reff, sres_ifile=sresfil, sres_fill=True, covar_ext=None)

    # Create a fake bin map
    bin_indx = numpy.arange(cube.nspec/4, dtype=int).reshape(cube.spatial_shape[0]//2,
                                                             cube.spatial_shape[0]//2)
    bin_indx = numpy.repeat(bin_indx, 2, axis=0)
    bin_indx = numpy.repeat(bin_indx, 2, axis=1)

    # Get the bin area
    bins, area = cube.binned_on_sky_area(bin_indx)

    assert numpy.array_equal(bins, numpy.arange(cube.nspec/4)), 'Bad bin list'


    assert numpy.allclose(area, 0.16), 'Bad area calculation'

    methods = available_reduction_assessments()
    i = numpy.where([m['key'] == 'SNRG' for m in methods])[0]
    assert len(i) == 1, 'Could not find correct reduction assessment definition.'

    # These numbers need to be adjusted for individual cubes
    # KHRR -- they kind of make sense for MUSE test cube, but I didn't check them in detail
    cen_wave = cube.central_wavelength(response_func=methods[i[0]]['response_func'],
                                       flag=cube.do_not_use_flags())
    #assert numpy.isclose(cen_wave, 4638.0), 'Central wavelength changed.'

    cen_wave = cube.central_wavelength(waverange=[4000,8000], flag=cube.do_not_use_flags(),
                                       fluxwgt=True)
    #assert numpy.isclose(cen_wave, 5895.7), 'Central wavelength changed.'

    cen_wave = cube.central_wavelength(waverange=[4000,8000], flag=cube.do_not_use_flags(),
                                       per_pixel=False)
    #assert numpy.isclose(cen_wave, 6044.9), 'Central wavelength changed.'

    sig, var, snr = cube.flux_stats(response_func=methods[i[0]]['response_func'])
    assert sig.shape == cube.spatial_shape, 'Should be shaped as a map.'
    assert isinstance(sig, numpy.ma.MaskedArray), 'Expected masked arrays'
    assert numpy.ma.amax(snr) > 60, 'S/N changed'

    # Try it with the linear cube
    #cube = MaNGADataCube.from_plateifu(7815, 3702, directory_path=remote_data_file(), log=False)
    #_sig, _var, _snr = cube.flux_stats(response_func=methods[i[0]]['response_func'])
    # TODO: Not sure why these are not closer.
    #assert numpy.absolute(numpy.ma.median((sig-_sig)/_sig)) < 0.01, \
    #        'Signal should be the same to better than 1%.'
    #assert numpy.absolute(numpy.ma.median((var-_var)/_var)) < 0.03, \
    #        'Variance should be the same to better than 3%.'
    #assert numpy.absolute(numpy.ma.median((snr-_snr)/_snr)) < 0.02, \
    #        'S/N should be the same to better than 2%.'





#def fit_one_cube_muse(cubefil, sresfil, directory_path=None, analysis_path=None,
#                      redshift=None, objra=None, objdec=None, vdisp=None,
#                      ellipticity=None, pa=None, reff=None, plate=None, ifudesign=None, ebvgal=None):
def fit_one_cube_muse(config_file, directory_path=None, analysis_path=None):


    # Grab the required input parameters
    #config_file = '{0}-{1}.cfg'.format(plt, ifu)
    #get_config(plt, ifu, config_file, drpall_file=drpall_file)

    #ifile = directory_path + cubefil

    # Read the datacube
    #cube = MUSEDataCube(ifile, z=redshift, objra=objra, objdec=objdec, vdisp=vdisp, ell=ellipticity, pa=pa,
    #                    reff=reff, sres_ifile=sresfil, sres_fill=True, covar_ext=None,
    #                    plate=plate, ifudesign=ifudesign, ebvgal=ebvgal)
    cube = MUSEDataCube.from_config(config_file)
    cube.mask[0,0,:] = True


    # Define how you want to analyze the data
    plan = AnalysisPlanSet([ AnalysisPlan(drpqa_key='SNRG',
                                          drpqa_clobber=False,
                                          bin_key='SQUARE2.0', #HYB10', #SQUARE2.0
                                          bin_clobber=False,
                                          continuum_key='MILESHCMPL11',
                                          continuum_clobber=False,
                                          elmom_key='EMOMMPL11',
                                          elmom_clobber=False,
                                          elfit_key='EFITMPL11SSP', #'EFITMPL9DB',
                                          elfit_clobber=False,
                                          spindex_key='INDXEN',
                                          spindex_clobber=False) ])

    # Run it!
    return manga_dap(cube, plan, verbose=2, directory_path=directory_path,
                     analysis_path=analysis_path)
#-----------------------------------------------------------------------------


if __name__ == '__main__':

    directory_path = '/Users/erickaguirre/mangadap/examples/'
    #directory_path = '/Users/erickaguirre/Desktop/SDSU_Research/Getting_used_to_MaNGA_DAP/'
    cubefil = 'NGC0000.fits'
    #redshift = 0.008764

    #sresfil = '/Users/rubin/Research/MUSE/gistwork/Tutorial/gistTutorial/configFiles/LSF-Config_MUSE_WFM'

    # inclination = 35deg
    #ell = 0.0
    #pa = 145.0
    #reff = 59.7 # this is 95% radius for HII regions; Crocker+96
    #vdisp = 200.0
    #objra = 334.1301507   # move this into meta
    #objdec = -21.43893547

    # Need to make up plate and ifudesign numbers
    plate = 100000
    ifudesign = 1
    #ebvgal = 0.0


    #config_file = '/Users/erickaguirre/Desktop/SDSU_Research/Getting_used_to_MaNGA_DAP/mangadap-100000-1-LINCUBE.ini'
    config_file = '/Users/erickaguirre/Desktop/SDSU_Research/Getting_used_to_MaNGA_DAP/mangadap-100000-1-LINCUBE.ini'


    #test_read_muse(cubefil, sresfil, directory_path=directory_path, analysis_path='./output',
    #                  redshift=redshift, objra=objra, objdec=objdec, vdisp=vdisp, ellipticity=ell, pa=pa, reff=reff)

    # Need to add OBJRA and OBJDEC to prihdr -- issues with centering
    #test_wcs_muse(cubefil, sresfil, directory_path=directory_path, analysis_path='./output',
    #               redshift=redshift, objra=objra, objdec=objdec, vdisp=vdisp, ellipticity=ell, pa=pa, reff=reff)

    #test_copyto_muse(cubefil, sresfil, directory_path=directory_path, analysis_path='./output',
    #                    redshift=redshift, objra=objra, objdec=objdec, vdisp=vdisp, ellipticity=ell, pa=pa, reff=reff)
    #test_stats_muse(cubefil, sresfil, directory_path=directory_path, analysis_path='./output',
    #                    redshift=redshift, objra=objra, objdec=objdec, vdisp=vdisp, ellipticity=ell, pa=pa, reff=reff)

    #fit_one_cube_muse(cubefil, sresfil, directory_path=directory_path, analysis_path='./output',
    #                  redshift=redshift, objra=objra, objdec=objdec, vdisp=vdisp,
    #                  ellipticity=ell, pa=pa, reff=reff, plate=plate, ifudesign=ifudesign,
    #                  ebvgal=ebvgal)

    fit_one_cube_muse(config_file, directory_path=directory_path, analysis_path='./output2.0_test_mycode')
    #fit_one_cube_muse(config_file, directory_path=directory_path, analysis_path='./output2.0_NGC4790')

    # Run QA -- note that these plate-ifu numbers do not match those in new config file

    daptype='SQUARE2.0-MILESHC-MASTARSSP'
    #daptype = 'SQUARE0.6-MILESHC-MASTARSSP'
    dapver='2.2.3dev'
    drpver ='v3_0_1'
    maps_file = '/Users/erickaguirre/mangadap/examples/output2.0_test/SQUARE2.0-MILESHC-MASTARSSP/100000/1/manga-100000-1-MAPS-SQUARE2.0-MILESHC-MASTARSSP.fits.gz'
    model_file = '/Users/erickaguirre/mangadap/examples/output2.0_test/SQUARE2.0-MILESHC-MASTARSSP/100000/1/manga-100000-1-LOGCUBE-SQUARE2.0-MILESHC-MASTARSSP.fits.gz'

    #spotcheck_images('./outputVOR_test_OGcode', daptype, plate, ifudesign, ofile=None, drpver=None, dapver=None)

    plan = AnalysisPlanSet([AnalysisPlan(drpqa_key='SNRG',
                                         drpqa_clobber=False,
                                         bin_key='VOR10',  # 'HYB10', # 'VOR10', # 'SQUARE2.0'
                                         bin_clobber=False,
                                         continuum_key='MILESHCMPL11',
                                         continuum_clobber=False,
                                         elmom_key='EMOMMPL11',
                                         elmom_clobber=False,
                                         elfit_key='EFITMPL11SSP',  # 'EFITMPL9DB',
                                         spindex_key='INDXEN')])

    #redux_path = '/Users/erickaguirre/mangadap/examples/test_cube/'
    #ppxffit_qa_plot(plate, ifudesign, plan, drpver=drpver, redux_path=redux_path, dapver=dapver, analysis_path='./output2.0_test',
                    #tpl_flux_renorm=None,cubefile=cubefil)

    # Need to make up plate and ifudesign numbers
    plate = 1
    ifudesign = 1

    #fit_residuals_muse(dapver, './output2.0_NGC4030', daptype, plate, ifudesign)
    #manga_dap_inspector(maps_file, model_file, ext=None, masked_spectra=True)