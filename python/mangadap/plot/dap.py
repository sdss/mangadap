from __future__ import division, print_function, absolute_import

import os
from os.path import join
import sys

import numpy as np
import pandas as pd
import seaborn as sns
from astropy.io import fits

from sdss.files import base_path

from mangadap import dapfile

import util


class DAP():
    """Container for information from DAP FITS file.

    Attributes:
        fits (DAPFile instance): Instance of the mangadap.dapfile.DAPFile
            class.
        drps (DataFrame): Contents of 'DRPS' extension.

        n_bins (int): Number of bins.
        n_spaxels (int): Number of spaxels.
        sqrt_n_spaxels (int): Square root of number of spaxels.

        header (astropy.io.fits.header.Header): Header of primary HDU.
        tplkey (str): Stellar templates used.
        mangaid (str): MaNGA ID.

        sinames (list): Spectral index names.
        siunits (Series): Spectral index units.
        indx (DataFrame): Spectral index measurements.
        indxerr (DataFrame): Spectral index measurement errors.

    """

    def __init__(self, path_data, file_kws):
        self.read_fits(path_data, file_kws)

    def read_fits(self, path_data, file_kws):
        """Read in DAP FITS file.
        """
        self.fits = dapfile.dapfile(directory_path=path_data, **file_kws)
        self.fits.open_hdu()
        self.fits.read_par()

    def get_drps(self):
        """Read in DRPS extension.
        """
        drps_in = self.fits.read_hdu_data('DRPS')
        self.drps = util.fitsrec_to_dataframe(drps_in)

    def calc_n(self):
        self.n_bins, self.n_pix = self.galaxy.shape
        self.n_spaxels = len(self.drps)
        if self.fits.mode == 'CUBE':
            self.sqrt_n_spaxels = np.sqrt(self.n_spaxels)
            if np.mod(self.sqrt_n_spaxels, 1) != 0:
                print('Sqrt of number of bins in cube is not an integer.')
            else:
                self.sqrt_n_spaxels = int(self.sqrt_n_spaxels)

    def get_header(self):
        """Read header info"""
        self.header = self.fits.hdu[0].header
        self.tplkey = self.header['TPLKEY']
        try:
            self.mangaid = self.header['MANGAID']
        except KeyError:
            self.mangaid = 'n/a'

    def get_spectral_index(self):
        """Get spectral index info."""
        # Read in spectral index names
        self.sinames = util.read_line_names(self.fits, ltype='specind')
        self.siunits = pd.Series(util.read_hdu_par(self.fits, ext='UNIT',
                                 ltype='specind'), index=self.sinames)
        
        # Read in spectral index measurements
        self.indx = util.read_vals(self.fits, hdu='SINDX', ext='INDX',
                                   columns=self.sinames)
        self.indxerr = util.read_vals(self.fits, hdu='SINDX',
                                      ext='INDXERR', columns=self.sinames)

        # calculate combination indices
        self.calc_Fe5270_5335()
        self.calc_CalII0p86()

    def calc_Fe5270_5335(self):
        """Combine Fe5270 and Fe5335 spectral indices."""
        columns = ['Fe5270', 'Fe5335']
        coeffs = np.array([0.72, 0.28])
        self.indx['Fe5270_5335'] = util.linear_combination(self.indx, columns,
                                                           coeffs)
        self.indxerr['Fe5270_5335'] = util.linear_combination_err(
                                        self.indxerr, columns, coeffs)

    def calc_CalII0p86(self):
        """Combine CaII0p86A, CaII0p86B, and CaII0p86C spectral indices."""
        columns = ['CaII0p86A', 'CaII0p86B', 'CaII0p86C']
        coeffs = np.array([1/3., 1/3., 1/3.])
        self.indx['CaII0p86'] = util.linear_combination(self.indx, columns,
                                                        coeffs)
        self.indxerr['CaII0p86'] = util.linear_combination_err(
                                        self.indxerr, columns, coeffs)

    def get_emline(self):
        """Get emission line info."""
        # Read in spectral index names
        self.elnames = util.read_line_names(self.fits, ltype='emission')
        self.elopar = self.fits.read_hdu_data('ELOPAR')

        self.get_elofit()

    def get_elofit(self):
        """Read results from emission line fits into DataFrames."""
        self.elofit = self.fits.read_hdu_data('ELOFIT')

        elovel_kws = dict(dapf=self.fits, hdu='ELOFIT',
                          columns=['vel', 'vdisp'])
        elonames_kws = dict(dapf=self.fits, hdu='ELOFIT',
                            columns=self.elnames)
        
        # EW: Enci Wang's fitting code
        self.kin_ew = util.read_vals(ext='KIN_EW', **elovel_kws)
        self.kinerr_ew = util.read_vals(ext='KINERR_EW', **elovel_kws)
        self.elomit_ew = util.read_vals(ext='ELOMIT_EW', **elonames_kws)
        self.ampl_ew = util.read_vals(ext='AMPL_EW', **elonames_kws)
        self.amplerr_ew = util.read_vals(ext='AMPLERR_EW', **elonames_kws)
        self.sinst_ew = util.read_vals(ext='SINST_EW', **elonames_kws)
        self.flux_ew = util.read_vals(ext='FLUX_EW', **elonames_kws)
        self.fluxerr_ew = util.read_vals(ext='FLUXERR_EW', **elonames_kws)
        self.ew_ew = util.read_vals(ext='EW_EW', **elonames_kws)
        self.ewerr_ew = util.read_vals(ext='EWERR_EW', **elonames_kws)
        
        # split 3D arrays into two 2D DataFrames
        ikin_ew_tmp = self.elofit['IKIN_EW'].byteswap().newbyteorder()
        self.ivel_ew = pd.DataFrame(ikin_ew_tmp[:, 0], columns=self.elnames)
        self.ivdisp_ew = pd.DataFrame(ikin_ew_tmp[:, 1], columns=self.elnames)
        
        ikinerr_ew_tmp = self.elofit['IKINERR_EW'].byteswap().newbyteorder()
        self.ivelerr_ew = pd.DataFrame(ikinerr_ew_tmp[:, 0], columns=self.elnames)
        self.ivdisperr_ew = pd.DataFrame(ikinerr_ew_tmp[:, 1], columns=self.elnames)
        
        # FB: Francesco Belfiore's fitting code
        self.kin_fb = util.read_vals(ext='KIN_FB', **elovel_kws)
        self.kinerr_fb = util.read_vals(ext='KINERR_FB', **elovel_kws)
        self.elomit_fb = util.read_vals(ext='ELOMIT_FB', **elonames_kws)
        self.ampl_fb = util.read_vals(ext='AMPL_FB', **elonames_kws)
        self.amplerr_fb = util.read_vals(ext='AMPLERR_FB', **elonames_kws)
        self.sinst_fb = util.read_vals(ext='SINST_FB', **elonames_kws)
        self.flux_fb = util.read_vals(ext='FLUX_FB', **elonames_kws)
        self.fluxerr_fb = util.read_vals(ext='FLUXERR_FB', **elonames_kws)
        self.ew_fb = util.read_vals(ext='EW_FB', **elonames_kws)
        self.ewerr_fb = util.read_vals(ext='EWERR_FB', **elonames_kws)
        
        # split 3D arrays into two 2D DataFrames
        ikin_fb_tmp = self.elofit['IKIN_FB'].byteswap().newbyteorder()
        self.ivel_fb = pd.DataFrame(ikin_fb_tmp[:, 0], columns=self.elnames)
        self.ivdisp_fb = pd.DataFrame(ikin_fb_tmp[:, 1], columns=self.elnames)
        
        ikinerr_fb_tmp = self.elofit['IKINERR_FB'].byteswap().newbyteorder()
        self.ivelerr_fb = pd.DataFrame(ikinerr_fb_tmp[:, 0], columns=self.elnames)
        self.ivdisperr_fb = pd.DataFrame(ikinerr_fb_tmp[:, 1], columns=self.elnames)

