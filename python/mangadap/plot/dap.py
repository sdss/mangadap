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
        header (astropy.io.fits.header.Header): Header of primary HDU.
        tplkey (str): Stellar templates used.
        mangaid (str): MaNGA ID.

        drps (DataFrame): Position, signal, noise, de-redshifting, bin
            assignment, and bin weighting of each input spectrum from the DRP
            product being analyzed.
        bins (DataFrame): Location, area, S/N, number of spectra per bin, and
            analysis flags of each bin. 

        wave_obs (array): Observed wavelengths (in Angstroms) of binned
            spectra.
        sres (array): Spectral resolution (lambda/delta lambda) in each
            wavelength channel.
        flux_obs (array): Flux (same as DRP units) of the binned spectra.
        ivar_obs (array): Inverse variance of the Nbins binned spectra.
        mask (array): Bad pixel mask of the binned spectra (0 = good,
            1 = bad).

        n_pix (int): Number of spectral channels.
        n_bins (int): Number of bins.
        n_spaxels (int): Number of spaxels.
        sqrt_n_spaxels (int): Square root of number of spaxels.

        elpar (DataFrame): Emission-line parameters used by the
            emission-line analysis (BHA 8/31/15: Is this accurate?).

        stfit_tplw (array): Optimal template weights in the stellar-continuum
            (pPXF) analysis
        stfit_addpoly (array): Coefficients of the additive legendre
            polynomials used in the stellar-continuum (pPXF) analysis.
        stfit_multpoly (array): Coefficients of the multiplicative legendre
            polynomials used in the stellar-continuum (pPXF) analysis.
        stfit_kin (DataFrame): Best-fitting parameters of the stellar LOSVD
            (velocity and velocity dispersion in km/s) determined by the
            stellar-continuum (pPXF) analysis.
        stfit_kinerr (DataFrame): Best-fitting parameters of the stellar LOSVD
            (velocity and velocity dispersion in km/s) determined by the
            stellar-continuum (pPXF) analysis.
        stfit_rchi2 (array): Chi-square per degree of freedom for the best-fit
            obtained by the stellar-continuum (pPXF) analysis. 

        smsk (array): Flag that pixel was/was not (0/1) included in the
            stellar continuum fit.
        smod (array): Best-fitting model returned by the stellar continuum
            fit.

        kin_ew (DataFrame): Gas kinematics (velocity and velocity dispersion
            in km/s) based on Enci Wang's code averaged over the valid
            individual measurements. Right now, the gas kinematics are only
            derived assuming a Gaussian LOSVD (no higher moments are fit).
        kinerr_ew (DataFrame): Gas kinematics errors (velocity and velocity
            dispersion in km/s) based on Enci Wang's code averaged over the
            valid individual measurements. Right now, the gas kinematics are
            only derived assuming a Gaussian LOSVD (no higher moments are
            fit).
        elomit_ew (DataFrame): Flag setting whether an emission line was/was
            not (1/0) omitted during the emission-line fitting using
            Enci Wang's code.
        ampl_ew (DataFrame): Best-fitting amplitudes of the emission lines
            obtained using Enci Wang's code in the same flux density units as
            the DRP spectra. Any omitted lines (ELOMIT_EW=1) should be
            ignored!
        amplerr_ew (DataFrame): Best-fitting amplitude errors of the emission
            lines obtained using Enci Wang's code in the same flux density
            units as the DRP spectra. Any omitted lines (ELOMIT_EW=1) should
            be ignored!
        ivel_ew (DataFrame): Best-fitting velocity (in km/s) of individual
            lines determined from the emission-line-only fits obtained using
            Enci Wang's code. Any omitted lines (ELOMIT_EW=1) should be
            ignored! 
        ivelerr_ew (DataFrame): Best-fitting velocity errors (in km/s) of
            individual lines determined from the emission-line-only fits
            obtained using Enci Wang's code. Any omitted lines (ELOMIT_EW=1)
            should be ignored! 
        ivdisp_ew (DataFrame): Best-fitting velocity dispersion (in km/s) of
            individual lines determined from the emission-line-only fits
            obtained using Enci Wang's code. Any omitted lines (ELOMIT_EW=1)
            should be ignored!
        ivdisperr_ew (DataFrame): Best-fitting velocity dispersion errors (in
            km/s) of individual lines determined from the emission-line-only
            fits obtained using Enci Wang's code. Any omitted lines
            (ELOMIT_EW=1) should be ignored! 
        sinst_ew (DataFrame): Interpolated spectral resolution (sigma in km/s)
            at the centroids of the fitted emission lines obtained using
            Enci Wang's code. Any omitted lines (ELOMIT_EW=1) should be
            ignored! 
        flux_ew (DataFrame): Best-fitting fluxes of the emission lines
            obtained by Enci Wang's code in the same flux units as the DRP
            spectra, integrated over wavelength. Any omitted lines
            (ELOMIT_EW=1) should be ignored! 
        fluxerr_ew (DataFrame): Best-fitting flux errors of the emission lines
            obtained by Enci Wang's code in the same flux units as the DRP
            spectra, integrated over wavelength. Any omitted lines
            (ELOMIT_EW=1) should be ignored! 
        ew_ew (DataFrame): Best-fitting equivalent widths (in angstroms) of
            the emission lines obtained by Enci Wang's code. Any omitted
            lines (ELOMIT_EW=1) should be ignored!
        ewerr_ew (DataFrame): Best-fitting equivalent width errors (in
            angstroms) of the emission lines obtained by Enci Wang's code. Any
            omitted lines (ELOMIT_EW=1) should be ignored!
        
        kin_fb (DataFrame): Gas kinematics (velocity and velocity dispersion
            in km/s) based on Francesco Belfiore's code averaged over the valid
            individual measurements. Right now, the gas kinematics are only
            derived assuming a Gaussian LOSVD (no higher moments are fit).
        kinerr_fb (DataFrame): Gas kinematics errors (velocity and velocity
            dispersion in km/s) based on Francesco Belfiore's code averaged
            over the valid individual measurements. Right now, the gas
            kinematics are only derived assuming a Gaussian LOSVD (no higher
            moments are fit).
        elomit_fb (DataFrame): Flag setting whether an emission line was/was
            not (1/0) omitted during the emission-line fitting using
            Francesco Belfiore's code.
        ampl_fb (DataFrame): Best-fitting amplitudes of the emission lines
            obtained using Francesco Belfiore's code in the same flux density
            units as the DRP spectra. Any omitted lines (ELOMIT_FB=1) should
            be ignored!
        amplerr_fb (DataFrame): Best-fitting amplitude errors of the emission
            lines obtained using Francesco Belfiore's code in the same flux
            density units as the DRP spectra. Any omitted lines (ELOMIT_FB=1)
            should be ignored!
        ivel_fb (DataFrame): Best-fitting velocity (in km/s) of individual
            lines determined from the emission-line-only fits obtained using
            Francesco Belfiore's code. Any omitted lines (ELOMIT_FB=1) should
            be ignored! 
        ivelerr_fb (DataFrame): Best-fitting velocity errors (in km/s) of
            individual lines determined from the emission-line-only fits
            obtained using Francesco Belfiore's code. Any omitted lines
            (ELOMIT_FB=1) should be ignored! 
        ivdisp_fb (DataFrame): Best-fitting velocity dispersion (in km/s) of
            individual lines determined from the emission-line-only fits
            obtained using Francesco Belfiore's code. Any omitted lines
            (ELOMIT_FB=1) should be ignored!
        ivdisperr_fb (DataFrame): Best-fitting velocity dispersion errors (in
            km/s) of individual lines determined from the emission-line-only
            fits obtained using Francesco Belfiore's code. Any omitted lines
            (ELOMIT_FB=1) should be ignored! 
        sinst_fb (DataFrame): Interpolated spectral resolution (sigma in km/s)
            at the centroids of the fitted emission lines obtained using
            Francesco Belfiore's code. Any omitted lines (ELOMIT_FB=1) should
            be ignored! 
        flux_fb (DataFrame): Best-fitting fluxes of the emission lines
            obtained by Francesco Belfiore's code in the same flux units as
            the DRP spectra, integrated over wavelength. Any omitted lines
            (ELOMIT_FB=1) should be ignored! 
        fluxerr_fb (DataFrame): Best-fitting flux errors of the emission lines
            obtained by Francesco Belfiore's code in the same flux units as
            the DRP spectra, integrated over wavelength. Any omitted lines
            (ELOMIT_FB=1) should be ignored! 
        ew_fb (DataFrame): Best-fitting equivalent widths (in angstroms) of
            the emission lines obtained by Francesco Belfiore's code. Any
            omitted lines (ELOMIT_FB=1) should be ignored!
        ewerr_fb (DataFrame): Best-fitting equivalent width errors (in
            angstroms) of the emission lines obtained by Francesco Belfiore's
            code. Any omitted lines (ELOMIT_FB=1) should be ignored!
        
        elomew (array): Best-fitting emission-line-only spectrum using
            Enci Wang's code. 
        elomfb (array): Best-fitting emission-line-only spectrum using
            Francesco Belfiore's code. 

        siwave (array): Wavelength in angstrom of D spectral channels in the
            spectra that have had their resolution matched to the
            spectral-index system.

        siflux (array): Flux (same as DRP units) of the object spectra (if
            fit, the emission-line model has been removed) with the resolution
            matched to the spectral-index system.

        siivar (array): Inverse variance in the SIFLUX image.

        simask (array): Bad pixel mask (0 = good, 1 = bad) for the SIFLUX
            image.

        siotpl (array): Optimal templates, ignoring any fitted polynomial and
            LOSVD, with the resolution matched to the spectral index system.

        siotplm (array): Mask for the unbroadened optimal templates.

        sibotpl (array): Optimal templates, ignoring any fitted polynomial,
            with the resolution matched to the spectral index system.

        sibotplm (array): Mask for the broadened optimal templates.

        sinames (list): Spectral index names.
        siomit (DataFrame): A flag giving if the index was (1) or was not (0)
            omitted from the analysis. The index will be omitted if any part
            of its unmasked wavelength range or that of the (unbroadened or
            broadened) optimal template (blueband start to redband end) falls
            outside of the observed spectral range. 
        indx (DataFrame): The "corrected" spectral indices, where the
            corrections are for the LOSVD based on the difference between the
            indices measured using the unbroadened and broadened optimal
            template spectra.
        indxerr (DataFrame): Spectral index measurement errors.
        indx_otpl (DataFrame): The measured spectral index using the optimal
            template with units given by the UNIT parameter in the SIPAR
            extension. These measurements are used to correct the indices
            measured for the object spectra for the effect of the LOSVD.
        indx_botpl (DataFrame): The measured spectral index using the
            broadened optimal template with units given by the UNIT
            parameter in the SIPAR extension. These measurements are used to
            correct the indices measured for the object spectra for the effect
            of the LOSVD.

    """

    def __init__(self, path_data, file_kws):
        self.read_fits(path_data, file_kws)

    def read_fits(self, path_data, file_kws):
        """Read in DAP FITS file."""
        self.fits = dapfile.dapfile(directory_path=path_data, **file_kws)
        self.fits.open_hdu()
        self.fits.read_par()

    def get_all_ext(self):
        """Read in all extensions from FITS file."""
        self.get_header()
        self.get_drps()
        self.get_bins()
        self.get_spectra()
        self.count_res_elements()
        self.get_elpar()
        self.get_stellar_cont_fit()
        # self.get_stellar_gas_fit()
        # self.get_elmod()

        self.get_emline_fit()
        self.get_elofit()

        self.get_elomew()
        self.get_elomfb()
        self.get_siwave()
        self.get_siflux()
        self.get_siivar()
        self.get_simask()
        self.get_siotpl()
        #self.get_siotplm()
        self.get_sibotpl()
        #self.get_sibotplm()
        self.get_sipar()
        self.get_sindx()

    def get_header(self):
        """Read header info"""
        self.header = self.fits.hdu[0].header
        self.tplkey = self.header['TPLKEY']
        try:
            self.mangaid = self.header['MANGAID']
        except KeyError:
            self.mangaid = 'n/a'

    def get_drps(self):
        """Read in DRPS extension."""
        drps_in = self.fits.read_hdu_data('DRPS')
        drps = util.fitsrec_to_dataframe(drps_in)
        self.drps = self.lowercase_colnames(drps)

    def get_bins(self):
        """Read in BINS extension"""
        bins_in = self.fits.read_hdu_data('BINS')
        bins = util.fitsrec_to_dataframe(bins_in)
        self.bins = self.lowercase_colnames(bins)

    def get_spectra(self):
        """Read in spectra.

        Read in observed wavelengths, fluxes, and inverse variances. Also read
        in spectral resolution and mask.
        """
        self.wave_obs = self.fits.read_hdu_data('WAVE')
        self.sres = self.fits.read_hdu_data('SRES')
        self.flux_obs = self.fits.read_hdu_data('FLUX')
        self.ivar_obs = self.fits.read_hdu_data('IVAR')
        self.mask = self.fits.read_hdu_data('MASK')

    def get_elpar(self):
        """Read in emission line parameters."""
        elpar_in = self.fits.read_hdu_data('ELPAR')
        elpar = util.fitsrec_to_dataframe(elpar_in)
        self.elpar = self.lowercase_colnames(elpar)
        self.elpar.elname = self.remove_hyphen(self.elpar.elname.values)

    def get_stellar_cont_fit(self):
        """Read in stellar continuum fits to the binned spectra."""
        kincols = ['vel', 'vdisp']
        stfit_in = self.fits.read_hdu_data('STFIT')
        self.stfit_tplw = stfit_in['TPLW']
        self.stfit_addpoly = stfit_in['ADDPOLY']
        self.stfit_multpoly = stfit_in['MULTPOLY']
        self.stfit_kin = pd.DataFrame(stfit_in['KIN'], columns=kincols)
        self.stfit_kinerr = pd.DataFrame(stfit_in['KINERR'], columns=kincols)
        self.stfit_rchi2 = stfit_in['RCHI2']
        self.smsk = self.fits.read_hdu_data('SMSK')
        self.smod = self.fits.read_hdu_data('SMOD')

    def get_stellar_gas_fit(self):
        """Read in star + gas fitting analysis."""
        pass

    def get_elmod(self):
        """Read in best-fitting emission-line only model."""
        pass

    def get_emline_fit(self):
        """Get emission line only fits."""
        # Read in spectral index names
        elopar_in = self.fits.read_hdu_data('ELOPAR')
        self.elopar = util.fitsrec_to_dataframe(elopar_in)
        self.elopar.ELNAME = self.remove_hyphen(self.elopar.ELNAME)
        
        self.get_elofit()
        self.get_elomew()
        self.get_elomfb()

    def get_elofit(self):
        """Read results from emission line fits into DataFrames."""
        elofit = self.fits.read_hdu_data('ELOFIT')
    
        elonames = self.elopar.ELNAME.values
        elovel_kws = dict(dapf=self.fits, hdu='ELOFIT',
                          columns=['vel', 'vdisp'])
        elonames_kws = dict(dapf=self.fits, hdu='ELOFIT', columns=elonames)
        
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
        ikin_ew = elofit['IKIN_EW'].byteswap().newbyteorder()
        self.ivel_ew = pd.DataFrame(ikin_ew[:, 0], columns=elonames)
        self.ivdisp_ew = pd.DataFrame(ikin_ew[:, 1], columns=elonames)
        
        ikinerr_ew = elofit['IKINERR_EW'].byteswap().newbyteorder()
        self.ivelerr_ew = pd.DataFrame(ikinerr_ew[:, 0], columns=elonames)
        self.ivdisperr_ew = pd.DataFrame(ikinerr_ew[:, 1], columns=elonames)
        
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
        ikin_fb = elofit['IKIN_FB'].byteswap().newbyteorder()
        self.ivel_fb = pd.DataFrame(ikin_fb[:, 0], columns=elonames)
        self.ivdisp_fb = pd.DataFrame(ikin_fb[:, 1], columns=elonames)
        
        ikinerr_fb = elofit['IKINERR_FB'].byteswap().newbyteorder()
        self.ivelerr_fb = pd.DataFrame(ikinerr_fb[:, 0], columns=elonames)
        self.ivdisperr_fb = pd.DataFrame(ikinerr_fb[:, 1], columns=elonames)

    def get_elomew(self):
        """Best fit emission-line-only spectrum (Enci Wang code)."""
        self.elomew = self.fits.read_hdu_data('ELOMEW')

    def get_elomfb(self):
        """Best fit emission-line-only spectrum (Francesco Belfiore code)."""
        self.elomfb = self.fits.read_hdu_data('ELOMFB')

    def get_siwave(self):
        """Resolution matched wavelengths to spectral-index system."""
        self.siwave = self.fits.read_hdu_data('SIWAVE')

    def get_siflux(self):
        """Resolution matched fluxes to spectral-index system."""
        self.siflux = self.fits.read_hdu_data('SIFLUX')

    def get_siivar(self):
        """Resolution matched inverse variances to spectral-index system."""
        self.siivar = self.fits.read_hdu_data('SIIVAR')

    def get_simask(self):
        """Bad pixel mask (0 = good, 1 = bad) for siflux."""
        self.simask = self.fits.read_hdu_data('SIMASK')

    def get_siotpl(self):
        """Optimal templates resolution matched to spectral index system."""
        self.siotpl = self.fits.read_hdu_data('SIOTPL')

    def get_siotplm(self):
        """Mask for the unbroadened optimal templates."""
        self.siotplm = self.fits.read_hdu_data('SIOTPLM')

    def get_sibotpl(self):
        """Broadened optimal templates resolution matched to spectral index
        system."""
        self.sibotpl = self.fits.read_hdu_data('SIBOTPL')

    def get_sibotplm(self):
        """Mask for the broadened optimal templates."""
        self.sibotplm = self.fits.read_hdu_data('SIBOTPLM')

    def get_sipar(self):
        """Spectral index parameters."""
        sipar_in = self.fits.read_hdu_data('SIPAR')
        self.sinames = list(sipar_in['SINAME'])
        sipar_tmp = dict(unit=sipar_in['UNIT'])
        for band in ['passband', 'blueband', 'redband']:
            for i, bedge in enumerate(['start', 'end']):
                sipar_tmp['_'.join((band, bedge))] = sipar_in[band][:, i]

        cols = ['passband_start', 'passband_end', 'blueband_start',
            'blueband_end', 'redband_start', 'redband_end', 'unit']
        self.sipar = pd.DataFrame(sipar_tmp, columns=cols, index=self.sinames)
        

    def get_sindx(self):
        """Read in spectral index info."""
        nm = self.sinames
        self.sindx = self.fits.read_hdu_data('SINDX')
        self.siomit = util.swap_byte(self.sindx['SIOMIT'], columns=nm)
        self.indx = util.swap_byte(self.sindx['INDX'], columns=nm)
        self.indxerr = util.swap_byte(self.sindx['INDXERR'], columns=nm)
        self.indx_otpl = util.swap_byte(self.sindx['INDX_OTPL'], columns=nm)
        self.indx_botpl = util.swap_byte(self.sindx['INDX_BOTPL'], columns=nm)
        
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

    def count_res_elements(self):
        """Count bins, spaxels, and pixels."""
        self.n_pix, self.n_bins = self.flux_obs.shape
        self.n_spaxels = len(self.drps)
        if self.fits.mode == 'CUBE':
            self.sqrt_n_spaxels = np.sqrt(self.n_spaxels)
            if np.mod(self.sqrt_n_spaxels, 1) != 0:
                print('Sqrt of number of bins in cube is not an integer.')
            else:
                self.sqrt_n_spaxels = int(self.sqrt_n_spaxels)

    def remove_hyphen(self, names):
        """Remove hyphens from emission line names."""
        return [name.replace('-', '').strip() for name in names]

    def lowercase_colnames(self, df):
        """Convert column names of a DataFrame to lowercase."""
        df.columns = [item.lower() for item in df.columns]
        return df

