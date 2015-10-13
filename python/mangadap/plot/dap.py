# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-
"""Read DAP FITS file into a python object."""

from __future__ import (division, print_function, absolute_import,
                        unicode_literals)

import copy

import numpy as np
import pandas as pd

from mangadap import dapfile
from mangadap.plot import util

class DAP():
    """Convenience class for information from DAP FITS file.

    Args:
        path_data (str): Path to DAP FITS file.
        paths_cfg (str): Full path to sdss_paths.ini file.
        file_kws (dict): DAP FITS file specifications.
        mpl4 (bool): True if version is MPL-4. Default is True.
        drpqa_file (bool): True if FITS file was produced for DRPQA. Default
            is False.

    Attributes:
        path_data (str): Path to DAP FITS file.
        paths_cfg (str): Full path to sdss_paths.ini file.
        plate (str): Plate ID.
        ifudesign (str): IFU design.
        plateifu (str): Plate ID and IFU design.
        mode (str): Binning mode (i.e., 'CUBE' or 'RSS').
        bintype (str): Binning type (e.g., 'NONE' or 'ALL')
        niter (str): Analysis iteration number.
        mpl4 (bool): True if version is MPL-4. Default is True.
        drpqa_file (bool): True if FITS file was produced for DRPQA. Default
            is False.

        fits (DAPFile instance): Instance of the mangadap.dapfile.DAPFile
            class.
        header (astropy.io.fits.header.Header): Header of primary HDU.
        tplkey (str): Stellar templates used.
        mangaid (str): MaNGA ID.
        nsa_redshift (float): NSA redshift.

        drps (DataFrame): Position, signal, noise, de-redshifting, bin
            assignment, and bin weighting of each input spectrum from the DRP
            product being analyzed.
        bins (DataFrame): Location, area, S/N, number of spectra per bin, and
            analysis flags of each bin.

        wave_obs (array): Observed wavelengths (in Angstroms) of binned
            spectra.
        sres (array): Spectral resolution (lambda/delta lambda) in each
            wavelength channel.
        flux_obs (array): Observed flux (same as DRP units) of the binned
            spectra.
        ivar_obs (array): Observed inverse variance of the binned spectra.
        mask (array): Bad pixel mask of the binned spectra (0 = good,
            1 = bad).
        wave_rest (array): Rest frame wavelengths (in Angstroms) of binned
            spectra.
        wave_rest_med (array): Median rest frame wavelengths (in Angstroms).
        dlam (array): Size of wavelength bins divided by wavelength (inverse
            of resolution).
        flux_rest (array): Rest frame flux of the binned spectra.
        ivar_rest (array): Rest frame inverse variance of the binned spectra.
        smod_rest (array): Rest frame stellar continuum fit.

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

        elband (DataFrame): Definition of the emission-line bandpasses used
            for the non-parametric moments of the flux.

        elmmnt (DataFrame): Non-parametric moments of the binned spectra over
            the defined emission-line bands.

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
        ikin_ew (DataFrame): Best-fitting velocity and velocity dispersion (in
            km/s) of individual lines determined from the emission-line-only
            fits obtained using Enci Wang's code. Any omitted lines
            (ELOMIT_EW=1) should be ignored!
        ikinerr_ew (DataFrame): Best-fitting velocity and velocity dispersion
            errors (in km/s) of individual lines determined from the
            emission-line-only fits obtained using Enci Wang's code. Any
            omitted lines (ELOMIT_EW=1) should be ignored!
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
        elomew (array): Best-fitting emission-line-only spectrum using
            Enci Wang's code.
        fullfit_ew (array): Emission line fit from Enci Wang's code plus
            stellar continuum fit to provide a full spectral fit.
        fullfit_ew_rest (array): Emission line fit from Enci Wang's code plus
            stellar continuum fit to provide a full spectral fit that has been
            shifted into the rest frame.

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
        ikin_fb (DataFrame): Best-fitting velocity and velocity dispersion (in
            km/s) of individual lines determined from the emission-line-only
            fits obtained using Francesco Belfiore's code. Any omitted lines
            (ELOMIT_EW=1) should be ignored!
        ikinerr_fb (DataFrame): Best-fitting velocity and velocity dispersion
            errors (in km/s) of individual lines determined from the
            emission-line-only fits obtained using Francesco Belfiore's code.
            Any omitted lines (ELOMIT_EW=1) should be ignored!
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
        fullfit_fb (array): Emission line fit from Francesco Belfiore's code
            plus stellar continuum fit to provide a full spectral fit.
        fullfit_fb_rest (array): Emission line fit from Francesco Belfiore's
            code plus stellar continuum fit to provide a full spectral fit
            that has been shifted into the rest frame.
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
        sindx (DataFrame): Results of the I spectral index measurements for
            each of the binned spectra.

        redshift (array): Redshift (stellar velocity / speed of light) for
            each bin.
        median_redshift (array): Median redshift for bins.

        ind_lam_good (list): indices of wavelength array within each good
            window for each bin

        stfit_chisq_pix (array): Chisq of stellar continuum fit at every
            pixel.
        stfit_chisq_bin (array): Average chisq of stellar continuum fit for
            each binned spectrum.
        stfit_resid_pix (array): Difference between observed flux and stellar
            continuum fit at every pixel.
        stfit_resid_bin (array): Average difference between observed flux and
            stellar continuum fit for each binned spectrum.
        stfit_resid_mod_pix: Difference between observed flux and stellar
            continuum fit normalized to stellar continnum fit at every pixel.
        stfit_resid_mod_bin: Average difference between observed flux and
            stellar continuum fit normalized to stellar continnum fit for each
            binned spectrum.
        stfit_resid_err_pix: Difference between observed flux and stellar
            continuum fit normalized to error in stellar continnum fit at
            every pixel.
        stfit_resid_err_bin: Average difference between observed flux and
            stellar continuum fit normalized to error in stellar continnum fit
            for each binned spectrum.

        stfit_resid_data_pix (array): Difference between observed flux and
            stellar continuum fit normalized to observed flux at every pixel.
        stfit_resid_data_bin (array): Average difference between observed flux
            and stellar continuum fit normalized to observed flux for each
            binned spectrum.
        stfit_resid_data_bin99 (array): 99th percentile of difference between
            observed flux and stellar continuum fit normalized to observed
            flux for each binned spectrum.
        stfit_resid_data_bin68 (array): 68th percentile of difference between
            observed flux and stellar continuum fit normalized to observed
            flux for each binned spectrum.

        signal (array): DRP signal but only for bins with signal.
        noise (array): DRP noise but only for bins with signal.
        snr (array): DRP signal-to-noise ratio for bins with non-zero signal.
        drpqa (DataFrame): Values to display in DRP QA plots.
        drpqa_err (DataFrame): Errors for DRP QA plots.
        basic_qa (DataFrame): Values to display in basic QA plots.
        basic_qa_err (DataFrame): Errors for basic QA plots.
        kinematics (DataFrame): Values to display in kinematics plots.
        kinematics_err (DataFrame): Errors for kinematics plots.
    """

    def __init__(self, path_data, paths_cfg, file_kws, mpl4=True,
                 drpqa_file=False):
        self.path_data = path_data
        self.paths_cfg = paths_cfg
        self.plate = file_kws['plate']
        self.ifudesign = file_kws['ifudesign']
        self.plateifu = '{plate}-{ifudesign}'.format(**file_kws)
        self.mode = file_kws['mode']
        self.bintype = file_kws['bintype']
        self.niter = file_kws['niter']
        self.mpl4 = mpl4
        self.drpqa_file = drpqa_file
        self.read_dap_fits(path_data, file_kws)

    def read_dap_fits(self, path_data, file_kws):
        """Read in DAP FITS file."""
        self.fits = dapfile.dapfile(directory_path=path_data, **file_kws)
        self.fits.open_hdu()
        self.fits.read_par()
        self.extnames = [hdu._summary()[0] for hdu in self.fits.hdu]

    def get_all_ext(self):
        """Read in all extensions from FITS file."""
        self.get_header()
        #for ext in self.extnames:
        self.get_drps()
        self.get_bins()
        self.get_spectra()
        self.count_res_elements()
        self.get_elpar()

        self.get_stellar_cont_fit()
        self.get_stellar_gas_fit()

        self.get_elmod()
        if self.mpl4:
            self.get_elband()
            self.get_elmmnt()

        self.get_elopar()

        self.get_elofit()
        self.get_elomew()
        self.get_elomfb()

        if not self.drpqa_file:
            self.get_siwave()
            self.get_siflux()
            self.get_siivar()
            self.get_simask()
            self.get_siotpl()
            if self.mpl4:
                self.get_siotplm()
            self.get_sibotpl()
            if self.mpl4:
                self.get_sibotplm()
            self.get_sipar()
            self.get_sindx()

        if not self.drpqa_file:
            self.calc_fullfit()
            self.get_nsa_redshift()
            self.deredshift_spectra()
            self.deredshift_velocities()

        if not self.drpqa_file:
            self.select_wave_range()
            self.calc_stfit_chisq()
            self.calc_stfit_resid_data()
            self.calc_stfit_resid_mod()
            self.calc_stfit_resid_err()
            self.calc_stfit_resid()

        self.calc_snr()
        self.make_drpqa()
        if not self.drpqa_file:
            self.make_kinematics()

    def read_hdu(self, extname):
        try:
            return self.fits.read_hdu_data(extname)
        except KeyError:
            print('Extension {} not found.'.format(extname))
            return None

    def get_header(self):
        """Read header info"""
        self.header = self.fits.hdu[0].header
        try:
            self.tplkey = self.header['TPLKEY']
        except KeyError:
            self.tplkey = None
        try:
            self.mangaid = self.header['MANGAID']
        except KeyError:
            self.mangaid = 'n/a'

    def get_drps(self):
        """Read in DRPS extension."""
        drps_in = self.read_hdu('DRPS')
        if drps_in is not None:
            drps = util.fitsrec_to_dataframe(drps_in)
            self.drps = util.lowercase_colnames(drps)

    def get_bins(self):
        """Read in BINS extension"""
        bins_in = self.read_hdu('BINS')
        if bins_in is not None:
            bins = util.fitsrec_to_dataframe(bins_in)
            self.bins = util.lowercase_colnames(bins)

    def get_spectra(self):
        """Read in spectra.

        Read in observed wavelengths, fluxes, and inverse variances. Also read
        in spectral resolution and mask.
        """
        self.wave_obs = self.read_hdu('WAVE')
        self.sres = self.read_hdu('SRES')

        self.flux_obs = self.read_hdu('FLUX')
        if self.flux_obs is not None:
            self.flux_obs = np.transpose(self.flux_obs)


        self.ivar_obs = self.read_hdu('IVAR')
        if self.ivar_obs is not None:
            self.ivar_obs = np.transpose(self.ivar_obs)

        self.mask = self.read_hdu('MASK')
        if self.mask is not None:
            self.mask = np.transpose(self.mask)

    def get_elpar(self):
        """Read in emission line parameters."""
        elpar_in = self.read_hdu('ELPAR')
        if elpar_in is not None:
            elpar = util.fitsrec_to_dataframe(elpar_in)
            self.elpar = util.lowercase_colnames(elpar)
            self.elpar.elname = util.remove_hyphen(self.elpar.elname.values)
        else:
            self.elpar = None

    def get_stellar_cont_fit(self):
        """Read in stellar continuum fits to the binned spectra."""
        kincols = ['vel', 'vdisp']
        stfit_in = self.read_hdu('STFIT')
        if stfit_in is not None:
            self.stfit_tplw = stfit_in['TPLW']
            self.stfit_addpoly = stfit_in['ADDPOLY']
            self.stfit_multpoly = stfit_in['MULTPOLY']
            self.stfit_kin = util.make_df(stfit_in['KIN'], columns=kincols)
            self.stfit_kinerr = util.make_df(stfit_in['KINERR'],
                                             columns=kincols)
            self.stfit_rchi2 = stfit_in['RCHI2']
        else:
            for it in ['stfit_tplw', 'stfit_addpoly', 'stfit_multpoly',
                       'stfit_kin', 'stfit_kinerr', 'stfit_rchi2']:
                self.__dict__[it] = None            

        self.smsk = self.read_hdu('SMSK')
        if self.smsk is not None:
            self.smsk = np.transpose(self.smsk)

        self.smod = self.read_hdu('SMOD')
        if self.smod is not None:
            self.smod = np.transpose(self.smod)

    def get_stellar_gas_fit(self):
        """Read in star + gas fitting analysis."""
        #self.sgmsk = np.transpose(self.fits.read_hdu_data('SGMSK'))
        #self.sgmod = np.transpose(self.fits.read_hdu_data('SGMOD'))
        pass

    def get_elmod(self):
        """Read in best-fitting emission-line only model."""
        pass

    def get_elband(self):
        """Get bandpasses used for non-parametric fits of emission lines."""
        elband_in = self.fits.read_hdu_data('ELBAND')
        elnames_in = list(elband_in['ELNAME'].byteswap().newbyteorder())
        elnames = util.remove_hyphen(elnames_in)
        elband_tmp = {}
        for band in ['bandpass', 'blueside', 'redside']:
            for i, bedge in enumerate(['start', 'end']):
                elband_tmp['_'.join((band, bedge))] = \
                    elband_in[band][:, i].byteswap().newbyteorder()

        cols = ['bandpass_start', 'bandpass_end', 'blueside_start',
                'blueside_end', 'redside_start', 'redside_end']
        self.elband = pd.DataFrame(elband_tmp, columns=cols, index=elnames)

    def get_elmmnt(self):
        """Get moments of non-parametric fits of emission lines."""
        elmmnt_in = self.fits.read_hdu_data('ELMMNT')
        cols = [it.lower() for it in elmmnt_in.columns.names]
        self.elmmnt = util.fitsrec_to_multiindex_df(elmmnt_in, cols,
                                                    self.elband.index)

    def get_elopar(self):
        """Get emission line only fit parameters."""
        elopar_in = self.fits.read_hdu_data('ELOPAR')
        elopar = util.fitsrec_to_dataframe(elopar_in)
        self.elopar = util.lowercase_colnames(elopar)
        self.elopar.elname = util.remove_hyphen(self.elopar.elname)

    def get_elofit(self):
        """Read results from emission line fits into DataFrames."""
        elofit = self.fits.read_hdu_data('ELOFIT')

        if elofit['KIN_EW'].shape[1] == 2:
            kincols = ['vel', 'vdisp']
        else:
            kincols = ['vel_flux_wt', 'vdisp_flux_wt', 'vel_verr_wt',
                       'vdisp_verr_wt', 'vel_flux_verr_wt',
                       'vdisp_flux_verr_wt', 'vel_uni_wt', 'vdisp_uni_wt']
        ikincols = ['vel', 'vdisp']
        elonames = self.elopar.elname.values
        elovel_kws = dict(dapf=self.fits, hdu='ELOFIT', columns=kincols)
        elonames_kws = dict(dapf=self.fits, hdu='ELOFIT', columns=elonames)

        # EW: Enci Wang's fitting code
        self.kin_ew = util.read_vals(ext='KIN_EW', **elovel_kws)
        self.kinerr_ew = util.read_vals(ext='KINERR_EW', **elovel_kws)
        if self.mpl4:
            self.kinstde_ew = util.read_vals(ext='KINSTDE_EW', **elovel_kws)
            self.nkin_ew = elofit['NKIN_EW']
        self.elomit_ew = util.read_vals(ext='ELOMIT_EW', **elonames_kws)
        self.ampl_ew = util.read_vals(ext='AMPL_EW', **elonames_kws)
        self.amplerr_ew = util.read_vals(ext='AMPLERR_EW', **elonames_kws)
        self.sinst_ew = util.read_vals(ext='SINST_EW', **elonames_kws)
        self.flux_ew = util.read_vals(ext='FLUX_EW', **elonames_kws)
        self.fluxerr_ew = util.read_vals(ext='FLUXERR_EW', **elonames_kws)
        self.ew_ew = util.read_vals(ext='EW_EW', **elonames_kws)
        self.ewerr_ew = util.read_vals(ext='EWERR_EW', **elonames_kws)

        ikin_ew_in = elofit['IKIN_EW'].byteswap().newbyteorder()
        ikin_ew = np.transpose(ikin_ew_in, (1, 2, 0))
        self.ikin_ew = util.arr_to_multiindex_df(ikin_ew, ikincols, elonames)

        ikinerr_ew_in = elofit['IKINERR_EW'].byteswap().newbyteorder()
        ikinerr_ew = np.transpose(ikinerr_ew_in, (1, 2, 0))
        self.ikinerr_ew = util.arr_to_multiindex_df(ikinerr_ew, ikincols,
                                                    elonames)

        # FB: Francesco Belfiore's fitting code
        self.kin_fb = util.read_vals(ext='KIN_FB', **elovel_kws)
        self.kinerr_fb = util.read_vals(ext='KINERR_FB', **elovel_kws)
        if self.mpl4:
            self.kinstde_fb = util.read_vals(ext='KINSTDE_FB', **elovel_kws)
            self.nkin_fb = elofit['NKIN_FB']
        self.elomit_fb = util.read_vals(ext='ELOMIT_FB', **elonames_kws)
        self.ampl_fb = util.read_vals(ext='AMPL_FB', **elonames_kws)
        self.amplerr_fb = util.read_vals(ext='AMPLERR_FB', **elonames_kws)
        self.sinst_fb = util.read_vals(ext='SINST_FB', **elonames_kws)
        self.flux_fb = util.read_vals(ext='FLUX_FB', **elonames_kws)
        self.fluxerr_fb = util.read_vals(ext='FLUXERR_FB', **elonames_kws)
        self.ew_fb = util.read_vals(ext='EW_FB', **elonames_kws)
        self.ewerr_fb = util.read_vals(ext='EWERR_FB', **elonames_kws)

        ikin_fb_in = elofit['IKIN_FB'].byteswap().newbyteorder()
        ikin_fb = np.transpose(ikin_fb_in, (1, 2, 0))
        self.ikin_fb = util.arr_to_multiindex_df(ikin_fb, ikincols, elonames)

        ikinerr_fb_in = elofit['IKINERR_FB'].byteswap().newbyteorder()
        ikinerr_fb = np.transpose(ikinerr_fb_in, (1, 2, 0))
        self.ikinerr_fb = util.arr_to_multiindex_df(ikinerr_fb, ikincols,
                                                    elonames)

    def get_elomew(self):
        """Best fit emission-line-only spectrum (Enci Wang code)."""
        self.elomew = np.transpose(self.fits.read_hdu_data('ELOMEW'))

    def get_elomfb(self):
        """Best fit emission-line-only spectrum (Francesco Belfiore code)."""
        self.elomfb = np.transpose(self.fits.read_hdu_data('ELOMFB'))

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
        self.sinames = list(sipar_in['SINAME'].byteswap().newbyteorder())
        sipar_tmp = dict(unit=sipar_in['UNIT'].byteswap().newbyteorder())
        for band in ['passband', 'blueband', 'redband']:
            for i, bedge in enumerate(['start', 'end']):
                sipar_tmp['_'.join((band, bedge))] = \
                    sipar_in[band][:, i].byteswap().newbyteorder()

        cols = ['passband_start', 'passband_end', 'blueband_start',
                'blueband_end', 'redband_start', 'redband_end', 'unit']
        self.sipar = pd.DataFrame(sipar_tmp, columns=cols, index=self.sinames)
        #self.siunits = pd.Series(unit_tmp_swap, index=self.sinames)

    def get_sindx(self):
        """Read in spectral index info."""
        sindx = self.fits.read_hdu_data('SINDX')
        cols = [it.lower() for it in sindx.columns.names]
        self.sindx = util.fitsrec_to_multiindex_df(sindx, cols, self.sinames)

        # calculate combination indices
        self.calc_Fe5270_5335()
        self.calc_CalII0p86()

    def calc_Fe5270_5335(self):
        """Combine Fe5270 and Fe5335 spectral indices."""
        columns = ['Fe5270', 'Fe5335']
        coeffs = np.array([0.72, 0.28])
        # Adding columns to a multiindex dataframe does not keep it cubic (i.e.,
        # it does not populate NaNs for other columns). This could be a problem
        # if I want to slice by sindx name, which will be different for columns
        # with derived specinds.
        # HOW SHOULD I FIX THIS?
        self.sindx['indx', 'Fe5270_5335'] = util.lin_comb(
            self.sindx.indx, columns, coeffs)
        self.sindx['indxerr', 'Fe5270_5335'] = util.lin_comb_err(
            self.sindx.indxerr, columns, coeffs)

    def calc_CalII0p86(self):
        """Combine CaII0p86A, CaII0p86B, and CaII0p86C spectral indices."""
        columns = ['CaII0p86A', 'CaII0p86B', 'CaII0p86C']
        coeffs = np.array([1/3., 1/3., 1/3.])
        self.sindx['indx', 'CaII0p86'] = util.lin_comb(
            self.sindx.indx, columns, coeffs)
        self.sindx['indxerr', 'CaII0p86'] = util.lin_comb_err(
            self.sindx.indxerr, columns, coeffs)

    def count_res_elements(self):
        """Count bins, spaxels, and pixels."""
        self.n_bins, self.n_pix = self.flux_obs.shape
        self.n_spaxels = len(self.drps)
        if self.fits.mode == 'CUBE':
            self.sqrt_n_spaxels = np.sqrt(self.n_spaxels)
            if np.mod(self.sqrt_n_spaxels, 1) != 0:
                print('Sqrt of number of bins in cube is not an integer.')
            else:
                self.sqrt_n_spaxels = int(self.sqrt_n_spaxels)







    def calc_fullfit(self):
        """Add emission line only and stellar continuum fits."""
        self.fullfit_ew = self.elomew + self.smod
        self.fullfit_fb = self.elomfb + self.smod

    def deredshift_spectra(self):
        """Deredshift spectra while conserving flux."""
        v_light = 299792.458
        self.redshift = self.stfit_kin['vel'].values / v_light
        self.median_redshift = np.median(self.redshift)

        # de-redshift spectra
        wave_obs_grid = np.ones((self.n_bins, self.n_pix)) * self.wave_obs
        self.wave_rest = (wave_obs_grid.T / (1. + self.redshift)).T
        # FIX: what do we need dlam for?
        self.wave_rest_med = self.wave_obs / (1. + self.median_redshift)
        self.dlam = ((self.wave_rest_med[1] - self.wave_rest_med[0]) /
                     self.wave_rest_med[0])

        # conserve flux by multiplying by (1+z)
        self.flux_rest = (self.flux_obs.T * (1. + self.redshift)).T
        self.ivar_rest = (self.ivar_obs.T * (1. + self.redshift)).T
        self.smod_rest = (self.smod.T * (1. + self.redshift)).T
        self.fullfit_ew_rest = (self.fullfit_ew.T * (1. + self.redshift)).T
        self.fullfit_fb_rest = (self.fullfit_fb.T * (1. + self.redshift)).T

    def get_nsa_redshift(self):
        """Get NSA redshift from drpall file."""
        try:
            drpall = util.read_drpall(self.paths_cfg)
            mask = (drpall['plateifu'] == self.plateifu)
            self.nsa_redshift = drpall['nsa_redshift'][mask][0]
        except ValueError:
            print('\nCould not read NSA redshift from drpall file.\n')


    def deredshift_velocities(self):
        """Deredshift stellar and emission line velocities."""
        stkincols = ['vel', 'vdisp']
        elkincols = copy.deepcopy(stkincols)
        # Emission line kinematics are calculated using several methods in
        # MPL-4+. Select flux-weighted kinematics ("flux_wt") by default. 
        if 'vel_flux_wt' in list(self.kin_ew.columns.values):
            elkincols = [it + '_flux_wt' for it in elkincols]

        if hasattr(self, 'nsa_redshift'):

            stvel, stvelerr = util.deredshift_velocities(self.nsa_redshift,
                vel=self.stfit_kin['vel'].values,
                velerr=self.stfit_kinerr['vel'].values)
            elvel, elvelerr = util.deredshift_velocities(self.nsa_redshift,
                vel=self.kin_ew[elkincols[0]].values,
                velerr=self.kinerr_ew[elkincols[0]].values)

            stkin = dict(vel=stvel, vdisp=self.stfit_kin['vdisp'].values)
            stkinerr = dict(vel=stvelerr,
                            vdisp=self.stfit_kinerr['vdisp'].values)
            elkin = {elkincols[0]: elvel,
                     elkincols[1]: self.kin_ew[elkincols[1]].values}
            elkinerr = {elkincols[0]: elvelerr,
                        elkincols[1]: self.kinerr_ew[elkincols[1]].values}

            self.stfit_kin_rest = pd.DataFrame(stkin, columns=stkincols)
            self.stfit_kinerr_rest = pd.DataFrame(stkinerr, columns=stkincols)
            self.kin_rest_ew = pd.DataFrame(elkin, columns=elkincols)
            self.kinerr_rest_ew = pd.DataFrame(elkinerr, columns=elkincols)

        else:
            print('Could not deredshift velocities.\n')


    # FIX: make lam_good a kwarg param to pass into dap.DAP
    def select_wave_range(self, lam_good=None):
        """
        Select wavelength range for calculating goodness-of-fit metrics.

        Args:
            lam_good (array, optional): Ranges of wavelengths over which to
                evaluate goodness-of-fit. Default is::
                    np.array([[3650, 4200],  # ends at 4230 sky line
                              [4250, 5400],  # ends at 5575 sky line
                              [5450, 6732]]) # ends at [SII]

        """
        if lam_good is None:
            lam_good = np.array([[3650., 4200.],   # ends at 4230 sky line
                                 [4250., 5400.],   # ends at 5575 sky line
                                 [5450., 6732.]])  # ends at [SII]
                                 #[8470., 8770.]]) # CaT

        n_regions = len(lam_good)
        self.ind_lam_good = []
        for b in range(self.n_bins):
            # find indices of wave_rest within each good window
            tmp = []
            for i in range(n_regions):
                tmp.append(np.where((self.wave_rest[b] > lam_good[i, 0]) &
                                    (self.wave_rest[b] < lam_good[i, 1]))[0])
            ind_wave_tmp = np.concatenate([tmp[i] for i in range(n_regions)])
            self.ind_lam_good.append(ind_wave_tmp)

    def calc_metric_per_bin(self, metric_pix):
        """Calculate average goodness-of-fit metric for binned spectra.

        Args:
            metric_pix (array): Goodness-of-fit metric at each spectral pixel.

        Returns:
            array: Average goodness-of-fit metric for each binned spectrum.
        """
        metric_bin = np.ones(self.n_bins) * np.nan
        for i in range(self.n_bins):
            metric_bin[i] = np.nanmean(metric_pix[i][self.ind_lam_good[i]])
        return metric_bin

    def calc_stfit_chisq(self):
        """Calculate chisq of stellar continuum fit.

        Calculate chisq at each spectral pixel and for each binned spectrum.
        """
        self.stfit_chisq_pix = (self.flux_obs - self.smod)**2. * self.ivar_obs
        self.stfit_chisq_pix[self.smsk == 1] = np.nan
        self.stfit_chisq_bin = self.calc_metric_per_bin(self.stfit_chisq_pix)

    def calc_stfit_resid(self):
        """Calculate residual of stellar continuum fit.

        Calculate difference between observed flux and stellar continuum fit
        at each spectral pixel and for each binned spectrum.
        """
        self.stfit_resid_pix = self.flux_obs - self.smod
        self.stfit_resid_pix[self.smsk == 1] = np.nan
        self.stfit_resid_bin = self.calc_metric_per_bin(self.stfit_resid_pix)

    def calc_stfit_resid_mod(self):
        """Calculate residual of stellar continuum fit normalized to fit.

        Calculate difference between observed flux and stellar continuum fit
        normalized to stellar continnum fit at each spectral pixel and for
        each binned spectrum.
        """
        self.stfit_resid_mod_pix = (self.flux_obs - self.smod) / self.smod
        self.stfit_resid_mod_pix[self.smsk == 1] = np.nan
        self.stfit_resid_mod_bin = self.calc_metric_per_bin(
                                        self.stfit_resid_mod_pix)

    def calc_stfit_resid_err(self):
        """Calculate residual of stellar continuum fit normalized to fit error.

        Calculate difference between observed flux and stellar continuum fit
        normalized to error in stellar continnum fit at each spectral pixel
        and for each binned spectrum.
        """
        self.stfit_resid_err_pix = ((self.flux_obs - self.smod) *
                                    np.sqrt(self.ivar_obs))
        self.stfit_resid_err_pix[self.smsk == 1] = np.nan
        self.stfit_resid_err_bin = self.calc_metric_per_bin(
                                        self.stfit_resid_err_pix)

    def calc_stfit_resid_data(self):
        """Calculate residual of stellar cont fit normalized to observed flux.

        Calculate difference between observed flux and stellar continuum fit
        normalized to eobserve flux at each spectral pixel and for each binned
        spectrum. Also calculate the 99th and 68th percentiles of the residual
        relative to the observed flux.
        """
        with np.errstate(divide='ignore', invalid='ignore'):
            self.stfit_resid_data_pix = (np.abs(self.flux_obs - self.smod) /
                                         self.flux_obs)
        self.stfit_resid_data_pix[self.smsk == 1] = np.nan
        self.stfit_resid_data_bin = self.calc_metric_per_bin(
                                        self.stfit_resid_data_pix)
        self.stfit_resid_data_bin99 = np.array(
            [np.percentile(self.stfit_resid_data_pix[i][self.smsk[i] == 0],
                           99.) for i in range(self.n_bins)])
        self.stfit_resid_data_bin68 = np.array(
            [np.percentile(self.stfit_resid_data_pix[i][self.smsk[i] == 0],
                           68.) for i in range(self.n_bins)])

    def calc_snr(self):
        """Calculate DRP signal-to-noise ratio."""
        self.signal = np.zeros(self.n_bins)
        self.noise = np.zeros(self.n_bins)
        for i, b in enumerate(self.drps.binid.values):
            if b > -1:
                self.signal[b] = self.drps.signal.values[i]
                self.noise[b] = self.drps.noise.values[i]
        self.snr = self.signal / self.noise

    def make_drpqa(self):
        """Create DRP QA DataFrame."""
        vals = dict(signal=self.signal, noise=self.noise,
                    Ha6564=self.flux_ew.Ha6564.values)
        errs = dict(signal=None, noise=None,
                    Ha6564=self.fluxerr_ew.Ha6564.values)
        columns = ['signal', 'noise', 'Ha6564']
        self.drpqa = pd.DataFrame(vals, columns=columns)
        self.drpqa_err = pd.DataFrame(errs, columns=columns)

    def make_basic_qa(self):
        """Create basic QA DataFrame."""
        vals = dict(signal=self.signal, noise=self.noise, snr=self.snr,
                    Ha6564=self.flux_ew.Ha6564.values,
                    resid_data_bin99=self.stfit_resid_data_bin99,
                    stfit_chisq=self.stfit_chisq_bin)
        errs = dict(signal=None, noise=None, snr=None,
                    Ha6564=self.fluxerr_ew.Ha6564.values,
                    resid_data_bin99=None, stfit_chisq=None)
        columns = ['signal', 'noise', 'snr', 'Ha6564', 'resid_data_bin99',
                   'stfit_chisq']
        self.basic_qa = pd.DataFrame(vals, columns=columns)
        self.basic_qa_err = pd.DataFrame(errs, columns=columns)

    def make_kinematics(self):
        """Create kinematics DataFrame."""
        elkincols = list(self.kin_ew.columns.values)
        vals = dict(stvel=self.stfit_kin_rest['vel'].values,
                    stvdisp=self.stfit_kin['vdisp'].values,
                    stfit_chisq=self.stfit_chisq_bin,
                    elvel=self.kin_rest_ew[elkincols[0]].values,
                    elvdisp=self.kin_ew[elkincols[1]].values,
                    stfit_resid99=self.stfit_resid_data_bin99)
        errs = dict(stvel=self.stfit_kinerr_rest['vel'].values,
                    stvdisp=self.stfit_kinerr['vdisp'].values,
                    stfit_chisq=None,
                    elvel=self.kinerr_rest_ew[elkincols[0]],
                    elvdisp=self.kinerr_ew[elkincols[1]],
                    stfit_resid99=None)
        columns = ['stvel', 'stvdisp', 'stfit_chisq', 'elvel', 'elvdisp',
                   'stfit_resid99']
        self.kinematics = pd.DataFrame(vals, columns=columns)
        self.kinematics_err = pd.DataFrame(errs, columns=columns)
