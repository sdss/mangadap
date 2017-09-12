#!/usr/bin/env python

import os
import warnings
import numpy as np
import matplotlib.pyplot as plt
import astropy.constants

from captools.ppxf import ppxf
from captools import capfit

# Function used to calculate 1 sigma errors of the residuals
# Currently the pixels in the same window are assigned with
# the same error
def calculate_noise(residuals):

    # Number of pixels in each window
    width = 140

    # Number of windows
    n_win = (residuals.shape[0]//width)

    # Assign each window with equal number of pixels
    cut = n_win*width
    resid_temp = residuals[:cut].reshape(n_win, -1)

    # Calculate noise (standard deviation) for each window
    stdd = np.std(resid_temp, axis=1, ddof=1)
    stdd = stdd.reshape(n_win, -1)

    # Create the noise array with the same size as galaxy array
    # There are some pixels uncovered by any window, which are 
    # assigned with the same sigma as the nearest window
    noise = np.repeat(stdd, width, axis=1).ravel()
    uncover = np.repeat(stdd[-1], residuals.shape[0]-cut)
    noise = np.append(noise, uncover)
    
    return noise

# Function: emline_fitter_with_ppxf
# Input is:
#   - wave0: wavelength vector (in vacuum); shape is (Nwave,)
#   - flux: observed, unbinned flux; masked array with shape
#     (Nspaxel,Nwave)
#   - noise: error in observed, unbinned flux; masked array with
#     shape (Nspaxel,Nwave)
#   - sres: spectral resolution (R=lambda/delta lambda) as a
#     function of wavelength for each unbinned spectrum; shape
#     is (Nspaxel,Nwave)
#   - flux_binned: binned flux; masked array with shape
#     (Nbin,Nwave)
#   - noise_binned: noise in binned flux; masked array with
#     shape (Nbin,Nwave)
#   - sres_binned: spectral resolution (R=lambda/delta lambda)
#     as a function of wavelength for each binned spectrum;
#     shape is (Nbin,Nwave)
#   - velscale: Velocity step per pixel
#   - velscale_ratio: Ratio of velocity step per pixel in the
#     observed data versus in the template data
#   - dv: Velocity offset between the galaxy and template data
#     due to the difference in the initial wavelength of the
#     spectra
#   - stars_vel: Velocity of the stellar component; shape is
#     (Nbins,)
#   - stars_sig: Velocity dispersion of the stellar component;
#     shape is (Nbins,)
#   - templates: Template flux; shape is (Ntpl,Ntplwave)
#   - gas_tpl: Boolean array for extracting gas templates; shape
#     is (Ntpl,)
#   - guess_vel: Initial guess velocity for the gas components;
#     shape is (Nbins,)
#   - guess_sig: Initial guess velocity dispersion for the gas
#     components; shape is (Nbins,)
#   - gas_names: Name of the gas templats; shape is (Ngastpl,)
#   - eml_wave: Rest wavelength of gas templates (in vacuum);
#     shape is (Ngastpl,)
#   - eml_component: Component vector for ppxf (each kinimatic
#     group has a unique index, and an index below zero will
#     automatically fix the moments); shape is (Nkin_group,)
#   - eml_tied: Tied vector for ppxf (to decide the tying 
#     relations among gas components); shape is (Nkin_group,
#     Nmoment)
#   - template_sres: spectral resolution (R=lambda/delta
#     lambda) as a function of wavelength for all the templates
#     templates; shape is (Ntplwave,)
#   - degree: Additive polynomial order
#   - mdegree: Multiplicative polynomial order
#   - x: On-sky spaxel x coordinates; shape is (Nspaxel,)
#   - y: On-sky spaxel y coordinates; shape is (Nspaxel,)
#   - x_binned: On-sky bin x coordinate; shape is (Nbin,)
#   - y_binned: On-sky bin y coordinate; shape is (Nbin,)
#   - mask_binned: Masks for binned spectra; shape is (Nbin,
#     Nwave)
#   - mask_drp: Masks for unbinned spectra; shape is (Nbin,
#     Nwave)
#   - nsa_z: NSA (or guess) redshift of the galaxy
#
# Output is:
#   - model_flux: stellar-continuum + emission-line model; shape
#     is (Nmod, Nwave); first axis is ordered by model ID number
#   - model_eml_flux: model emission-line flux only; shape is
#     (Nmod, Nwave); first axis is ordered by model ID number
#   - model_mask: boolean or bit mask for fitted models; shape
#     is (Nmod, Nwave); first axis is ordered by model ID number
#   - model_binid: ID numbers assigned to each spaxel with a
#     fitted model; any spaxel without a model should have
#     model_binid = -1; the number of >-1 IDs must be Nmod;
#     shape is (Nspaxel,)
#   - weights: weights of the templates (the best fit template 
#     can be obtained by setting: bestemp = templates @ weights);
#     shape is (Nmod, Ntpl)
#   - weights_err: 1 sigma formal error of weights; shape is
#     (Nmod, Ntpl)
#   - mweights: coefficients of multiplicative Legendre polynomials
#     of order > 0; shape is (Nmod, mdegree)
#   - eml_kin: Kinematics (velocity and velocity dispersion) of
#     each emission line; shape is (Nmod,Neml,Nkin)
#   - eml_kinerr: Error in the kinematics of each emission line
#   - eml_sigmacorr: Quadrature corrections required to obtain
#     the astrophysical velocity dispersion; shape is
#     (Nmod,Neml); corrections are expected to be applied as
#     follows:
#       sigma = numpy.ma.sqrt( numpy.square(eml_kin[:,:,1])
#                               - eml_sigmacorr)
#     NOTE THAT: eml_sigmacorr here are actually in unit of 
#     (km/s)^2 and could have negative values
# --------------------------------------------------------------

def emline_fitter_with_ppxf(wave0, flux, noise, sres, flux_binned, noise_bin,
                            velscale, velscale_ratio, dv, vel, sig, templates, gas_tpl, 
                            vel_gas, sig_gas, gas_names, eml_wave, eml_component, 
                            eml_tied, templates_sres, degree, mdegree, x, y, xbin, 
                            ybin, mask_binned, mask_drp, nsa_z, debug=False,
                            mode='bin_to_spaxel'):

    if debug:
        warnings.warn('JUST DEBUGGING.  NO EMISSION-LINE FITS PERFORMED!!')
        n_spaxels, nwave = flux.shape
        n_eml = np.sum(gas_tpl)
        masked_fraction = np.sum(mask_drp, axis=1)/nwave
        indx = masked_fraction > 1/3
        model_binid = np.full(n_spaxels, -1, dtype=int)
        model_binid[indx] = np.arange(np.sum(indx))

        model_flux = np.zeros(flux.shape, dtype=float)
        model_eml_flux = np.zeros(flux.shape, dtype=float)
        model_mask = np.zeros(flux.shape, dtype=float)
        
        eml_flux = np.zeros([n_spaxels,n_eml], dtype=float)
        eml_fluxerr = np.zeros([n_spaxels,n_eml], dtype=float)
        eml_kin = np.zeros([n_spaxels,n_eml,2], dtype=float)
        eml_kinerr = np.zeros([n_spaxels,n_eml,2], dtype=float)
    
        eml_sigmacorr = np.zeros([n_spaxels,n_eml], dtype=float)
        return model_flux, model_eml_flux, model_mask, model_binid, eml_flux, eml_fluxerr, \
                    eml_kin, eml_kinerr, eml_sigmacorr

#    plate = '8258'
#    ifu = '6102'
#    npz_file = '/Users/mac/Documents/New_test_results/first_'+plate+'_'+ifu+'.npz'
#    npz_file = 'test.npz'
    
    # templates should have the shape (Nwave, Ntpl) for ppxf
    templates = templates.T
    
    # Total number of bins
    n_binnum = flux_binned.shape[0]
    
    # Total number of spaxels
    n_spaxels = flux.shape[0]
    
    # Total number of emission lines
    n_eml = np.sum(gas_tpl)
    
    # Total number of stellar templates
    n_temps = np.sum(~gas_tpl)

    gas_vel = np.zeros(n_binnum)
    gas_sig = np.zeros(n_binnum)
    gas_flux_binned = np.zeros([n_binnum,n_eml])
    star_weights = np.zeros([n_binnum, n_temps])
    star_weights_err = np.zeros([n_binnum, n_temps])
    bestfit_template = np.zeros([n_binnum, templates.shape[0]])
    
    # Testing use
    #beg = n_binnum if os.path.isfile(npz_file) else 0
    
    # 'bin' mode contains only one iteration, fitting the binned data
    # 'bin_to_spaxel' mode contains two iterations, fitting both the binned and
    # unbinned data
    if isinstance(mode, str) and mode == 'bin':
        beg1 = 0
        beg2 = n_spaxels
        n = n_binnum
        model_shape = flux_binned.shape
        model_binid = None
    elif isinstance(mode, str) and mode == 'bin_to_spaxel':
        beg1 = 0
        beg2 = 0
        n = n_spaxels
        model_shape = flux.shape
        model_binid = np.zeros(n_spaxels)
        id_num = 0
        # Get the index of the nearest bin for every spaxel
        nearest_bin = np.argmin((x[:, None] - xbin)**2 + (y[:, None] - ybin)**2, axis=1)
    
    # Initialize output arrays
    eml_flux = np.zeros([n,n_eml])
    eml_fluxerr = np.zeros([n,n_eml])
    eml_kin = np.zeros([n,n_eml,2])
    eml_kinerr = np.zeros([n,n_eml,2])
    weights = np.zeros([n,templates.shape[1]])
    weights_err = np.zeros([n,templates.shape[1]])
    mweights = np.zeros([n,mdegree])

    eml_sigmacorr = np.zeros([n,n_eml])

    model_flux = np.zeros(model_shape)
    model_eml_flux = np.zeros(model_shape)
    model_mask = np.zeros(model_shape)
    
    # First iteration
    for i in range(beg1, n_binnum):
        print('Fitting bin: {0}/{1}'.format(i+1,n_binnum))
        
        # Normalize spectrum to avoid numerical issues
        galaxy = flux_binned[i,:]/np.median(flux_binned[i,:])
        
        noise_binned = noise_bin[i,:]/np.median(flux_binned[i,:])
        
        # Use the first moment of Halpha line as the initial velocity
        # guess, and use the stellar velocity dispersion as initial
        # sigma guess for gas components
        #kin = [vel_gas[i], 100.]
        kin = [vel_gas[i], sig[i]]
        
        # Assign all emission lines with single component to tie their
        # kinematics (and the same for stellar templates)
        component = np.asarray([0]*n_temps + [1]*n_eml) if mode == 'bin_to_spaxel'\
        else np.append(np.zeros(n_temps,dtype=int), eml_component+1)
        
        # The negative moments automatically fix the stellar kinematics
        moments = [-2, 2] if mode == 'bin_to_spaxel'\
        else [-2] + [2]*np.max(component+1)
        
        # Tie emission lines if one selects the 'bin' mode
        tied = None if mode == 'bin_to_spaxel' else eml_tied
        
        # Initial guess for components (for stellar components, these
        # are fixed values)
        start = [[vel[i], sig[i]]]
        start += [kin for j in range(len(moments)-1)]
        
        # Tie kinematics of all emission lines while fixing the stellar kinematics in the first
        # fit
        pp = ppxf(templates, galaxy, noise_binned, velscale, start, velscale_ratio=velscale_ratio,
                  plot=False, moments=moments, degree=degree, mdegree=mdegree, tied=tied,
                  mask=mask_binned[i,:], vsyst=dv, lam=wave0, component=component, quiet=True,
                  gas_component=component>0, gas_names=gas_names)
        
        # Calculate the residuals
        resid = galaxy - pp.bestfit
        
        # Calculate the noises
        NOISE = calculate_noise(resid)
        
        # Mask out pixels with residuals > 3 sigmas
        mask = (abs(resid) < (3*NOISE)) & mask_binned[i,:]
        
        # Add a three-sigma cut in the second fit while tying the kinematics as the first one
        pp = ppxf(templates, galaxy, noise_binned, velscale, start, velscale_ratio=velscale_ratio,
                  mask=mask, plot=False, moments=moments, degree=degree, mdegree=mdegree,
                  vsyst=dv, lam=wave0, component=component, quiet=True, tied=tied,
                  gas_component=component>0, gas_names=gas_names)
        
        gas_vel[i] = pp.sol[1][0]
        gas_sig[i] = pp.sol[1][1]
        # Scale back the fluxes; note that the fluxes here are caculated by direct summing rather
        # than integration! If one uses ppxf_util to construct input templates, the ouput fluxes
        # should be corrected using pixel step near the lines
        gas_flux_binned[i] = pp.gas_flux*np.median(flux_binned[i,:])
        star_weights[i] = pp.weights[0:n_temps]
        bestfit_template[i,:] = templates[:,~gas_tpl] @ star_weights[i]
        
        design_matrix = pp.matrix/pp.noise[:, None]
        star_weights_err[i] = capfit.cov_err(design_matrix)[1][0:n_temps]
        
        # The ouput of the fitting are parsed in the following lines
        # (for 'bin' mode)
        if mode == 'bin':
            eml_kin[i][:,0] = np.full(eml_kin.shape[1], pp.sol[1][0])
            eml_kin[i][:,1] = pp.sol[:,1][eml_component+1]

            # The errors are calculated according to docstring in ppxf.py
            eml_kinerr[i][:,0] = np.full(eml_kin.shape[1], pp.error[1][0])*np.sqrt(pp.chi2)
            eml_kinerr[i][:,1] = pp.error[:,1][eml_component+1]*np.sqrt(pp.chi2)

            # Full models of flux and models of emission lines
            model_flux[i,:] = pp.bestfit*np.median(flux_binned[i,:])

            spectra = pp.matrix[:,degree+1:]
            gas_spectrum = spectra[:,pp.gas_component].dot(pp.weights[pp.gas_component])
            model_eml_flux[i,:] = gas_spectrum*np.median(flux_binned[i,:])

            # Masks for models
            model_mask[i,:] = mask

            # Calculate the sigma corrections
            FWHM = wave0/np.max(sres, axis=0)
            eml_wave_new = eml_wave*np.exp(pp.sol[1][0]/astropy.constants.c.to('km/s').value)
            index_new = np.argmin(abs(wave0-eml_wave_new[:,None]),axis=1)
            FWHM_diff2 = FWHM[index_new]**2 - FWHM[index]**2
            eml_sigmacorr[i] = FWHM_diff2/((2.355*step)**2)*(velscale**2)

            # Weights and weight errors of gas tempaltes
            # Note that the star weights errors should never be used
            weights[i] = np.append(star_weights[i], pp.weights[1:])
            weights_err[i] = np.append(star_weights_err[i], pp.gas_flux_error)

            # Coefficients of multiplicative Legendre polynomials of order > 0
            mweights[i] = pp.mpolyweights
        
        #plt.clf()
        #pp.plot()
        #plt.show()
    
#    # Save the results (only for testing)
#    np.savez_compressed('/Users/mac/Documents/New_test_results/first_'+plate+'_'+ifu+'.npz', 
#                        gas_vel = gas_vel, gas_sig = gas_sig,
#                        gas_flux_binned = gas_flux_binned,
#                        bestfit_template = bestfit_template)

    # Second iteration
    # Load data from the first iteration (only for testing)
#    first_iter = np.load('/Users/mac/Documents/New_test_results/first_'+plate+'_'+ifu+'.npz')
#    gas_vel = first_iter['gas_vel']
#    gas_sig = first_iter['gas_sig']
#    bestfit_template = first_iter['bestfit_template']
    
     # Testing use
#    if beg == 0:
#        print('writing {0}'.format(npz_file))
#        np.savez_compressed(npz_file, gas_vel = gas_vel, gas_sig = gas_sig,
#                            gas_flux_binned = gas_flux_binned, bestfit_template = bestfit_template)
#    else:
#        print('reading {0}'.format(npz_file))
#        inp = np.load(npz_file)
#        gas_vel = inp['gas_vel']
#        gas_sig = inp['gas_sig']
#        bestfit_template = inp['bestfit_template']
    
#    for i in range(1084,n_spaxels):
    for i in range(beg2, n_spaxels):
        print('Fitting spaxel: {0}/{1}'.format(i+1,n_spaxels))
        
        # Here start the first half of the second iteration
        
        mask_DRP = mask_drp[i,:]
        
        # Continue the loop if the current spaxel has no spectrum
        # TODO: decide if 1/3 is appropriate to indicate insufficient
        # pixels?
        if sum(~mask_DRP) > flux[i,:].shape[0]/3:
            model_binid[i] = -1
            continue
        else:
            model_binid[i] = id_num
            id_num += 1
        
        # Normalize spectrum to avoid numerical issues
        galaxy = flux[i,:]/np.median(flux[i,:])
                
        noise_spaxel = noise[i,:]/np.median(flux[i,:])
        noise_spaxel[(noise_spaxel<=0.)|(~np.isfinite(noise_spaxel))] = 1.
        
        # This is the index of the nearest bin to the spaxel being used
        k = nearest_bin[i]
        
        gas_templates = templates[gas_tpl]
        templates_sec = np.column_stack([bestfit_template[k], gas_templates])
        
        # Only one stellar template is used in the second iteration
        n_temps = 1
        
        # This is the gas kinematics [V, sigma] of the nearest bin which
        # will be used as the initial guess for emission-line fitting
        kin = [gas_vel[k], gas_sig[k]]
        
        # Combine the input gas components and stellar ones
        component = np.asarray([0]*n_temps + [1]*n_eml)
        
        moments = [-2, 2]
        
        # Adopt the same value in the first iteration for star component
        start = [[vel[k], sig[k]]]
        # Adopt the same starting value for all gas component
        start += [kin for j in range(len(moments)-1)]
        
        # In the third fit, the tying of kinematics is the same as the first two while the initial
        # guess for gas components, input noise and mask are different
        pp = ppxf(templates_sec, galaxy, noise_spaxel, velscale, start, velscale_ratio=velscale_ratio,
                  mask=mask_DRP, plot=False, moments=moments, degree=degree, mdegree=mdegree,
                  vsyst=dv, lam=wave0, component=component, quiet=True,
                  gas_component=component>0, gas_names=gas_names)
        
        # Calculate the residuals
        resid = galaxy - pp.bestfit
        
        # Calculate the noises
        NOISE = calculate_noise(resid)
        
        # Mask out pixels with residuals > 3 sigmas
        mask = (abs(resid) < (3*NOISE))&mask_DRP
        #print('Outliers: '+str(3109 - sum(mask)))
        
        # The forth fit adds a three-sigma cut after the third one
        pp = ppxf(templates_sec, galaxy, noise_spaxel, velscale, start, velscale_ratio=velscale_ratio,
                  mask=mask, plot=False, moments=moments, degree=degree, mdegree=mdegree,
                  vsyst=dv, lam=wave0, component=component, quiet=True,
                  gas_component=component>0, gas_names=gas_names)
        
        # Here start the second half of the second iteration
        
        # This is the gas kinematics [V, sigma] of the current spaxel from the last 
        # call of ppxf which will be used as the initial guess for emission-line fitting
        kin = [pp.sol[1][0], pp.sol[1][1]]
        
        # Combine the input gas components and stellar ones
        # Now each kinematic group has a unique component index (>0)
        component = np.append(np.zeros(n_temps,dtype=int), eml_component+1)
        
        moments = [-2] + [2]*np.max(component+1)
        
        # Adopt the same value in the first iteration for star component
        start = [[vel[k], sig[k]]]
        # Adopt the same starting value for all gas component
        start += [kin for j in range(len(moments)-1)]
        
        # In the fifth fit, the velocities and sigmas of emission lines are tied as desired
        pp = ppxf(templates_sec, galaxy, noise_spaxel, velscale, start, velscale_ratio=velscale_ratio,
                  mask=mask_DRP, quiet=True, plot=False, moments=moments, degree=degree,
                  mdegree=mdegree, vsyst=dv, lam=wave0, component=component, tied=eml_tied,
                  gas_component=component>0, gas_names=gas_names)
        
        # Calculate the residuals
        resid = galaxy - pp.bestfit
        
        # Calculate the noises
        NOISE = calculate_noise(resid)
        
        # Mask out pixels with residuals > 3 sigmas
        mask = (abs(resid) < (3*NOISE))&mask_DRP
        #print('Outliers: '+str(3109 - sum(mask)))
        
        # Add a three-sigma cut in the sixth fit
        pp = ppxf(templates_sec, galaxy, noise_spaxel, velscale, start, velscale_ratio=velscale_ratio,
                  mask=mask, quiet=True, plot=False, moments=moments, degree=degree,
                  mdegree=mdegree, vsyst=dv, lam=wave0, component=component, tied=eml_tied,
                  gas_component=component>0, gas_names=gas_names)
        
        # The ouput of the fitting are parsed in the following lines
        
        eml_kin[i][:,0] = np.full(eml_kin.shape[1], pp.sol[1][0])
        eml_kin[i][:,1] = pp.sol[:,1][eml_component+1]
        
        # The errors are calculated according to docstring in ppxf.py
        eml_kinerr[i][:,0] = np.full(eml_kin.shape[1], pp.error[1][0])*np.sqrt(pp.chi2)
        eml_kinerr[i][:,1] = pp.error[:,1][eml_component+1]*np.sqrt(pp.chi2)
        
        # Integral fluxes and their errors of all emission lines (not in the final output)
        # Reset pixel-step for each emline (since ppxf returns flux caculated by direct summing)
        eml_wave0 = eml_wave*(1+nsa_z)
        index = np.argmin(abs(wave0-eml_wave0[:,None]),axis=1)
        step = wave0[index+1] - wave0[index]
        # Remain to be examined if /(1+nsa_z) correction is needed...
        eml_flux[i] = pp.gas_flux*np.median(flux[i,:])*step
        eml_fluxerr[i] = pp.gas_flux_error*np.sqrt(pp.chi2)*np.median(flux[i,:])*step
        
        # Full models of flux and models of emission lines
        model_flux[i,:] = pp.bestfit*np.median(flux[i,:])
        
        spectra = pp.matrix[:,degree+1:]
        gas_spectrum = spectra[:,pp.gas_component].dot(pp.weights[pp.gas_component])
        model_eml_flux[i,:] = gas_spectrum*np.median(flux[i,:])
        
        # Masks for models
        model_mask[i,:] = mask
        
        # Calculate the sigma corrections
        FWHM = wave0/np.max(sres, axis=0)
        eml_wave_new = eml_wave*np.exp(pp.sol[1][0]/astropy.constants.c.to('km/s').value)
        index_new = np.argmin(abs(wave0-eml_wave_new[:,None]),axis=1)
        FWHM_diff2 = FWHM[index_new]**2 - FWHM[index]**2
        eml_sigmacorr[i] = FWHM_diff2/((2.355*step)**2)*(velscale**2)
        
        # Weights and weight errors of gas tempaltes
        # Note that the star weights errors should never be used
        weights[i] = np.append(star_weights[k]*pp.weights[0], pp.weights[1:])
        weights_err[i] = np.append(star_weights_err[k]*pp.weights[0], pp.gas_flux_error)
        
        # Coefficients of multiplicative Legendre polynomials of order > 0
        mweights[i] = pp.mpolyweights
        
        #plt.clf()
        #pp.plot()
        #plt.show()
    
#    # Save the results (only for testing)
#    np.savez_compressed('/Users/mac/Documents/New_test_results/second_new'+plate+'_'+ifu+'.npz', 
#                        model_flux = model_flux, model_eml_flux = model_eml_flux, model_mask = model_mask,
#                        model_binid = model_binid, eml_flux = eml_flux, eml_fluxerr = eml_fluxerr,
#                        eml_kin = eml_kin, eml_kinerr = eml_kinerr, eml_sigmacorr = eml_sigmacorr)
#    
#    second_iter = np.load('/Users/mac/Documents/New_test_results/second_new'+plate+'_'+ifu+'.npz')
#    model_flux = second_iter['model_flux']
#    model_eml_flux = second_iter['model_eml_flux']
#    model_mask = second_iter['model_mask']
#    model_binid = second_iter['model_binid']
#    eml_flux = second_iter['eml_flux']
#    eml_fluxerr = second_iter['eml_fluxerr']
#    eml_kin = second_iter['eml_kin']
#    eml_kinerr = second_iter['eml_kinerr']
#    eml_sigmacorr = second_iter['eml_sigmacorr']
    
    return model_flux, model_eml_flux, model_mask, model_binid, weights, weights_err, mweights,\
           eml_kin, eml_kinerr, eml_sigmacorr


