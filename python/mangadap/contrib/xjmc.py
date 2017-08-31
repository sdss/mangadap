#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
import astropy.constants

from ppxf import ppxf

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
#   - wave: wavelength vector; shape is (Nwave,)
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
#   - stars_templates: Stellar templates; shape is (Nstartpl,
#     Ntplwave)
#   - guess_vel: Initial guess velocity for the gas components;
#     shape is (Nbins,)
#   - guess_sig: Initial guess velocity dispersion for the gas
#     components; shape is (Nbins,)
#   - gas_templates: Gas template flux; shape is (Ngastpl,
#     Ntplwave)
#   - gas_names: Name of the gas templats; shape is (Ngastpl,)
#   - template_sres: spectral resolution (R=lambda/delta
#     lambda) as a function of wavelength for all the templates
#     templates; shape is (Ntplwave,)
#   - degree: Additive polynomial order
#   - mdegree: Multiplicative polynomial order
#   - x: On-sky spaxel x coordinates; shape is (Nspaxel,)
#   - y: On-sky spaxel y coordinates; shape is (Nspaxel,)
#   - x_binned: On-sky bin x coordinate; shape is (Nbin,)
#   - y_binned: On-sky bin y coordinate; shape is (Nbin,)

def emline_fitter_with_ppxf(wave0, flux, noise, sres, flux_binned, noise_bin,
                            velscale, velscale_ratio, dv, vel, sig, stars_templates,
                            vel_gas, sig_gas, gas_templates, gas_templates_binned, gas_names, 
                            templates_sres, degree, mdegree, x, y, xbin, ybin, 
                            mask_binned, mask_drp, nsa_z):
    
    
    # First iteration
    plate = '8258'
    ifu = '6102'
    
    # stellar templates should have the shape (Nwave, Ntpl) for ppxf
    stars_templates = stars_templates.T    
    
    # Total number of bins
    n_binnum = flux_binned.shape[0]
    
    gas_vel = np.zeros(n_binnum)
    gas_sig = np.zeros(n_binnum)
    gas_flux_binned = np.zeros([n_binnum,11])
    bestfit_template = np.zeros([n_binnum, stars_templates.shape[0]])
    
    for i in range(n_binnum):
        
        # Normalize spectrum to avoid numerical issues
        galaxy = flux_binned[i,:]/np.median(flux_binned[i,:])
        
        noise_binned = noise_bin[i,:]/np.median(flux_binned[i,:])
        
        # Combines the stellar and gaseous templates into a single array
        # during the PPXF fit they will be assigned a different kinematic
        # COMPONENT value
        #
        templates = np.column_stack([stars_templates, gas_templates_binned])

        n_temps = stars_templates.shape[1]
        
        # This are the moments of Halpha [V, sigma] which will be used as
        # the initial guess for emission-line fitting
        #kin = [vel_gas[i], 100.]
        kin = [vel_gas[i], sig[i]]
        
        # Single stellar kinematic component=0 for all templates
        component = [0]*n_temps

        # 'Hdelta' 'Hgamma' 'Hbeta' 'Halpha'
        # The Balmer lines all share the same LOSVD, so they
        # are assigned the same kinematic component=1
        component += [1, 1, 1, 1]

        # ''[OII]3726' '[OII]3729' '[SII]6716' '[SII]6731' '[OIII]5007d' '[OI]6300d' '[NII]6583d'
        # Seven lines, all on the same kinematic component on the first iteration
        component += [1, 1, 1, 1, 1, 1, 1]
        component = np.asarray(component)
        
        # The negative moments automatically fix the stellar kinematics
        moments = [-2, 2]
        
        start = [[vel[i], sig[i]]]
        start += [kin for j in range(len(moments)-1)]
        
        # Using constant noise to make the first fit
        pp = ppxf(templates, galaxy, noise_binned, velscale, start, velscale_ratio=velscale_ratio,
                  plot=False, moments=moments, degree=degree, mdegree=mdegree, mask=mask_binned[i,:],
                  vsyst=dv, lam=wave0, component=component, quiet=True,
                  gas_component=component>0, gas_names=gas_names)
        
        # Calculate the residuals
        resid = galaxy - pp.bestfit
        
        # Calculate the noises
        NOISE = calculate_noise(resid)
        
        # Mask out pixels with residuals > 3 sigmas
        mask = (abs(resid) < (3*NOISE)) & mask_binned[i,:]
        
        # Do the second fit with new noise parameter and masked pixels
        pp = ppxf(templates, galaxy, noise_binned, velscale, start, velscale_ratio=velscale_ratio, mask=mask, 
                  plot=False, moments=moments, degree=degree, mdegree=mdegree,
                  vsyst=dv, lam=wave0, component=component, quiet=True,
                  gas_component=component>0, gas_names=gas_names)
        
        gas_vel[i] = pp.sol[1][0]
        gas_sig[i] = pp.sol[1][1]
        gas_flux_binned[i] = pp.gas_flux*np.median(flux_binned[i,:])
        bestfit_template[i,:] = stars_templates @ pp.weights[0:n_temps]        
        
        plt.clf()
        pp.plot()
        plt.show()
        #plt.savefig('/Users/mac/Downloads/ppxf_work/fits_first_iter_7815-6101/fits_iter_one_'+str(i).zfill(4))    
    
    # Save the results (only for testing)
    np.savez_compressed('/Users/mac/Documents/New_test_results/first_'+plate+'_'+ifu+'.npz', 
                        gas_vel = gas_vel, gas_sig = gas_sig,
                        gas_flux_binned = gas_flux_binned,
                        bestfit_template = bestfit_template)
    
     
    # Second iteration
    # Load data from the first iteration (only for testing)
    first_iter = np.load('/Users/mac/Documents/New_test_results/first_'+plate+'_'+ifu+'.npz')
    gas_vel = first_iter['gas_vel']
    gas_sig = first_iter['gas_sig']
    bestfit_template = first_iter['bestfit_template']
    
    # Get the index of the nearest bin for every spaxel
    nearest_bin = np.argmin((x[:, None] - xbin)**2 + (y[:, None] - ybin)**2, axis=1)
    
    # Total number of spaxels
    n_spaxels = nearest_bin.shape[0]
    
    eml_flux = np.zeros([n_spaxels,11])
    eml_fluxerr = np.zeros([n_spaxels,11])
    eml_kin = np.zeros([n_spaxels,11,2])
    eml_kinerr = np.zeros([n_spaxels,11,2])
    
    eml_sigmacorr = np.zeros([n_spaxels,11])
    
    model_flux = np.zeros(flux.shape)
    model_eml_flux = np.zeros(flux.shape)
    model_mask = np.zeros(flux.shape)
    model_binid = np.zeros(n_spaxels)
    id_num = 0
    
    # Need to fix here    
    for i in range(n_spaxels):
        
        # Here start the first half of the second iteration
        
        mask_DRP = mask_drp[i,:]
        
        # Continue the loop if the current spaxel has no spectra
        if sum(~mask_DRP) > flux[i,:].shape[0]/3:
            model_binid[i] = -1
            continue
        else:
            model_binid[i] = id_num
            id_num += 1
        
        galaxy = flux[i,:]/np.median(flux[i,:])   # Normalize spectrum to avoid numerical issues
                
        noise_spaxel = noise[i,:]/np.median(flux[i,:])
        noise_spaxel[(noise_spaxel<=0.)|(~np.isfinite(noise_spaxel))] = 1.
        
        # This is the index of the nearest bin to the spaxel being used
        k = nearest_bin[i]
        
        templates = np.column_stack([bestfit_template[k], gas_templates])

        n_temps = 1
        
        # This is the gas kinematics [V, sigma] of the nearest bin which
        # will be used as the initial guess for emission-line fitting
        kin = [gas_vel[k], gas_sig[k]]
        
        component = [0]*n_temps  # Single stellar kinematic component=0 for all templates

        # 'Hdelta' 'Hgamma' 'Hbeta' 'Halpha'
        # The Balmer lines all share the same LOSVD, so they
        # are assigned the same kinematic component=1
        component += [1, 1, 1, 1]

        # ''[OII]3726' '[OII]3729' '[SII]6716' '[SII]6731' '[OIII]5007d' '[OI]6300d' '[NII]6583d'
        # Seven lines, all on the same component to tie the kinematics 
        component += [1, 1, 1, 1, 1, 1, 1]
        component = np.asarray(component)
        
        moments = [-2, 2]
        
        start = [[vel[k], sig[k]]]  # adopt the same value in the first iteration for star component
        start += [kin for j in range(len(moments)-1)]  # adopt the same starting value for all gas component
        
        pp = ppxf(templates, galaxy, noise_spaxel, velscale, start, velscale_ratio=velscale_ratio, mask=mask_DRP,
                  plot=False, moments=moments, degree=degree, mdegree=mdegree,
                  vsyst=dv, lam=wave0, component=component, quiet=True,
                  gas_component=component>0, gas_names=gas_names)
        
        # Calculate the residuals
        resid = galaxy - pp.bestfit
        
        # Calculate the noises
        NOISE = calculate_noise(resid)
        
        # Mask out pixels with residuals > 3 sigmas
        mask = (abs(resid) < (3*NOISE))&mask_DRP
        #print('Outliers: '+str(3109 - sum(mask)))
        
        pp = ppxf(templates, galaxy, noise_spaxel, velscale, start, velscale_ratio=velscale_ratio, mask=mask,
                  plot=False, moments=moments, degree=degree, mdegree=mdegree,
                  vsyst=dv, lam=wave0, component=component, quiet=True,
                  gas_component=component>0, gas_names=gas_names)
        
        #plt.clf()
        #pp.plot()
        #plt.show()
        
        # Here start the second half of the second iteration
        
        # This is the gas kinematics [V, sigma] of the current spaxel from the last 
        # call of ppxf which will be used as the initial guess for emission-line fitting
        kin = [pp.sol[1][0], pp.sol[1][1]]
        
        # These templates are used only when the stellar weights are set to be free
        #templates = np.column_stack([stars_templates, gas_templates])
        #n_temps = stars_templates.shape[1]
        
        component = [0]*n_temps  # Single stellar kinematic component=0 for all templates

        # 'Hdelta' 'Hgamma' 'Hbeta' 'Halpha'
        # The Balmer lines all share the same LOSVD, so they
        # are assigned the same kinematic component=1
        component += [1, 1, 1, 1]

        # ''[OII]3726' '[OII]3729' '[SII]6716' '[SII]6731' '[OIII]5007d' '[OI]6300d' '[NII]6583d'
        # Seven lines, all on different kinematic components
        # to allow for different sigmas
        component += [2, 2, 3, 3, 4, 5, 6]
        component = np.asarray(component)
        
        moments = [-2, 2, 2, 2, 2, 2, 2]
        
        start = [[vel[k], sig[k]]]  # adopt the same value in the first iteration for star component
        start += [kin for j in range(len(moments)-1)]  # adopt the same starting value for all gas component
        
        # Tie the velocity of all gas components while leaving the sigma free
        tied = [['', ''] for j in range(len(moments))]
        for j in range(2, len(moments)):
            tied[j][0] = 'p[2]'
        
        pp = ppxf(templates, galaxy, noise_spaxel, velscale, start, velscale_ratio=velscale_ratio, mask=mask_DRP,
                  quiet=True, plot=False, moments=moments, degree=degree, mdegree=mdegree,
                  vsyst=dv, lam=wave0, component=component, tied=tied,
                  gas_component=component>0, gas_names=gas_names)
        
        # Calculate the residuals
        resid = galaxy - pp.bestfit
        
        # Calculate the noises
        NOISE = calculate_noise(resid)
        
        # Mask out pixels with residuals > 3 sigmas
        mask = (abs(resid) < (3*NOISE))&mask_DRP
        #print('Outliers: '+str(3109 - sum(mask)))
        
        pp = ppxf(templates, galaxy, noise_spaxel, velscale, start, velscale_ratio=velscale_ratio, mask=mask, quiet=True,
                  plot=False, moments=moments, degree=degree, mdegree=mdegree,
                  vsyst=dv, lam=wave0, component=component, tied=tied,
                  gas_component=component>0, gas_names=gas_names)
                        
        eml_kin[i][:,0] = np.full(eml_kin.shape[1], pp.sol[1][0])
        eml_component_num = [4,2,2,1,1,1]
        emlsig = np.array([j[1] for j in pp.sol[1:7]])
        eml_kin[i][:,1] = np.repeat(emlsig, eml_component_num)
        
        eml_kinerr[i][:,0] = np.full(eml_kin.shape[1], pp.error[1][0])*np.sqrt(pp.chi2)
        emlsigerr = np.array([j[1] for j in pp.error[1:7]])*np.sqrt(pp.chi2)
        eml_kinerr[i][:,1] = np.repeat(emlsigerr, eml_component_num)
        
        # Reset pixel-step for each emline in pPXF...
        eml_wave = np.array([4102.,4341.,4861.,6563.,3726.,3729.,6716.,6731.,5007.,6300.,6583.])
        eml_wave0 = eml_wave*(1+nsa_z)
        index = np.argmin(abs(wave0-eml_wave0[:,None]),axis=1)
        step = wave0[index+1] - wave0[index]
        
        # Remain to be examined if /(1+nsa_z) correction is needed...
        eml_flux[i] = pp.gas_flux*np.median(flux[i,:])*step
        eml_fluxerr[i] = pp.gas_flux_error*np.sqrt(pp.chi2)*np.median(flux[i,:])*step
        
        model_flux[i,:] = pp.bestfit*np.median(flux[i,:])
        
        spectra = pp.matrix[:,degree+1:]
        gas_spectrum = spectra[:,pp.gas_component].dot(pp.weights[pp.gas_component])
        model_eml_flux[i,:] = gas_spectrum*np.median(flux[i,:])
        
        model_mask[i,:] = mask
        
        # Calculate the sigma corrections
        FWHM = wave0/sres[0,:]
        eml_wave_new = eml_wave*np.exp(pp.sol[1][0]/astropy.constants.c.to('km/s').value)
        index_new = np.argmin(abs(wave0-eml_wave_new[:,None]),axis=1)
        FWHM_diff2 = FWHM[index_new]**2 - FWHM[index]**2
        eml_sigmacorr[i] = FWHM_diff2/((2.355*step)**2)*(velscale**2)
        
        #plt.clf()
        #pp.plot()
        #plt.show()
        #plt.savefig('/Users/mac/Downloads/ppxf_work/fits_second_iter_7815-6101/fits_iter_two_'+str(i).zfill(4))
    
    # Save the results (only for testing)
    np.savez_compressed('/Users/mac/Documents/New_test_results/second_new'+plate+'_'+ifu+'.npz', 
                        model_flux = model_flux, model_eml_flux = model_eml_flux, model_mask = model_mask,
                        model_binid = model_binid, eml_flux = eml_flux, eml_fluxerr = eml_fluxerr,
                        eml_kin = eml_kin, eml_kinerr = eml_kinerr, eml_sigmacorr = eml_sigmacorr)
    
    '''
    second_iter = np.load('/Users/mac/Documents/New_test_results/second_new'+plate+'_'+ifu+'.npz')
    model_flux = second_iter['model_flux']
    model_eml_flux = second_iter['model_eml_flux']
    model_mask = second_iter['model_mask']
    model_binid = second_iter['model_binid']
    eml_flux = second_iter['eml_flux']
    eml_fluxerr = second_iter['eml_fluxerr']
    eml_kin = second_iter['eml_kin']
    eml_kinerr = second_iter['eml_kinerr']
    eml_sigmacorr = second_iter['eml_sigmacorr']
    '''
    
    return model_flux, model_eml_flux, model_mask, model_binid, eml_flux, eml_fluxerr, eml_kin, eml_kinerr, eml_sigmacorr