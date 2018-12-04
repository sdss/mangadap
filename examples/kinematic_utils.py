from speclite import redshift
import numpy as np
from scipy import ndimage


def deredshift(flux, wave, ivar, stellar_vel, stellar_sig, c, z): # Deredshift voronoi-binned data
    #Deredshift AND Dispersion-broaden cube. Output should be cube of deredshifted fluxes and wavelenghts
    #TODO: Does the ivar need to be broadened, too?
    velscale = np.log(wave[1]/wave[0])*c #Sould be 69 km/s
    sig_max = np.max(stellar_sig)
    #First, broaden
    FWHM_diff, sigma_corr = [],[]
    flux_smooth = np.empty((flux.shape[0], flux.shape[1]))#, flux.shape[2]))
    for i in range(0,flux.shape[0]): #for every bin
        if stellar_sig[i] <= sig_max:
            FWHM_diff = np.sqrt((((2.355*sig_max)/velscale)**2)-((2.355*stellar_sig[i])/velscale)**2)
            sigma_corr = FWHM_diff/2.355 #?z
            flux_smooth[i,:] = ndimage.filters.gaussian_filter1d(flux[i,:], sigma_corr) #may have to switch this line
        else:
            FWHM_diff = 0
            sigma_corr =0
            flux_smooth[i,:] = flux[i,:]


    z_map = (stellar_vel/c) +z #Mask should come into this somewhere.

    waves, fluxes, ivars = np.empty((flux.shape[0], flux.shape[1])),np.empty((flux.shape[0], flux.shape[1])), np.empty((flux.shape[0], flux.shape[1]))
    #waves, fluxes = np.empty((flux.shape[0], flux.shape[1])), np.empty((flux.shape[0], flux.shape[1]))
    for i in range(0, flux.shape[0]): #for every bin
        rules = [dict(name='wave', exponent=+1, array_in = wave), dict(name='flux', exponent=-1, array_in=flux_smooth[i,:]), dict(name='ivar', exponent=2, array_in=ivar[i,:])]#, dict(name='ivar', exponent=+2)]
        #rules = [dict(name='wave', exponent=+1, array_in = wave), dict(name='flux', exponent=-1, array_in=flux_smooth[i])]#, dict(name='ivar', exponent=+2)]
        result = redshift(z_in = z_map[i], z_out=0, rules=rules)
        waves[i,:] = result['wave']
        fluxes[i,:] = result['flux']
        ivars[i,:] = result['ivar']
    return waves, fluxes, ivars

def create_single_spec(fluxes, ivars, bin_disk, signal, Rb):
    #Function creates a single bulge and disk spectrum from deredshifted cube and weights by 'signal' (proxy for m_flux?)
    bulge_specs, disk_specs, bulge_ivars, disk_ivars, = np.zeros(fluxes.shape[1]), np.zeros(fluxes.shape[1]), np.zeros(fluxes.shape[1]), np.zeros(fluxes.shape[1])
    bulge_weights, disk_weights = [0], [0]
    for i in range(0, len(bin_disk)):
        if bin_disk[i] < Rb: #if the centre of the bin is within a bulge effective radius.
            bulge_specs = np.vstack((bulge_specs, fluxes[i,:]))
            bulge_ivars = np.vstack((bulge_ivars, ivars[i,:]))
            bulge_weights.append(signal[i])
        elif bin_disk[i] > (2*Rb):
            disk_specs = np.vstack((disk_specs, fluxes[i,:]))
            disk_ivars = np.vstack((disk_ivars, ivars[i,:]))
            disk_weights.append(signal[i])
    #Calculate averages
    bulge_spec = np.average(bulge_specs, axis=0, weights=bulge_weights)
    disk_spec = np.average(disk_specs, axis=0, weights=disk_weights)
    bulge_ivar = np.average(bulge_ivars, axis=0, weights = bulge_weights)  #NOTE: Is averaging the right thing to do with the IVARS?
    disk_ivar = np.average(disk_ivars, axis=0, weights = disk_weights)
    return bulge_spec, bulge_ivar, disk_spec, disk_ivar
