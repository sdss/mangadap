# Force Python 3 behavior
from __future__ import division
from __future__ import print_function
from __future__ import absolute_import

# Check long/int definition
import sys
if sys.version > '3':
    long = int

# Imports
import numpy
import os.path
from matplotlib import cm
from matplotlib import mlab
from matplotlib import pyplot
from matplotlib import image
from matplotlib import ticker.NullFormatter

from scipy.interpolate import interp1d

# DAP imports
from mangadap.util.parser import arginp_to_list
from mangadap.util.exception_tools import print_frame
from mangadap.drpfile import drpfile
from mangadap.drpcomplete import drpcomplete
from mangadap.mangampl import mangampl
from mangadap.dapfile import dapfile

class qaplot:
        """
        Class used to produce quality assurance plots upon the
        completion of the DAP.
        """

        def __init__(self, dapver, plate, ifudesign, mode, bintype, niter, mplver=None):

            """
            Initializes the qaplot object.  Currently, only allowed to
            use with a single DRP/DAP file combination.

            TODO: Allow console input?

            TODO: Provide bintype and niter, or just do all available files?

            ARGUMENTS:
                -    dapver: DAP version used
                -     plate: specified plate
                - ifudesign: specified ifudesign
                -      mode: specified mode ('CUBE' or 'RSS')
                -   bintype: binning type uses (defines DAP file name)
                -     niter: plan iteration number (defines DAP file name)
                -    mplver: MPL version used, also sets the DRP
                             version; defaults to 'MPL-2'
            """

            # Set the mplversion
            try:
                self.mpl = mangampl(mplver)
            except Exception as e:
                print_frame('Exception')
                raise Exception('Undefined MPL:'+e)

            # Set the drpfile
            self.drpf = drpfile(self.mpl.drpver, plate, ifudesign, mode)

            # Set the dapfile and check that it exists
            try:
                self.dapf = dapfile(dapver, plate, ifudesign, mode, bintype, niter)
            except Exception as e:
                print_frame('Exception')
                raise Exception('Undefined DAP file:'+e)

            # Set the drpcomplete file
            self.drpc = drpcomplete(dapver=self.dapf.dapver, drpver=self.mpl.drpver)

            # Check it has the data for this plate/ifudesign
            try:
                self.drpc_indx = self.drpc.entry_index(self.dapf.plate, self.dapf.ifudesign)
            except Exception as e:
                print_frame('Exception')
                print(e)
                print("Ignoring DRPComplete data.")
                self.drpc_indx = None


        def DRP_inputs(self, ofile=None):
            """
            Produce a plot that summarizes some of the DRP inputs.

            2 X 3 layout

            tl: finding chart (if available)
            tm: gri composite from DRP file (if available; CUBE only)
            tr: "white light" integration of DRP data (if available; good pixels only)

            ml: fraction of good pixels per spectrum (fiducial positions for RSS)
            mm: good spectrum mask (fiducial positions for RSS)
            mr: growth of good pixel fraction; mark good fraction
                threshold, resulting percentage of used spectra &
                percentage of used good pixels.

            """

            # Create the data to plot
            # TODO: All images are expected to be square

            # Finding chart png file
            fc_file = self.drpf.finding_chart_path()
            if not os.path.exists(fc):
                fc_file = None

            # gri DRP image data
            # Returns None, None, None if DRP file is an RSS image
            gri_x, gri_y, gri_z = self.drpf.gri_composite()

            # "white light" integration of DRP data
            wl_x, wl_y, wl_z = self.drpf.white_light(self.dapf.mask_drppix_flags())

            # Fraction of good pixels per spectrum
            gp_x, gp_y, gp_z = self.dapf.good_pixel_fraction()
            
            # Good spectrum mask
            gs_z = self.dapf.good_spectra()

            # Good pixel growth
            good_pixel_fraction = numpy.sort(gp_z.reshape)
            n = len(good_pixel_fraction)
            good_pixel_fraction_growth = range(1,n+1)/numpy.float64(n)

            gp_threshold = self.dapf.good_pixel_fraction_threshold()
            spectrum_fraction = \
                        (interp1d(good_pixel_fraction, good_pixel_fraction_growth))(gp_threshold)

            # Create the plot

            w,h = pyplot.figaspect(1.0)
            fig = pyplot.figure(1,figsize=(w,h))

            j

            height = width = 0.45
            left = 0.05
            bottom = 0.275

            ax_rc.xaxis.set_major_formatter(ticker.NullFormatter())
            
            ax_fc = 

ax_stink = plt.axes([left, bottom, width, height])
ax_12701 = plt.axes([left+width, bottom, width, height])

#delta = 0.5
## Pixel center coordinates
#x = y = np.arange(-3.0, 3.5, delta)
#
#print(x)
#print(y)
#
#X, Y = np.meshgrid(x, y)
#Z1 = mlab.bivariate_normal(X, Y, 1.0, 1.0, 0.0, 0.0)
#Z2 = mlab.bivariate_normal(X, Y, 1.5, 0.5, 1, 1)
#Z = Z2-Z1  # difference of Gaussians

#im = plt.imshow(Z, interpolation='nearest', cmap=cm.RdYlGn,
#                origin='lower', vmax=abs(Z).max(), vmin=-abs(Z).max(),
#                #extent=[-3,3,-3,3])
## from the lower edge of the first pixel to the upper edge of the last pixel
#                extent=[-3.25,3.25,-3.25,3.25])
#
#plt.scatter(X, Y)


img=mpimg.imread('stinkbug.png')
lum_img = img[:,:,0]
im = ax_stink.imshow(lum_img, interpolation='nearest', cmap=cm.hot, aspect='auto')

img=mpimg.imread('12701.png')
#lum_img = img[:,:,0]
im = ax_12701.imshow(img, origin='upper', aspect='auto')


plt.show()
            




            
            


            



#            ml: signal image (from DAP file)
#            mm: noise image (from DAP file)
#            mr: S/N image (from DAP file)
#            bl-bm: S/N vs. on-sky r (From DAP file)
#            br: Table with input quantites:
#                    plate-ifudesign-mode
#                    mangaID
#                    redshift, guess dispersion
#                    ell, PA, Reff
#                    wavelength range used for S/N
#                    threshold for inclusion in S/N calculation
            

            

            

            

                



