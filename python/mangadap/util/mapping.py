"""

Defines some utility routines used to map a provided set of quantities.
These routines are largely adapted from code provided by Michele
Cappellari.

*Source location*:
    $MANGADAP_DIR/python/mangadap/util/mapping.py

*Python2/3 compliance*::

    from __future__ import division
    from __future__ import print_function
    from __future__ import absolute_import
    from __future__ import unicode_literals
    
    import sys
    if sys.version > '3':
        long = int

*Imports*::

    import numpy
    from matplotlib import pyplot
    from mpl_toolkits.axes_grid1 import make_axes_locatable
    from matplotlib.ticker import MaxNLocator

*Revision history*:
    | **10 Jun 2015**: Pulled out functions by Michele Cappellari into a
        utility file (K. Westfall; KBW)

.. _`numpy.ma.MaskedArray`: http://docs.scipy.org/doc/numpy/reference/maskedarray.baseclass.html#numpy.ma.MaskedArray

.. _`pyplot.imshow`: http://matplotlib.org/api/pyplot_api.html#matplotlib.pyplot.imshow

"""

from __future__ import print_function
from __future__ import division
from __future__ import absolute_import
from __future__ import unicode_literals

import sys
if sys.version > '3':
    long = int

import numpy
from matplotlib import pyplot
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.ticker import MaxNLocator

##############################################################################

def masked_pixelized_image(x, y, z, pixelscale=1.0, zmin=None, zmax=None, imshow_prep=False,
                           fill_value=0.0):
    r"""

    Provided a set of pixelized data, return a masked image with a type
    `numpy.ma.MaskedArray`_.  The :math:`x` coordinates are organized
    along rows and the :math:`y` coordinates are organized along
    columns.  I.e., i.e., img[0,1] = :math:`z(x_0,y_1)`, unless
    *imshow_prep* is requested (see below).
    
    Args:
        x (numpy.array): X coordinates of the pixels
        y (numpy.array): Y coordinates of the pixels
        z (numpy.array): Image values at :math:`x,y`.

        pixelscale (float): (Optional) Pixelscale of the image in
            arcsec/pixel.

        zmin (float): (Optional) Minimum z value to include in the
            output image.  Default is to allow all pixel values.

        zmax (float): (Optional) Maximum z value to include in the
            output image.  Default is to allow all pixel values.

        imshow_prep (bool): (Optional) Prepare the matrix for use with
            `pyplot.imshow`_.  If *imshow_prep* is True, before output,
            the matrix is reordered such that increasing :math:`x`
            values are along columns and increasing :math:`y` values are
            along rows; i.e., the output is the transpose of the default
            behavior.  The appropriate call to imshow that will then put
            increasing x values along the abcissa and increasing y
            values along the ordinate is, e.g.::

                ext, img = masked_pixelized_image(x, y, z, imshow_prep=True)
                pyplot.imshow(img, interpolation='nearest', extent=ext, origin='lower')

            Note that origin **must** be "lower" when calling
            `pyplot.imshow`_ on the array produced by this routine for
            the :math:`x,y` ordering to be as expected!

        fill_value (float): (Optional) The default value to use for
            pixels without any data.

    Returns:
        numpy.ndarray: A four-element array with the extent of the image
            from the bottom edge to the top edge of the x and y pixels,
            respectively.
        numpy.ma.MaskedArray: The image data and associated (boolean)
            mask.

    """

    # Image extent
    ext = numpy.empty(4, dtype=numpy.float64)
    ext[0], ext[1] = numpy.amin(x)-pixelscale/2., numpy.amax(x)+pixelscale/2.
    ext[2], ext[3] = numpy.amin(y)-pixelscale/2., numpy.amax(y)+pixelscale/2.

    # Image size
    nx, ny = numpy.floor((ext[1]-ext[0])/pixelscale), numpy.floor((ext[3]-ext[2])/pixelscale)
#    print('MAP:  nx,ny: {0},{1}'.format(nx, ny))
#    print('MAP: extent: {0}'.format(ext))

    # Image and mask values
    img = numpy.full((nx,ny), fill_value, dtype=numpy.float64)
    i = numpy.floor((x.ravel()-ext[0])/pixelscale).astype(numpy.int)
    j = numpy.floor((y.ravel()-ext[2])/pixelscale).astype(numpy.int)
    i[i == nx] = nx-1
    j[j == ny] = ny-1
#    print(i)
#    print(j)
    img[i,j] = z
    mask = numpy.full((nx,ny), True, dtype=bool)
    mask[i,j] = False

    # Mask values
    if zmin is not None:
        mask[(img < zmin)] = True
    if zmax is not None:
        mask[(img > zmax)] = True

#    # Flip the x axis (sky-right coordinates?)
#    if xflip:
#        img = numpy.fliplr(img)
#        mask = numpy.fliplr(mask)
#        ext[0], ext[1] = ext[1], ext[0]
#        print('EXTENT: {0}'.format(ext))

    # Prep the image for direct display using matplotlib.pyplot.imshow
    if imshow_prep:
        img = numpy.transpose(img)
        mask = numpy.transpose(mask)

    return ext, numpy.ma.masked_array(img, mask, fill_value=fill_value)


def map_quantity(x, y, z, zmin=None, zmax=None, ncolors=64, dots=False, cmap='coolwarm',
                 colorbar=False, nticks=7, clabel=None, flux=None, fixpdf=False, xlim=None,
                 xlabel=None, ylim=None, ylabel=None, pixelscale=None, fill_value=0.0, **kwargs):

    """
    Copyright (C) 2013-2014, Michele Cappellari
    E-mail: cappellari_at_astro.ox.ac.uk
    http://purl.org/cappellari/software
    
    *Revision history*:
        | V1.0: Michele Cappellari, Paranal, 11 November 2013
        | V1.0.1: Clip values before contouring. MC, Oxford, 26 February
            2014
        | V1.0.2: Include SAURON colormap. MC, Oxford, 29 January 2014
        | V1.0.3: Call set_aspect(1). MC, Oxford, 22 February 2014
        | V1.0.4: Call autoscale_view(tight=True). Overplot small dots
            by default.  MC, Oxford, 25 February 2014
        | V1.0.5: Use axis('image'). MC, Oxford, 29 March 2014
        | V1.0.6: Allow changing colormap. MC, Oxford, 29 July 2014
        | V1.0.7: Added optional fixpdf keyword to remove PDF artifacts
            like `stackoverflow_15822159`_ .  Make nice tick levels for
            colorbar. Added nticks keyword for colorbar.  MC, Oxford, 16
            October 2014

    .. _stackoverflow_15822159: http://stackoverflow.com/questions/15822159/aliasing-when-saving-matplotlib-filled-contour-plot-to-pdf-or-eps
 
    """
    if zmin is None:
        zmin = numpy.min(z)
    if zmax is None:
        zmax = numpy.max(z)

#    print('raveling')

    x, y, z = map(numpy.ravel, [x, y, z])
#    print('getting axes')
    ax = pyplot.gca()

    # Provide a pixelized map
    if pixelscale is not None:
#        print('pixelizing')
        ext, img = masked_pixelized_image(x, y, z, pixelscale, zmin=zmin, zmax=zmax,
                                          fill_value=fill_value, imshow_prep=True)

#        print('made pixelized image')
        cs = ax.imshow(img, interpolation='nearest', cmap='coolwarm', vmin=zmin,
                       vmax=zmax, extent=ext, origin='lower')
#        print('output imshow')
       
    # Provide an interpolated map
    else:
        levels = numpy.linspace(zmin, zmax, ncolors) 
        cs = ax.tricontourf(x, y, z.clip(zmin, zmax), levels=levels, cmap='coolwarm')

    ax.axis('image')  # Equal axes and no rescaling
    ax.minorticks_on()
    ax.tick_params(length=10, which='major')
    ax.tick_params(length=5, which='minor')
#    print('axes limits: {0}'.format(ax.get_xlim()))
    if xlim is not None:
        ax.set_xlim(xlim)
    if ylim is not None:
        ax.set_ylim(ylim)
    if xlabel is not None:
        ax.set_xlabel(xlabel)
    if ylabel is not None:
        ax.set_ylabel(ylabel)

    if flux is not None:
        ax.tricontour(x, y, -2.5*numpy.log10(flux/numpy.max(flux).ravel()),
                      levels=numpy.arange(20), colors='k') # 1 mag contours

    if fixpdf:  # remove white contour lines in PDF at expense of larger file size
        ax.tricontour(x, y, z.clip(zmin, zmax), levels=levels, zorder=0,
                      cmap='coolwarm')

    if dots:
        ax.plot(x, y, '.k', markersize=kwargs.get("markersize", 3))

    if colorbar:
        divider = make_axes_locatable(ax)
        cax = divider.append_axes("right", size="5%", pad=0.1)
        ticks = MaxNLocator(nticks).tick_values(zmin, zmax)
        cbar = pyplot.colorbar(cs, cax=cax, ticks=ticks)
        if clabel:
            cbar.set_label(clabel)

#    print('done map')

    return cs

