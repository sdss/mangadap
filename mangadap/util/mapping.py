# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-
"""
Defines some utility routines used to map a provided set of quantities.

----

.. include license and copyright
.. include:: ../include/copy.rst

----

.. include common links, assuming primary doc root is up one directory
.. include:: ../include/links.rst
"""
import numpy
import warnings
from scipy import ndimage

from astropy.wcs import WCS

from matplotlib import pyplot, patches
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
        x (numpy.ndarray): X coordinates of the pixels
        y (numpy.ndarray): Y coordinates of the pixels
        z (numpy.ndarray): Image values at :math:`x,y`.
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
            `pyplot.imshow`_ on the array produced by this
            routine for the :math:`x,y` ordering to be as expected!
        fill_value (float): (Optional) The default value to use for
            pixels without any data.

    Returns:
        numpy.ndarray, `numpy.ma.MaskedArray`_: The first object
        returned is a four-element array with the extent of the image
        from the bottom edge to the top edge of the x and y pixels,
        respectively.  The second object returned is the image data and
        associated (boolean) mask.
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

    # Prep the image for direct display using pyplot.imshow
    if imshow_prep:
        img = numpy.transpose(img)
        mask = numpy.transpose(mask)

    return ext, numpy.ma.MaskedArray(img, mask=mask, fill_value=fill_value)


def map_quantity(x, y, z, zmin=None, zmax=None, ncolors=64, dots=False, cmap='coolwarm',
                 colorbar=False, nticks=7, clabel=None, flux=None, fixpdf=False, xlim=None,
                 xlabel=None, ylim=None, ylabel=None, pixelscale=None, fill_value=0.0,
                 contour_levels=None, **kwargs):

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

    if cmap is None:
        cmap = 'coolwarm'

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
        cs = ax.imshow(img, interpolation='nearest', cmap=cmap, vmin=zmin,
                       vmax=zmax, extent=ext, origin='lower')
#        print('output imshow')

        if contour_levels is not None:
            nx = img.shape[0]
            dx = (ext[1]-ext[0])/nx
            ny = img.shape[1]
            dy = (ext[3]-ext[2])/ny
            ximg, yimg = numpy.meshgrid(ext[0]+dx*(numpy.arange(nx)+0.5),
                                        ext[2]+dy*(numpy.arange(ny)+0.5), indexing='xy')
            cs1 = ax.contour(ximg, yimg, img, levels=contour_levels, colors='k', linewidths=1.5)
       
    # Provide an interpolated map
    else:
        levels = numpy.linspace(zmin, zmax, ncolors) 
        cs = ax.tricontourf(x, y, z.clip(zmin, zmax), levels=levels, cmap=cmap)

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
                      cmap=cmap)

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


def map_extent(hdu, ext, offset=True):
    """
    Get the on-sky extent of a map using the provided WCS coordinates.

    Args:
        hdu (`astropy.io.fits.HDUList`_): Fits HDU list.
        ext (str): Extension of the file with the map.
        offset (bool): (**Optional**)  Return the extent as the offset
            from the coordinates of the object in arcseconds.  The WCS
            is assumed to provide RA and declination in units of
            degrees, and the header must contain the coordinates of the
            object with keywords `OBJRA` and `OBJDEC`.  If offset is
            True and these keywords do not exist, the function will
            throw a warning and proceed *without* applying any offset.
            Default is just to return the extent of the WCS coordinates.

    Returns:
        list: List of four floats: minimum and maximum x and minimum and
        maximum y.

    """
    if offset:
        try:
            objra = hdu['PRIMARY'].header['OBJRA']
            objdec = hdu['PRIMARY'].header['OBJDEC']
        except Exception as e:
            print(e)
            warnings.warn('Cannot find OBJRA and/or OBJDEC in the header of extension PRIMARY.'
                          '  No offset applied.')
            offset = False

    # Does the extension contain a cube or a single map image
    cube = len(hdu[ext].data.shape) > 2

    # Get the shape of the image
    if cube:
        nc,ny,nx = hdu[ext].data.shape
    else:
        ny,nx = hdu[ext].data.shape

    # Create a grid with the pixel coordinates
    x,y = numpy.meshgrid( numpy.arange(nx)+1, numpy.arange(ny)+1,
                          indexing='ij')

    # Reshape them into an array with X,Y coordinates
    xy = numpy.array([x.reshape(nx*ny), y.reshape(nx*ny)]).T
    if cube:
        xy = numpy.append(xy, numpy.ones(nx*ny).reshape(-1,1), axis=1)

    # Read the World Coordinate System information from the fits header
    wcs = WCS(header=hdu[ext].header, fix=False)

    # Convert the pixel coordinates into world coordinates
    XY = wcs.all_pix2world(xy, 1)
    x = XY[:,0].reshape(nx,ny)
    y = XY[:,1].reshape(nx,ny)

    # So now the x coordinate of:
    #  hdu[ext].data[:,j,i]
    # is
    #  x[j,i]
    if not offset:
        return [ numpy.amax(x), numpy.amin(x), numpy.amin(y), numpy.amax(y) ]

    # Convert the world coordinates to offset from center in arcsec
    objra = hdu['PRIMARY'].header['OBJRA']
    objdec = hdu['PRIMARY'].header['OBJDEC']
    x = (x-objra)*numpy.cos(numpy.radians(objdec))*3600.
    y = (y-objdec)*3600.

    # Return the extent
    return [numpy.amax(x), numpy.amin(x), numpy.amin(y), numpy.amax(y)]


def map_center_pixel_offset(hdu_1, hdu_2, ext, quiet=False):
    """
    Determine the offset in pixels between the object center in two the
    maps of two fits files.

    Args:
        hdu_1 (`astropy.io.fits.HDUList`_): Fits HDU list
            for the reference map.
        hdu_2 (`astropy.io.fits.HDUList`_): Fits HDU list
            for the comparison map.
        ext (str): Extension *in both files* with the maps.
        quiet (bool): (**Optional**)  Suppress terminal output.
    
    Returns:
        float: Two floats with the offset in x and y between the two maps.
    """
    cube = len(hdu_1[ext].data.shape) > 2

    objra = hdu_1['PRIMARY'].header['OBJRA']
    objdec = hdu_1['PRIMARY'].header['OBJDEC']

    if not quiet:
        print('Coordinates of 1st object: {0:10.6f} {1:9.5f}'.format(objra, objdec))
        print('Coordinates of 2nd object: {0:10.6f} {1:9.5f}'.format(
                                                                hdu_2['PRIMARY'].header['OBJRA'],
                                                                hdu_2['PRIMARY'].header['OBJDEC']))

    wcs_1 = WCS(header=hdu_1[ext].header, fix=False)
    wcs_2 = WCS(header=hdu_2[ext].header, fix=False)

    XY = numpy.zeros(3 if cube else 2).reshape(1,-1)
    XY[0,0] = objra
    XY[0,1] = objdec

    xy_1 = wcs_1.all_world2pix(XY, 1)
    xy_2 = wcs_2.all_world2pix(XY, 1)
    if not quiet:
        print('Object pixel coordinates in 1st HDU: {0:6.2f} {1:6.2f}'.format(xy_1[0,0], xy_1[0,1]))
        print('Object pixel coordinates in 2nd HDU: {0:6.2f} {1:6.2f}'.format(xy_2[0,0], xy_2[0,1]))
        print('Pixel shifts: {0} {1}'.format(xy_2[0,0]-xy_1[0,0], xy_2[0,1]-xy_1[0,1]))

    return (xy_2[0,0]-xy_1[0,0]), (xy_2[0,1]-xy_1[0,1])


def match_map_arrays(arr1, ext1, arr2, ext2, dx, dy, tol=1e-6, truncate=True, pixelscale=0.5):
    """
    Match the centers and extent of two map arrays by applying the
    previously calculated pixel offsets.  See
    :func:`map_center_pixel_offset`.

    This is just a wrapper for the two supporting functions below.

    Both maps *must* be two dimensional; i.e., a single map, not a set
    of maps arranged in different channels.

    Args:
        arr1 (`numpy.ma.MaskedArray`_): Reference map data array.
        ext1 (list): Reference map extent in arcseconds.
        arr2 (`numpy.ma.MaskedArray`_): Comparison map data array.
        ext1 (list): Comparison map extent in arcseconds.
        dx (float): **Pixel** offset in the *first* dimension (x).
        dy (float): **Pixel** offset in the *second* dimension (y).
        tol (float): (**Optional**) Tolerance for non-integer pixel
            shifts.
        truncate (float): (**Optional**) Truncate the shape of the
            shifted array to match the unshifted array.  If False, the
            unshifted array is padded to match the shifted array.
        pixelscale (float): (**Optional**) The pixel scale (extent units
            per array pixel) that *must be* common to both input arrays.

    Returns:
        numpy.ma.MaskedArray, list: The two matched arrays and the
        extent that is common to both.

    """
    # Must apply a sub-pixel shift
    if abs(dx-numpy.around(dx))>tol or abs(dy-numpy.around(dy))>tol:
        return _match_map_arrays_sub_pixel_shift(arr1, ext1, arr2, ext2, -dx, -dy,
                                                 truncate=truncate, pixelscale=pixelscale) \
                            if arr1.size > arr2.size else \
                    _match_map_arrays_sub_pixel_shift(arr2, ext2, arr1, ext1, dx, dy, swap=True,
                                                      truncate=truncate, pixelscale=pixelscale)
    # Can just apply a whole pixel shift
    return _match_map_arrays_integer_pixel_shift(arr1, ext1, arr2, ext2, -int(numpy.around(dx)),
                                                 -int(numpy.around(dy))) \
                            if arr1.size > arr2.size else \
                    _match_map_arrays_integer_pixel_shift(arr2, ext2, arr1, ext1,
                                                          int(numpy.around(dx)),
                                                          int(numpy.around(dy)), swap=True)


def _match_map_arrays_sub_pixel_shift(arr1, ext1, arr2, ext2, dx, dy, swap=False, truncate=False,
                                      pixelscale=0.5):

    # Pad the second array with an even number of pixels
    add = int(abs(dx))+1 if abs(dx) > abs(dy) else int(abs(dy))+1
    add *= 2
    if arr2.shape[0]+add - arr1.shape[0] > 0 and (arr2.shape[0]+add - arr1.shape[0]) % 2 > 1:
        add += 2
    newshape = tuple(numpy.array(arr2.shape) + add)

    # Pad and copy the array and its mask
    _arr2 = numpy.ma.masked_all(newshape, dtype=numpy.float)
    _arr2[add//2:add//2+arr2.shape[0],add//2:add//2+arr2.shape[1]] = arr2[:,:]
    _gpm = numpy.invert(numpy.ma.getmaskarray(_arr2)).astype(float)
    # Shifts are in pixels, so adjust given the new coordinates
    dx -= add//2
    dy -= add//2
    newext2 = list(numpy.array(ext2) + numpy.array([add/2, -add/2, -add/2, add/2])*pixelscale)

    # Get the input and output coordinates
    newx, newy = numpy.meshgrid(numpy.arange(newshape[0]), numpy.arange(newshape[1]), indexing='ij')
    outcoo = numpy.array([newx.ravel(), newy.ravel()])
    incoo = numpy.array([newx.ravel()-dx, newy.ravel()-dy])
    indx = ~(incoo[0,:] < 0) & ~(incoo[1,:] < 0)

    # Shift the second array and set the mask
    shifted_arr2 = numpy.ma.masked_all(newshape, dtype=numpy.float)
    shifted_arr2[outcoo[0,indx], outcoo[1,indx]] = \
        ndimage.map_coordinates(_arr2.filled(fill_value=0.0), incoo[:,indx], order=1)
    shifted_gpm = numpy.zeros(newshape, dtype=numpy.float)
    shifted_gpm[outcoo[0,indx], outcoo[1,indx]] = \
        ndimage.map_coordinates(_gpm, incoo[:,indx], order=1)
    shifted_arr2[ shifted_gpm < 0.8 ] = numpy.ma.masked

    # Return the two arrays
    sub = arr1.shape[0]-shifted_arr2.shape[0]
    if sub == 0:
        return (shifted_arr2, arr1, ext1) if swap else (arr1, shifted_arr2, ext1)
    # Pad the shifted array up to the shape of the first array
    if sub > 0:
        _shifted_arr2 = numpy.ma.masked_all(arr1.shape, dtype=numpy.float)
        _shifted_arr2[:shifted_arr2.shape[0],:shifted_arr2.shape[1]] = shifted_arr2[:,:]
        return (_shifted_arr2, arr1, ext1) if swap else (arr1, _shifted_arr2, ext1)

    sub *= -1

    # If truncation is requested, set the shifted 2nd array to the same
    # size as the first input array
    if truncate:
        _shifted_arr2 = numpy.ma.masked_all(arr1.shape, dtype=numpy.float)
        _shifted_arr2[:,:] = shifted_arr2[:arr1.shape[0],:arr1.shape[1]]
        return (_shifted_arr2, arr1, ext1) if swap else (arr1, _shifted_arr2, ext1)

    # Pad the first array to the size of the shifted array
    _arr1 = numpy.ma.masked_all(shifted_arr2.shape, dtype=numpy.float)
    _arr1[:arr1.shape[0],:arr1.shape[1]] = arr1[:,:]
    return (shifted_arr2, _arr1, newext2) if swap else (_arr1, shifted_arr2, newext2)


# Arr1 must be bigger or the same size as Arr2
def _match_map_arrays_integer_pixel_shift(arr1, ext1, arr2, ext2, dx, dy, swap=False):
    if arr1.size < arr2.size:
        raise ValueError('First array must be larger than the second.')
    if abs(dx) > (arr1.shape[0] - arr2.shape[0])/2 or abs(dy) > (arr1.shape[1] - arr2.shape[1])/2:
        raise ValueError('First array is not large enough to accommodate second array and shift.')
    _arr = numpy.ma.MaskedArray(numpy.zeros(arr1.shape, dtype=float),
                                mask=numpy.ones(arr1.shape, dtype=bool))
    _arr[dx:dx+arr2.shape[0],dy:dy+arr2.shape[1]] = arr2[:,:]
    return (_arr, arr1, ext1) if swap else (arr1, _arr, ext1) 


def map_beam_patch(extent, ax, pos=(0.1,0.1), **kwargs):
    width = extent[0]-extent[1]
    return patches.Circle(pos, 2.5/width/2, transform=ax.transAxes, **kwargs)



def permute_wcs_axes(wcs, axes):
    r"""
    Permute the axes in an `astropy.wcs.WCS`_ object.

    The number of axes must match the number axes defined by the
    `astropy.wcs.WCS`_ object. Axes are iteratively swapped until
    they're in the requested order. For example, to transpose the
    axes of a 3D WCS, set ``axes=[2,1,0]``.

    .. note::

        Method uses the :func:`deepcopy` method of the
        `astropy.wcs.WCS`_ object when altering and returning a new
        object. If the axes are already sorted, this function returns
        a deepcopy of the input ``wcs``.

    Args:
        wcs (`astropy.wcs.WCS`_):
            Object with the world-coordinate system.
        axes (array-like):
            New locations for the current axes; i.e., the axis
            currently at index ``i`` (0-indexed) is moved to be at
            index ``axis[i]``. For example, to swap the axes for a 2D
            WCS, set ``axes=[1,0]``.

    Returns:
        `astropy.wcs.WCS`_: The world-coordinate system with
        re-ordered axes.

    Raises:
        ValueError:
            Raised if the length of the ``axes`` vector is not the
            same as the number of axes defined by the
            `astropy.wcs.WCS`_ object, or if ``axes`` includes
            undefined axes (i.e., a value less than 0 or >= the
            number of axes in the WCS.)
    """
    if len(axes) != wcs.wcs.naxis:
        raise ValueError('Number of specified axes must match the WCS object.')
    if not numpy.array_equal(numpy.arange(len(axes)), numpy.sort(axes)):
        raise ValueError('Must select defined axes; 0 ... wcs.wcs.naxis-1.')
    _wcs = wcs.deepcopy()
    if numpy.array_equal(axes, numpy.sort(axes)):
        return _wcs
    _axes = numpy.atleast_1d(axes).copy()
    for i,a in enumerate(_axes):
        if a == i:
            continue
        _wcs = _wcs.swapaxes(a,i)
        _axes[i], _axes[a] = i, a
    return _wcs

