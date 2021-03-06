import os
import time
import argparse

import numpy

from astropy.io import fits

from matplotlib import pyplot
from matplotlib.ticker import NullFormatter
from matplotlib.widgets import Button, RadioButtons, AxesWidget

from mangadap.dapfits import DAPCubeBitMask
from mangadap.par.analysisplan import AnalysisPlanSet
from mangadap.util.bitmask import BitMask
from mangadap.config import defaults
from mangadap.proc.spatiallybinnedspectra import SpatiallyBinnedSpectra
from mangadap.proc.stellarcontinuummodel import StellarContinuumModel
from mangadap.proc.emissionlinemodel import EmissionLineModel

#-----------------------------------------------------------------------
# Slider class
class Slider(AxesWidget):
    """
    Copied from matplotlib.  Just wanted no val to be plotted because I
    wanted to defined it explicitly.

    A slider representing a floating point range.

    For the slider to remain responsive you must maintain a
    reference to it.

    The following attributes are defined
      *ax*        : the slider :class:`matplotlib.axes.Axes` instance

      *val*       : the current slider value

      *vline*     : a :class:`matplotlib.lines.Line2D` instance
                     representing the initial value of the slider

      *poly*      : A :class:`matplotlib.patches.Polygon` instance
                     which is the slider knob

      *valfmt*    : the format string for formatting the slider text

      *label*     : a :class:`matplotlib.text.Text` instance
                     for the slider label

      *closedmin* : whether the slider is closed on the minimum

      *closedmax* : whether the slider is closed on the maximum

      *slidermin* : another slider - if not *None*, this slider must be
                     greater than *slidermin*

      *slidermax* : another slider - if not *None*, this slider must be
                     less than *slidermax*

      *dragging*  : allow for mouse dragging on slider

    Call :meth:`on_changed` to connect to the slider event
    """
    def __init__(self, ax, label, valmin, valmax, valinit=0.5, valfmt='%1.2f',
                 closedmin=True, closedmax=True, slidermin=None,
                 slidermax=None, dragging=True, **kwargs):
        """
        Create a slider from *valmin* to *valmax* in axes *ax*.

        Additional kwargs are passed on to ``self.poly`` which is the
        :class:`matplotlib.patches.Rectangle` that draws the slider
        knob.  See the :class:`matplotlib.patches.Rectangle` documentation for
        valid property names (e.g., *facecolor*, *edgecolor*, *alpha*, ...).

        Parameters
        ----------
        ax : Axes
            The Axes to put the slider in

        label : str
            Slider label

        valmin : float
            The minimum value of the slider

        valmax : float
            The maximum value of the slider

        valinit : float
            The slider initial position

        label : str
            The slider label

        valfmt : str
            Used to format the slider value, fprint format string

        closedmin : bool
            Indicate whether the slider interval is closed on the bottom

        closedmax : bool
            Indicate whether the slider interval is closed on the top

        slidermin : Slider or None
            Do not allow the current slider to have a value less than
            `slidermin`

        slidermax : Slider or None
            Do not allow the current slider to have a value greater than
            `slidermax`


        dragging : bool
            if the slider can be dragged by the mouse

        """
        AxesWidget.__init__(self, ax)

        self.valmin = valmin
        self.valmax = valmax
        self.val = valinit
        self.valinit = valinit
        self.poly = ax.axvspan(valmin, valinit, 0, 1, **kwargs)

        self.vline = ax.axvline(valinit, 0, 1, color='r', lw=1)

        self.valfmt = valfmt
        ax.set_yticks([])
        ax.set_xlim((valmin, valmax))
        ax.set_xticks([])
        ax.set_navigate(False)

        self.connect_event('button_press_event', self._update)
        self.connect_event('button_release_event', self._update)
        if dragging:
            self.connect_event('motion_notify_event', self._update)
        self.label = ax.text(-0.02, 0.5, label, transform=ax.transAxes,
                             verticalalignment='center',
                             horizontalalignment='right')

        if valfmt is not None:
            self.valtext = ax.text(1.02, 0.5, valfmt % valinit, transform=ax.transAxes,
                                   verticalalignment='center', horizontalalignment='left')

        self.cnt = 0
        self.observers = {}

        self.closedmin = closedmin
        self.closedmax = closedmax
        self.slidermin = slidermin
        self.slidermax = slidermax
        self.drag_active = False

    def _update(self, event):
        """update the slider position"""
        if self.ignore(event):
            return

        if event.button != 1:
            return

        if event.name == 'button_press_event' and event.inaxes == self.ax:
            self.drag_active = True
            event.canvas.grab_mouse(self.ax)

        if not self.drag_active:
            return

        elif ((event.name == 'button_release_event') or
              (event.name == 'button_press_event' and
               event.inaxes != self.ax)):
            self.drag_active = False
            event.canvas.release_mouse(self.ax)
            return

        val = event.xdata
        if val <= self.valmin:
            if not self.closedmin:
                return
            val = self.valmin
        elif val >= self.valmax:
            if not self.closedmax:
                return
            val = self.valmax

        if self.slidermin is not None and val <= self.slidermin.val:
            if not self.closedmin:
                return
            val = self.slidermin.val

        if self.slidermax is not None and val >= self.slidermax.val:
            if not self.closedmax:
                return
            val = self.slidermax.val

        self.set_val(val)


    def set_val(self, val):
        xy = self.poly.xy
        xy[2] = val, 1
        xy[3] = val, 0
        self.poly.xy = xy
        if self.valfmt is not None:
            self.valtext.set_text(self.valfmt % val)
        if self.drawon:
            self.ax.figure.canvas.draw_idle()
        self.val = val
        if not self.eventson:
            return
        for cid, func in self.observers.items():
            func(val)

    def on_changed(self, func):
        """
        When the slider value is changed, call *func* with the new
        slider position

        A connection id is returned which can be used to disconnect
        """
        cid = self.cnt
        self.observers[cid] = func
        self.cnt += 1
        return cid

    def disconnect(self, cid):
        """remove the observer with connection id *cid*"""
        try:
            del self.observers[cid]
        except KeyError:
            pass

    def reset(self):
        """reset the slider to the initial value if needed"""
        if (self.val != self.valinit):
            self.set_val(self.valinit)


#-----------------------------------------------------------------------
# Vertical Slider class
class VertSlider(AxesWidget):
    """
    Pulled from: http://stackoverflow.com/questions/25934279/add-a-vertical-slider-with-matplotlib

    A slider representing a floating point range

    The following attributes are defined
      *ax*        : the slider :class:`matplotlib.axes.Axes` instance

      *val*       : the current slider value

      *vline*     : a :class:`matplotlib.lines.Line2D` instance
                     representing the initial value of the slider

      *poly*      : A :class:`matplotlib.patches.Polygon` instance
                     which is the slider knob

      *valfmt*    : the format string for formatting the slider text

      *label*     : a :class:`matplotlib.text.Text` instance
                     for the slider label

      *closedmin* : whether the slider is closed on the minimum

      *closedmax* : whether the slider is closed on the maximum

      *slidermin* : another slider - if not *None*, this slider must be
                     greater than *slidermin*

      *slidermax* : another slider - if not *None*, this slider must be
                     less than *slidermax*

      *dragging*  : allow for mouse dragging on slider

    Call :meth:`on_changed` to connect to the slider event
    """
    def __init__(self, ax, label, valmin, valmax, valinit=0.5, valfmt='%1.2f',
                 closedmin=True, closedmax=True, slidermin=None,
                 slidermax=None, dragging=True, **kwargs):
        """
        Create a slider from *valmin* to *valmax* in axes *ax*

        *valinit*
            The slider initial position

        *label*
            The slider label

        *valfmt*
            Used to format the slider value

        *closedmin* and *closedmax*
            Indicate whether the slider interval is closed

        *slidermin* and *slidermax*
            Used to constrain the value of this slider to the values
            of other sliders.

        additional kwargs are passed on to ``self.poly`` which is the
        :class:`matplotlib.patches.Rectangle` which draws the slider
        knob.  See the :class:`matplotlib.patches.Rectangle` documentation
        valid property names (e.g., *facecolor*, *edgecolor*, *alpha*, ...)
        """
        AxesWidget.__init__(self, ax)

        self.valmin = valmin
        self.valmax = valmax
        self.val = valinit
        self.valinit = valinit
        self.poly = ax.axhspan(valmin, valinit, 0, 1, **kwargs)

        self.vline = ax.axhline(valinit, 0, 1, color='r', lw=1)

        self.valfmt = valfmt
        ax.set_xticks([])
        ax.set_ylim((valmin, valmax))
        ax.set_yticks([])
        ax.set_navigate(False)

        self.connect_event('button_press_event', self._update)
        self.connect_event('button_release_event', self._update)
        if dragging:
            self.connect_event('motion_notify_event', self._update)
        self.label = ax.text(0.5, 1.03, label, transform=ax.transAxes,
                             verticalalignment='center',
                             horizontalalignment='center')

        if self.valfmt is not None:
            self.valtext = ax.text(0.5, -0.03, valfmt % valinit, transform=ax.transAxes,
                                   verticalalignment='center', horizontalalignment='center')

        self.cnt = 0
        self.observers = {}

        self.closedmin = closedmin
        self.closedmax = closedmax
        self.slidermin = slidermin
        self.slidermax = slidermax
        self.drag_active = False

    def _update(self, event):
        """update the slider position"""
        if self.ignore(event):
            return

        if event.button != 1:
            return

        if event.name == 'button_press_event' and event.inaxes == self.ax:
            self.drag_active = True
            event.canvas.grab_mouse(self.ax)

        if not self.drag_active:
            return

        elif ((event.name == 'button_release_event') or
              (event.name == 'button_press_event' and
               event.inaxes != self.ax)):
            self.drag_active = False
            event.canvas.release_mouse(self.ax)
            return

        val = event.ydata
        if val <= self.valmin:
            if not self.closedmin:
                return
            val = self.valmin
        elif val >= self.valmax:
            if not self.closedmax:
                return
            val = self.valmax

        if self.slidermin is not None and val <= self.slidermin.val:
            if not self.closedmin:
                return
            val = self.slidermin.val

        if self.slidermax is not None and val >= self.slidermax.val:
            if not self.closedmax:
                return
            val = self.slidermax.val

        self.set_val(val)

    def set_val(self, val):
        xy = self.poly.xy
        xy[1] = 0, val
        xy[2] = 1, val
        self.poly.xy = xy
        if self.valfmt is not None:
            self.valtext.set_text(self.valfmt % val)
        if self.drawon:
            self.ax.figure.canvas.draw()
        self.val = val
        if not self.eventson:
            return
        for cid, func in self.observers.items():
            func(val)

    def on_changed(self, func):
        """
        When the slider value is changed, call *func* with the new
        slider position

        A connection id is returned which can be used to disconnect
        """
        cid = self.cnt
        self.observers[cid] = func
        self.cnt += 1
        return cid

    def disconnect(self, cid):
        """remove the observer with connection id *cid*"""
        try:
            del self.observers[cid]
        except KeyError:
            pass

    def reset(self):
        """reset the slider to the initial value if needed"""
        if (self.val != self.valinit):
            self.set_val(self.valinit)


#-----------------------------------------------------------------------------

class MaNGATargetBitMask(BitMask):
    r"""
    Structure with the DRP mask bits.
    """
    def __init__(self, sdss_maskbits=None, trgt='MANGA_TARGET1'):
        _sdss_maskbits = defaults.sdss_maskbits_file() if sdss_maskbits is None else sdss_maskbits
        tmp = BitMask.from_par_file(_sdss_maskbits, trgt)
        keys, descr = tmp._init_objs()
        super(MaNGATargetBitMask, self).__init__(keys, descr=descr)


def channel_bitmask():
    return BitMask(['ERRORS', 'MASKED'])


def pull_maps(hdu, ext, toerr=False):

    # Copy the map data
    maps = hdu[ext].data.copy().T
    if len(maps.shape) not in [ 2, 3 ]:
        raise ValueError('WTF')

    bits = numpy.zeros(maps.shape[-1] if len(maps.shape) == 3 else 1, dtype=numpy.uint8)
    bm = channel_bitmask()

    try:
        err_ext = hdu[ext].header['ERRDATA']
        maps_err = numpy.ma.power(hdu[err_ext].data.copy().T, -0.5) if toerr else \
                        hdu[err_ext].data.copy().T
        bits[:] = bm.turn_on(bits[:], 'ERRORS')
    except:
        maps_err = numpy.zeros(maps.shape, dtype=numpy.float)

    try:
        msk_ext = hdu[ext].header['QUALDATA']
        maps_msk = hdu[msk_ext].data.copy().T
        bits[:] = bm.turn_on(bits[:], 'MASKED')
    except:
        maps_msk = numpy.zeros(maps.shape, dtype=numpy.uint8)

    if len(maps.shape) == 2:
        return numpy.array([ ext ]), bits, maps.reshape(*maps.shape,1), \
                    maps_err.reshape(*maps.shape,1), maps_msk.reshape(*maps.shape,1)

    nmaps = maps.shape[2]
    names = numpy.full(nmaps, ext, dtype=str)
    ndig = int(numpy.log10(nmaps))+1

    names = numpy.array([ '{0}:{1}'.format(ext, hdu[ext].header['C'+'{0}'.format(i+1).zfill(ndig)])
                            for i in range(nmaps) ])

    return names, bits, maps, maps_err, maps_msk


def build_maps(hdu, ext=['BIN_MFLUX', 'BIN_SNR', 'STELLAR_VEL', 'STELLAR_SIGMA']):

    # TODO: Deal with the transpose

    # Always get the bin ID extension
    binid = hdu['BINID'].data.copy().T

    # Build the maps cube
    _ext = numpy.array([ 'BIN_MFLUX', 'BIN_SNR', 'STELLAR_VEL', 'STELLAR_SIGMA']) \
                if ext is None else numpy.asarray(ext)

    print('Building maps from extensions: {0}'.format(_ext))

    ndim = numpy.array([ len(hdu[e].data.shape) for e in _ext ])
    spatial_shape = hdu[_ext[ndim==2][0]].shape if numpy.sum(ndim == 2) > 0 \
                            else hdu[_ext[ndim==3][0]].shape[1:3]

    channel_name = numpy.empty(0, dtype=numpy.str)
    channel_bits = numpy.empty(0, dtype=numpy.uint8)
    map_cube = numpy.empty((*spatial_shape,0), dtype=numpy.float)
    map_cube_ivar = numpy.empty((*spatial_shape,0), dtype=numpy.float)
    map_cube_mask = numpy.empty((*spatial_shape,0), dtype=numpy.uint64)
    names = numpy.empty(0, dtype=numpy.str)
    zmin = numpy.empty(0, dtype=numpy.float)
    zmax = numpy.empty(0, dtype=numpy.float)
    for i in range(len(_ext)):
        _names, _bits, _maps, _maps_ivar, _maps_mask = pull_maps(hdu, _ext[i])
        channel_name = numpy.append( channel_name, _names )
        channel_bits = numpy.append( channel_bits, _bits )
        map_cube = numpy.append( map_cube, _maps, axis=2 )
        map_cube_ivar = numpy.append( map_cube_ivar, _maps_ivar, axis=2 )
        map_cube_mask = numpy.append( map_cube_mask, _maps_mask, axis=2 )
        _maps = numpy.ma.MaskedArray(_maps, mask=_maps_mask>0) \
                    if numpy.sum(_maps_mask>0) < _maps_mask.size else numpy.ma.MaskedArray(_maps)
        mins = numpy.ma.amin(_maps) if len(_maps.shape) == 2 else \
                    numpy.ma.amin(_maps.reshape(-1,_maps.shape[-1]), axis=0)
        maxs = numpy.ma.amax(_maps) if len(_maps.shape) == 2 else \
                    numpy.ma.amax(_maps.reshape(-1,_maps.shape[-1]), axis=0)
        zmin = numpy.append( zmin, mins)
        zmax = numpy.append( zmax, maxs)

    return binid, channel_name, channel_bits, map_cube, map_cube_ivar, map_cube_mask, zmin, zmax


def build_model_spectra(hdu, masked=True):

    bm = DAPCubeBitMask()
    
    wave = hdu['WAVE'].data.copy()
    if masked:
        masked_spec = numpy.ma.MaskedArray(hdu['FLUX'].data.copy().T,
                                           mask=bm.flagged(hdu['MASK'].data.T,
                                                           flag=['IGNORED', 'FORESTAR',
                                                                 'FLUXINVALID', 'IVARINVALID']))

        model_spec = numpy.ma.MaskedArray(hdu['MODEL'].data.copy().T,
                                          mask=bm.flagged(hdu['MASK'].data.T,
                                                          flag=['FITIGNORED', 'FITFAILED']))

        model_emission = numpy.ma.MaskedArray(hdu['MODEL'].data.T - hdu['STELLAR'].data.T,
                                              mask=bm.flagged(hdu['MASK'].data.T,
                                                              flag=['FITIGNORED', 'FITFAILED',
                                                                    'ELIGNORED', 'ELFAILED']))

        model_continuum = numpy.ma.MaskedArray(hdu['MODEL'].data.T - hdu['EMLINE'].data.T,
                                               mask=bm.flagged(hdu['MASK'].data.T,
                                                               flag=['FITIGNORED', 'FITFAILED']))

#        model_emission = numpy.ma.MaskedArray(hdu['EMLINE'].data.T + hdu['EMLINE_BASE'].data.T,
#                                              mask=bm.flagged(hdu['MASK'].data.T,
#                                                              flag=['FITIGNORED', 'FITFAILED',
#                                                                    'ELIGNORED', 'ELFAILED']))
#
#        model_continuum = numpy.ma.MaskedArray(hdu['MODEL'].data.T - hdu['EMLINE'].data.T
#                                                    - hdu['EMLINE_BASE'].data.T,
#                                               mask=bm.flagged(hdu['MASK'].data.T,
#                                                               flag=['FITIGNORED', 'FITFAILED']))
    else:
        masked_spec = hdu['FLUX'].data.copy().T
        model_spec = hdu['MODEL'].data.copy().T
        model_emission = hdu['MODEL'].data.T - hdu['STELLAR'].data.T
        model_continuum = hdu['MODEL'].data.T - hdu['EMLINE'].data.T
#        model_emission = hdu['EMLINE'].data.T + hdu['EMLINE_BASE'].data.T
#        model_continuum = hdu['MODEL'].data.T - hdu['EMLINE'].data.T - hdu['EMLINE_BASE'].data.T

    return wave, masked_spec, model_spec, model_continuum, model_emission

                                      
                                        

#def restructure_hdu(hdu):
#    ext = numpy.array([ h.name for h in hdu ])
#    ndim = numpy.array([ 0 if h.data is None else len(h.data.shape) for h in hdu ])
#    MaNGAFits.restructure_cube(hdu, ext=ext[ndim==3])
#    MaNGAFits.restructure_map(hdu, ext=ext[ndim==2])


class ImageMarker:
    def __init__(self, s, x_slider, y_slider, pos=None):
        self.s = s
        self.x = x_slider
        self.y = y_slider
        self.reset_pos(pos)
    def __call__(self, val=None):
        self.apply_offsets()
        pyplot.draw()
    def reset_pos(self, pos=None):
        self.pos = (int(self.x.val), int(self.y.val)) if pos is None else pos
    def apply_offsets(self):
        self.reset_pos()
        self.s.set_offsets(self.pos)


class UpdateableMap:
    def __init__(self, fig, im, cbar, channel_names, channel_bits, map_cube, map_cube_ivar,
                 map_cube_mask, zmin, zmax, current=0):
        self.fig = fig
        self.im = im
        self.cbar = cbar

        # Colorbar sliders
        zmin_slider_ax = fig.add_axes([0.055, 0.1, 0.6/2.5, 0.04], facecolor='1')
        self.zmin_slider = Slider(zmin_slider_ax, 'Min', 0, 1, valinit=0.0, valfmt=None,
                                  facecolor='k')
        self.zmin_slider.on_changed(self.change_zmin)

        zmax_slider_ax = fig.add_axes([0.055, 0.14, 0.6/2.5, 0.04], facecolor='0')
        self.zmax_slider = Slider(zmax_slider_ax, 'Max', 0, 1, valinit=1.0, valfmt=None,
                                  facecolor='w')
        self.zmax_slider.on_changed(self.change_zmax)

        self.bitmask = channel_bitmask()

        self.current = current
        self.nmap = map_cube.shape[-1]
        self.channel_names = channel_names
        self.channel_bits = channel_bits
        self.map_cube = map_cube
        self.map_cube_ivar = map_cube_ivar
        self.map_cube_mask = map_cube_mask
        self.zmin = zmin
        self.zmax = zmax

        self.msk_zmin = zmin.copy()
        self.msk_zmax = zmax.copy()
        indx = self.bitmask.flagged(self.channel_bits, flag='MASKED')
        self.map_data = numpy.ma.MaskedArray(self.map_cube,
                                             mask=map_cube_mask>0).reshape(-1,len(indx))[:,indx]
#        self.msk_zmin[indx] = numpy.ma.amin(self.map_data, axis=0)
#        self.msk_zmax[indx] = numpy.ma.amax(self.map_data, axis=0)

        self.snr_zmin = zmin.copy()
        self.snr_zmax = zmax.copy()
        indx = self.bitmask.flagged(self.channel_bits, flag='ERRORS')
        self.map_data = numpy.ma.MaskedArray(self.map_cube*numpy.sqrt(self.map_cube_ivar),
                                             mask=map_cube_mask>0).reshape(-1,len(indx))[:,indx]
        self.snr_zmin[indx] = numpy.ma.amin(self.map_data, axis=0)
        self.snr_zmax[indx] = numpy.ma.amax(self.map_data, axis=0)
#        print(self.snr_zmin)
#        print(self.snr_zmax)

        self.map_data = map_cube.copy()

        self.rng = [self.msk_zmin[self.current], self.msk_zmax[self.current]]
        self.clim = self.rng.copy()
        self.min_label = zmin_slider_ax.text(1.05, 0.5, '{0:.2f}'.format(self.clim[0]))
        self.max_label = zmax_slider_ax.text(1.05, 0.5, '{0:.2f}'.format(self.clim[1]))
        
        self.zlabel = zmax_slider_ax.text(0.5, 1.4, '{0}'.format(self.channel_names[self.current]),
                                          horizontalalignment='center', verticalalignment='center',
                                          transform=zmax_slider_ax.transAxes)
        self.masked = False
        self.snr = False

        mask_button_ax = fig.add_axes([0.015+0*0.13/2.5, 0.02, 0.13/2.5, 0.05])
        self.mask_button = Button(mask_button_ax, 'Mask', color='0.7', hovercolor='0.9')
        self.mask_button.on_clicked(self.toggle_mask)

        snr_button_ax = fig.add_axes([0.015+1*0.13/2.5, 0.02, 0.13/2.5, 0.05])
        self.snr_button = Button(snr_button_ax, 'SNR', color='0.7', hovercolor='0.9')
        self.snr_button.on_clicked(self.toggle_snr)

        frst_button_ax = fig.add_axes([0.015+2*0.13/2.5+0.01, 0.02, 0.13/2.5, 0.05])
        self.frst_button = Button(frst_button_ax, 'First', color='0.7', hovercolor='0.8')
        self.frst_button.on_clicked(self.first)

        prev_button_ax = fig.add_axes([0.015+3*0.13/2.5+0.01, 0.02, 0.13/2.5, 0.05])
        self.prev_button = Button(prev_button_ax, 'Prev', color='0.7', hovercolor='0.8')
        self.prev_button.on_clicked(self.prev)

        next_button_ax = fig.add_axes([0.015+4*0.13/2.5+0.01, 0.02, 0.13/2.5, 0.05])
        self.next_button = Button(next_button_ax, 'Next', color='0.7', hovercolor='0.8')
        self.next_button.on_clicked(self.next)

        last_button_ax = fig.add_axes([0.015+5*0.13/2.5+0.01, 0.02, 0.13/2.5, 0.05])
        self.last_button = Button(last_button_ax, 'Last', color='0.7', hovercolor='0.8')
        self.last_button.on_clicked(self.last)

        pyplot.draw()


    def get_masked_data(self):
        return numpy.ma.MaskedArray(self.map_cube.copy(), mask=self.map_cube_mask>0)
    def get_snr_data(self):
        return self.map_cube*numpy.sqrt(self.map_cube_ivar)
    def change_zmin(self, val=None):
        self.clim[0] = numpy.diff(self.rng)[0]*self.zmin_slider.val+self.rng[0]
        self.min_label.set_text('{0:.2f}'.format(self.clim[0]))
        self.im.set_clim(self.clim)
        pyplot.draw()
    def change_zmax(self, val=None):
        self.clim[1] = numpy.diff(self.rng)[0]*self.zmax_slider.val+self.rng[0]
        self.max_label.set_text('{0:.2f}'.format(self.clim[1]))
        self.im.set_clim(self.clim)
        pyplot.draw()
    def toggle_mask(self, val=None):
        if self.masked:
            self.masked = False
            self.map_data = self.get_snr_data() if self.snr else self.map_cube.copy()
#            self.map_data = self.map_cube.copy()
        else:
            self.masked = True
            self.map_data = self.get_snr_data() if self.snr else self.map_cube.copy()
            self.map_data = numpy.ma.MaskedArray(self.map_data, mask=self.map_cube_mask>0)
#            self.map_data = self.map_cube_mask.copy()
        t = self.mask_button.hovercolor
        self.mask_button.hovercolor = self.mask_button.color
        self.mask_button.color = t
        self.update_map()
    def toggle_snr(self, val=None):
        if self.snr:
            self.snr = False
            self.map_data = self.get_masked_data() if self.masked else self.map_cube.copy()
        else:
            self.snr = True
            self.map_data = self.get_masked_data() if self.masked else self.map_cube.copy()
            self.map_data = self.map_data*numpy.sqrt(self.map_cube_ivar)
        t = self.snr_button.hovercolor
        self.snr_button.hovercolor = self.snr_button.color
        self.snr_button.color = t
        self.update_map()
#    def update_map_data(self):
#        self.im.set_data(self.map_data[:,:,self.current].T)
#        pyplot.draw()
    def first(self, val=None):
        self.current = 0
        self.update_map()
    def prev(self, val=None):
        self.current = self.nmap-1 if self.current == 0 else self.current - 1
        self.update_map()
    def next(self, val=None):
        self.current = 0 if self.current == self.nmap-1 else self.current + 1
        self.update_map()
    def last(self, val=None):
        self.current = self.nmap-1
        self.update_map()
    def update_map(self):
        self.im.set_data(self.map_data[:,:,self.current].T)
        self.zmin_slider.reset()
        self.zmax_slider.reset()
#        if self.masked:
#            self.rng = [self.msk_zmin[self.current], self.msk_zmax[self.current]]
        if self.snr:
            self.rng = [self.snr_zmin[self.current], self.snr_zmax[self.current]]
        else:
            self.rng = [self.msk_zmin[self.current], self.msk_zmax[self.current]]
        self.clim = self.rng.copy()
        self.im.set_clim(*self.clim)
        self.min_label.set_text('{0:.2f}'.format(self.clim[0]))
        self.max_label.set_text('{0:.2f}'.format(self.clim[1]))
        self.zlabel.set_text('{0}'.format(self.channel_names[self.current]))
        pyplot.draw()


class UpdateSpectrumPlot:
    def __init__(self, binid, masked_spec, model_continuum, model_emission, model_spec, specl,
                 contl, fitl, continuum_residl, em_fitl, residl, image_marker, spx_coo_txt,
                 binid_txt):

        self.binid = binid
        self.masked_spec = masked_spec
        self.model_continuum = model_continuum
        self.model_emission = model_emission
        self.model_spec = model_spec
        self.specl = specl
        self.contl = contl
        self.fitl = fitl
        self.cresidl = continuum_residl
        self.em_fitl = em_fitl
        self.residl = residl
        self.marker = image_marker
        self.cootxt = spx_coo_txt
        self.bintxt = binid_txt

    def __call__(self, val=None):
        self.marker.apply_offsets()
        self.specl.set_ydata(self.masked_spec[self.marker.pos[0],self.marker.pos[1],:])
        self.contl.set_ydata(self.model_continuum[self.marker.pos[0],self.marker.pos[1],:])
        self.fitl.set_ydata(self.model_spec[self.marker.pos[0],self.marker.pos[1],:])
        self.cresidl.set_ydata(self.masked_spec[self.marker.pos[0],self.marker.pos[1],:]
                                - self.model_continuum[self.marker.pos[0],self.marker.pos[1],:])
        self.em_fitl.set_ydata(self.model_emission[self.marker.pos[0],self.marker.pos[1],:])
        self.residl.set_ydata(self.masked_spec[self.marker.pos[0],self.marker.pos[1],:]
                                - self.model_spec[self.marker.pos[0],self.marker.pos[1],:])

        self.cootxt.set_text('Spaxel = {0}'.format(self.marker.pos))
        self.bintxt.set_text('Bin IDs = {0}'.format(self.binid[self.marker.pos[0], 
                                                               self.marker.pos[1],:]))

        pyplot.draw()
        

def generate_file_paths(obs, analysisplan, drpver=None, redux_path=None, dapver=None,
                        directory_path=None, analysis_path=None):

    # Generate the paths to the MAPS and LOGCUBE files:
    bin_method = SpatiallyBinnedSpectra.define_method(analysisplan['bin_key'])
    sc_method = StellarContinuumModel.define_method(analysisplan['continuum_key'])
    el_method = EmissionLineModel.define_method(analysisplan['elfit_key'])
    method = defaults.dap_method(bin_method['key'], sc_method['fitpar']['template_library_key'],
                                 el_method['continuum_tpl_key'])
    directory_path = defaults.dap_method_path(method, plate=obs['plate'],
                                              ifudesign=obs['ifudesign'], drpver=drpver,
                                              dapver=dapver, analysis_path=analysis_path) \
                                if directory_path is None else str(directory_path)

    #-------------------------------------------------------------------
    # Read the maps to plot
    maps_file = os.path.join(directory_path,
                             defaults.dap_file_name(obs['plate'], obs['ifudesign'], method,
                                                    mode='MAPS'))
    print('MAP file: {0}'.format(maps_file))
    model_file = os.path.join(directory_path,
                              defaults.dap_file_name(obs['plate'], obs['ifudesign'], method,
                                                     mode='LOGCUBE'))
    print('LOGCUBE file: {0}'.format(model_file))

    return maps_file, model_file


def parse_plt_ifu(maps_file):
    file_name = maps_file.split('/')[-1]
    return tuple([ int(n) for n in file_name.split('-')[1:3] ])


def manga_dap_inspector(maps_file, model_file, ext=None, masked_spectra=True):

    # Get the plate and ifu numbers
    plt, ifu = parse_plt_ifu(maps_file)

    #-------------------------------------------------------------------
    # Read the maps to plot
    print('Reading MAP file: {0}'.format(maps_file))
    hdu = fits.open(maps_file)
    # Get the map data
    binid, channel_names, channel_bits, map_cube, map_cube_ivar, map_cube_mask, \
            map_zmin, map_zmax = build_maps(hdu, ext=ext)

    # Get the MaNGA ID
    mangaid = hdu['PRIMARY'].header['MANGAID']

    # Determine the sample flags
    trg_flags = MaNGATargetBitMask().flagged_bits(hdu['PRIMARY'].header['MNGTARG1'])

    target_type='None'
    if 'SECONDARY_v1_2_0' in trg_flags:
        target_type = 'S'
    if 'PRIMARY_v1_2_0' in trg_flags:
        target_type = 'P'
    if 'COLOR_ENHANCED_v1_2_0' in trg_flags:
        target_type = 'P+'

    # ELSE?
    hdu.close()
    del hdu

    #-------------------------------------------------------------------
    # Read the spectra and models to plot
    print('Reading MODEL file: {0}'.format(model_file))
    hdu = fits.open(model_file)
    wave, masked_spec, model_spec, model_continuum, model_emission \
                    = build_model_spectra(hdu, masked=masked_spectra)
    hdu.close()
    del hdu

    # Start the plot
    w,h = pyplot.figaspect(1)
    fig = pyplot.figure(figsize=(2.5*w,h))

    # ------------------------------------------------------------------
    # Create the channel image plot
    on_sky_xlim = numpy.array([0, map_cube.shape[0]]) - 0.5
    on_sky_ylim = numpy.array([0, map_cube.shape[1]]) - 0.5
    image_ax = fig.add_axes([0.055, 0.35, 0.6/2.5, 0.6], facecolor='0.95')
    image_ax.tick_params(axis='both', which='both', direction='in')
    image_cax = fig.add_axes([0.055+0.6/2.5+0.01, 0.35, 0.01, 0.6])
    image_ax.set_xlim(on_sky_xlim)
    image_ax.set_ylim(on_sky_ylim)
    image_ax.grid(True, which='major', color='0.8', zorder=0, linestyle='-')

    # Create the initial plot, color bar, and labels
    im = image_ax.imshow(map_cube[:,:,0].T, origin='lower', interpolation='nearest',
                         cmap='RdBu_r', vmin=map_zmin[0], vmax=map_zmax[0], zorder=3)
    cbar = pyplot.colorbar(im, cax=image_cax)
    image_ax.text(0.5, 1.03, '{0}-{1}; MaNGA ID={2}; {3}'.format(plt, ifu, mangaid, target_type),
            horizontalalignment='center', verticalalignment='center', transform=image_ax.transAxes)
    image_ax.text(-0.18, 0.5, 'Y (pix)', horizontalalignment='center', verticalalignment='center',
            transform=image_ax.transAxes, rotation='vertical')
    image_ax.text(0.5, -0.17, 'X (pix)', horizontalalignment='center', verticalalignment='center',
            transform=image_ax.transAxes)

    # Start the marker at the center of the map
    pos = tuple(*(numpy.array([map_cube.shape[0:2]])/2).astype(int))
    marker = image_ax.scatter(*pos,marker='s', s=20, color='k', alpha=0.5, zorder=4)
    # ------------------------------------------------------------------


    # ------------------------------------------------------------------
    # Add the buttons and sliders for changing the map properties
    map_selector = UpdateableMap(fig, im, cbar, channel_names, channel_bits, map_cube,
                                 map_cube_ivar, map_cube_mask, map_zmin, map_zmax)
    # ------------------------------------------------------------------


    # ------------------------------------------------------------------
    # Add the sliders for the map position
    x_slider_ax = fig.add_axes([0.055, 0.27, 0.6/2.5, 0.04], facecolor='0.9')
    x_slider = Slider(x_slider_ax, '', *on_sky_xlim, valinit=pos[0], valfmt='%d', facecolor='k',
                      alpha=0.5)
    y_slider_ax = fig.add_axes([0.02, 0.35, 0.04/2.5, 0.6], facecolor='0.9')
    y_slider = VertSlider(y_slider_ax, '', *on_sky_ylim, valinit=pos[1], valfmt='%d', facecolor='k',
                          alpha=0.5)

    coo_selector = ImageMarker(marker, x_slider, y_slider)
    # ------------------------------------------------------------------


    # ------------------------------------------------------------------
    # Add the spectrum plots
    wavelim = [ wave[0], wave[-1] ]
    fluxlim = [ numpy.ma.amin(masked_spec[pos[0],pos[1],:])-1,
                numpy.ma.amax(masked_spec[pos[0],pos[1],:])+1 ]
    continuum_residlim = [ 
                numpy.ma.amin(masked_spec[pos[0],pos[1],:]-model_continuum[pos[0],pos[1],:])-0.1,
                numpy.ma.amax(masked_spec[pos[0],pos[1],:]-model_continuum[pos[0],pos[1],:])+0.1 ]
    residlim = [ numpy.ma.amin(masked_spec[pos[0],pos[1],:]-model_spec[pos[0],pos[1],:])-0.1,
                 numpy.ma.amax(masked_spec[pos[0],pos[1],:]-model_spec[pos[0],pos[1],:])+0.1 ]

    spectrum_ax = fig.add_axes([ 0.40, 0.31, 0.58, 0.6 ], facecolor='0.95')
    spectrum_ax.set_xlim(wavelim)
    spectrum_ax.set_ylim(fluxlim)
    spectrum_ax.minorticks_on()
    spectrum_ax.tick_params(axis='both', which='both', direction='in')
    spectrum_ax.tick_params(which='major', length=6)
    spectrum_ax.tick_params(which='minor', length=3)
    spectrum_ax.xaxis.set_major_formatter(NullFormatter())
    spectrum_ax.grid(True, which='major', color='0.8', zorder=0, linestyle='-')
    spectrum_ax.grid(True, which='minor', color='0.85', zorder=0, linestyle=':')

    wave_ax = spectrum_ax.twiny()
    wave_ax.set_xlim(wavelim)
    wave_ax.minorticks_on()
    wave_ax.tick_params(axis='both', which='both', direction='in')
    wave_ax.tick_params(which='major', length=6)
    wave_ax.tick_params(which='minor', length=3)

    continuum_resid_ax = fig.add_axes([ 0.40, 0.11, 0.58, 0.20 ], facecolor='0.95',
                                        sharex=spectrum_ax)
    continuum_resid_ax.set_xlim(wavelim)
    continuum_resid_ax.set_ylim(continuum_residlim)
    continuum_resid_ax.minorticks_on()
    continuum_resid_ax.tick_params(axis='both', which='both', direction='in')
    continuum_resid_ax.tick_params(which='major', length=6)
    continuum_resid_ax.tick_params(which='minor', length=3)
    continuum_resid_ax.grid(True, which='major', color='0.8', zorder=0, linestyle='-')
    continuum_resid_ax.grid(True, which='minor', color='0.85', zorder=0, linestyle=':')

    
    resid_ax = fig.add_axes([ 0.40, 0.01, 0.58, 0.10 ], facecolor='0.95', sharex=spectrum_ax)
    resid_ax.set_xlim(wavelim)
    resid_ax.set_ylim(residlim)
    resid_ax.minorticks_on()
    resid_ax.tick_params(axis='both', which='both', direction='in')
    resid_ax.tick_params(which='major', length=6)
    resid_ax.tick_params(which='minor', length=3)
    resid_ax.grid(True, which='major', color='0.8', zorder=0, linestyle='-')
    resid_ax.grid(True, which='minor', color='0.85', zorder=0, linestyle=':')

    
    spec, = spectrum_ax.step(wave, masked_spec[pos[0],pos[1],:],
                             where='mid', lw=0.5, color='k', zorder=3)
    cont, = spectrum_ax.plot(wave, model_continuum[pos[0],pos[1],:],
                             lw=0.5, color='b', zorder=4)
    fit, = spectrum_ax.plot(wave, model_spec[pos[0],pos[1],:], lw=0.5,
                            color='g', zorder=5)

    continuum_resid, = continuum_resid_ax.step(wave,
                           masked_spec[pos[0],pos[1],:] - model_continuum[pos[0],pos[1],:],
                           where='mid', lw=0.5, color='k', zorder=3)
    em_fit, = continuum_resid_ax.plot(wave, model_emission[pos[0],pos[1],:],
                                      lw=0.5, color='r', zorder=4)

    resid, = resid_ax.step(wave, masked_spec[pos[0],pos[1],:] - model_spec[pos[0],pos[1],:],
                           where='mid', lw=0.5, color='k', zorder=3)

    spectrum_ax.text(0.5, 1.1, r'Wavelength ($\AA$)', horizontalalignment='center',
                     verticalalignment='center', transform=spectrum_ax.transAxes)
    spectrum_ax.text(1.02, 0.5, 'Flux', horizontalalignment='center', verticalalignment='center',
                     transform=spectrum_ax.transAxes, rotation='vertical')
    continuum_resid_ax.text(1.02, 0.5, r'Cont. $\Delta$', horizontalalignment='center',
                            verticalalignment='center', transform=continuum_resid_ax.transAxes,
                            rotation='vertical')
    resid_ax.text(1.02, 0.5, r'$\Delta$', horizontalalignment='center',
                  verticalalignment='center', transform=resid_ax.transAxes, rotation='vertical')

    spx_coo_txt = spectrum_ax.text(0.98, 0.92, 'Spaxel = {0}'.format(coo_selector.pos),
                                   horizontalalignment='right', verticalalignment='center',
                                   transform=spectrum_ax.transAxes, fontsize=8)
    binid_txt = spectrum_ax.text(0.98, 0.87, 'Bin IDs = {0}'.format(
                                        binid[coo_selector.pos[0],coo_selector.pos[1],:]),
                                 horizontalalignment='right', verticalalignment='center',
                                 transform=spectrum_ax.transAxes, fontsize=8)

    spectrum_selector = UpdateSpectrumPlot(binid, masked_spec, model_continuum, model_emission,
                                           model_spec, spec, cont, fit, continuum_resid, em_fit,
                                           resid, coo_selector, spx_coo_txt, binid_txt)
    x_slider.on_changed(spectrum_selector)
    y_slider.on_changed(spectrum_selector)
    # ------------------------------------------------------------------

    pyplot.show()

def parse_args(options=None, return_parser=False):
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('--maps_file', type=str, default=None, help='Output MAPS file')
    parser.add_argument('--model_file', type=str, default=None, help='Output LOGCUBE model file')
#    parser.add_argument('--drp_file', type=str, default=None, help='DRP LOGCUBE file')

    parser.add_argument('--ext', type=str, nargs='*', help='List of map extensions to include',
                        default=None)

    parser.add_argument('--maskspec', action='store_true', default=False,
                        help='Use masks for spectra.')

    # - Add "from scratch" keyword that runs main DAP classes instead of
    # just reading MAPS and LOGCUBE files

    # - Use Marvin instead of just reading fits files; specify auto
    # mode...

    if return_parser:
        return parser

    return parser.parse_args() if options is None else parser.parse_args(options)


def main(args):

    # Can it proceed?
    maps_file, model_file = args.maps_file, args.model_file
    if maps_file is None or model_file is None:
        raise ValueError('Must provide both a MAPS file and a model LOGCUBE file.')

    manga_dap_inspector(maps_file, model_file, ext=args.ext, masked_spectra=args.maskspec)
