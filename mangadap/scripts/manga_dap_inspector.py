
import os
import time
import argparse

from IPython import embed

import numpy

from astropy.io import fits

from matplotlib import pyplot
from matplotlib.ticker import NullFormatter
from matplotlib.widgets import Button, AxesWidget, Cursor, RangeSlider

from mangadap.dapfits import DAPCubeBitMask
from mangadap.util.bitmask import BitMask
from mangadap.config import defaults

from mangadap.scripts import scriptbase


#-----------------------------------------------------------------------
# Pointer class
class Pointer(AxesWidget):
    """
    """
    def __init__(self, ax, pos, **kwargs):
        AxesWidget.__init__(self, ax)

        self.cursor = Cursor(ax, useblit=True, **kwargs)

        self.connect_event('button_press_event', self._update)
        self.connect_event('button_release_event', self._update)
        self.cnt = 0
        self.observers = {}
        self.drag_active = False
        self.pos = pos

    def _update(self, event):
        """update the pointer position"""
        if self.ignore(event):
            return

        if event.button != 1:
            return

        if event.name == 'button_press_event' and event.inaxes is self.ax:
            self.drag_active = True
            event.canvas.grab_mouse(self.ax)

        if not self.drag_active:
            return

        elif event.name == 'button_release_event' \
                or (event.name == 'button_press_event' and event.inaxes is not self.ax):
            self.drag_active = False
            event.canvas.release_mouse(self.ax)
            return

        self.set_val((event.xdata, event.ydata))

    def set_val(self, val):
        self.pos = val
        if not self.eventson:
            return
        for cid, func in self.observers.items():
            func(self.pos)

    def on_changed(self, func):
        """
        When the pointer position is changed, call *func* with the new
        position.

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
    maps = hdu[ext].data
    if maps.ndim not in [ 2, 3 ]:
        raise ValueError('WTF')

    bm = channel_bitmask()
    bits = numpy.zeros(maps.shape[0] if maps.ndim == 3 else 1, dtype=bm.minimum_dtype())

    try:
        err_ext = hdu[ext].header['ERRDATA']
        maps_err = numpy.ma.power(hdu[err_ext].data, -0.5) if toerr else hdu[err_ext].data
        bits = bm.turn_on(bits, 'ERRORS')
    except:
        maps_err = numpy.zeros(maps.shape, dtype=numpy.float)

    try:
        msk_ext = hdu[ext].header['QUALDATA']
        maps_msk = hdu[msk_ext].data
        bits = bm.turn_on(bits, 'MASKED')
    except:
        maps_msk = numpy.zeros(maps.shape, dtype=numpy.uint8)

    if len(maps.shape) == 2:
        return numpy.array([ ext ]), bits, maps.reshape(1,*maps.shape), \
                    maps_err.reshape(1,*maps.shape), maps_msk.reshape(1,*maps.shape)

    nmaps = maps.shape[0]
    names = numpy.full(nmaps, ext, dtype=str)
    ndig = int(numpy.log10(nmaps))+1

    names = numpy.array([ '{0}:{1}'.format(ext, hdu[ext].header['C'+'{0}'.format(i+1).zfill(ndig)])
                            for i in range(nmaps) ])

    return names, bits, maps, maps_err, maps_msk


def build_maps(hdu, ext=['BIN_MFLUX', 'BIN_SNR', 'STELLAR_VEL', 'STELLAR_SIGMA']):

    # TODO: Deal with the transpose

    # Always get the bin ID extension
    binid = hdu['BINID'].data

    # Build the maps cube
    _ext = numpy.array([ 'BIN_MFLUX', 'BIN_SNR', 'STELLAR_VEL', 'STELLAR_SIGMA']) \
                if ext is None else numpy.asarray(ext)

    print('Building maps from extensions: {0}'.format(_ext))

    ndim = numpy.array([ hdu[e].data.ndim for e in _ext ])
    spatial_shape = hdu[_ext[ndim==2][0]].shape if numpy.sum(ndim == 2) > 0 \
                            else hdu[_ext[ndim==3][0]].shape[1:3]

    channel_name = numpy.empty(0, dtype=object)
    channel_bits = numpy.empty(0, dtype=numpy.uint8)
    map_cube = numpy.empty((0,*spatial_shape), dtype=float)
    map_cube_ivar = numpy.empty((0,*spatial_shape), dtype=float)
    map_cube_mask = numpy.empty((0,*spatial_shape), dtype=numpy.uint64)
    names = numpy.empty(0, dtype=object)
    zmin = numpy.empty(0, dtype=float)
    zmax = numpy.empty(0, dtype=float)
    for i in range(len(_ext)):
        _names, _bits, _maps, _maps_ivar, _maps_mask = pull_maps(hdu, _ext[i])
        channel_name = numpy.append( channel_name, _names )
        channel_bits = numpy.append( channel_bits, _bits )
        map_cube = numpy.append(map_cube, _maps, axis=0)
        map_cube_ivar = numpy.append(map_cube_ivar, _maps_ivar, axis=0)
        map_cube_mask = numpy.append(map_cube_mask, _maps_mask, axis=0)
        _maps = numpy.ma.MaskedArray(_maps, mask=_maps_mask>0) \
                    if numpy.sum(_maps_mask>0) < _maps_mask.size else numpy.ma.MaskedArray(_maps)
        mins = numpy.ma.amin(_maps) if _maps.ndim == 2 else numpy.ma.amin(_maps, axis=(1,2))
        maxs = numpy.ma.amax(_maps) if _maps.ndim == 2 else numpy.ma.amax(_maps, axis=(1,2))
        zmin = numpy.append(zmin, mins)
        zmax = numpy.append(zmax, maxs)

    return binid, channel_name, channel_bits, map_cube, map_cube_ivar, map_cube_mask, zmin, zmax


def build_model_spectra(hdu, masked=True):

    if masked:
        bm = DAPCubeBitMask()
        masked_spec = numpy.ma.MaskedArray(hdu['FLUX'].data,
                                           mask=bm.flagged(hdu['MASK'].data,
                                                           flag=['IGNORED', 'FORESTAR',
                                                                 'FLUXINVALID', 'IVARINVALID']))

        model_spec = numpy.ma.MaskedArray(hdu['MODEL'].data,
                                          mask=bm.flagged(hdu['MASK'].data,
                                                          flag=['FITIGNORED', 'FITFAILED']))

        model_emission = numpy.ma.MaskedArray(hdu['MODEL'].data - hdu['STELLAR'].data,
                                              mask=bm.flagged(hdu['MASK'].data,
                                                              flag=['FITIGNORED', 'FITFAILED',
                                                                    'ELIGNORED', 'ELFAILED']))

        model_continuum = numpy.ma.MaskedArray(hdu['MODEL'].data - hdu['EMLINE'].data,
                                               mask=bm.flagged(hdu['MASK'].data,
                                                               flag=['FITIGNORED', 'FITFAILED']))

    else:
        masked_spec = hdu['FLUX'].data
        model_spec = hdu['MODEL'].data
        model_emission = hdu['MODEL'].data - hdu['STELLAR'].data
        model_continuum = hdu['MODEL'].data - hdu['EMLINE'].data

    return hdu['WAVE'].data, masked_spec, model_spec, model_continuum, model_emission

                                      
def select_spectrum(spec, pos):
    return spec[:,pos[1],pos[0]]


class ImageMarker:
    def __init__(self, s, pointer, pos=None):
        self.s = s
        self.pointer = pointer
        self.reset_pos(pos)
    def __call__(self, val=None):
        self.apply_offsets()
        pyplot.draw()
    def reset_pos(self, pos=None):
        self.pos = tuple([int(numpy.round(p)) for p in self.pointer.pos]) if pos is None else pos
    def apply_offsets(self):
        self.reset_pos()
        self.s.set_offsets(self.pos)


class UpdateableRangeSlider(RangeSlider):
    def __init__(self, ax, label, valmin, valmax, valinit=None, valfmt=None, closedmin=True,
                 closedmax=True, dragging=True, valstep=None, orientation="horizontal",
                 track_color='lightgrey', handle_style=None, **kwargs):

        super().__init__(ax, label, valmin, valmax, valinit=valinit, valfmt=valfmt,
                         closedmin=closedmin, closedmax=closedmax, dragging=dragging,
                         valstep=valstep, orientation=orientation, track_color=track_color,
                         handle_style=handle_style, **kwargs)
        # Save the axes reference internally
        self.ax = ax
        self.label.set_position((0.5, 1.0))
        self.label.set_verticalalignment('bottom')
        self.label.set_horizontalalignment('center')
        self.label.set_weight('bold')

    def update_range(self, rng, label):
        self._active_handle = None
        xy = self.poly.get_xy()
        xy[:,0] = numpy.roll(numpy.repeat(rng, 3)[:-1],-1)
        self.poly.set_xy(xy)
        self.valmin, self.valmax = rng
        self.valinit = numpy.array(rng)
        self._handles[0].set_xdata(numpy.array([rng[0]]))
        self._handles[1].set_xdata(numpy.array([rng[1]]))
        self.ax.set_xlim(rng)
        self.label.set_text(label)
        self.set_val(rng)


class UpdateableMap:
    def __init__(self, fig, im, cbar, channel_names, channel_bits, map_cube, map_cube_ivar,
                 map_cube_mask, zmin, zmax, current=0):
        self.fig = fig
        self.im = im
        self.cbar = cbar
        self.bitmask = channel_bitmask()

        self.current = current
        self.nmap = map_cube.shape[0]
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

        self.snr_zmin = zmin.copy()
        self.snr_zmax = zmax.copy()
        indx = self.bitmask.flagged(self.channel_bits, flag='ERRORS')
        self.map_data = numpy.ma.MaskedArray(self.map_cube*numpy.sqrt(self.map_cube_ivar),
                                             mask=map_cube_mask>0).reshape(-1,len(indx))[:,indx]
        self.snr_zmin[indx] = numpy.ma.amin(self.map_data, axis=0)
        self.snr_zmax[indx] = numpy.ma.amax(self.map_data, axis=0)

        self.map_data = map_cube.copy()
        self.rng = (self.msk_zmin[self.current], self.msk_zmax[self.current])

        # Colorbar range slider
        self.slider_ax = pyplot.axes([0.02, 0.08, 0.7/2.5, 0.04])
        self.zslider = UpdateableRangeSlider(self.slider_ax, self.channel_names[self.current],
                                             self.rng[0], self.rng[1], valinit=self.rng)
        self.zslider.on_changed(self.change_range)

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
    def change_range(self, val):
        self.im.set_clim(*val)
        pyplot.draw()
    def toggle_mask(self, val=None):
        if self.masked:
            self.masked = False
            self.map_data = self.get_snr_data() if self.snr else self.map_cube.copy()
        else:
            self.masked = True
            self.map_data = self.get_snr_data() if self.snr else self.map_cube.copy()
            self.map_data = numpy.ma.MaskedArray(self.map_data, mask=self.map_cube_mask>0)
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
        self.im.set_data(self.map_data[self.current])
        if self.snr:
            self.rng = (self.snr_zmin[self.current], self.snr_zmax[self.current])
        else:
            self.rng = (self.msk_zmin[self.current], self.msk_zmax[self.current])
        self.im.set_clim(*self.rng)
        self.zslider.update_range(self.rng, self.channel_names[self.current])
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
        _masked_spec = select_spectrum(self.masked_spec, self.marker.pos)
        _model_spec = select_spectrum(self.model_spec, self.marker.pos)
        _model_continuum = select_spectrum(self.model_continuum, self.marker.pos)
        _model_emission = select_spectrum(self.model_emission, self.marker.pos)
        self.specl.set_ydata(_masked_spec)
        self.contl.set_ydata(_model_continuum)
        self.fitl.set_ydata(_model_spec)
        self.cresidl.set_ydata(_masked_spec - _model_continuum)
        self.em_fitl.set_ydata(_model_emission)
        self.residl.set_ydata(_masked_spec - _model_spec)

        self.cootxt.set_text(f'Spaxel = {self.marker.pos}')
        self.bintxt.set_text(f'Bin IDs = {self.binid[:,self.marker.pos[1], self.marker.pos[0]]}')

        pyplot.draw()
        

def parse_plt_ifu(maps_file):
    file_name = maps_file.split('/')[-1]
    return tuple([ int(n) for n in file_name.split('-')[1:3] ])


def manga_dap_inspector(maps_file, model_file, ext=None, masked_spectra=True):

    # Get the plate and ifu numbers
    plt, ifu = parse_plt_ifu(maps_file)

    #-------------------------------------------------------------------
    # Read the maps to plot
    print('Reading MAP file: {0}'.format(maps_file))
    with fits.open(maps_file) as hdu:
        # Get the map data
        binid, channel_names, channel_bits, map_cube, map_cube_ivar, map_cube_mask, \
                map_zmin, map_zmax = build_maps(hdu, ext=ext)

        # Get the MaNGA ID
        mangaid = hdu['PRIMARY'].header['MANGAID'] if 'MANGAID' in hdu['PRIMARY'].header else None

        # Determine the sample flags
        if 'MNGTARG1' in hdu['PRIMARY'].header:
            trg_flags = MaNGATargetBitMask().flagged_bits(hdu['PRIMARY'].header['MNGTARG1'])
        else:
            trg_flags = None

    target_type = 'None'
    if trg_flags is not None:
        if 'SECONDARY_v1_2_0' in trg_flags:
            target_type = 'S'
        if 'PRIMARY_v1_2_0' in trg_flags:
            target_type = 'P'
        if 'COLOR_ENHANCED_v1_2_0' in trg_flags:
            target_type = 'P+'

    #-------------------------------------------------------------------
    # Read the spectra and models to plot
    print('Reading MODEL file: {0}'.format(model_file))
    with fits.open(model_file) as hdu:
        wave, masked_spec, model_spec, model_continuum, model_emission \
                = build_model_spectra(hdu, masked=masked_spectra)

    # Start the plot
    w,h = pyplot.figaspect(1)
    fig = pyplot.figure(figsize=(2.5*w,h))

    # ------------------------------------------------------------------
    # Create the channel image plot
    on_sky_xlim = numpy.array([0, map_cube.shape[2]]) - 0.5
    on_sky_ylim = numpy.array([0, map_cube.shape[1]]) - 0.5
    image_ax = fig.add_axes([0.02, 0.25, 0.7/2.5, 0.7], facecolor='0.95')
    image_ax.tick_params(axis='both', which='both', direction='in')
    image_cax = fig.add_axes([0.02+0.7/2.5+0.01, 0.25, 0.01, 0.7])
    image_ax.set_xlim(on_sky_xlim)
    image_ax.set_ylim(on_sky_ylim)
    image_ax.grid(True, which='major', color='0.8', zorder=0, linestyle='-')

    # Create the initial plot, color bar, and labels
    im = image_ax.imshow(map_cube[0], origin='lower', interpolation='nearest',
                         cmap='RdBu_r', vmin=map_zmin[0], vmax=map_zmax[0], zorder=3)
    cbar = pyplot.colorbar(im, cax=image_cax)
    image_ax.text(0.5, 1.03, '{0}-{1}; MaNGA ID={2}; {3}'.format(plt, ifu, mangaid, target_type),
            horizontalalignment='center', verticalalignment='center', transform=image_ax.transAxes)
    image_ax.set_ylabel('Y [pix]')
    image_ax.set_xlabel('X [pix]')

    # Start the marker at the center of the map
    pos = tuple(*(numpy.array([map_cube.shape[1:][::-1]])/2).astype(int))
    marker = image_ax.scatter(*pos, marker='s', s=20, color='k', alpha=0.5, zorder=4)
    # ------------------------------------------------------------------


    # ------------------------------------------------------------------
    # Add the buttons and sliders for changing the map properties
    map_selector = UpdateableMap(fig, im, cbar, channel_names, channel_bits, map_cube,
                                 map_cube_ivar, map_cube_mask, map_zmin, map_zmax)
    # ------------------------------------------------------------------


    # ------------------------------------------------------------------
    # Add the image marker
    coo_pointer = Pointer(im.axes, pos, color='C3', lw=0.5, alpha=0.5)
    coo_selector = ImageMarker(marker, coo_pointer)
    # ------------------------------------------------------------------

    # ------------------------------------------------------------------
    # Add the spectrum plots
    _masked_spec = select_spectrum(masked_spec, pos)
    _model_spec = select_spectrum(model_spec, pos)
    _model_continuum = select_spectrum(model_continuum, pos)
    _model_emission = select_spectrum(model_emission, pos)

    wavelim = [wave[0], wave[-1]]
    fluxlim = [numpy.ma.amin(_masked_spec)-1, numpy.ma.amax(_masked_spec)+1]
    continuum_residlim = [numpy.ma.amin(_masked_spec-_model_continuum)-0.1,
                          numpy.ma.amax(_masked_spec-_model_continuum)+0.1]
    residlim = [numpy.ma.amin(_masked_spec-_model_spec)-0.1,
                numpy.ma.amax(_masked_spec-_model_spec)+0.1]

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

    
    spec, = spectrum_ax.step(wave, _masked_spec, where='mid', lw=0.5, color='k', zorder=3)
    cont, = spectrum_ax.plot(wave, _model_continuum, lw=0.5, color='b', zorder=4)
    fit, = spectrum_ax.plot(wave, _model_spec, lw=0.5, color='g', zorder=5)

    continuum_resid, = continuum_resid_ax.step(wave, _masked_spec - _model_continuum,
                                               where='mid', lw=0.5, color='k', zorder=3)
    em_fit, = continuum_resid_ax.plot(wave, _model_emission, lw=0.5, color='r', zorder=4)

    resid, = resid_ax.step(wave, _masked_spec - _model_spec,
                           where='mid', lw=0.5, color='k',zorder=3)

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
                                        binid[:,coo_selector.pos[1],coo_selector.pos[0]]),
                                 horizontalalignment='right', verticalalignment='center',
                                 transform=spectrum_ax.transAxes, fontsize=8)

    spectrum_selector = UpdateSpectrumPlot(binid, masked_spec, model_continuum, model_emission,
                                           model_spec, spec, cont, fit, continuum_resid, em_fit,
                                           resid, coo_selector, spx_coo_txt, binid_txt)
    coo_pointer.on_changed(spectrum_selector)
    # ------------------------------------------------------------------

    pyplot.show()


class MangaDapInspector(scriptbase.ScriptBase):

    @classmethod
    def name(cls):
        """
        Return the name of the executable.
        """
        return 'manga_dap_inspector'

    @classmethod
    def get_parser(cls, width=None):

        parser = super().get_parser(description='Inspect results of the MaNGA DAP', width=width)

        parser.add_argument('--maps_file', type=str, default=None, help='Output MAPS file')
        parser.add_argument('--model_file', type=str, default=None, help='Output LOGCUBE model file')

        parser.add_argument('--ext', type=str, nargs='*', help='List of map extensions to include',
                            default=None)

        parser.add_argument('--maskspec', action='store_true', default=False,
                            help='Use masks for spectra.')
        return parser

    @staticmethod
    def main(args):

        # Can it proceed?
        maps_file, model_file = args.maps_file, args.model_file
        if maps_file is None or model_file is None:
            raise ValueError('Must provide both a MAPS file and a model LOGCUBE file.')

        manga_dap_inspector(maps_file, model_file, ext=args.ext, masked_spectra=args.maskspec)


