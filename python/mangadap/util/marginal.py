# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-
"""
Provides a set of functions to handle resampling.

*License*:
    Copyright (c) 2019, SDSS-IV/MaNGA Pipeline Group
        Licensed under BSD 3-clause license - see LICENSE.rst

*Source location*:
    $MANGADAP_DIR/python/mangadap/util/marginal.py

*Revision history*:
    | **13 Sep 2019**: Original implementation by K.  Westfall (KBW)
"""

import numpy
from matplotlib import pyplot

class MarginalizedDistribution:
    r"""
    Provide a utility to construct 1D and 2D marginalized distributions
    of a set of samples in an N-dimensional feature space.

    Args:
        samples (array-like):
            The samples of the full N-dimensional distribution.  The
            shape of the array must be :math:`(N_{\rm samples}, N_{\rm
            feature})`.
        bins (:obj:`int`, array-like, optional):
            The number of bins for each feature.  If a single integer,
            the same number of bins are used for all features;
            otherwise, must have shape :math:`(N_{\rm feature},)`.
        span (array-like, optional):
            The minimum and maximum span for the binned data.  Provide a
            2-element list for a single span for all features, or
            provide an array with shape :math:`(N_{\rm feature},2)` to
            define the span of each feature.  If `None`, the span is set
            to the minimum and maximum of each feature.
        log (:obj:`int`, :obj:`bool`, array-like, optional):
            Bin the features geometrically.  Single booleans are applied
            to all features.  Boolean arrays provide the matching list
            of features to bin; the shape must be :math:`(N_{\rm
            feature},)`.  Single or arrays of integers select the
            indices of the feature(s) to geometrically bin.

    """
    def __init__(self, samples, bins=10, span=None, log=None):
        self.samples = numpy.atleast_2d(samples)
        if self.samples.ndim != 2:
            raise ValueError('Samples must be provided in a 2D array: (Nsamples, Ndim).')
        self.nsamp, self.ndim = self.samples.shape
        if self.ndim < 2:
            raise ValueError('Samples must be for 2 or more dimensions.')
        self._set_sampling(bins, span, log)
        self._pdf = numpy.empty(self.ndim, dtype=object)
        self._joint = numpy.empty((self.ndim*self.ndim - self.ndim)//2, dtype=object)

    def _set_sampling(self, bins, span, log):
        self.bins = numpy.atleast_1d(bins) \
                        if hasattr(bins, '__len__') else numpy.full(self.ndim, bins, dtype=int)
        if self.bins.ndim > 1 or self.bins.size != self.ndim:
            raise ValueError('Must provide a single number of bins or one per sample dimension.')

        self.span = numpy.array([[numpy.amin(self.samples[:,i]),numpy.amax(self.samples[:,i])]
                                        for i in range(self.ndim)], dtype=float) \
                        if span is None else numpy.atleast_2d(span)
        if self.span.ndim != 2:
            raise ValueError('Range must be 2D.')
        if self.span.shape[0] == 1:
            self.span = numpy.tile(self.span, (self.ndim,1))
        if self.span.shape[1] != 2:
            raise ValueError('Second dim of input span must be 2; start and end.')

        self.log = numpy.zeros(self.ndim, dtype=bool)
        if log is not None:
            _log = numpy.atleast_1d(log)
            if _log.dtype == 'int':
                self.log[_log] = True
            elif _log.dtype == 'bool':
                if _log.size == 1:
                    _log = numpy.repeat(_log, self.ndim)
                if _log.size != self.ndim:
                    raise ValueError('Log entry has incorrect length.')
                self.log = _log
            else:
                raise ValueError('Log entry has incorrect dtype.')
        if numpy.any(self.log):
            _samples = self.samples.copy()
            _samples[:,self.log] = numpy.log10(self.samples[:,self.log])
            self.samples = _samples
            self.span[self.log,:] = numpy.log10(sel.span[self.log,:])

        self.start = self.span[:,0]
        self.step = numpy.squeeze(numpy.diff(self.span, axis=1))/self.bins

        self.grid = numpy.array([s + d*numpy.arange(b+1) 
                                    for s, b, d in zip(self.start, self.bins, self.step)],
                                dtype = object)

    def _bin_oned(self, i):
        bini = numpy.ma.MaskedArray(
                    numpy.floor((self.samples[:,i]-self.start[i])/self.step[i]).astype(int))
        bini[(bini < 0) | (bini >= self.bins[i])] = numpy.ma.masked
        uniq, cnt = numpy.unique(bini.compressed(), return_counts=True)
        n = numpy.zeros(self.bins[i], dtype=int)
        n[uniq] = cnt
        return n

    def _bin_twod(self, i, j):
        bini = numpy.floor((self.samples[:,i] - self.start[i])/self.step[i]).astype(int)
        bini[ (bini < 0) | (bini >= self.bins[i]) ] = -1
        binj = numpy.floor((self.samples[:,j] - self.start[j])/self.step[j]).astype(int)
        binj[ (binj < 0) | (binj >= self.bins[j]) ] = -1

        binij = numpy.ma.MaskedArray(bini*self.bins[j] + binj, mask=(bini < 0) | (binj < 0))

        uniq, cnt = numpy.unique(binij.compressed(), return_counts=True)
        n = numpy.zeros((self.bins[i],self.bins[j]),dtype=int)
        _i = uniq//self.bins[j]
        _j = uniq - _i*self.bins[j]
        n[_i,_j] = cnt
        return n

    def cdf(self, i):
        return numpy.sort(self.samples[:,i]), (numpy.arange(self.nsamp)+1)/self.nsamp, 

    def pdf(self, i):
        if self._pdf[i] is None:
            self._pdf[i] = self._bin_oned(i)
        g = (self.grid[i][:-1] + self.grid[i][1:])/2
        return (numpy.power(10, g) if self.log[i] else g), self._pdf[i]

    def _joint_indx(self, i, j):
        if i > j:
            return self._joint_indx(j,i)[0], True
        return (j if i == 0 else numpy.cumsum(numpy.arange(self.ndim-1)[::-1])[i-1]+j)-1, False
    
    def joint(self, i, j):
        if i == j:
            return self.pdf(i)
        ii, transpose = self._joint_indx(i,j)
        if self._joint[ii] is None:
            self._joint[ii] = self._bin_twod(j,i) if transpose else self._bin_twod(i,j)
        x = numpy.power(10, self.grid[i]) if self.log[i] else self.grid[i]
        y = numpy.power(10, self.grid[j]) if self.log[j] else self.grid[j]
        return x, y, (self._joint[ii].T if transpose else self._joint[ii])


class CornerPlot(MarginalizedDistribution):
    def __init__(self, samples, bins=10, span=None, log=None, label=None, hist_scl=1.0):
        super(CornerPlot, self).__init__(samples, bins=bins, span=span, log=log)
        if label is None:
            self.label = [ 'Feature {0}'.format(i+1) for i in range(self.ndim) ]
        self.hist_scl = hist_scl
        self.build_axes()
        self.show()

    @classmethod
    def from_dist(cls, dist, label=None, hist_scl=1.0):
        return cls(dist.samples, bins=dist.bins, span=dist.span, log=dist.log, label=label,
                   hist_scl=hist_scl)

    def _axes_bounds(self):
        bounds = [0.1, 0.1, 0.88, 0.88]
        sep = 0.01
        dx = (bounds[2] - (self.ndim-1)*sep)/(self.ndim + self.hist_scl - 1)
        dy = (bounds[3] - (self.ndim-1)*sep)/(self.ndim + self.hist_scl - 1)
        # The x limits are always used in reverse order
        x = bounds[0]+numpy.arange(self.ndim)[::-1]*(dx+sep)
        y = bounds[1]+numpy.arange(self.ndim)*(dy+sep)
        self.axes_bounds = numpy.empty((self.ndim, self.ndim), dtype=object)
        # Diagonal has histograms in reverse order
        self.axes_bounds[numpy.arange(self.ndim),numpy.arange(self.ndim)] \
            = [[_x, _y, dx, dy] for _x,_y in zip(x,y)]
        # Rescale the height of the histogram plots
        self.axes_bounds[0,0][2] *= self.hist_scl
        for i in range(1,self.ndim):
            self.axes_bounds[i,i][3] *= self.hist_scl

        for i in range(self.ndim-1):
            self.axes_bounds[i+1:,i] = [[_x, y[i], dx, dy] for _x in x[i+1:]]

        print(self.axes_bounds)
        return

        # Bottom row of scatter plots show all quantities against the
        # first one

        self.axes_bounds[1:,0] = [ [_x, y[0], dx, dy] for _x in x[1:]]
        # Remainder of rows start with the first quantity along the y
        # axes and all but the last one along the x axis, in reverse
        # order
        for i in range(self.ndim-1):
            self.axes_bounds[i,i+1:-1] = [ (_x, y[i+1], dx, dy) for _x in x[i+2:]]

    def build_axes(self):
        w,h = pyplot.figaspect(1)
        self.fig = pyplot.figure(figsize=(1.5*w,1.5*h))

        # Define each plot boundary
        self._axes_bounds()

        # Build the group of axes
        self.axes = numpy.empty((self.ndim, self.ndim), dtype=object)
        for i in range(self.ndim):
            for j in range(self.ndim):
                if self.axes_bounds[i,j] is None:
                    continue
                self.axes[i,j] = self.fig.add_axes(self.axes_bounds[i,j])

    def show(self):
        for i in range(self.ndim):
            x, y = self.pdf(i)
            if i == 0:
                y, x = x, y
                self.axes[i,i].set_ylabel(self.label[i])
            else:
                self.axes[i,i].set_xlabel(self.label[i])
            self.axes[i,i].step(x, y, where='mid')
        for i in range(self.ndim):
            for j in range(i+1,self.ndim):
                if i == j or self.axes[j,i] is None:
                    continue
                x, y, img = self.joint(j,i)
                img = numpy.ma.MaskedArray(img, mask=img < 1)
                X, Y = numpy.meshgrid(x, y, indexing='ij')
                self.axes[j,i].pcolormesh(X,Y,img)
                self.axes[j,i].set_xlabel(self.label[j])
                self.axes[j,i].set_ylabel(self.label[i])
        pyplot.show()


