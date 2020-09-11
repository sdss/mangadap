"""
Implements a base class that performs linear least-squares fitting of
1D basis functions constructed via a pseudo-Vandermonde matrix.

----

.. include license and copyright
.. include:: ../include/copy.rst

----

.. include common links, assuming primary doc root is up one directory
.. include:: ../include/links.rst
"""
import os
import time

from IPython import embed

import numpy

from . import modeling

# TODO: Generalize this even further so that a user can provide their
# own "basis" functions. And change the name of this module to
# basis.py?

class Vander1D:
    r"""
    Base class for least-squares fitting using linear combinations of 1D
    basis functions.

    Args:
        vander (callable):
            Function used to generate the pseudo-Vandermonde matrix
            of the basis functions; e.g.,
            `numpy.polynomial.legendre.legvander`_. Calling sequence
            must be ``vander(x,order)``.
        x (array-like):
            Ordinate for fitting.
        order (:obj:`int`):
            Order of the fit.
        rcond (:obj:`float`, optional):
            Relative condition number of the fit. Singular values
            smaller than this relative to the largest singular value
            will be ignored. The default value is len(x)*eps, where eps
            is the relative precision of the float type, about 2e-16 in
            most cases.

    Attributes:
        x (`numpy.ndarray`_):
            See the instantiation argument.  Shape is :math:`(N_{\rm
            x},)`.
        v (`numpy.ndarray`_):
            Pseudo-Vandermonde matrix of the request fitting order,
            :math:`o`.  Shape is :math:`(o+1, N_{\rm x})`.
        rcond (:obj:`float`):
            See the instantiation argument.
    """

    def __init__(self, vander, x, order, rcond=None):
        self.x = numpy.atleast_1d(x)
        if self.x.ndim != 1:
            raise ValueError('Ordinate must be a vector.')

        self.v = vander(self.x, order).T
        self.rcond = self.x.size * numpy.finfo(self.x.dtype).eps if rcond is None else rcond

    @property
    def size(self):
        """The size of the coordinate vector."""
        return self.x.size

    @property
    def shape(self):
        """The shape of the coordinate vector."""
        return self.x.shape

    @property
    def order(self):
        """The order of the function."""
        return self.v.shape[0]-1

    def fit(self, y, w=None, err=None, mask=None, rej_iter=0, rej_lo=3., rej_hi=3., rej_win=None,
            model=False):
        r"""
        Fit the provided data.

        The ordinate is provided on the class instantiation.

        If provided, the errors, :math:`\epsilon`, are converted to
        weights to use during the fitting, where :math:`\omega_i =
        1/\epsilon_i`. Both weights, :math:`w`, and errors can be
        provided; however, in this case, the weights used during the
        fit are are :math:`\omega_i = w_i / \epsilon_i`.

        .. todo::
            
            Allow for the nominal calculation of the errors in the
            coefficients.

        Args:
            y (array-like):
                A 1D or 2D array with the 1D functions to fit. Shape
                is :math:`(N_x,)` or :math:`(N_{\rm vec},N_x)`, where
                :math:`N_x` is the length of the ordinate
                (:attr:`x`). Input can be a masked array.
            w (array-like, optional):
                A 1D or 2D array with the weights for each pixel.
                Shape is :math:`(N_x,)` or :math:`(N_{\rm vec},N_x)`,
                where :math:`N_x` is the length of the ordinate
                (:attr:`x`). If ``w`` is 1D, but ``y`` is 2D, the
                same weights are applied to all vectors.
            err (array-like, optional):
                Error in ``y``, used for weighting. Shape must match
                ``y``.
            mask (array-like, optional):
                A boolean mask; values to ignore should have ``mask
                == True``. Shape must match ``y``. Treated
                identically as a mask provided by passing in a
                `numpy.ma.MaskedArray`_ for ``y``. Note that
                providing a weight of 0 is identical to masking the
                value.
            rej_iter (:obj:`int`, optional):
                Number of fit rejection iterations.  If 0, no rejection
                iterations are performed.  If -1, rejections are
                performed until no points are rejected in a given
                iteration.
            rej_lo (scalar-like, optional):
                The lower sigma rejection threshold.
            rej_hi (scalar-like, optional):
                The higher sigma rejection threshold.
            rej_win (:obj:`int`, optional):
                The length of a boxcar window to use in a calculation of
                a *local* standard deviation used during the rejection
                iterations.  This is useful if the S/N varies strongly
                within the provided vector being fit.  If None, no
                windowing is applied, meaning a *global* standard
                deviation is used for rejection.
            model (:obj:`bool`, optional):
                If True, return the model as well as the best fitting
                coefficients.

        Returns:
            `numpy.ndarray`, :obj:`tuple`: As many as three objects
            can be returned. Regardless of the input, the function
            always returns an array with the best-fitting
            coefficients, with shape :math:`(o+1,)` or :math:`(N_{\rm
            vec},o+1)`, depending on the input shape of ``y``. If
            ``model`` is True, an array with the best-fitting models
            are also returned (same shape as ``y``). If rejection
            operations are performed, a boolean array flagging the
            pixels excluded from the fit is also returned. Note that
            any data that will have been ignored based on the input
            (the element is masked or has 0 weight) is *not* flagged
            as having been rejected.
        """
        # Check the input vectors to fit
        _y = numpy.atleast_1d(y)
        if _y.ndim != 1 and _y.ndim != 2:
            raise ValueError('Input y coordinates must be a single vector or 2D array.')
        if _y.shape[-1] != self.x.size:
            raise ValueError('The last dimension of y must have the same length as the ordinate.')

        # Setup the weights
        if w is None:
            # Weights are needed if rejecting
            _w = None if rej_iter == 0 else numpy.ones(_y.shape, dtype=float) 
        else:
            _w = numpy.atleast_1d(w)
            if _w.ndim != 1 and _w.ndim != 2:
                raise ValueError('Input weights must be a single vector or 2D array.')
            if _w.ndim == 1 and _w.shape != self.x.shape:
                raise ValueError('If a single vector, the weights must have the same shape as '
                                 'the ordinate vector.')
            if rej_iter != 0 and _w.ndim != _y.ndim:
                # If rejecting, weights must be applied to each vector
                # separately
                _w = numpy.tile(_w, (_y.shape[0],1))
            if _w.ndim == 2 and _w.shape != _y.shape:
                raise ValueError('If an array, the weights must have the same shape as the data '
                                 'to fit.')
            if numpy.any(_w < 0):
                raise ValueError('Weights must be positive!')
        if isinstance(_y, numpy.ma.MaskedArray):
            _m = numpy.logical_not(numpy.ma.getmaskarray(_y)).astype(float)
            _w = _m if _w is None else _m * _w
            _y = _y.filled(0.0)
        if mask is not None:
            _m = numpy.logical_not(numpy.atleast_1d(mask)).astype(float)
            if _m.shape != _y.shape:
                raise ValueError('Mask must have the same shape as the data to fit.')
            _w = _m if _w is None else _m * _w

        # Check the input error and add it to the weights, if provided
        if err is not None:
            # Below is faster than numpy.ma.power(err, -1).filled(0.0)
            # and avoids errors when err is 0.
            _e = err.filled(0.0) if isinstance(err, numpy.ma.MaskedArray) else err
            _e = (_e > 0)/(_e + (_e == 0.))
            if _e.shape != _y.shape:
                raise ValueError('Input errors must have the same shape as the data to fit.')
            _w = _e if _w is None else _e * _w

        # Fit the data
        c = self._single_fit(_y, _w)

        # If no rejection iterations, we're done
        if rej_iter == 0:
            return (c, self.model(c)) if model else c

        # Initial mask
        init_good = _w > 0

        # Iteratively reject
        i = 0
        while i < rej_iter or rej_iter < 0:
            # Setup a masked array to hold the weighted residuals
            resid = _w * numpy.ma.MaskedArray(_y - self.model(c), mask=numpy.logical_not(_w > 0))
            # Reject pixels
            rej = modeling.reject_residuals_1d(resid, lo=rej_lo, hi=rej_hi, boxcar=rej_win)
            # Determine if any were rejected
            refit = numpy.any(rej, axis=-1)
            if not numpy.any(refit):
                # None were rejected so we're done
                break
            # Set the weight to 0 for the rejected data
            _w *= numpy.logical_not(rej).astype(float)
            # Refit
            c = self._single_fit(_y, _w)
            # Increment iteration
            i += 1

        rejected = numpy.logical_not(_w > 0) & init_good
        return (c, self.model(c), rejected) if model else (c, rejected)

    def model(self, c):
        r"""
        Return the model function constructed using the provided
        coefficients.

        Args:
            c (array-like):
                The coefficients for the model with shape shape
                :math:`(o+1,)` or :math:`(N_{\rm vec},o+1)`, where
                :math:`o` is the order of the function.

        Returns:
            `numpy.ndarray`_: Array with the evaluated functions. The
            output shape is :math:`(N_x,)` or :math:`(N_{\rm
            vec},N_x)` depending on the input shape of ``c``.
        """
        return numpy.dot(c, self.v)

    def _single_fit(self, y, w):
        r"""
        Perform a single fitting iteration.

        Args:
            y (`numpy.ndarray`_):
                The data to fit.  Shape is :math:`(N_{\rm x},)` or
                :math:`(N_{\rm vec}, N_{\rm x})`.  If the latter, a fit
                is done to each vector.
            w (`numpy.ndarray`_):
                Weights.  Shape is :math:`(N_{\rm x},)` or
                :math:`(N_{\rm vec}, N_{\rm x})`.  If the latter, this
                call is equivalent to looping over :math:`N_{\rm vec}`
                and fitting each vector in ``y`` individually.  Can be
                None, and if so the fit is unweighted.
        
        Returns:
            `numpy.ndarray`_: The best-fitting coefficients for the
            basis functions with order :math:`o`.  Shape is
            :math:`(o+1,)` or :math:`(N_{\rm vec},o+1)`.
        """
        if w is None:
            # Unweighted fit. For multiple y vectors, make sure
            # returned array has the shape described by the return
            # statement above.
            return numpy.linalg.lstsq(self.v.T, y.T, rcond=self.rcond)[0].T
        
        if y.ndim > 1 and w.shape == y.shape:
            # If weighting, can only solve all fits simultaneously if
            # the same weights are used for all vectors. Therefore,
            # need to loop through all y vectors here to accommodate
            # weights that may be different fro each y vector.
            c = numpy.zeros((y.shape[0], self.v.shape[0]), dtype=float)
            for i in range(y.shape[0]):
                c[i] = self._single_fit(y[i], w[i])
            return c

        # w is a single vector, used when fitting one or more y
        # vectors.
        lhs = self.v * w
        rhs = y * w
        scl = numpy.sqrt(numpy.square(lhs).sum(axis=1))
        scl[scl == 0] == 1
        c = numpy.linalg.lstsq(lhs.T/scl, rhs.T, rcond=self.rcond)[0].T
        return c/scl


class Legendre1D(Vander1D):
    """
    Least-squares fitting of 1D Legendre functions.

    Args:
        x (array-like):
            Ordinate for fitting.
        order (:obj:`int`):
            Order of the fit.
        rcond (:obj:`float`, optional):
            Relative condition number of the fit. Singular values
            smaller than this relative to the largest singular value
            will be ignored. The default value is len(x)*eps, where eps
            is the relative precision of the float type, about 2e-16 in
            most cases.
    """
    def __init__(self, x, order, rcond=None):
        super(Legendre1D, self).__init__(numpy.polynomial.legendre.legvander, x, order,
                                         rcond=rcond)

