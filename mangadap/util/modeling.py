"""
Module with generic functions used when modeling data.

----

.. include license and copyright
.. include:: ../include/copy.rst

----

.. include common links, assuming primary doc root is up one directory
.. include:: ../include/links.rst
"""

import numpy
from .filter import BoxcarFilter

def reject_residuals_1d(resid, lo=3.0, hi=3.0, boxcar=None):
    r"""
    Reject fit residuals.

    Args:
        resid (`numpy.ndarray`_, numpy.ma.MaskedArray`_):
            Weighted fit residuals. If input as a `numpy.ndarray`_,
            all data is included in the calculation of the local or
            global standard deviation. To exclude values from this
            calculation, the input *must* be a masked array. Shape
            must be either :math:`(N_x,)` or :math:`(N_{\rm
            vec},N_x})`. If 2D, the rejection is performed separately
            for each *row*.
        lo (:obj:`float`, optional):
            Sigma rejection for low values.
        hi (:obj:`float`, optional):
            Sigma rejection for high values.
        boxcar (:obj:`int`, optional):
            Size of a boxcar to use if using local sigma
            rejection. If None, rejection is based on the global
            standard deviation.

    Returns:
        `numpy.ndarray`_: Boolean array flagging values to reject.
    """
    if boxcar is None:
        mean = numpy.ma.mean(resid, axis=-1)
        std = numpy.ma.std(resid, axis=-1)
        # NOTE: Transposes avoid the need for numpy.newaxis usage.
        return ((numpy.ma.asarray(resid).data.T < mean - lo * std) \
                | (numpy.ma.asarray(resid).data.T > mean + hi * std)).T \
                    & numpy.logical_not(numpy.ma.getmaskarray(resid))

    bf = BoxcarFilter(boxcar, lo=lo, hi=hi, niter=1, y=resid, local_sigma=True)
    return numpy.squeeze(bf.output_mask) & numpy.logical_not(numpy.ma.getmaskarray(resid))


# TODO: Should maybe change the behavior of this so that if rng=None,
# the range is set to x[[0,-1]].
def scaled_coordinates(x, rng=None):
    r"""
    Scale the provided coordinates.

    The coordinates are scaled such that over the provided range the
    coordinates are rescaled to the range :math:`[-1, 1]`.

    .. warning::

        If ``rng`` is None, the returned object is identical to ``x``
        (i.e., it's not a copy of it).

    Args:
        x (`numpy.ndarray`_):
            Coordinate array
        rng (array-like, optional):
            The range over which to rescale to [-1,1]. If None, the
            input coordinates are simply returned.

    Returns:
        `numpy.ndarray`_: The rescaled coordinates.

    """
    return x if rng is None else 2 * (x - rng[0]) / (rng[1] - rng[0]) - 1 


class FitResiduals:
    def __init__(self, y, model, gpm=None):
        self.y = y
        self.model = model
        self.gpm = gpm
    def __call__(self, a):
        resid = (self.y-self.model(a))
        return resid if self.gpm is None else resid[self.gpm]


class FitChiSquare:
    def __init__(self, y, e, model, gpm=None):
        self.y = y
        self.e = e
        self.model = model
        _gpm = self.e > 0
        self.gpm = _gpm if gpm is None else gpm & _gpm
    def __call__(self, a):
        return (self.y[self.gpm]-self.model(a)[self.gpm])/self.e[self.gpm]


