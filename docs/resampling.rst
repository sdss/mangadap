
.. include:: include/links.rst

.. _resampling:

Resampling
==========

The DAP is picky about the spectral sampling required to fit the spectra; the
spectra must be geometrically sampled in wavelength (constant steps in, e.g.,
:math:`\log\lambda`.  Ideally, the data-reduction procedures would provide this
sampling directly for your data (even if this is a rather unnatural sampling of
data from a spectrograph).  However, if the spectra are, e.g., sampled linearly
in wavelength, the DAP provides a utility to robustly resample the spectra.

.. warning::

    The DAP is, in fact, *very* picky about the sampling.  Even if you expect
    the sampling is geometric, the DAP may throw an error if the wavelength
    vector does not yield a constant geometric step to (nearly) numerical
    accuracy.  If that's the case, please `Submit an issue`_.  The likely
    solution will be to use :class:`~mangadap.util.sampling.Resample` to
    resample the data anyway, or you can follow the riskier approach to fool the
    DAP by remaking the wavelength vector so the step is constant to numerical
    accuracy (e.g., using ``numpy.logspace``).

The example script `test_resample.py
<https://github.com/sdss/mangadap/blob/master/examples/test_resample.py>`__
illustrates how to use :class:`~mangadap.util.sampling.Resample` to resample a
spectrum.  This class is a generalization of the ``log_rebin`` function provided
by Michele Cappellari in the `pPXF`_ package.  The script uses the resampling
code to de-redshift example spectra from MaNGA, and compares the results to a
brute force approach and, if it is installed in your environment, `spectres`_.
The script should result in identical results for the brute force, `spectres`_,
and :class:`~mangadap.util.sampling.Resample`; however,
:class:`~mangadap.util.sampling.Resample` is ~100 times faster than the
brute-force approach and ~10 times faster than `spectres`_.

For another practical example of how to use
:class:`~mangadap.util.sampling.Resample`, specifically for a full datacube, see
the resampling of a KCWI datacube here: :ref:`fitdatacube-kcwi`.


