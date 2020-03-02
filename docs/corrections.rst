
.. _DAP Overview Paper: <https://ui.adsabs.harvard.edu/abs/2019AJ....158..231W/abstract>

.. _corrections:

MAPS Corrections
================

The :ref:`datamodel-maps` provide values that *must be corrected by
the user* using provided corrections. It's recommended that you take
advantage of the convenience methods in `Marvin
<http://sdss-marvin.readthedocs.io/en/stable/>`_ to apply these
corrections.

Velocity-Dispersion Measurements
--------------------------------

Note that the stellar and gas velocity dispersions *must be corrected
for instrumental resolution effects* to obtain the astrophysical Doppler
broadening.

The corrected gas velocity dispersion is:

.. math::

    \sigma_{\rm gas}^2 = \sigma_{\rm gas,obs}^2 - \sigma_{\rm inst}^2

where :math:`\sigma_{\rm gas}` and :math:`\sigma_{\rm inst}` are
provided in, respectively, the ``EMLINE_GSIGMA`` and
``EMLINE_INSTSIGMA`` extensions of the :ref:`datamodel-maps`.

The corrected stellar velocity dispersion is:

.. math::

    \sigma_\ast^2 = \sigma_{\rm obs}^2 - \delta\sigma_{\rm inst}^2

where :math:`\sigma_{\rm obs}` and :math:`\delta\sigma_{\rm inst}` are
provided in, respectively, the ``STELLAR_SIGMA`` and
``STELLAR_SIGMACORR`` extensions of the :ref:`datamodel-maps`.

**In both cases**, beware of imaginary numbers.  That is, when the
correction is larger than the provided value, the above equations result
in taking the sqrt of a negative number.  Specifically for the stellar
velocity dispersions, we recommend you consult Section 7.7 of the `DAP
Overview Paper`_ for some usage guidelines and discussion.

In particular, we have found that it is important to understand the
error-convolved *distribution* of stellar-velocity-dispersion
measurements when analyzing the data.  For example, ignoring anything
that has a converged pPXF fit, even at low S/N and low dispersion, will
yield a biased determination of the mean (or median) dispersion as a
function of radius.  Further assessments of the reliability of the data
is an ongoing process:  Any additional assessments of the data along
these lines from the collaboration is more than welcome.

Stellar velocity dispersions are currently provided for two approaches
to the calculation:  A nominal correction is calculated using the
quadrature difference between the instrumental dispersion of the
template and galaxy spectra over the fitted wavelength range.  This is
the correction provided in MPL-5, MPL-7/DR15.  In MPL-8 and later, we
also provide a correction based on a fit of the optimal template with
and without the resolution matched to the MaNGA data.  **For now**, use
the correction in the first channel of the ``STELLAR_SIGMACORR``
extension until the data in the second channel can be vetted.

Spectral-Index Measurements
---------------------------

Corrections that account for the effect of the velocity dispersion on
the spectral indices are provided, as discussed in Section 10.1 of the
`DAP Overview Paper`_.  Unlike the Firefly VAC, these corrections *must
be applied by the user*.  To apply the corrections, you have to know the
unit of each index.  For angstrom units:

.. math::

    \mathcal{I}^c_a k= \mathcal{I}_a \delta\mathcal{I}_a

and for magnitude units:

.. math::

    \mathcal{I}^c_a k= \mathcal{I}_a + \delta\mathcal{I}_a

where the raw index measurements, :math:`\mathcal{I}_a`, and the
correction, :math:`\delta\mathcal{I}_a` are provided in, respectively,
the ``SPECINDEX`` and ``SPECINDEX_CORR`` extensions of the
:ref:`datamodel-maps`.


Usage Example
-------------

A :ref:`gettingstarted-maps-example` is included as our discussion of
:ref:`gettingstarted`.

