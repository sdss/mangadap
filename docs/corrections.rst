
.. include:: include/links.rst

.. _corrections:

MAPS Corrections
================

The :ref:`datamodel-maps` provide values that *must be corrected by
the user* using provided corrections. It's recommended that you take
advantage of the convenience methods in `Marvin`_ to apply these
corrections.

----

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

**In both cases**, beware of imaginary numbers. That is, when the
correction is larger than the provided value, the above equations
result in taking the sqrt of a negative number. Specifically for the
stellar velocity dispersions, we recommend you consult Section 7.7 of
`Westfall et al. (2019, AJ, 158, 231)`_ for some usage guidelines and
discussion.

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

The `Marvin`_ method that applies these corrections is
`marvin.tools.quantities.map.Map.inst_sigma_correction
<https://sdss-marvin.readthedocs.io/en/latest/reference/quantities.html#marvin.tools.quantities.map.Map.inst_sigma_correction>`_.

----

Spectral-Index Measurements
---------------------------

Corrections that account for the effect of the velocity dispersion on
the spectral indices are provided, as discussed in Section 10.1 of
the `Westfall et al. (2019, AJ, 158, 231)`_. Unlike, e.g., the
`Firefly VAC`_, these corrections *must be applied by the user*. To
apply the corrections, you have to know the unit of each index. For
angstrom units (or for unitless bandhead/color indices):

.. math::

    \mathcal{I}^c = \mathcal{I}\ \delta\mathcal{I}

and for magnitude units:

.. math::

    \mathcal{I}^c = \mathcal{I} + \delta\mathcal{I}

where the raw index measurements, :math:`\mathcal{I}`, and the
correction, :math:`\delta\mathcal{I}` are provided in, respectively,
the ``SPECINDEX`` and ``SPECINDEX_CORR`` extensions of the
:ref:`datamodel-maps`. Correction are identical for both index
definitions, :math:`{\mathcal I}_{\rm WT}` and :math:`{\mathcal
I}_{\rm BF}`; see :ref:`spectralindices`. Corrections for the weights
should only be applied when aggregating *corrected* indices, and the
weight corrections are multiplicative:

.. math::

    w^c = w\ \delta w .

The `Marvin`_ method that applies these corrections is
`marvin.tools.quantities.map.Map.specindex_correction
<https://sdss-marvin.readthedocs.io/en/latest/reference/quantities.html#marvin.tools.quantities.map.Map.specindex_correction>`_.

----

Usage Example
-------------

A :ref:`gettingstarted-maps-example` is included as our discussion of
:ref:`gettingstarted`.

