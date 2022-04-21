The absorption-line index calculations are performed as defined/used
by Worthey (1994) and Trager et al. (1998) (:math:`{\mathcal
I}_{\rm WT}`), as well as defined/used by Burstein et al. (1984) and
Faber et al. (1985) (:math:`{\mathcal I}_{\rm BF}`).

Specifically, let
    
.. math::

    S(y) \equiv \int_{\lambda_1}^{\lambda_2} y\ d\lambda 
            \approx \sum_i y_i\ {\rm d}p_i\ {\rm d}\lambda_i,

where :math:`{\rm d}p_i` is the fraction of pixel :math:`i` (with
width :math:`{\rm d}\lambda_i`) in the passband defined by
:math:`\lambda_1 < \lambda < \lambda_2`; the discrete sum is
performed by :func:`~mangadap.proc.bandpassfilter.passband_integral`.
Also, let :math:`f` be the spectrum flux density and define a linear
continuum using two sidebands ("blue" and "red"):

.. math::

    C(\lambda) = (\langle f\rangle_{\rm red} - \langle f\rangle_{\rm blue})\
        \frac{\lambda - \lambda_{\rm blue}}{\lambda_{\rm red}-\lambda_{\rm blue}}
         + \langle f\rangle_{\rm blue},


where :math:`\lambda_{\rm blue}` and :math:`\lambda_{\rm red}` are
the wavelengths at the center of the two sidebands --- e.g.,
:math:`\lambda_{\rm red} = (\lambda_{1,{\rm red}} + \lambda_{2,{\rm
red}})/2` --- and

.. math::

    \langle y\rangle = S(y)/S(1).

When **no pixels are masked** in the spectrum, :math:`S(1) =
(\lambda_2 - \lambda_1) \equiv \Delta\lambda`.

Following the Worthey (1994) definition, we can then calculate
the absorption-line spectral indices as follows:

.. math::

    {\mathcal I}_{\rm WT} = \left\{
            \begin{array}{ll}
                S(1 - f/C), & \mbox{for angstrom units} \\[3pt]
                -2.5\log\left[\langle f/C\rangle\right], & \mbox{for magnitude units}
            \end{array}\right..

The difficulty with the Worthey et al. definitions is that it makes
it difficult to construct an aggregate index measurement based on
individual index measurements from multiple spectra; however, this is
straight-forward under definitions closer to those provided by
Burstein et al. (1984) and Faber et al. (1985).

Let:

.. math::

    {\mathcal I}_{\rm BF} = \left\{
            \begin{array}{ll}
                S(1) - S(f)/C_0, & \mbox{for angstrom units} \\[3pt]
                -2.5\log\left[\langle f\rangle/C_0\right], & \mbox{for magnitude units}
            \end{array}\right.,

where :math:`C_0` is the value of the linear (not necessarily "flat")
continuum at the center of the main passband. Given that the
continuum is a linear function, :math:`S(C) = C_0 \Delta\lambda`;
i.e., the integral of the continuum over the passband is
mathematically identical to the continuum sampled at the center of
the bandpass times the bandpass width.

Now, we can calculate a weighted sum of indices using the value of
:math:`C_0` for each index as the weight and assume that no pixels
are masked such that we can replace :math:`S(1)` with
:math:`\Delta\lambda` to find:

.. math::

    \begin{eqnarray}
    \frac{\sum_i C_{0,i} {\mathcal I}_{\rm BF}}{\sum_i C_{0,i}}
        & = & \frac{\sum_i C_{0,i} (\Delta\lambda - S(f)_i / C_{0,i})}{\sum_i C_{0,i}} \\
        & = & \Delta\lambda - \frac{\sum_i S(f)_i}{\sum_i C_{0,i}}
    \end{eqnarray}.

That is, this weighted sum of the individual indices is
mathematically identical (to within the limits of how error affects
the construction of the linear continuum) to the index measured for
the sum (or mean) of the individual spectra. Similarly for the
indices in magnitude units:

.. math::

    \begin{eqnarray}
    -2.5\log\left[\frac{\sum_i C_{0,i} 10^{-0.4 {\mathcal I}_{\rm BF}}}{\sum_i C_{0,i}}\right]
        & = & -2.5\log\left[\frac{\sum_i C_{0,i} (S(f)_i / C_{0,i})}
                {\Delta\lambda \sum_i C_{0,i}}\right] \\
        & = & -2.5\log\left[\frac{\sum_i S(f)_i}{\Delta\lambda \sum_i C_{0,i}}\right]
    \end{eqnarray}.

Given the ease with which one can combine indices in the latter
definition, the DAP calculates both :math:`{\mathcal I}_{\rm WT}` and
:math:`{\mathcal I}_{\rm BF}`.

