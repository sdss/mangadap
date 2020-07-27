
Bandhead, or color, indices simply measure the ratio of fluxes in two
sidebands. Following the nomenclature defined for the
:ref:`spectralindices-absorption`, a color index is:

.. math::

    {\mathcal I} = \frac{\langle f\rangle_0}{\langle f\rangle_1},

where :math:`\langle y\rangle = S(y)/S(1)` and the two sidebands are
denoted with subscripts 0 and 1. The "order" of the index selects if
the index is calculated as the red-to-blue flux ratio or the
blue-to-red flux ratio (e.g., D4000 is defined as a red-to-blue
index, whereas TiOCvD is defined as a blue-to-red index).

Given the simple definition of color indices, a combined index for
the sum (or mean) of multiple spectra can be calculated by
constructing the weighted-mean index, where the continuum in the
denominator is used as the weight:

.. math::

    \frac{\sum_i \langle f_i\rangle_1 {\mathcal I}_i}{\sum_i  \langle f_i\rangle_1}
        = \frac{\sum_i \langle f_i\rangle_0} {\sum_i \langle f_i\rangle_1}
        = \frac{\langle\sum_i  f_i\rangle_0} {\langle\sum_i f_i\rangle_1}
