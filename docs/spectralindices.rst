
.. |ang|   unicode:: U+212B

.. include:: include/links.rst

.. _spectralindices:

Spectral Indices
================

The spectral indices to be measured by the DAP are divided into two
groups: (1) absorption-line indices that measure the equivalent width
of an absorption feature and (2) bandhead or color indices that
measure the ratio of fluxes is in two passbands. Both sets of indices
are measured using
:class:`~mangadap.proc.spectralindices.SpectralIndices`; see
:ref:`spectral-index-measurements`.

----

.. _spectralindices-absorption:

Absorption-line Indices
-----------------------

Calculation
~~~~~~~~~~~

.. include:: include/absindices.rst

See :class:`~mangadap.proc.spectralindices.AbsorptionLineIndices`.

Index Parameters
~~~~~~~~~~~~~~~~

The table below provides a compilation of absorption-line indices.
Recent survey-level runs of the DAP have included all of these
measurements.

+-----------+-----------------------+-----------------------+----------------------+-------+-------+------+
| Name      | Main Passband (|ang|) | Blue Sideband (|ang|) | Red Sideband (|ang|) | Frame | Units |  Ref |
+===========+=======================+=======================+======================+=======+=======+======+
| CN1       |  4142.125 -- 4177.125 |  4080.125 -- 4117.625 | 4244.125 -- 4284.125 |   air |   mag | [1]_ |
+-----------+-----------------------+-----------------------+----------------------+-------+-------+------+
| CN2       |  4142.125 -- 4177.125 |  4083.875 -- 4096.375 | 4244.125 -- 4284.125 |   air |   mag | [1]_ |
+-----------+-----------------------+-----------------------+----------------------+-------+-------+------+
| Ca4227    |  4222.250 -- 4234.750 |  4211.000 -- 4219.750 | 4241.000 -- 4251.000 |   air | |ang| | [1]_ |
+-----------+-----------------------+-----------------------+----------------------+-------+-------+------+
| G4300     |  4281.375 -- 4316.375 |  4266.375 -- 4282.625 | 4318.875 -- 4335.125 |   air | |ang| | [1]_ |
+-----------+-----------------------+-----------------------+----------------------+-------+-------+------+
| Fe4383    |  4369.125 -- 4420.375 |  4359.125 -- 4370.375 | 4442.875 -- 4455.375 |   air | |ang| | [1]_ |
+-----------+-----------------------+-----------------------+----------------------+-------+-------+------+
| Ca4455    |  4452.125 -- 4474.625 |  4445.875 -- 4454.625 | 4477.125 -- 4492.125 |   air | |ang| | [1]_ |
+-----------+-----------------------+-----------------------+----------------------+-------+-------+------+
| Fe4531    |  4514.250 -- 4559.250 |  4504.250 -- 4514.250 | 4560.500 -- 4579.250 |   air | |ang| | [1]_ |
+-----------+-----------------------+-----------------------+----------------------+-------+-------+------+
| C24668    |  4634.000 -- 4720.250 |  4611.500 -- 4630.250 | 4742.750 -- 4756.500 |   air | |ang| | [1]_ |
+-----------+-----------------------+-----------------------+----------------------+-------+-------+------+
| Hb        |  4847.875 -- 4876.625 |  4827.875 -- 4847.875 | 4876.625 -- 4891.625 |   air | |ang| | [1]_ |
+-----------+-----------------------+-----------------------+----------------------+-------+-------+------+
| Fe5015    |  4977.750 -- 5054.000 |  4946.500 -- 4977.750 | 5054.000 -- 5065.250 |   air | |ang| | [1]_ |
+-----------+-----------------------+-----------------------+----------------------+-------+-------+------+
| Mg1       |  5069.125 -- 5134.125 |  4895.125 -- 4957.625 | 5301.125 -- 5366.125 |   air |   mag | [1]_ |
+-----------+-----------------------+-----------------------+----------------------+-------+-------+------+
| Mg2       |  5154.125 -- 5196.625 |  4895.125 -- 4957.625 | 5301.125 -- 5366.125 |   air |   mag | [1]_ |
+-----------+-----------------------+-----------------------+----------------------+-------+-------+------+
| Mgb       |  5160.125 -- 5192.625 |  5142.625 -- 5161.375 | 5191.375 -- 5206.375 |   air | |ang| | [1]_ |
+-----------+-----------------------+-----------------------+----------------------+-------+-------+------+
| Fe5270    |  5245.650 -- 5285.650 |  5233.150 -- 5248.150 | 5285.650 -- 5318.150 |   air | |ang| | [1]_ |
+-----------+-----------------------+-----------------------+----------------------+-------+-------+------+
| Fe5335    |  5312.125 -- 5352.125 |  5304.625 -- 5315.875 | 5353.375 -- 5363.375 |   air | |ang| | [1]_ |
+-----------+-----------------------+-----------------------+----------------------+-------+-------+------+
| Fe5406    |  5387.500 -- 5415.000 |  5376.250 -- 5387.500 | 5415.000 -- 5425.000 |   air | |ang| | [1]_ |
+-----------+-----------------------+-----------------------+----------------------+-------+-------+------+
| Fe5709    |  5696.625 -- 5720.375 |  5672.875 -- 5696.625 | 5722.875 -- 5736.625 |   air | |ang| | [1]_ |
+-----------+-----------------------+-----------------------+----------------------+-------+-------+------+
| Fe5782    |  5776.625 -- 5796.625 |  5765.375 -- 5775.375 | 5797.875 -- 5811.625 |   air | |ang| | [1]_ |
+-----------+-----------------------+-----------------------+----------------------+-------+-------+------+
| NaD       |  5876.875 -- 5909.375 |  5860.625 -- 5875.625 | 5922.125 -- 5948.125 |   air | |ang| | [1]_ |
+-----------+-----------------------+-----------------------+----------------------+-------+-------+------+
| TiO1      |  5936.625 -- 5994.125 |  5816.625 -- 5849.125 | 6038.625 -- 6103.625 |   air |   mag | [1]_ |
+-----------+-----------------------+-----------------------+----------------------+-------+-------+------+
| TiO2      |  6189.625 -- 6272.125 |  6066.625 -- 6141.625 | 6372.625 -- 6415.125 |   air |   mag | [1]_ |
+-----------+-----------------------+-----------------------+----------------------+-------+-------+------+
| HDeltaA   |  4083.500 -- 4122.250 |  4041.600 -- 4079.750 | 4128.500 -- 4161.000 |   air | |ang| | [2]_ |
+-----------+-----------------------+-----------------------+----------------------+-------+-------+------+
| HGammaA   |  4319.750 -- 4363.500 |  4283.500 -- 4319.750 | 4367.250 -- 4419.750 |   air | |ang| | [2]_ |
+-----------+-----------------------+-----------------------+----------------------+-------+-------+------+
| HDeltaF   |  4091.000 -- 4112.250 |  4057.250 -- 4088.500 | 4114.750 -- 4137.250 |   air | |ang| | [2]_ |
+-----------+-----------------------+-----------------------+----------------------+-------+-------+------+
| HGammaF   |  4331.250 -- 4352.250 |  4283.500 -- 4319.750 | 4354.750 -- 4384.750 |   air | |ang| | [2]_ |
+-----------+-----------------------+-----------------------+----------------------+-------+-------+------+
| CaHK      |  3899.5   -- 4003.5   |  3806.5   -- 3833.8   | 4020.7   -- 4052.4   |   air | |ang| | [3]_ |
+-----------+-----------------------+-----------------------+----------------------+-------+-------+------+
| CaII1     |  8484.0   -- 8513.0   |  8474.0   -- 8484.0   | 8563.0   -- 8577.0   |   air | |ang| | [4]_ |
+-----------+-----------------------+-----------------------+----------------------+-------+-------+------+
| CaII2     |  8522.0   -- 8562.0   |  8474.0   -- 8484.0   | 8563.0   -- 8577.0   |   air | |ang| | [4]_ |
+-----------+-----------------------+-----------------------+----------------------+-------+-------+------+
| CaII3     |  8642.0   -- 8682.0   |  8619.0   -- 8642.0   | 8700.0   -- 8725.0   |   air | |ang| | [4]_ |
+-----------+-----------------------+-----------------------+----------------------+-------+-------+------+
| Pa17      |  8461.0   -- 8474.0   |  8474.0   -- 8484.0   | 8563.0   -- 8577.0   |   air | |ang| | [4]_ |
+-----------+-----------------------+-----------------------+----------------------+-------+-------+------+
| Pa14      |  8577.0   -- 8619.0   |  8563.0   -- 8577.0   | 8619.0   -- 8642.0   |   air | |ang| | [4]_ |
+-----------+-----------------------+-----------------------+----------------------+-------+-------+------+
| Pa12      |  8730.0   -- 8772.0   |  8700.0   -- 8725.0   | 8776.0   -- 8792.0   |   air | |ang| | [4]_ |
+-----------+-----------------------+-----------------------+----------------------+-------+-------+------+
| MgICvD    |  5165.0   -- 5220.0   |  5125.0   -- 5165.0   | 5220.0   -- 5260.0   |   vac | |ang| | [5]_ |
+-----------+-----------------------+-----------------------+----------------------+-------+-------+------+
| NaICvD    |  8177.0   -- 8205.0   |  8170.0   -- 8177.0   | 8205.0   -- 8215.0   |   vac | |ang| | [5]_ |
+-----------+-----------------------+-----------------------+----------------------+-------+-------+------+
| MgIIR     |  8801.9   -- 8816.9   |  8777.4   -- 8789.4   | 8847.4   -- 8857.4   |   vac | |ang| | [5]_ |
+-----------+-----------------------+-----------------------+----------------------+-------+-------+------+
| FeHCvD    |  9905.0   -- 9935.0   |  9855.0   -- 9880.0   | 9940.0   -- 9970.0   |   vac | |ang| | [5]_ |
+-----------+-----------------------+-----------------------+----------------------+-------+-------+------+
| NaI       |  8168.500 -- 8234.125 |  8150.000 -- 8168.400 | 8235.250 -- 8250.000 |   air | |ang| | [6]_ |
+-----------+-----------------------+-----------------------+----------------------+-------+-------+------+
| bTiO      |  4758.500 -- 4800.000 |  4742.750 -- 4756.500 | 4827.875 -- 4847.875 |   air |   mag | [7]_ |
+-----------+-----------------------+-----------------------+----------------------+-------+-------+------+
| aTiO      |  5445.000 -- 5600.000 |  5420.000 -- 5442.000 | 5630.000 -- 5655.000 |   air |   mag | [7]_ |
+-----------+-----------------------+-----------------------+----------------------+-------+-------+------+
| CaH1      |  6357.500 -- 6401.750 |  6342.125 -- 6356.500 | 6408.500 -- 6429.750 |   air |   mag | [7]_ |
+-----------+-----------------------+-----------------------+----------------------+-------+-------+------+
| CaH2      |  6775.000 -- 6900.000 |  6510.000 -- 6539.250 | 7017.000 -- 7064.000 |   air |   mag | [7]_ |
+-----------+-----------------------+-----------------------+----------------------+-------+-------+------+
| NaISDSS   |  8180.0   -- 8200.0   |  8143.0   -- 8153.0   | 8233.0   -- 8244.0   |   air | |ang| | [8]_ |
+-----------+-----------------------+-----------------------+----------------------+-------+-------+------+
| TiO2SDSS  |  6189.625 -- 6272.125 |  6066.625 -- 6141.625 | 6422.0   -- 6455.0   |   air |   mag | [8]_ |
+-----------+-----------------------+-----------------------+----------------------+-------+-------+------+

Input Data Format
~~~~~~~~~~~~~~~~~

The parameters that define the absorption-line-index calculations are
provided via the
:class:`~mangadap.par.absorptionindexdb.AbsorptionIndexDB` object,
which is built using an `SDSS-style parameter file`_. The core level
class that calculates the raw absorption-line indices is
:class:`~mangadap.proc.spectralindices.AbsorptionLineIndices`.

The columns of the parameter file are:

+---------------+-----------+-----------------------------------------------------------------------+
| Parameter     | Format    | Description                                                           |
+===============+===========+=======================================================================+
| ``index``     | int       | Unique integer identifier of the absorption-line index. **Must** be   |
|               |           | unique.                                                               |
+---------------+-----------+-----------------------------------------------------------------------+
| ``name``      | str       | Name of the index.  **Must** be unique.                               |
+---------------+-----------+-----------------------------------------------------------------------+
| ``primary``   | float[2]  | A two-element vector with the starting and ending wavelength for the  |
|               |           | primary passband surrounding the absorption feature(s).               |
+---------------+-----------+-----------------------------------------------------------------------+
| ``blueside``  | float[2]  | A two-element vector with the starting and ending wavelength for a    |
|               |           | passband to the blue of the primary band.                             |
+---------------+-----------+-----------------------------------------------------------------------+
| ``redside``   | float[2]  | A two-element vector with the starting and ending wavelength for a    |
|               |           | passband to the red of the primary band.                              |
+---------------+-----------+-----------------------------------------------------------------------+
| ``waveref``   | str       | The reference frame of the wavelengths; must be either 'air' for air  |
|               |           | or 'vac' for vacuum.                                                  |
+---------------+-----------+-----------------------------------------------------------------------+
| ``units``     | str       | Units for the absorption index, which must be either 'ang' or 'mag'.  |
+---------------+-----------+-----------------------------------------------------------------------+
| ``component`` | int       | **Never used**: Binary flag (0-false,1-true) that the index is a      |
|               |           | component of a composite index.  If true (1), all components with the |
|               |           | same NAME are added together to form the composite index.             |
+---------------+-----------+-----------------------------------------------------------------------+

and an example file might look like this:

.. code-block:: c

    typedef struct {
        int index;
        char name[9];
        double primary[2];
        double blueside[2];
        double redside[2];
        char waveref[3];
        char units[3];
        int component;
    } DAPABI;

    DAPABI  1  CN1        { 4142.125  4177.125 }  { 4080.125  4117.625 }  { 4244.125  4284.125 }  air  mag  0
    DAPABI  2  CN2        { 4142.125  4177.125 }  { 4083.875  4096.375 }  { 4244.125  4284.125 }  air  mag  0

Note that the functionality implied by the ``component`` parameter has
been a notional future development for the module, but has never been
implemented.  However, unfortunately, it's still a required element of
the database file.

Changing the absorption-line index parameters
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The absorption-line indices are measured by
:class:`~mangadap.proc.spectralindices.SpectralIndices`; see
:ref:`spectral-index-measurements`.  A set of parameter files that
define a list of absorption-line index sets are provided with the DAP
source distribution and located at
``$MANGADAP_DIR/mangadap/data/absorption_indices``.  There are a few
methods that you can use to change the set of absorption-line index parameters
used by :class:`~mangadap.proc.spectralindices.SpectralIndices`:

    #. To use one of the existing parameter databases, you can change
       the ``absorption_indices`` keyword in the
       :class:`~mangadap.proc.spectralindices.SpectralIndices`
       configuration file. The keyword should be the capitalized root
       of the parameter filename. E.g., to use
       ``$MANGADAP_DIR/mangadap/data/absorption_indices/lickindx.par``,
       set the keyword to ``LICKINDX``.

    #. To use a *new* parameter database, write the file and save it
       in the ``$MANGADAP_DIR/mangadap/data/absorption_indices/``
       directory, and then change the relevant configuration file in
       the same way as described above.

----

.. _spectralindices-bandhead:

Bandhead or Color Indices
-------------------------

Calculation
~~~~~~~~~~~

.. include:: include/bhdindices.rst

See :class:`~mangadap.proc.spectralindices.BandheadIndices`.

Index Parameters
~~~~~~~~~~~~~~~~

The table below provides a compilation of bandhead and color indices.
Recent survey-level runs of the DAP have included all of these
measurements.

+--------+-----------------------+----------------------+-------+-------------------+-------+-------+
| Name   | Blue Sideband (|ang|) | Red Sideband (|ang|) | Frame |         Integrand | Order |   Ref |
+========+=======================+======================+=======+===================+=======+=======+
| D4000  |          3750 -- 3950 |         4050 -- 4250 |   air | :math:`F_\nu`     |   R/B |  [9]_ |
+--------+-----------------------+----------------------+-------+-------------------+-------+-------+
| Dn4000 |          3850 -- 3950 |         4000 -- 4100 |   air | :math:`F_\nu`     |   R/B | [10]_ |
+--------+-----------------------+----------------------+-------+-------------------+-------+-------+
| TiOCvD |          8835 -- 8855 |         8870 -- 8890 |   vac | :math:`F_\lambda` |   B/R |  [5]_ |
+--------+-----------------------+----------------------+-------+-------------------+-------+-------+

Input Data Format
~~~~~~~~~~~~~~~~~

The parameters that define the bandhead index calculations are provided
via the :class:`~mangadap.par.bandheadindexdb.BandheadIndexDB` object,
which is built using an `SDSS-style parameter file`_.  The core level
class that calculates the raw bandhead indices is
:class:`~mangadap.proc.spectralindices.BandheadIndices`.

The columns of the parameter file are:

+---------------+-----------+-----------------------------------------------------------------------+
| Parameter     | Format    | Description                                                           |
+===============+===========+=======================================================================+
| ``index``     | int       | Unique integer identifier of the absorption-line index. **Must** be   |
|               |           | unique.                                                               |
+---------------+-----------+-----------------------------------------------------------------------+
| ``name``      | str       | Name of the index.  **Must** be unique.                               |
+---------------+-----------+-----------------------------------------------------------------------+
| ``blueside``  | float[2]  | A two-element vector with the starting and ending wavelength for a    |
|               |           | passband to the blue of the primary band.                             |
+---------------+-----------+-----------------------------------------------------------------------+
| ``redside``   | float[2]  | A two-element vector with the starting and ending wavelength for a    |
|               |           | passband to the red of the primary band.                              |
+---------------+-----------+-----------------------------------------------------------------------+
| ``waveref``   | str       | The reference frame of the wavelengths; must be either 'air' for air  |
|               |           | or 'vac' for vacuum.                                                  |
+---------------+-----------+-----------------------------------------------------------------------+
| ``integrand`` | str       | Integrand within the passband for the construction of the index,      |
|               |           | which must be either 'fnu' or 'flambda'.                              |
+---------------+-----------+-----------------------------------------------------------------------+
| ``order``     | str       | Define the order to use when constructing the index.  The options are |
|               |           | either a ratio of red-to-blue or blue-to-red, which are respectively  |
|               |           | selected using 'r_b' or 'b_r'.                                        |
+---------------+-----------+-----------------------------------------------------------------------+

and an example file might look like this:

.. code-block:: c

    typedef struct {
        int index;
        char name[9];
        double blueside[2];
        double redside[2];
        char waveref[3];
        char integrand[7];
        char order[3];
    } DAPBHI;

    DAPBHI  1  D4000      { 3750.000  3950.000 }  { 4050.000  4250.000 }  air      fnu  r_b
    DAPBHI  2  Dn4000     { 3850.000  3950.000 }  { 4000.000  4100.000 }  air      fnu  r_b
    DAPBHI  3  TiOCvD     { 8835.000  8855.000 }  { 8870.000  8890.000 }  vac  flambda  b_r

Changing the bandhead index parameters
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The bandhead and color indices are measured by
:class:`~mangadap.proc.spectralindices.SpectralIndices`; see
:ref:`spectral-index-measurements`. A set of parameter files that
define a list of bandhead index sets are provided with the DAP source
distribution and located at
``$MANGADAP_DIR/mangadap/data/bandhead_indices``. There are a few
methods that you can use to change the set of bandhead-index
parameters used by
:class:`~mangadap.proc.spectralindices.SpectralIndices`:

    #. To use one of the existing parameter databases, you can change
       the ``bandhead_indices`` keyword in the
       :class:`~mangadap.proc.spectralindices.SpectralIndices`
       configuration file. The keyword should be the capitalized root
       of the parameter filename. E.g., to use
       ``$MANGADAP_DIR/mangadap/data/bandhead_indices/bhbasic.par``,
       set the keyword to ``BHBASIC``.

    #. To use a *new* parameter database, write the file and save it
       in the ``$MANGADAP_DIR/mangadap/data/bandhead_indices/``
       directory, and then change the relevant configuration file in
       the same way as described above.

----

.. [1] `Trager et al. (1998, ApJS, 116, 1) <https://ui.adsabs.harvard.edu/abs/1998ApJS..116....1T/abstract>`_
.. [2] `Worthey & Ottaviani (1997, ApJS, 111, 377) <https://ui.adsabs.harvard.edu/abs/1997ApJS..111..377W/abstract>`_
.. [3] `Serven et al. (2005, ApJ, 627, 754) <https://ui.adsabs.harvard.edu/abs/2005ApJ...627..754S/abstract>`_
.. [4] `Cenarro et al. (2001, MNRAS, 326, 959) <https://ui.adsabs.harvard.edu/abs/2001MNRAS.326..959C/abstract>`_; however, note that each index is considered seperately, which is *not* the definition provided in Table 4 of their paper.
.. [5] `Conroy & van Dokkum (2012, ApJ, 747, 69) <https://ui.adsabs.harvard.edu/abs/2012ApJ...747...69C/abstract>`_
.. [6] `Spiniello et al. (2012, ApJL, 753, 32) <https://ui.adsabs.harvard.edu/abs/2012ApJ...753L..32S/abstract>`_
.. [7] `Spiniello et al. (2014, MNRAS, 438, 1483) <https://ui.adsabs.harvard.edu/abs/2014MNRAS.438.1483S/abstract>`_
.. [8] `La Barbera et al. (2013, MNRAS, 433, 3017) <https://ui.adsabs.harvard.edu/abs/2013MNRAS.433.3017L/abstract>`_
.. [9] `Bruzual (1983, ApJ, 273, 105) <https://ui.adsabs.harvard.edu/abs/1983ApJ...273..105B/abstract>`_
.. [10] `Balogh et al. (1999, ApJ, 527, 54) <https://ui.adsabs.harvard.edu/abs/1999ApJ...527...54B/abstract>`_


