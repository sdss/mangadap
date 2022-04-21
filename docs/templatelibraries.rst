
.. include:: include/links.rst

.. _templatelibraries:

Spectral Template Libraries
===========================

The two main full-spectrum-fitting modules of the DAP, the stellar
kinematics module
(:class:`~mangadap.proc.stellarcontinuummodel.StellarContinuumModel`)
and the emission-line module
(:class:`~mangadap.proc.emissionlinemodel.EmissionLineModel`), use
template libraries to fit the stellar continuum. We describe below
the template libraries included in the DAP repository and the format
needed to add new libraries.

Included Libraries
------------------

The following spectral template libraries are included with the DAP distribution:

+-----------------+------------+---------+-------------+-------------------------------------------------------------------------------------------------------------------+
| Key             |  Reference |    Type |    Pedigree | Comments                                                                                                          |
+=================+============+=========+=============+===================================================================================================================+
| BC03            |  [1]_ [2]_ |     SPS |   Empirical | `README <https://github.com/sdss/mangadap/blob/master/mangadap/data/spectral_templates/bc03_mpajhu/README>`__     |
+-----------------+------------+---------+-------------+-------------------------------------------------------------------------------------------------------------------+
| BPASS           |       [3]_ |     SPS | Theoretical | `README <https://github.com/sdss/mangadap/blob/master/mangadap/data/spectral_templates/bpass/README>`__           |
+-----------------+------------+---------+-------------+-------------------------------------------------------------------------------------------------------------------+
| M11ELODIE       |       [4]_ |     SPS |   Empirical | `README <https://github.com/sdss/mangadap/blob/master/mangadap/data/spectral_templates/m11_elodie/README>`__      |
+-----------------+------------+---------+-------------+-------------------------------------------------------------------------------------------------------------------+
| M11MARCS        |       [4]_ |     SPS | Theoretical | `README <https://github.com/sdss/mangadap/blob/master/mangadap/data/spectral_templates/m11_marcs/README>`__       |
+-----------------+------------+---------+-------------+-------------------------------------------------------------------------------------------------------------------+
| M11MILES        |       [4]_ |     SPS |   Empirical | `README <https://github.com/sdss/mangadap/blob/master/mangadap/data/spectral_templates/m11_miles/README>`__       |
+-----------------+------------+---------+-------------+-------------------------------------------------------------------------------------------------------------------+
| M11STELIB       |       [4]_ |     SPS |   Empirical | `README <https://github.com/sdss/mangadap/blob/master/mangadap/data/spectral_templates/m11_stelib/README>`__      |
+-----------------+------------+---------+-------------+-------------------------------------------------------------------------------------------------------------------+
| M11STELIBZSOL   |       [4]_ |     SPS |   Empirical | `README <https://github.com/sdss/mangadap/blob/master/mangadap/data/spectral_templates/m11_stelib_zsol/README>`__ |
+-----------------+------------+---------+-------------+-------------------------------------------------------------------------------------------------------------------+
| MASTARHC        |       [5]_ | Stellar |   Empirical | `README <https://github.com/sdss/mangadap/blob/master/mangadap/data/spectral_templates/mastarhc/README>`__        |
+-----------------+------------+---------+-------------+-------------------------------------------------------------------------------------------------------------------+
| MASTARHC2       |       [5]_ | Stellar |   Empirical | `README <https://github.com/sdss/mangadap/blob/master/mangadap/data/spectral_templates/mastarhc_v2/README>`__     |
+-----------------+------------+---------+-------------+-------------------------------------------------------------------------------------------------------------------+
| MASTARSSP       |      [10]_ |     SPS |   Empirical | `README <https://github.com/sdss/mangadap/blob/master/mangadap/data/spectral_templates/mastar_ssp_v1.0/README>`__ |
+-----------------+------------+---------+-------------+-------------------------------------------------------------------------------------------------------------------+
| MILES           |       [6]_ | Stellar |   Empirical | `README <https://github.com/sdss/mangadap/blob/master/mangadap/data/spectral_templates/miles/README>`__           |
+-----------------+------------+---------+-------------+-------------------------------------------------------------------------------------------------------------------+
| MILESAVG        |       [6]_ | Stellar |   Empirical | `README <https://github.com/sdss/mangadap/blob/master/mangadap/data/spectral_templates/miles_avg/README>`__       |
+-----------------+------------+---------+-------------+-------------------------------------------------------------------------------------------------------------------+
| MILESHC         |       [7]_ | Stellar |   Empirical | `README <https://github.com/sdss/mangadap/blob/master/mangadap/data/spectral_templates/miles_cluster/README>`__   |
+-----------------+------------+---------+-------------+-------------------------------------------------------------------------------------------------------------------+
| MILESTHIN       |       [6]_ | Stellar |   Empirical | `README <https://github.com/sdss/mangadap/blob/master/mangadap/data/spectral_templates/miles_thin/README>`__      |
+-----------------+------------+---------+-------------+-------------------------------------------------------------------------------------------------------------------+
| MIUSCAT         |       [8]_ |     SPS |   Empirical | `README <https://github.com/sdss/mangadap/blob/master/mangadap/data/spectral_templates/miuscat/README>`__         |
+-----------------+------------+---------+-------------+-------------------------------------------------------------------------------------------------------------------+
| MIUSCATTHIN     |       [8]_ |     SPS |   Empirical | `README <https://github.com/sdss/mangadap/blob/master/mangadap/data/spectral_templates/miuscat_thin/README>`__    |
+-----------------+------------+---------+-------------+-------------------------------------------------------------------------------------------------------------------+
| STELIB          |       [9]_ | Stellar |   Empirical | `README <https://github.com/sdss/mangadap/blob/master/mangadap/data/spectral_templates/stelib/README>`__          |
+-----------------+------------+---------+-------------+-------------------------------------------------------------------------------------------------------------------+

Template Library Datamodel
--------------------------

.. _templatelibraries-input:

Input Data Format
~~~~~~~~~~~~~~~~~

The primary constraint on the format of the spectra to be used as
templates is that they are read using
:func:`~mangadap.util.fileio.read_template_spectrum`.  That is:

    * Each spectrum must be in its own fits file.
    * The flux data must be in the first extension.
    * The wavelength vector is constructed from the header WCS (see
      :func:`~mangadap.util.fileio.wavelength_vector`). The wavelength
      range of each spectrum need not be the same.
    * The inverse variance or wavelength dependent spectral
      resolution (:math:`R=\lambda/\delta\lambda`) can be optionally
      provided. If the spectral resolution is not provided, it should
      be defined by the :ref:`templatelibraries-definition`.

.. _templatelibraries-definition:

Template Library Definition
~~~~~~~~~~~~~~~~~~~~~~~~~~~

Template libraries are defined for use in the DAP using the
configuration files in
``$MANGADAP_DIR/mangadap/config/spectral_templates``. These
configuration files are parsed by
:func:`~mangadap.proc.templatelibrary.available_template_libraries`,
which provides a list of
:class:`~mangadap.proc.templatelibrary.TemplateLibraryDef` instances
that can be selected using a keyword.

The critical components to the definition of the template library are:

    * ``file_search``: the search string used by `glob.glob`_ to find
      the fits files with the template library spectra and
    * ``fwhm`` or ``sres_ext``: the FWHM of the (Gaussian)
      line-spread function in angstroms or the name of the extension
      in each template fits file with the spectral resolution,
      :math:`R=\lambda/\delta\lambda`.

Output
~~~~~~

When instantiating a template library, the processed library will be
written to disk, depending on the ``hardcopy`` argument; by default
this argument is True. If the template library has already been
processed and written to disk, the instantiation of the object will
skip processing the library and just read the result of the previous
instantiation.

The path (``directory_path``) and name of the file
(``processed_file``) can be defined upon instantiating the object. If
not provided, the default path is set by
:func:`~mangadap.config.defaults.dap_common_path`, and the output
file name is set by :func:`~mangadap.config.defaults.dap_file_name`.

The format of the output file is:

+-----+-------------+----------------------------------------------------------------+
| HDU | Name        | Description                                                    |
+=====+=============+================================================================+
|   0 | ``PRIMARY`` | Empty                                                          |
+-----+-------------+----------------------------------------------------------------+
|   1 | ``WAVE``    | Wavelength vector                                              |
+-----+-------------+----------------------------------------------------------------+
|   2 | ``FLUX``    | Flux array                                                     |
+-----+-------------+----------------------------------------------------------------+
|   3 | ``MASK``    | Bitmask values; see                                            |
|     |             | :class:`~mangadap.proc.templatelibrary.TemplateLibraryBitMask` |
+-----+-------------+----------------------------------------------------------------+
|   4 | ``SPECRES`` | Spectral resolution                                            |
+-----+-------------+----------------------------------------------------------------+
|   5 | ``SIGOFF``  | If needed, the offset in km/s between the target               |
|     |             | resolution of the spectra and the actual resolution            |
|     |             | acheived.  See                                                 |
|     |             | :func:`~mangadap.util.resolution.match_spectral_resolution`.   |
+-----+-------------+----------------------------------------------------------------+

Adding new template libraries
-----------------------------

Adding new template libraries is relatively straight-forward. First,
make sure the templates adhere to the :ref:`templatelibraries-input`.
Then, you can either use the library by defining it programmatically,
or by adding a new configuration file to the
``$MANGADAP_DIR/mangadap/config/spectral_templates`` directory. See
the TemplateLibrary :ref:`templatelibrary-usage` and the
:ref:`templatelibrary-dev-example` provided by the DAP
:ref:`development`.

----

.. [1] `Bruzual & Charlot (2003, MNRAS, 344, 1000) <https://ui.adsabs.harvard.edu/abs/2003MNRAS.344.1000B/abstract>`_
.. [2] https://wwwmpa.mpa-garching.mpg.de/SDSS/DR4/
.. [3] `Eldridge et al. (2017, PASA, 34, 58) <https://ui.adsabs.harvard.edu/abs/2017PASA...34...58E/abstract>`_
.. [4] `Maraston & Strömbäch (2011, MNRAS, 418, 2785) <https://ui.adsabs.harvard.edu/abs/2011MNRAS.418.2785M/abstract>`_
.. [5] `Yan et al. (2019, ApJ, 883, 175) <https://ui.adsabs.harvard.edu/abs/2019ApJ...883..175Y/abstract>`_
.. [6] `Falcón-Barroso et al. (2011, A&A, 532, 95) <https://ui.adsabs.harvard.edu/abs/2011A%26A...532A..95F/abstract>`_
.. [7] `Westfall et al. (2019, AJ, 158, 231) <https://ui.adsabs.harvard.edu/abs/2019AJ....158..231W/abstract>`_
.. [8] `Vazdekis et al. (2012, MNRAS, 424, 157) <https://ui.adsabs.harvard.edu/abs/2012MNRAS.424..157V/abstract>`_
.. [9] `Le Borgne et al. (2003, A&A, 402, 433) <https://ui.adsabs.harvard.edu/abs/2003A%26A...402..433L/abstract>`_
.. [10] `Maraston et al. (2020, MNRAS, 496, 2962)`_

