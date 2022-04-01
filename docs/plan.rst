
.. include:: include/links.rst

.. _plan:

Analysis Plans
==============

The DAP uses a `toml`_ file (parsed using `tomli`_) to define one or more
"analysis plans," which define how the DAP will analyze the provided datacube.
The `toml`_ file is used to instantiate an
:class:`~mangadap.config.analysisplan.AnalysisPlan` object and parsed using 
:func:`~mangadap.config.analysisplan.AnalysisPlan.from_toml`.

If no `toml`_ file is provided, the DAP uses the default parameters for all
analysis, and the name of the plan is set to ``default``; see
:func:`~mangadap.config.analysisplan.AnalysisPlan.default`.

When providing a `toml`_ file with an altered analysis plan, each plan is
identified by a top-level dictionary keyword.  For example, you could name the
plan ``default`` and provide a keyword identifier for that plan:

.. code-block:: toml

    [default]
     key = 'HYB10-MILESHC-MASTARHC2'

The optional ``key`` keyword specifies the identifier for the plan separately
from the dictionary key.  If no ``key`` entry is provided, the name of the plan
is the same as the top-level dictionary keyword.

If you want to apply more than one analysis plan, you can have multiple plans in
a single `toml`_ file and the DAP will execute them in series.  For example, the
following defines two plans with separate keywords:

.. code-block:: toml

    [plan1]
     key = 'HYB10-MILESHC-MASTARHC2'

    [plan2]
     key = 'HYB10-MILESHC-MASTARSSP'

The nested components of the `toml`_ file need only provide **alterations** to
the default parameters for each analysis module.  For example, if you simply
want to turn off the treatment of covariance in the reduction assessments
module, you would do the following:

.. code-block:: toml

    [default]
     key = 'myplan'

    [default.rdxqa]
     covariance = false

This only alters the single provided parameter and leaves the other parameters
at their default values.

.. warning::

    The parameter sets for each analysis module have an associated keyword.
    These keywords are used in the name of the DAP reference files.  In the
    example above, the user will have changed the ``covariance`` value, but the
    keyword of the ``rdxqa`` plan will have remained unchanged.  Because the DAP
    will not redo analysis steps that it thinks are already complete (indicated
    by the presence of the reference file), this means that the change may not
    take effect if the reference file is present.  Make sure to either change
    the name of the reference file when changing parameters on-the-fly, or set
    the ``overwrite`` keyword to ``true``.  See :ref:`execution`.

Below, we specify how to alter the parameters for all of the main DAP modules
and provide their default values.  In addition to these details, make sure you
understand how the DAP uses the ``key`` for each module and how the modules
understand whether or not they've already been completed; see :ref:`execution`.

.. note::

    For examples of fully defined plan files, see:
        * `mangadap/config/default_plan.toml
          <https://github.com/sdss/mangadap/blob/master/mangadap/config/default_plan.toml>`__:
          The default plan, when not overridden by the user.
        * `mangadap/data/tests/global_bin.toml
          <https://github.com/sdss/mangadap/blob/master/mangadap/data/tests/global_bin.toml>`__:
          A plan used by a DAP unit test to analyze a datacube by binning all
          spaxels into a single spectrum.
        * `mangadap/data/tests/dr17.toml
          <https://github.com/sdss/mangadap/blob/master/mangadap/data/tests/dr17.toml>`__:
          A plan that mimics the survey-level plan used to produce the MaNGA DAP
          products released in DR17.

----

.. contents:: Main Module Parameters
    :depth: 1
    :local:

Basic Reduction Assessments
---------------------------

Parameters guiding the basic assessments of the datacube to be analyzed are set
by the :class:`~mangadap.proc.reductionassessments.ReductionAssessmentDef`
object.  Its parameters and default values are:

.. include:: tables/reductionassessmentdef.rst

To `toml`_ file section used to change these parameters is ``.rdxqa``.  For
example, to change the keyword identifier for the parameter set for the
``default`` plan, your `toml`_ file would include

.. code-block:: toml

    [default.rdxqa]
     key = 'test'

This module has no nested parameter dictionaries.

Once parsed, the paramaters for the reduction assessments are kept in the
``rdxqa`` attribute of the :class:`~mangadap.config.analysisplan.AnalysisPlan`
object; see :ref:`plan-inspect`.

Spatial Binning
---------------

Parameters guiding the spatial binning of the datacube are set
by the :class:`~mangadap.proc.spatiallybinnedspectra.SpatiallyBinnedSpectraDef`
object.  Its parameters and default values are:

.. include:: tables/spatiallybinnedspectradef.rst

The `toml`_ file section used to change these parameters is ``.binning``.  For
example, to change the keyword identifier for the parameter set for the
``default`` plan, your `toml`_ file would include

.. code-block:: toml

    [default.binning]
     key = 'test'

This module has two nested parameter dictionaries, one that defines the method
and parameters used for the spatial binning and one defining the parameters used
for the spectral stacking.

Once parsed, the paramaters for the spatial binning are kept in the ``binning``
attribute of the :class:`~mangadap.config.analysisplan.AnalysisPlan` object; see
:ref:`plan-inspect`.

Spatial Binning Methods
+++++++++++++++++++++++

Selection of the spatial binning method is done using the ``spatial_method``
keyword in the ``.binning`` dictionary on the provided `toml`_ file.  This is
not listed in the table above because it is only used to set the relevant
``binpar``, ``binclass``, and ``binfunc`` components of the
:class:`~mangadap.proc.spatiallybinnedspectra.SpatiallyBinnedSpectraDef` object.
To change any of the spatial binning parameters, change the values in the
``.binning.spatial`` nested dictionary.  For example, to change the S/N
threshold for the Voronoi binning, your `toml`_ file would include:

.. code-block:: toml

    [default.binning.spatial]
     target_snr = 10

The currently available spatial binning methods are:

Global
~~~~~~

Selected using ``spatial_method = 'global'``, and puts all the datacube spectra
into a single bin.  No other parameters are defined.

Radial
~~~~~~

Selected using ``spatial_method = 'radial'``, and bins spectra radially
according to the isophotal properties provided by the datacube metadata; see
:func:`~mangadap.datacube.datacube.DataCube.populate_metadata`.  The available
parameters are:

.. include:: tables/radialbinningpar.rst

Voronoi
~~~~~~~

Selected using ``spatial_method = 'voronoi'``, and bins spectra using `vorbin`_.
The available parameters are:

.. include:: tables/voronoibinningpar.rst

Square
~~~~~~

Selected using ``spatial_method = 'square'``, and rebins spaxels to a more
coarse binning with square spaxels.  The available parameters are:

.. include:: tables/squarebinningpar.rst

Spectral Stacking Methods
+++++++++++++++++++++++++

There is currently only one spectral stacking method, meaning there is no
``spectral_method`` keyword expected by the main ``.binning`` dictionary.  To change
any of the spectral stacking parameters, change the values in the
``.binning.spectral`` nested dictionary.  For example, to change the stacking operation,
your `toml`_ file would include:

.. code-block:: toml

    [default.binning.spectral]
     operation = 'sum'

The available spectral stacking parameters are:

.. include:: tables/spectralstackpar.rst

Stellar Continuum Model (Stellar Kinematics)
--------------------------------------------

Parameters guiding the stellar continuum modeling of the datacube, which is
currently only used to model the stellar kinematics, are set by the
:class:`~mangadap.proc.stellarcontinuummodel.StellarContinuumModelDef` object.
Its parameters and default values are:

.. include:: tables/stellarcontinuummodeldef.rst

The `toml`_ file section used to change these parameters is ``.continuum``.  For
example, to change the keyword identifier for the parameter set for the
``default`` plan, your `toml`_ file would include

.. code-block:: toml

    [default.continuum]
     key = 'test'

Selection of the stellar continuum modeling method is done using the
``fit_method`` keyword in the ``.continuum`` dictionary on the provided `toml`_
file.  Currently, this can only be ``'ppxf'``, which selects the
:class:`~mangadap.proc.ppxffit.PPXFFit` class to perform the modeling.  The
parameters that govern the modeling are set by the
:class:`~mangadap.proc.ppxffit.PPXFFitPar` object.  Its parameters and default
values are:

.. include:: tables/ppxffitpar.rst

Changes to these parameters must be made in the ``.continuum.fit`` section of
the `toml`_ file.  For example, to change the degree of the additive polynomial,
your `toml`_ file should include:

.. code-block:: toml

    [default.continuum.fit]
     degree = 10

Finally, the :class:`~mangadap.proc.ppxffit.PPXFFit` modeling approach requires
a set of spectral templates.  These are set using the
``.continuum.fit.templates`` section of the `toml`_ file; see
:ref:`templatelibraries`.  One can select a template library distributed with
the DAP by specifying *only* the library keyword.  For example, to specify the
MILESHC library, your `toml`_ file would include:

.. code-block:: toml

    [default.continuum.fit.templates]
     key = 'MILESHC'

However, you could also include an user-defined template library by defining all
or most of the required parameters held by the
:class:`~mangadap.proc.templatelibrary.TemplateLibraryDef` object:

.. include:: tables/templatelibrarydef.rst

.. note::

    The templates uses by the stellar-continuum and the emission-line fitting
    modules do *not* need to be the same.

Once parsed, the parameters for the continuum fitting are kept in the
``continuum`` attribute of the :class:`~mangadap.config.analysisplan.AnalysisPlan`
object; see :ref:`plan-inspect`.

Non-parametric Emission-line Properties
---------------------------------------

Parameters guiding non-parametric measurements of the emission-line properties
(using a moment analysis) are set by the
:class:`~mangadap.proc.emissionlinemoments.EmissionLineMomentsDef` object.  Its
parameters and default values are:

.. include:: tables/emissionlinemomentsdef.rst

The `toml`_ file section used to change these parameters is ``.eline_moments``.  For
example, to change the keyword identifier for the parameter set for the
``default`` plan, your `toml`_ file would include

.. code-block:: toml

    [default.eline_moments]
     key = 'test'

Particular attention should be given to the ``passbands`` parameter.  This can
either be a keyword selecting one of the emission-line bandpass databases
distributed with the DAP or the filename of a new SDSS parameter file with a
user-defined database; see :ref:`emissionlines-moments`.

Once parsed, the parameters for the emission-line moment measurements are kept
in the ``elmom`` attribute of the
:class:`~mangadap.config.analysisplan.AnalysisPlan` object; see
:ref:`plan-inspect`.

Gaussian Emission-line Modeling
-------------------------------

Parameters guiding the Gaussian emission-line modeling are set by the
:class:`~mangadap.proc.emissionlinemodel.EmissionLineModelDef` object.
Its parameters and default values are:

.. include:: tables/emissionlinemodeldef.rst

The `toml`_ file section used to change these parameters is ``.eline_fits``.  For
example, to change the keyword identifier for the parameter set for the
``default`` plan, your `toml`_ file would include

.. code-block:: toml

    [default.eline_fits]
     key = 'test'

Selection of the emission-line modeling method is done using the ``fit_method``
keyword in the ``.eline_fits`` dictionary on the provided `toml`_ file.
Currently, this can only be ``'sasuke'``, which selects the
:class:`~mangadap.proc.sasuke.Sasuke` class to perform the modeling.  The
parameters that govern the modeling are set by the
:class:`~mangadap.proc.sasuke.SasukePar` object.  Its parameters and default
values are:

.. include:: tables/sasukepar.rst

Changes to these parameters must be made in the ``.eline_fits.fit`` section of
the `toml`_ file.  For example, to change the degree of the multiplicative polynomial,
your `toml`_ file should include:

.. code-block:: toml

    [default.eline_fits.fit]
     mdegree = 10

Particular attention should be given to the ``emission_lines`` parameter.  This can
either be a keyword selecting one of the emission-line databases
distributed with the DAP or the filename of a new SDSS parameter file with a
user-defined database; see :ref:`emissionlines-modeling`.

Once parsed, the parameters for the emission-line modeling are kept in the
``elfit`` attribute of the :class:`~mangadap.config.analysisplan.AnalysisPlan`
object; see :ref:`plan-inspect`.

Finally, the :class:`~mangadap.proc.sasuke.Sasuke` modeling approach requires
a set of spectral templates.  These are set using the
``.emline_fits.fit.templates`` section of the `toml`_ file; see
:ref:`templatelibraries`.  One can select a template library distributed with
the DAP by specifying *only* the library keyword.  For example, to specify the
MASTARSSP library, your `toml`_ file would include:

.. code-block:: toml

    [default.eline_fits.fit.templates]
     key = 'MASTARSSP'

However, you could also include an user-defined template library by defining all
or most of the required parameters held by the
:class:`~mangadap.proc.templatelibrary.TemplateLibraryDef` object:

.. include:: tables/templatelibrarydef.rst

.. note::

    The templates uses by the stellar-continuum and the emission-line fitting
    modules do *not* need to be the same.

Spectral Indices
----------------

Parameters guiding spectral index measurements are set by the
:class:`~mangadap.proc.spectralindices.SpectralIndicesDef` object.  Its
parameters and default values are:

.. include:: tables/spectralindicesdef.rst

The `toml`_ file section used to change these parameters is ``.indices``.  For
example, to change the keyword identifier for the parameter set for the
``default`` plan, your `toml`_ file would include

.. code-block:: toml

    [default.indices]
     key = 'test'

Particular attention should be given to the ``absindex`` and ``bandhead``
parameters.  These can either be keywords selecting one of the bandpass
databases distributed with the DAP or the filename of a new SDSS parameter file
with a user-defined database; see :ref:`spectralindices`.

Once parsed, the parameters for the spectral index measurements are kept in the
``sindx`` attribute of the :class:`~mangadap.config.analysisplan.AnalysisPlan`
object; see :ref:`plan-inspect`.

----

.. _plan-inspect:

Inspecting the parameters
-------------------------

All adjustments to the parameters used by the DAP should be made via the input
`toml`_ file.  However, if you'd like to directly inspect the parameters being
used for each plan, you can get these by direct interaction with the
:class:`~mangadap.config.analysisplan.AnalysisPlan` object.

For example, you can read and inspect the DR17 analysis plan as follows.
Assuming you've installed the MaNGA DAP source code is in a `mangadap`
directory:

.. code-block:: python

    from mangadap.tests.util import data_test_file
    from mangadap.config.analysisplan import AnalysisPlan
    plan_file = data_test_file('dr17.toml')
    plan = AnalysisPlan.from_toml(plan_file)

You then have access to (most of) the parameters that will be used throughout
the analysis.  For example:

.. code-block:: python

    # How many times will the datacube be analyzed?
    >>> plan.nplans
    4
    # What are the plan keys?
    >>> list(plan.keys())
    dict_keys(['plan1', 'plan2', 'plan3', 'plan4'])
    # What is the keyword of the reduction assessment module for the first plan
    >>> plan.rdxqa['plan1']['key']
    'SNRG'
    # Or you can view the full parameter set
    >>> plan.rdxqa['plan1']
    Parameter           Value                    Default                  Type           Callable
    ---------------------------------------------------------------------------------------------
    key                 SNRG                     SNRG                     str            False
    waverange           None                     None                     ndarray, list  False
    response_func_file  gunn_2001_g_response.db  gunn_2001_g_response.db  str            False
    in_vacuum           True                     True                     bool           False
    covariance          True                     True                     bool           False
    minimum_frac        0.8                      0.8                      int, float     False
    overwrite           False                    False                    bool           False
    # What are the template-library parameters used for the emission-line
    # fitting in the 2nd plan?
    >>> plan.elfit['plan2']['fitpar']['template_library']
    Parameter         Value                      Default               Type           Callable
    ------------------------------------------------------------------------------------------
    key               MASTARSSP                  MILESHC               str            False
    file_search       mastar_ssp_v1.0/*.fits.gz  miles_cluster/*.fits  str            False
    fwhm              2.5                        2.5                   int, float     False
    sres_ext          SRES                       None                  str            False
    in_vacuum         True                       False                 bool           False
    wave_limit        None                       None                  ndarray, list  False
    lower_flux_limit  None                       None                  int, float     False
    log10             True                       False                 bool           False
    # What database is being used to define the absorption line indices?
    >>> plan.sindx['plan1']['absindex']
    'EXTINDX'

.. note::

    Importantly, not all "parameters" are initially available, and some will
    change over the course of the DAP execution.  Some parameters are "filled"
    as the code progresses using any ``fill`` methods (e.g.,
    :func:`~mangadap.proc.sasuke.SasukePar.fill`) associated with the relevant 
    parameter set.  These methods use the provided data to fill in parameters
    that are specific to each datacube or spaxel.  For example, the
    ``'guess_redshift'`` for the emission-line parameters is initially set to
    the default (``plan.elfit['plan1']['fitpar']['guess_redshift']`` is 0.);
    however, the guess redshift is replaced by the first moment estimates from
    the emission-line moment analysis.


