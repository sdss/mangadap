
.. include:: include/links.rst

.. _execution:

Execution
=========

The MaNGA DAP is, of course, geared toward working at the
survey-level for data produced by the MaNGA data-reduction pipeline
(DRP). However, we have begun to generalize the usage of the DAP for
datacubes from other instruments.

What follows below describes both how to use the current MaNGA DAP to
analyze a datacube from the MaNGA survey, as well as how to execute
the survey-level batch mode that analyzes all the MaNGA datacubes
within a given directory structure.

Examples are available (or will be soon) that demonstrate
:ref:`fitonespec` and :ref:`fitonecube` using the DAP software, for
both MaNGA and other integral-field instruments.

Input files
-----------

.. _execution-analysis-plan:

DAP AnalysisPlan
~~~~~~~~~~~~~~~~

The DAP uses an `SDSS-style parameter file`_ to define one or more
methods to use when analyzing any given MaNGA datacube. Each method,
or "analysis plan", is defined by a set of six keywords that identify
the method to use for each of the DAP's six main
:ref:`workflow-analysis-modules`. The ``AnalysisPlan`` parameter file
allows you to analyze the same datacube multiple ways in a single
execution of the DAP; however, note that this is no different than
executing the DAP once per analysis method. The ``AnalysisPlan`` is
required to execute the DAP; however, the DAP provides a
:ref:`execution-analysis-plan-default` that will be used if no file
is provided.

An example ``AnalysisPlan`` parameter file looks like this (this is
exactly the file used for MPL-10):

.. code-block:: c

    typedef struct {
        char drpqa_key[8];
        int drpqa_clobber;
        char bin_key[8];
        int bin_clobber;
        char continuum_key[8];
        int continuum_clobber;
        char elmom_key[8];
        int elmom_clobber;
        char elfit_key[8];
        int elfit_clobber;
        char spindex_key[8];
        int spindex_clobber;
    } DAPPLAN;

    #           DRP QA    BINNING        CONTINUUM      MOMENTS       LINEPROF   SPECINDEX
    #        ---------  ---------  ---------------  -----------  -------------  ----------
    DAPPLAN    SNRG  0     SPX  0   MILESHCMPL10 0  EMOMMPL10 0    EFITMPL10 0   INDXEN  0
    DAPPLAN    SNRG  0   VOR10  0   MILESHCMPL10 0  EMOMMPL10 0    EFITMPL10 0   INDXEN  0
    DAPPLAN    SNRG  0   HYB10  0   MILESHCMPL10 0  EMOMMPL10 0  EFITMPL10DB 0   INDXEN  0

The configuration of each module is set by two values: a configuration
key and a flag indicating if any existing results should be overwritten
(0 = False, 1 = True).  Each keyword points to a unique configuration
file, and each of these keywords produce a unique instance of a
parameter set that drives the relevant analysis module.  The link
between the relevant structure element root (e.g., ``drpqa``), the
location of its configuration files, and the object used to read the
configurations are as follows:

+-----------+-------------------------------------+--------------------------------------------------------------------------+
|    struct |  ``$MANGADAP_DIR/mangadap/config/`` |                                                               DAP object |
+===========+=====================================+==========================================================================+
|     drpqa |           ``reduction_assessments`` |      :class:`~mangadap.proc.reductionassessments.ReductionAssessmentDef` |
+-----------+-------------------------------------+--------------------------------------------------------------------------+
|       bin |                 ``spatial_binning`` | :class:`~mangadap.proc.spatiallybinnedspectra.SpatiallyBinnedSpectraDef` |
+-----------+-------------------------------------+--------------------------------------------------------------------------+
| continuum |      ``stellar_continuum_modeling`` |   :class:`~mangadap.proc.stellarcontinuummodel.StellarContinuumModelDef` |
+-----------+-------------------------------------+--------------------------------------------------------------------------+
|     elmom |           ``emission_line_moments`` |       :class:`~mangadap.proc.emissionlinemoments.EmissionLineMomentsDef` |
+-----------+-------------------------------------+--------------------------------------------------------------------------+
|     elfit |          ``emission_line_modeling`` |           :class:`~mangadap.proc.emissionlinemodel.EmissionLineModelDef` |
+-----------+-------------------------------------+--------------------------------------------------------------------------+
|   spindex |                ``spectral_indices`` |               :class:`~mangadap.proc.spectralindices.SpectralIndicesDef` |
+-----------+-------------------------------------+--------------------------------------------------------------------------+

When setting the keyword values in the analysis-plan file, they must be
recognized as a method defined in the relevant configuration directory.
The DAP will only execute properly if at least the first three steps
have valid keywords.  The remaining three modules can be skipped (i.e.,
the emission-line moments, emission-line-model parameters, and spectral
indices are not measured) by setting their keyword to ``None``; the
primary DAP output will still be produced but with empty arrays for
those extensions/channels that would normally be populated by the
skipped analysis steps.

.. _execution-analysis-plan-default:

Default AnalysisPlan
++++++++++++++++++++

If executed without an ``AnalysisPlan`` parameter file, the command-line
execution of the DAP will use a default plan; see
:func:`mangadap.par.analysisplan.AnalysisPlanSet.default`.

The current default plan uses the following keys:

.. include:: tables/default_analysisplan.rst

.. _execution-config:

The DAP Datacube Configuration File
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The DAP uses a configuration (``ini``) file to set the datacube to be
analyzed and provides some relevant metadata. These configuration
files are generated at the survey-level by
:func:`~mangadap.survey.drpcomplete.write_config`. However, we also
provide the ``$MANGADAP_DIR/bin/write_dap_config`` script that will
generate the relevant configuration files if you have the DRPall or
DRPComplete file. As with all the DAP scripts, you can use the ``-h``
command-line option to get the usage:

.. include:: help/write_dap_config.rst

To construct the configuration file for datacube 7815-3702, run:

.. code-block:: console

    write_dap_config 7815 3702 mangadap-7815-3702.ini -a drpall-v3_0_1.fits

to produce:

.. code-block:: ini

    # Auto-generated configuration file
    # Fri 28 Feb 2020 16:57:19

    [default]
    drpver
    redux_path
    directory_path
    plate = 7815
    ifu = 3702
    log = True
    sres_ext
    sres_fill
    covar_ext
    z = 2.9382300e-02
    vdisp
    ell = 1.1084400e-01
    pa = 1.6324500e+02
    reff = 3.7749500e+00

Use the relevant keywords to change the paths or the extensions used
for the spectral resolution and spatial correlation matrix (e.g.,
``GCORREL``).

.. note::

    If the isophotal ellipticity, ``ell``, is such that
    :math:`\epsilon < 0` or :math:`\epsilon > 1`, the DAP will adopt
    a default value of 0. The DAP accepts any value for the position
    angle, ``pa``, but imposes periodic limits (e.g., 380 deg is
    converted to 20 deg). If the effective radius, ``reff``, is
    :math:`R_{\rm eff} < 0` or undefined, the DAP uses :math:`R_{\rm
    eff} = 1`.

.. _execution-mangadap:

DAP command-line script
-----------------------

The main DAP script is ``$MANGADAP_DIR/bin/manga_dap``, which is a
simple wrapper of :func:`~mangadap.survey.manga_dap.manga_dap`.  With the
DAP installed, you can call the script directly from the command line:

.. include:: help/manga_dap.rst

The DAP allows you to define your own datacube class, as long as it
is derived from :class:`~mangadap.datacube.datacube.DataCube`. You can
then specify that your data should be instantiated with that derived
class using the ``cube_module`` and ``cube_object`` arguments; these
default to ``mangadap.datacube`` and ``MaNGADataCube``, respectively.

When running the DAP on a MaNGA datacube, you have to provide a
configuration file; however, for derived classes, you may be able to
fully instantiate the relevant data using the datacube file, which is
why we've provided the ``-f`` option.

Note that the analysis plan file is an *optional* argument. If it is
not given, the DAP will use the
:ref:`execution-analysis-plan-default`.

To run the DAP on a single datacube using the default analysis plan,
and assuming you have the DRPall file and the relevant LOGCUBE and
LOGRSS files (see the warning below) in the current directory, you
could execute the DAP as follows:

.. code-block:: console

    write_dap_config 7815 3702 mangadap-7815-3702.ini -a drpall-v3_0_1.fits
    manga_dap -c mangadap-7815-3702.ini -vv -log mangadap-7815-3702.log -d . -a dap_output

This will analyze the datacube for observation 7815-3702 using the
default analysis plan, with verbose output and a log written to
``mangadap-7815-3702.log``, and with the root directory for all the
DAP output (except for the log) set to ``dap_output``.

.. warning::

    When running the DAP on a MaNGA datacube, you should have both
    the DRP ``LOGRSS`` and ``LOGCUBE`` files in the
    ``directory_path`` if you want to account for the
    :ref:`spatialcovariance`! If the ``LOGRSS`` files are not
    present, the DAP will throw a warning and continue, which means
    that the warning can get buried among all the other messages and
    likely missed.

Programmatic execution
----------------------

Alternatively, ``$MANGADAP_DIR/examples/fit_one_cube.py`` (see
:ref:`fitonecube`) provides a programmatic approach to running the
exact same script that is executed by the ``manga_dap`` command-line
script. The code provides a way to generate :ref:`execution-config`
directly from the DRPall file, instead of from a file, and it
directly defines the ``AnalysisPlan`` object with a hard-coded set of
keywords. Using this script as an example, one could construct a
script that programmatically analyzes a large set of MaNGA datacubes.

.. _execution-rundap:

Batch execution using automatically generated scripts
-----------------------------------------------------

The survey-level execution of the DAP uses the
``$MANGADAP_DIR/bin/rundap`` script, which is a simple wrapper of
:class:`~mangadap.scripts.rundap.rundap`.  This script

 * sets up the DAP output directory structure
 * either confirms that a provided list of datacubes to analyze exist on
   disk or trolls the DRP directory structure to find all or some subset
   of available datacubes to analyze
 * creates :ref:`execution-config` for each ``plateifu`` to be analyzed,
 * creates a script file for each ``plateifu`` that can be sourced to
   execute the DAP and the associated QA plots,
 * creates scripts that execute the plate-level QA plots,
 * creates scripts that build the DAPall file and its QA plots, and
 * submits the scripts to the Utah cluster.

The last step uses an SDSS python package called ``pbs``, which isn't
required for the more general-purpose use of the ``rundap`` script
discussed here.  With the DAP installed, you can call the script
directly from the command line:

.. include:: help/rundap.rst

If a :ref:`execution-analysis-plan` is not provided, the scripts will
use the :ref:`execution-analysis-plan-default`.

An example call of this script that will only construct scripts for
the analysis of observation 7443-12701 using the default AnalysisPlan
is:

.. code-block:: console
    
    rundap --platelist 7443 --ifudesignlist 12701 --redux_path /path/with/drp/output/ --analysis_path /path/for/dap/output/ -vv --log

In this call, I've specified that the DRP data is in
``/path/with/drp/output/`` and that the DAP output should be placed
in ``/path/for/dap/output/`` instead of using the default
:ref:`datamodel-directory-structure`.

The script file this call produces is written to
``/path/for/dap/output/log/[time]/7495/12704/mangadap-7495-12704``,
where ``[time]`` is a time stamp of when ``rundap`` was executed. (If
you execute ``rundap`` multiple times, it will create new directories
using new time stamps each time.) The lines of the script file for
each plate-ifu:

 - touches the ``*.started`` file
 - executes ``manga_dap``
 - executes a series of QA plotting scripts
 - touches the ``*.done`` file 

The example script generated by the above command would look something like this:

.. code-block:: bash

    # Auto-generated batch file
    # Wed 27 May 2020 11:50:45

    touch /uufs/chpc.utah.edu/common/home/sdss/mangawork/manga/spectro/analysis/v3_0_1/3.0.1/log/27May2020T17.50.29UTC/7443/12701/mangadap-7443-12701.started

    OMP_NUM_THREADS=1 manga_dap -c /uufs/chpc.utah.edu/common/home/sdss/mangawork/manga/spectro/analysis/v3_0_1/3.0.1/common/7443/12701/mangadap-7443-12701-LOGCUBE.ini -r /uufs/chpc.utah.edu/common/home/sdss/mangawork/manga/spectro/redux/v3_0_1 -a /uufs/chpc.utah.edu/common/home/sdss/mangawork/manga/spectro/analysis/v3_0_1/3.0.1 -p /uufs/chpc.utah.edu/common/home/sdss/mangawork/manga/spectro/analysis/v3_0_1/3.0.1/log/27May2020T17.50.29UTC/mpl10_plan.par --log /uufs/chpc.utah.edu/common/home/sdss/mangawork/manga/spectro/analysis/v3_0_1/3.0.1/log/27May2020T17.50.29UTC/7443/12701/mangadap-7443-12701.log -vv

    OMP_NUM_THREADS=1 dap_ppxffit_qa 7443 12701 --analysis_path /uufs/chpc.utah.edu/common/home/sdss/mangawork/manga/spectro/analysis/v3_0_1/3.0.1 --plan_file /uufs/chpc.utah.edu/common/home/sdss/mangawork/manga/spectro/analysis/v3_0_1/3.0.1/log/27May2020T17.50.29UTC/mpl10_plan.par

    OMP_NUM_THREADS=1 spotcheck_dap_maps 7443 12701 --analysis_path /uufs/chpc.utah.edu/common/home/sdss/mangawork/manga/spectro/analysis/v3_0_1/3.0.1 --plan_file /uufs/chpc.utah.edu/common/home/sdss/mangawork/manga/spectro/analysis/v3_0_1/3.0.1/log/27May2020T17.50.29UTC/mpl10_plan.par

    OMP_NUM_THREADS=1 dap_fit_residuals 7443 12701 --analysis_path /uufs/chpc.utah.edu/common/home/sdss/mangawork/manga/spectro/analysis/v3_0_1/3.0.1 --plan_file /uufs/chpc.utah.edu/common/home/sdss/mangawork/manga/spectro/analysis/v3_0_1/3.0.1/log/27May2020T17.50.29UTC/mpl10_plan.par

    touch /uufs/chpc.utah.edu/common/home/sdss/mangawork/manga/spectro/analysis/v3_0_1/3.0.1/log/27May2020T17.50.29UTC/7443/12701/mangadap-7443-12701.done


To execute the script, you would then run:

.. code-block:: bash

    source /path/for/dap/output/v3_0_1/3.0.1/log/27May2020T17.50.29UTC/7443/12701/mangadap-7443-12701

The ``rundap`` script allows you to construct scripts for all datacubes
it can find on disk, all IFUs on a given plate, all combinations of a
set of plate and IFU numbers, or for a specified list of ``plateifu``
IDs.

.. note::

    The ``rundap`` script constructs the
    :class:`~mangadap.survey.drpcomplete.DRPComplete` object and
    writes its associated fits file (see
    :ref:`metadatamodel-drpcomplete`). The data compiled into this
    database can either be drawn from the DRPall file or from the
    plateTargets data in ``mangacore``; the latter is the only reason
    the DAP has ``mangacore`` as a dependency. For general use, you
    should have ``rundap`` use the DRPall file. The use of the
    plateTargets data is only necessary in the rare case when the DAP
    is executed before the relevant DRPall file has been constructed.

To write the post-processing scripts, execute ``rundap`` with the
``--post`` and ``--post_plots`` options.  This produces two additional
types of scripts:

.. todo::

    update the scripts below!


 - Scripts to produce QA plots for all IFUs on a given plate.  This file
   is written to, e.g.,
   ``/path/for/dap/output/log/01Nov2019T16.58.40UTC/7443/7443_fitqa``
   and looks like this:

   .. code-block:: bash

        # Auto-generated batch file
        # Fri 01 Nov 2019 10:58:52

        touch /path/for/dap/output/v2_7_1/2.4.1/log/01Nov2019T16.58.40UTC/7443/7443_fitqa.started

        dap_plate_fit_qa 7443 --analysis_path /path/for/dap/output/v2_7_1/2.4.1

        touch /path/for/dap/output/v2_7_1/2.4.1/log/01Nov2019T16.58.40UTC/7443/7443_fitqa.done

 - A script that builds the :ref:`metadatamodel-dapall` and writes its QA plots.  This
   file is written to, e.g.,
   ``/path/for/dap/output/v2_7_1/2.4.1/log/01Nov2019T16.58.40UTC/build_dapall``
   and looks like this:

   .. code-block:: bash

        # Auto-generated batch file
        # Fri 01 Nov 2019 10:58:52

        touch /path/for/dap/output/v2_7_1/2.4.1/log/01Nov2019T16.58.40UTC/build_dapall.started

        construct_dapall --drpver v2_7_1 -r /path/with/drp/output/v2_7_1 --dapver 2.4.1 -a /path/for/dap/output/v2_7_1/2.4.1 -vv

        dap_dapall_qa --drpver v2_7_1 --redux_path /path/with/drp/output/v2_7_1 --dapver 2.4.1 --analysis_path /path/for/dap/output/v2_7_1/2.4.1

        touch /path/for/dap/output/v2_7_1/2.4.1/log/01Nov2019T16.58.40UTC/build_dapall.done


.. todo::

    Check the details here.

In the automated run of the DAP, any entry in the
:ref:`metadatamodel-drpcomplete` file that meets the criteria set by
:func:`mangadap.survey.drpcomplete.DRPComplete.can_analyze` will be
analyzed. Currently the relevant criteria are:
 
 - ``MANGAID != NULL``
 - ``MANGA_TARGET1 > 0 | MANGA_TARGET3 > 0``
 - ``VEL > -500``

An important consequence of the selection above is that *any targets
without a provided redshift will not be analyzed by the DAP*, unless
it has replacement redshift in the :ref:`metadatamodel-redshift-fix`.
Ancillary targets not analyzed by the DAP are likely because a
redshift was not available.

