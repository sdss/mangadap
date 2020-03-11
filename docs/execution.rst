
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

Test cases for how to use the DAP to analyze datacubes from other
instruments/surveys will be added soon.

Input files
-----------

.. _execution-analysis-plan:

The DAP AnalysisPlan
~~~~~~~~~~~~~~~~~~~~

The DAP uses an `SDSS parameter file
<https://www.sdss.org/dr15/software/par/>`_ to define one or more
methods to use when analyzing any given MaNGA datacube. Each method,
or "analysis plan", is defined by a set of six keywords that identify
the method to use for each of the DAP's six main
:ref:`workflow-analysis-modules`. The 
``AnalysisPlan`` parameter file
allows you to analyze the same datacube multiple ways in a single
execution of the DAP; however, note that this is no different than
executing the DAP once per analysis method. The ``AnalysisPlan`` is
required to execute the DAP; however, the DAP provides a
`execution-analysis-plan-default`_ that will be used if no file is
provided.

An example ``AnalysisPlan`` parameter file looks like this (this is
exactly the file used for MPL-9):

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
    DAPPLAN    SNRG  0     SPX  0   MILESHCMPL9  0  EMOMMPL9  0    EFITMPL9  0   INDXEN  0
    DAPPLAN    SNRG  0   VOR10  0   MILESHCMPL9  0  EMOMMPL9  0    EFITMPL9  0   INDXEN  0
    DAPPLAN    SNRG  0   HYB10  0   MILESHCMPL9  0  EMOMMPL9  0  EFITMPL9DB  0   INDXEN  0

The configuration of each module is set by two values: a configuration
key and a flag indicating if any existing results should be overwritten
(0 = False, 1 = True).  Each keyword points to a unique configuration
file, and each of these keywords produce a unique instance of a
parameter set that drives the relevant analysis module.  The link
between the relevant structure element root (e.g., ``drpqa``), the
location of its configuration files, and the object used to read the
configurations are as follows:

+-----------+-------------------------------------+-------------------------------------------------------------------------+
|    struct |  ``$MANGADAP_DIR/mangadap/config/`` |                                                              DAP object |
+===========+=====================================+=========================================================================+
|     drpqa |           ``reduction_assessments`` |      :class:`mangadap.proc.reductionassessments.ReductionAssessmentDef` |
+-----------+-------------------------------------+-------------------------------------------------------------------------+
|       bin |                 ``spatial_binning`` | :class:`mangadap.proc.spatiallybinnedspectra.SpatiallyBinnedSpectraDef` |
+-----------+-------------------------------------+-------------------------------------------------------------------------+
| continuum |      ``stellar_continuum_modeling`` |   :class:`mangadap.proc.stellarcontinuummodel.StellarContinuumModelDef` |
+-----------+-------------------------------------+-------------------------------------------------------------------------+
|     elmom |           ``emission_line_moments`` |       :class:`mangadap.proc.emissionlinemoments.EmissionLineMomentsDef` |
+-----------+-------------------------------------+-------------------------------------------------------------------------+
|     elfit |          ``emission_line_modeling`` |           :class:`mangadap.proc.emissionlinemodel.EmissionLineModelDef` |
+-----------+-------------------------------------+-------------------------------------------------------------------------+
|   spindex |                ``spectral_indices`` |               :class:`mangadap.proc.spectralindices.SpectralIndicesDef` |
+-----------+-------------------------------------+-------------------------------------------------------------------------+

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

.. somehow generate this automatically

If executed without an ``AnalysisPlan`` parameter file, the command-line
execution of the DAP will use a default plan; see
:func:`mangadap.par.analysisplan.AnalysisPlanSet.default`.

The current default plan uses the following keys:

+---------------+----------------+
|     drpqa_key |           SNRG |
+---------------+----------------+
|       bin_key |          HYB10 |
+---------------+----------------+
| continuum_key |    MILESHCMPL9 |
+---------------+----------------+
|     elmom_key |       EMOMMPL9 |
+---------------+----------------+
|     elfit_key |     EFITMPL9DB |
+---------------+----------------+
|   spindex_key |         INDXEN |
+---------------+----------------+

.. _execution-config:

The DAP Datacube Configuration File
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The DAP uses a configuration (``ini``) file to set the datacube to be
analyzed and provide some relevant metadata. These configuration
files are generated at the survey-level by a combination of
:func:`mangadap.survey.drpcomplete.write_config` and
:func:`mangadap.datacube.MaNGADataCube.write_config`. However, we
also provide the ``$MANGADAP_DIR/bin/write_dap_config`` script that
will generate the relevant configuration files if you have the DRPall
or DRPComplete file. As with all the DAP scripts, you can use the
``-h`` command-line option to get the usage:

.. code-block::

    $ write_dap_config -h
    usage: write_dap_config [-h] (-c DRPCOMPLETE | -a DRPALL)
                            [--sres_ext SRES_EXT] [--sres_fill SRES_FILL]
                            [--covar_ext COVAR_EXT] [--drpver DRPVER]
                            [--redux_path REDUX_PATH]
                            [--directory_path DIRECTORY_PATH] [-o]
                            plate ifudesign ofile

    positional arguments:
      plate                 Plate number
      ifudesign             IFU design number
      ofile                 Output file name

    optional arguments:
      -h, --help            show this help message and exit
      -c DRPCOMPLETE, --drpcomplete DRPCOMPLETE
                            DRP complete fits file
      -a DRPALL, --drpall DRPALL
                            DRPall fits file
      --sres_ext SRES_EXT   Spectral resolution extension to use. Default set by
                            MaNGADataCube class.
      --sres_fill SRES_FILL
                            If present, use interpolation to fill any masked
                            pixels in the spectral resolution vectors. Default set
                            by MaNGADataCube class.
      --covar_ext COVAR_EXT
                            Use this extension to define the spatial correlation
                            matrix. Default set by MaNGADataCube class.
      --drpver DRPVER       DRP version. Default set by MaNGADataCube class.
      --redux_path REDUX_PATH
                            Path to the top-level DRP reduction directory. Default
                            set by MaNGADataCube class.
      --directory_path DIRECTORY_PATH
                            Exact path to the directory with the MaNGA DRP
                            datacube. The name of the file itself must match the
                            nominal MaNGA DRP naming convention. Default set by
                            MaNGADataCube class.
      -o, --overwrite       Overwrite any existing files.

To construct the configuration file for datacube 7815-3702, executing:

.. code-block::

    write_dap_config 7815 3702 mangadap-7815-3702.cfg -a drpall-v2_7_1.fits

produces the following file:

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
    angle, ``pa``, but imposes periodic limits (i.e., 380 deg is
    converted set to 20 deg). If the effective radius, ``reff``, is
    :math:`R_{\rm eff} < 0` or undefined, the DAP uses :math:`R_{\rm
    eff} = 1`.

.. _execution-mangadap:

DAP command-line script
-----------------------

The main DAP script is ``$MANGADAP_DIR/bin/manga_dap``, which is a
simple wrapper of :func:`mangadap.survey.manga_dap.manga_dap`.  With the
DAP installed, you can call the script directly from the command line:

.. code-block::

    $ manga_dap -h
    usage: manga_dap [-h] (-c CONFIG | -f CUBEFILE) [-p PLAN] [-m CUBE_MODULE]
                     [-o CUBE_OBJECT] [--dbg] [--log LOG] [-v] [--drpver DRPVER]
                     [-r REDUX_PATH] [-d DIRECTORY_PATH] [--dapver DAPVER]
                     [-a ANALYSIS_PATH]

    optional arguments:
      -h, --help            show this help message and exit
      -c CONFIG, --config CONFIG
                            Configuration file used to instantiate the relevant
                            DataCube derived class. (default: None)
      -f CUBEFILE, --cubefile CUBEFILE
                            Name of the file with the datacube data. Must be
                            possible to instantiate the relevant DataCube derived
                            class directly from the file only. (default: None)
      -p PLAN, --plan PLAN  SDSS parameter file with analysis plan. If not
                            provided, a default plan is used. (default: None)
      -m CUBE_MODULE, --cube_module CUBE_MODULE
                            The name of the module that contains the DataCube
                            derived class. (default: mangadap.datacube)
      -o CUBE_OBJECT, --cube_object CUBE_OBJECT
                            The name of the DataCube derived class object.
                            (default: MaNGADataCube)
      --dbg                 Run manga_dap in debug mode (default: False)
      --log LOG             File name for runtime log (default: None)
      -v, --verbose         Set verbosity level; can be omitted and set up to -vv
                            (default: 0)
      --drpver DRPVER       DRP version (default: None)
      -r REDUX_PATH, --redux_path REDUX_PATH
                            Top-level directory with the DRP products; defaults to
                            $MANGA_SPECTRO_REDUX/$MANGADRP_VER (default: None)
      -d DIRECTORY_PATH, --directory_path DIRECTORY_PATH
                            Path directly to directory with DRP file to analyze
                            (default: None)
      --dapver DAPVER       DAP version (default: None)
      -a ANALYSIS_PATH, --analysis_path ANALYSIS_PATH
                            Top-level output directory for the DAP results;
                            defaults to
                            $MANGA_SPECTRO_ANALYSIS/$MANGADRP_VER/$MANGADAP_VER
                            (default: None)

The DAP allows you to define your own datacube class, as long as it
is derived from :class:`mangadap.datacube.datacube.DataCube`. You can
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

.. code-block::

    write_dap_config 7815 3702 mangadap-7815-3702.cfg -a drpall-v2_7_1.fits
    manga_dap -c mangadap-7815-3702.cfg -vv -log mangadap-7815-3702.log -d . -a dap_output

This will analyze the datacube for observation 7815-3702 using the
default analysis plan, with verbose output and a log written to
``mangadap-7815-3702.log`` and with the root directory for all the
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
:class:`mangadap.scripts.rundap.rundap`.  This script

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

.. code-block:: bash

    $ rundap -h
    usage: rundap [-h] [--clobber] [-v] [--quiet] [--print_version] [--loose]
                  [--mplver MPLVER] [--redux_path REDUX_PATH] [--dapver DAPVER]
                  [--analysis_path ANALYSIS_PATH] [--plan_file PLAN_FILE]
                  [--platelist PLATELIST] [--ifudesignlist IFUDESIGNLIST]
                  [--list_file LIST_FILE] [--combinatorics] [--sres_ext SRES_EXT]
                  [--sres_fill SRES_FILL] [--covar_ext COVAR_EXT]
                  [--use_plttargets] [--plttargets PLTTARGETS] [--on_disk] [--log]
                  [--no_proc] [--no_plots] [--post] [--post_plots] [--dapall]
                  [--label LABEL] [--nodes NODES] [--cpus CPUS] [--fast QOS]
                  [--umask UMASK] [--walltime WALLTIME] [--toughness] [--create]
                  [--submit] [--progress] [--queue QUEUE]

    optional arguments:
      -h, --help            show this help message and exit
      --clobber             if all selected, will run dap for all
                            plates/ifudesigns/modes regardless of state (default:
                            False)
      -v, --verbose         Set verbosity level for manga_dap; can be omitted and
                            set up to -vv (default: 0)
      --quiet               suppress screen output (default: False)
      --print_version       print DAP version and stop (default: False)
      --loose               Only throw warnings if the versioning is not
                            identically as it should be for the designated MPL
                            (default: False)
      --mplver MPLVER       select MPL version to analyze (default: None)
      --redux_path REDUX_PATH
                            main DRP output path (default: None)
      --dapver DAPVER       optional output version, different from product
                            version (default: None)
      --analysis_path ANALYSIS_PATH
                            main DAP output path (default: None)
      --plan_file PLAN_FILE
                            parameter file with the MaNGA DAP execution plan to
                            use instead of the default (default: None)
      --platelist PLATELIST
                            set list of plates to reduce (default: None)
      --ifudesignlist IFUDESIGNLIST
                            set list of ifus to reduce (default: None)
      --list_file LIST_FILE
                            A file with the list of plates and ifudesigns to
                            analyze (default: None)
      --combinatorics       force execution of all permutations of the provided
                            lists (default: False)
      --sres_ext SRES_EXT   Spectral resolution extension to use. Default set by
                            MaNGADataCube class. (default: None)
      --sres_fill SRES_FILL
                            If present, use interpolation to fill any masked
                            pixels in the spectral resolution vectors. Default set
                            by MaNGADataCube class. (default: None)
      --covar_ext COVAR_EXT
                            Use this extension to define the spatial correlation
                            matrix. Default set by MaNGADataCube class. (default:
                            None)
      --use_plttargets      Use platetargets files instead of the DRPall file to
                            generate the DRP complete database (default: False)
      --plttargets PLTTARGETS
                            path to plateTargets file(s); if provided will force
                            update to drpcomplete fits file (default: None)
      --on_disk             When using the DRPall file to collate the data for
                            input to the DAP, search for available DRP files on
                            disk instead of using the DRPall file content.
                            (default: False)
      --log                 Have the main DAP executable produce a log file
                            (default: False)
      --no_proc             Do NOT perform the main DAP processing steps (default:
                            False)
      --no_plots            Do NOT create QA plots (default: False)
      --post                Create/Submit the post-processing scripts (default:
                            False)
      --post_plots          Create/Submit the post-processing plotting scripts
                            (default: False)
      --dapall              Wait for any individual plate-ifu processes to finish
                            and thenupdate the DAPall file (default: False)
      --label LABEL         label for cluster job (default: None)
      --nodes NODES         number of nodes to use in cluster (default: 1)
      --cpus CPUS           number of cpus to use per node. Default is to use all
                            available; otherwise, set to minimum of provided
                            number and number of processors per node (default:
                            None)
      --fast QOS            qos state (default: None)
      --umask UMASK         umask bit for cluster job (default: 0027)
      --walltime WALLTIME   walltime for cluster job (default: 240:00:00)
      --toughness           turn off hard keyword for cluster submission (default:
                            True)
      --create              use the pbs package to create the cluster scripts
                            (default: False)
      --submit              submit the scripts to the cluster (default: False)
      --progress            instead of closing the script, report the progress of
                            the analysis on the cluster; this is required if you
                            want to submit the DAPall script immediately after
                            completing the individual cube analysis (default:
                            False)
      --queue QUEUE         set the destination queue (default: None)

If :ref:`execution-analysis-plan` is not provided, the scripts will
use the :ref:`execution-analysis-plan-default`.

An example call of this script that will only construct scripts for
the analysis of observation 7443-12701 using the default AnalysisPlan
is:

.. code-block:: bash
    
    rundap --platelist 7443 --ifudesignlist 12701 --redux_path /path/with/drp/output/ --analysis_path /path/for/dap/output/ -vv --log

In this call, I've specified that the DRP data is in
``/path/with/drp/output/`` and that the DAP output should be placed
in ``/path/for/dap/output/`` instead of using the default
:ref:`datamodel-directory-structure`.

The **script file** this call produces is written to
``/path/for/dap/output/log/[time]/7495/12704/mangadap-7495-12704``,
where ``[time]`` is a time stamp of when ``rundap`` was executed.  (If
you execute ``rundap`` multiple times, it will create new directories
using new time stamps each time.)  The lines of the script file for each
plate-ifu:

 - touches the ``*.started`` file
 - executes ``manga_dap``
 - executes a series of QA plotting scripts
 - touches the ``*.done`` file 

.. todo::

    Update the scripts produced by ``rundap``.

The example script generated by the above command would look something like this:

.. code-block:: bash

    # Auto-generated batch file
    # Fri 01 Nov 2019 10:58:52

    touch /path/for/dap/output/v2_7_1/2.4.1/log/01Nov2019T16.58.40UTC/7443/12701/mangadap-7443-12701.started

    manga_dap /path/for/dap/output/v2_7_1/2.4.1/common/7443/12701/mangadap-7443-12701-LOGCUBE-input.par -r /path/with/drp/output/v2_7_1 -a /path/for/dap/output/v2_7_1/2.4.1 --log /path/for/dap/output/v2_7_1/2.4.1/log/01Nov2019T16.58.40UTC/7443/12701/mangadap-7443-12701.log -vv

    dap_ppxffit_qa 7443 12701 --analysis_path /path/for/dap/output/v2_7_1/2.4.1

    spotcheck_dap_maps 7443 12701 --analysis_path /path/for/dap/output/v2_7_1/2.4.1

    dap_fit_residuals 7443 12701 --analysis_path /path/for/dap/output/v2_7_1/2.4.1

    touch /path/for/dap/output/v2_7_1/2.4.1/log/01Nov2019T16.58.40UTC/7443/12701/mangadap-7443-12701.done

To execute the script, you would then run:

.. code-block:: bash

    source /path/for/dap/output/v2_7_1/2.4.1/log/01Nov2019T16.58.40UTC/7443/12701/mangadap-7443-12701

The ``rundap`` script allows you to construct scripts for all datacubes
it can find on disk, all IFUs on a given plate, all combinations of a
set of plate and IFU numbers, or for a specified list of ``plateifu``
IDs.

.. note::

    The ``rundap`` script constructs the
    :class:`mangadap.survey.drpcomplete.DRPComplete` object and
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


In the automated run of the DAP, any entry in the
:ref:`metadatamodel-drpcomplete` file that meets the criteria set by
:func:`mangadap.survey.drpcomplete.DRPComplete.can_analyze` will be
analyzed. Currently the relevant criteria are (check against the
documentation of the relevant function):
 
 - ``MANGAID != NULL``
 - ``MANGA_TARGET1 > 0 | MANGA_TARGET3 > 0``
 - ``VEL > -500``

An important consequence of the selection above is that *any
ancillary targets without a provided redshift will not be analyzed by
the DAP*, unless it has replacement redshift in the
:ref:`metadatamodel-redshift-fix`. Ancillary targets not analyzed by
the DAP are likely because a redshift was not available.

