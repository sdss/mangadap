Execution
=========

At the moment, the DAP execution is tailored toward its execution at the
survey level.  On a best effort basis, we are trying to simplify the
inputs to the main executable, as well as abstracting the DAP front-end
to more easily work with non-MaNGA data.

Input files
-----------

1. The DAP AnalysisPlan
~~~~~~~~~~~~~~~~~~~~~~~

The DAP uses an `SDSS parameter file
<https://www.sdss.org/dr15/software/par/>`_ to define one or more
methods to use when analyzing any given MaNGA datacube.  Each method, or
"analysis plan", is defined by a set of six keywords that identify the
configuration to use for each of the DAP's six main analysis modules.
To execute the DAP, you must construct this parameter file with at least
one analysis plan.  For example, the MPL-9 analysis-plan parameter file
is:

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

+-----------+--------------------------------------------+-------------------------------------------------------------------------+
|    struct |  ``$MANGADAP_DIR/python/mangadap/config/`` |                                                              DAP object |
+===========+============================================+=========================================================================+
|     drpqa |                  ``reduction_assessments`` |      :class:`mangadap.proc.reductionassessments.ReductionAssessmentDef` |
+-----------+--------------------------------------------+-------------------------------------------------------------------------+
|       bin |                        ``spatial_binning`` | :class:`mangadap.proc.spatiallybinnedspectra.SpatiallyBinnedSpectraDef` |
+-----------+--------------------------------------------+-------------------------------------------------------------------------+
| continuum |             ``stellar_continuum_modeling`` |   :class:`mangadap.proc.stellarcontinuummodel.StellarContinuumModelDef` |
+-----------+--------------------------------------------+-------------------------------------------------------------------------+
|     elmom |                  ``emission_line_moments`` |       :class:`mangadap.proc.emissionlinemoments.EmissionLineMomentsDef` |
+-----------+--------------------------------------------+-------------------------------------------------------------------------+
|     elfit |                 ``emission_line_modeling`` |           :class:`mangadap.proc.emissionlinemodel.EmissionLineModelDef` |
+-----------+--------------------------------------------+-------------------------------------------------------------------------+
|   spindex |                       ``spectral_indices`` |               :class:`mangadap.proc.spectralindices.SpectralIndicesDef` |
+-----------+--------------------------------------------+-------------------------------------------------------------------------+

When setting the keyword values in the plan file, they must be
recognized as a method defined in the relevant configuration directory.
The DAP will only execute properly if at least the first three steps
have valid keywords.  The remaining three modules can be skipped (i.e.,
the emission-line moments, emission-line-model parameters, and spectral
indices are not measured) by setting their keyword to ``None``; the
primary DAP output will still be produced but with empty arrays for the
skipped analysis steps.

2. The DAP ObsInputPar
~~~~~~~~~~~~~~~~~~~~~~

The DAP uses another `SDSS parameter file
<https://www.sdss.org/dr15/software/par/>`_ to set the observation it is
to analyze, as well as some input parameters that it uses during its
execution.  Specifically, it requires:

+---------------+----------------------------------------------------------------+
+     Parameter | Description                                                    |
+===============+================================================================+
|     ``plate`` | The MaNGA plate number                                         |
+---------------+----------------------------------------------------------------+
| ``ifudesign`` | The MaNGA IFU number                                           |
+---------------+----------------------------------------------------------------+
|      ``mode`` | The MaNGA DRP output mode (should always be ``CUBE``)          |
+---------------+----------------------------------------------------------------+
|       ``vel`` | Guess systemic velocity, :math:`cz`, in km/s                   |
+---------------+----------------------------------------------------------------+
|     ``vdisp`` | Guess stellar velocity dispersion; if less than 0, the DAP     |
|               | assumes 100 km/s.                                              |
+---------------+----------------------------------------------------------------+
|       ``ell`` | Isophotal ellipticity (:math:`\epsilon = 1-b/a`).              |
+---------------+----------------------------------------------------------------+
|        ``pa`` | Position angle, :math:`0 \leq \phi_0 < 360`, of the            |
|               | isophotal ellipse defined as the angle (deg) from N through E. |
+---------------+----------------------------------------------------------------+
|      ``reff`` | Effective (half-light) radius, :math:`R_{\rm eff}`.            |
+---------------+----------------------------------------------------------------+

**Notes:** If :math:`\epsilon < 0` or :math:`\epsilon > 1`, the DAP will
adopt a default value of 0.  The DAP accepts any value for
:math:`\phi_0` but imposes the periodic limits (i.e., 380 deg is
converted set to 20 deg).  If :math:`R_{\rm eff} < 0`, the DAP uses
:math:`R_{\rm eff} = 1`. 

An example file looks like this:

.. code-block:: c

    typedef struct {
        long plate;
        long ifudesign;
        char mode[4];
        double vel;
        double vdisp;
        double ell;
        double pa;
        double reff;
    } DAPPAR;

    DAPPAR 7443 12701 CUBE  6.1391504e+03  1.0000000e+02  3.4162802e-01  1.5024400e+02  5.7011299e+00

DAP command-line script
-----------------------

The main DAP script is ``$MANGADAP_DIR/bin/manga_dap``, which is a
simple wrapper of :func:`mangadap.survey.manga_dap.manga_dap`.  With the
DAP installed, you can call the script directly from the command line:

.. code-block:: bash

    $ manga_dap -h
    usage: manga_dap [-h] [--dbg] [--log LOG] [-v] [--drpver DRPVER]
                     [-r REDUX_PATH] [-d DIRECTORY_PATH] [--dapver DAPVER]
                     [-s DAP_SRC] [-a ANALYSIS_PATH]
                     obs plan

    positional arguments:
      obs                   SDSS parameter file with observational input
      plan                  SDSS parameter file with analysis plan

    optional arguments:
      -h, --help            show this help message and exit
      --dbg                 Run manga_dap in debug mode
      --log LOG             File name for runtime log
      -v, --verbose         Set verbosity level; can be omitted and set up to -vv
      --drpver DRPVER       DRP version
      -r REDUX_PATH, --redux_path REDUX_PATH
                            Top-level directory with the DRP products; defaults to
                            $MANGA_SPECTRO_REDUX/$MANGADRP_VER
      -d DIRECTORY_PATH, --directory_path DIRECTORY_PATH
                            Path directly to directory with DRP file to analyze
      --dapver DAPVER       DAP version
      -s DAP_SRC, --dap_src DAP_SRC
                            Top-level directory with the DAP source code; defaults
                            to $MANGADAP_DIR
      -a ANALYSIS_PATH, --analysis_path ANALYSIS_PATH
                            Top-level output directory for the DAP results;
                            defaults to
                            $MANGA_SPECTRO_ANALYSIS/$MANGADRP_VER/$MANGADAP_VER

.. warning::

    When running the DAP, you should have both the DRP ``LOGRSS`` and
    ``LOGCUBE`` files if you want to account for covariance!  If the
    ``LOGRSS`` files are not present, the DAP will throw a warning and
    continue, and the warnings can get buried in all the other messages.

An example execution of the DAP might look like this:

.. code-block:: bash

    manga_dap mangadap-7495-12704-LOGCUBE-input.par plan.par --log mangadap-7495-12704.log -vv

where ``mangadap-7495-12704-LOGCUBE-input.par`` and ``plan.par`` are,
respectively, `The DAP ObsInputPar`_ and `The DAP AnalysisPlan`_ files.

Programmatic excution
---------------------

Alternatively, ``$MANGADAP_DIR/examples/fit_one_cube.py`` provides a
programmatic approach to running the exact same script that is executed
by the ``manga_dap`` command-line script.  The code provides a way to
generate the `The DAP ObsInputPar`_ object directly from the DRPall
file, instead of from a file, and it directly defines the
``AnalysisPlan`` object with a hard-coded set of keywords.  Using this
script as an example, one could construct a script that programmatically
analyzes a large set of datacubes.

Batch execution using automatically generated scripts
-----------------------------------------------------

The survey-level execution of the DAP uses the
``$MANGADAP_DIR/bin/rundap`` script, which is a simple wrapper of
:class:`mangadap.survey.rundap.rundap`.  This script

 * sets up the DAP output directory structure
 * either confirms that a provided list of datacubes to analyze exist on
   disk or trolls the DRP directory structure to find all or some subset
   of available datacubes to analyze
 * creates `The DAP ObsInputPar`_ parameter file for each ``plateifu``
   to be analyzed,
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
                  [--list_file LIST_FILE] [--combinatorics] [--use_plttargets]
                  [--plttargets PLTTARGETS] [--on_disk] [--log] [--no_proc]
                  [--no_plots] [--post] [--post_plots] [--dapall] [--label LABEL]
                  [--nodes NODES] [--cpus CPUS] [--fast QOS] [--umask UMASK]
                  [--walltime WALLTIME] [--toughness] [--create] [--submit]
                  [--progress] [--queue QUEUE]
    
    optional arguments:
      -h, --help            show this help message and exit
      --clobber             if all selected, will run dap for all
                            plates/ifudesigns/modes regardless of state
      -v, --verbose         Set verbosity level for manga_dap; can be omitted and
                            set up to -vv
      --quiet               suppress screen output
      --print_version       print DAP version and stop
      --loose               Only throw warnings if the versioning is not
                            identically as it should be for the designated MPL
      --mplver MPLVER       select MPL version to analyze
      --redux_path REDUX_PATH
                            main DRP output path
      --dapver DAPVER       optional output version, different from product
                            version
      --analysis_path ANALYSIS_PATH
                            main DAP output path
      --plan_file PLAN_FILE
                            parameter file with the MaNGA DAP execution plan to
                            use instead of the default
      --platelist PLATELIST
                            set list of plates to reduce
      --ifudesignlist IFUDESIGNLIST
                            set list of ifus to reduce
      --list_file LIST_FILE
                            a file with the list of plates, ifudesigns, and modes
                            to analyze
      --combinatorics       force execution of all permutations of the provided
                            lists
      --use_plttargets      Use platetargets files instead of the DRPall file to
                            generate the DRP complete database
      --plttargets PLTTARGETS
                            path to plateTargets file(s); if provided will force
                            update to drpcomplete fits file
      --on_disk             When using the DRPall file to collate the data for
                            input to the DAP, search for available DRP files on
                            disk instead of using the DRPall file content.
      --log                 Have the main DAP executable produce a log file
      --no_proc             Do NOT perform the main DAP processing steps
      --no_plots            Do NOT create QA plots
      --post                Create/Submit the post-processing scripts
      --post_plots          Create/Submit the post-processing plotting scripts
      --dapall              Wait for any individual plate-ifu processes to finish
                            and then update the DAPall file
      --label LABEL         label for cluster job
      --nodes NODES         number of nodes to use in cluster
      --cpus CPUS           number of cpus to use per node. Default is to use all
                            available; otherwise, set to minimum of provided
                            number and number of processors per node
      --fast QOS            qos state
      --umask UMASK         umask bit for cluster job
      --walltime WALLTIME   walltime for cluster job
      --toughness           turn off hard keyword for cluster submission
      --create              use the pbs package to create the cluster scripts
      --submit              submit the scripts to the cluster
      --progress            instead of closing the script, report the progress of
                            the analysis on the cluster; this is required if you
                            want to submit the DAPall script immediately after
                            completing the individual cube analysis
      --queue QUEUE         set the destination queue
    

Note that you still need a parameter file with `The DAP AnalysisPlan`_
details; all datacubes selected by the ``rundap`` execution will be
analyzed with the same ``AnalysisPlan``.  An example call of this script that will only construct scripts for the analysis of observation 7443-12701 is:

.. code-block:: bash
    
    rundap --platelist 7443 --ifudesignlist 12701 --redux_path /path/with/drp/output/ --analysis_path /path/for/dap/output/ --plan_file /path/to/plan/file/plan.par -vv --log

In this call, I've specified that the DRP data is in
``/path/with/drp/output/`` and that the DAP output should be placed in
``/path/for/dap/output/`` instead of using the default directory
structure based on the `Local Environment Setup`_.  The script file this
call produces is written to
``/path/for/dap/output/log/[time]/7495/12704/mangadap-7495-12704``,
where ``[time]`` is a time stamp of when ``rundap`` was executed.  (If
you execute ``rundap`` multiple times, it will create new directories
using new time stamps.)  The lines of the script file for each plate-ifu:

 - touches the ``*.started`` file
 - executes manga_dap
 - executes a series of QA plotting scripts
 - touches the ``*.done`` file 

The example script generated by the above command would look something like this:

.. code-block:: bash

    # Auto-generated batch file
    # Fri 01 Nov 2019 10:58:52

    touch /path/for/dap/output/v2_7_1/2.4.1/log/01Nov2019T16.58.40UTC/7443/12701/mangadap-7443-12701.started

    manga_dap /path/for/dap/output/v2_7_1/2.4.1/common/7443/12701/mangadap-7443-12701-LOGCUBE-input.par /path/for/dap/output/v2_7_1/2.4.1/log/01Nov2019T16.58.40UTC/plan.par -r /path/with/drp/output/v2_7_1 -a /path/for/dap/output/v2_7_1/2.4.1 --log /path/for/dap/output/v2_7_1/2.4.1/log/01Nov2019T16.58.40UTC/7443/12701/mangadap-7443-12701.log -vv

    ppxffit_qa 7443 12701 --analysis_path /path/for/dap/output/v2_7_1/2.4.1 --plan_file /path/for/dap/output/v2_7_1/2.4.1/log/01Nov2019T16.58.40UTC/plan.par

    spotcheck_dap_maps 7443 12701 --analysis_path /path/for/dap/output/v2_7_1/2.4.1 --plan_file /path/for/dap/output/v2_7_1/2.4.1/log/01Nov2019T16.58.40UTC/plan.par

    dap_fit_residuals 7443 12701 --analysis_path /path/for/dap/output/v2_7_1/2.4.1 --plan_file /path/for/dap/output/v2_7_1/2.4.1/log/01Nov2019T16.58.40UTC/plan.par

    touch /path/for/dap/output/v2_7_1/2.4.1/log/01Nov2019T16.58.40UTC/7443/12701/mangadap-7443-12701.done

To execute the script, you would then run:

.. code-block:: bash

    source /path/for/dap/output/log/01Nov2019T16.58.40UTC/7443/12701/mangadap-7443-12701

The ``rundap`` script allows you to construct scripts for all datacubes
it can find on disk, all IFUs on a given plate, all combinations of a
set of plate and IFU numbers, or for a specified list of ``plateifu``
IDs.

.. note:: 
    
    The ``rundap`` script constructs the
    :class:`mangadap.survey.drpcomplete.DRPComplete` object and writes
    its associated fits file.  The data compiled into this database can
    either be drawn from the DRPall file or from the plateTargets data
    in ``mangacore``; the latter is the only reason the DAP has
    ``mangacore`` as a dependency.  For general use, you should have
    ``rundap`` use the DRPall file.  The use of the plateTargets data is
    only necessary in the rare case when the DAP is executed before the
    relevant DRPall file has been constructed.

To write the post-processing scripts, execute ``rundap`` with the
``--post`` and ``--post_plots`` options.  This produces two additional
types of scripts:

    1. Scripts to produce QA plots for all IFUs on a given plate.  This
    file is written to, e.g.,
    ``/path/for/dap/output/log/01Nov2019T16.58.40UTC/7443/7443_fitqa``
    and looks like this:

    .. code-block:: bash

        # Auto-generated batch file
        # Fri 01 Nov 2019 10:58:52

        touch /path/for/dap/output/v2_7_1/2.4.1/log/01Nov2019T16.58.40UTC/7443/7443_fitqa.started

        dap_plate_fit_qa 7443 --analysis_path /path/for/dap/output/v2_7_1/2.4.1 --plan_file /path/for/dap/output/v2_7_1/2.4.1/log/01Nov2019T16.58.40UTC/plan.par

        touch /path/for/dap/output/v2_7_1/2.4.1/log/01Nov2019T16.58.40UTC/7443/7443_fitqa.done

    2. A script that builds the DAPall file and writes its QA plots.
    This file is written to, e.g.,
    ``/path/for/dap/output/v2_7_1/2.4.1/log/01Nov2019T16.58.40UTC/build_dapall``
    and looks like this:

    .. code-block:: bash

        # Auto-generated batch file
        # Fri 01 Nov 2019 10:58:52

        touch /path/for/dap/output/v2_7_1/2.4.1/log/01Nov2019T16.58.40UTC/build_dapall.started

        construct_dapall /path/for/dap/output/v2_7_1/2.4.1/log/01Nov2019T16.58.40UTC/plan.par --drpver v2_7_1 -r /path/with/drp/output/v2_7_1 --dapver 2.4.1 -a /path/for/dap/output/v2_7_1/2.4.1 -vv

        dap_dapall_qa --drpver v2_7_1 --redux_path /path/with/drp/output/v2_7_1 --dapver 2.4.1 --analysis_path /path/for/dap/output/v2_7_1/2.4.1 --plan_file /path/for/dap/output/v2_7_1/2.4.1/log/01Nov2019T16.58.40UTC/plan.par

        touch /path/for/dap/output/v2_7_1/2.4.1/log/01Nov2019T16.58.40UTC/build_dapall.done


