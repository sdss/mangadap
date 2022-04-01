
.. include:: include/links.rst

.. _execution:

Execution
=========

The MaNGA DAP was originally developed to analyze data from the SDSS-IV/MaNGA
Survey, provided reduction products from the MaNGA data-reduction pipeline
(DRP).  *However, it is now easy to apply the MaNGA DAP to any IFU datacube.*

Below, we describe how to use the MaNGA DAP to analyze MaNGA
datacubes, as well as how to execute the survey-level batch mode that analyzes
all the MaNGA datacubes within a given directory structure.

.. warning::

    Much of the survey-level execution of the code needs to be updated for the
    front end changes in version>4.x.  If not removed outright, out of date
    documentation is flagged in what follows.

To analyze non-MaNGA datacubes, see both :ref:`fitdatacube` and the description
of the MaNGA DAP :ref:`plan`.  Also, for an example of how to fit single spectra
using MaNGA DAP core algorithms, see :ref:`fitonespec`.

Input files
-----------

.. _execution-analysis-plan:

DAP AnalysisPlan
~~~~~~~~~~~~~~~~

The DAP uses a `toml`_ file to define parameters used by its modules to analyze
the provided datacube.  This is described in detail by :ref:`plan`.

.. _execution-analysis-plan-default:

Default AnalysisPlan
++++++++++++++++++++

If executed without an ``AnalysisPlan`` parameter file, the command-line
execution of the DAP will use a default plan; see
:func:`mangadap.config.analysisplan.AnalysisPlan.default` and `default_plan.toml
<https://github.com/sdss/mangadap/blob/plan/mangadap/config/default_plan.toml>`__.

.. _execution-config:

The DAP Datacube Configuration File
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The DAP uses a configuration (``ini``) file to set the datacube to be analyzed
and provides some relevant metadata. These configuration files are generated at
the survey-level by :func:`~mangadap.survey.drpcomplete.write_config`. However,
we also provide the ``write_dap_config`` script that will generate the relevant
configuration files if you have the DRPall or DRPComplete file. As with all the
DAP scripts, you can use the ``-h`` command-line option to get the usage:

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

The DAP allows you to define your own datacube class, as long as it is derived
from :class:`~mangadap.datacube.datacube.DataCube`. You can then specify that
your data should be instantiated with that derived class using the
``cube_module`` argument; this defaults to
``mangadap.datacube.datacube.MaNGADataCube``, expected a MaNGA datacube.

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

    write_dap_config 7815 3702 mangadap-7815-3702.ini -a drpall-v3_1_1.fits --directory_path .
    manga_dap -c mangadap-7815-3702.ini -vv --log mangadap-7815-3702.log -d . -o dap_output

This will analyze the datacube for observation 7815-3702 using the default
analysis plan, with verbose output and a log written to
``mangadap-7815-3702.log``, and with the root directory for all the DAP output
(except for the log) set to ``dap_output``.  These commands can be successfully
run in the ``$MANGADAP_DIR/mangadap/data/remote`` directory if you've run the
``download_test_data.py`` script.

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

The ``manga_dap`` executable script is simple, in that it only (1) reads the
datacube, (2) sets the analysis plan, and (3) executes the main DAP analysis
wrapper function, :func:`mangadap.survey.manga_dap.manga_dap`.  It is
straight-forward then to construct a script that executes the DAP
programmatically for many datacubes.

.. _execution-rundap:

Batch execution using automatically generated scripts
-----------------------------------------------------

.. warning::

    **The batch execution of the DAP (even for MaNGA) is currently out-of-date.
    These will be fixed soon.**

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

Execution Recovery
------------------

If the code faults, all of the modules completed up until the fault occurred
should have created a "reference" file that will effectively allow the code to
pick up where it left off.  The mechanism used to determine if it can do this
relies simply on the existence of the expected output file.  This means that if
you change on of the parameters in your input :ref:`plan` file *without*
changing the keyword identifier for that module parameter set, the code may read
in the existing file and keep going without incorporating your parameter change.

If you're testing the performance with different parameter values, either
perform the testing with different output directories, always remember to change
the keyword for the relevant module parameter set, set the overwrite keyword for
the relevant module (and each subsequent module!) to True, or just nuke the
directory and start again.

Quality Assessment Plots
------------------------

The survey-level execution of the MaNGA DAP constructs automatically generated
plots that provide spotchecks of the performance of the DAP; see
:ref:`qualityassessment`.  The main QA plots for the analysis of a single
datacube have a similar calling sequence as the main :ref:`execution-mangadap`.
The three main scripts are ``spotcheck_dap_maps``, ``dap_ppxffit_qa``, and
``dap_fit_residuals``.  Following the example above, these can be created after
analysing the test datacube as follows:

.. code-block:: console

    spotcheck_dap_maps -c mangadap-7815-3702.ini -vv --log mangadap-7815-3702.log -d . -o dap_output
    dap_ppxffit_qa -c mangadap-7815-3702.ini -vv --log mangadap-7815-3702.log -d . -o dap_output
    dap_fit_residuals -c mangadap-7815-3702.ini -vv --log mangadap-7815-3702.log -d . -o dap_output

Local Environment Setup for Survey-Level MaNGA Analysis
-------------------------------------------------------

The DAP uses environmental variables to define the paths to specific data and
other repositories, when executed for MaNGA data. If these are not defined,
default values will be used; see `the initialization of the mangadap.config
module
<https://github.com/sdss/mangadap/blob/plan/mangadap/config/__init__.py>`__.
The relevant environmental variables, their default, and their usage are
provided below.

+----------------------------+-------------------------------------+------------------------------------------------+
|                   Variable |                             Default |                                       Comments |
+============================+=====================================+================================================+
| ``MANGADRP_VER``           | ``v3_1_1`` (i.e., MPL-11)           | Version of the DRP, used for path construction |
+----------------------------+-------------------------------------+------------------------------------------------+
| ``MANGA_SPECTRO_REDUX``    | ``$HOME/MaNGA/redux``               | Root path for the reduced data                 |
+----------------------------+-------------------------------------+------------------------------------------------+
| ``MANGADAP_VER``           | ``mangadap.__version__``            | Version of the DAP, used for path construction |
+----------------------------+-------------------------------------+------------------------------------------------+
| ``MANGA_SPECTRO_ANALYSIS`` | ``$HOME/MaNGA/analysis``            | Root path for the analysis data                |
+----------------------------+-------------------------------------+------------------------------------------------+

These environmental variables can be added to, e.g., your
``.bash_profile`` file in your home directory or be included in a script
that is sourced when you want to run the DAP.  The lines added to your
``.bash_profile`` file could look something like this:

.. code-block:: bash

    export MANGA_SPECTRO_REDUX=/Volumes/MaNGA/redux
    export MANGADRP_VER=v3_1_1

    export MANGA_SPECTRO_ANALYSIS=/Volumes/MaNGA/analysis
    export MANGADAP_VER=3.1.0

.. note::

 * Importantly, note that ``$MANGADAP_VER`` is **only** used to set the
   path names, not to select the specific version of the DAP that
   should be used. The version of the DAP used is always the one
   installed by your python environment.
 * The DAP checks that these variables are defined *every time it is
   imported*.
 * Some of these same variables are defined by `Marvin`_. It is
   possible to have both Marvin and the DAP point to the same
   directory, but beware that this may mean that some of the files
   get overwritten!
 * Two additional variables (``$MANGACORE_VER`` and
   ``$MANGACORE_DIR``) are used in a specific mode of survey-level
   execution of the DAP. However, this is a niche usage mode and is
   effectively never used. See :ref:`execution-rundap`.
 * The DAP expects to find the DRP ``LOGCUBE`` *and* ``LOGRSS`` files
   in the directory
   ``$MANGA_SPECTRO_REDUX/$MANGADRP_VER/[PLATE]/stack``, where
   ``[PLATE]`` is the desired plate number. The ``LOGRSS`` files are
   required if you want to properly account for
   :ref:`spatialcovariance`. This path can be altered when executing
   the DAP.
 * The DAP expects to find/write data to
   ``$MANGA_SPECTRO_ANALYSIS/$MANGADRP_VER/$MANGADAP_VER``. This path
   can be altered when executing the DAP, but the subdirectory
   structure used by the DAP to organize its outputs within this root
   directory cannot currently be changed.


