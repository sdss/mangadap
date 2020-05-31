# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-
"""
Defines the class used to automate the execution of the MaNGA DAP at
Utah.

Revision history
----------------

    | **?? Nov 2014**: Architecture setup by J. Brownstein
        <joelbrownstein@astro.utah.edu>
    | **18 Nov 2014**: Handover to K. Westfall (KBW) for development.
        (Attempt to) conform style to `PEP 8`_, except that lines are 99
        characters long (docstrings still limited to 72 characters).
        Also (attempt to) conform to docstrings style according to `PEP
        257`_.
    | **01 Dec 2014**: (KBW) Committed to SVN
    | **02 Dec 2014**: (KBW) Changed drpver and idlutilsver to mangaver
    | **13 Jan 2015**: (KBW) Changed it back, added platetargets
    | **05 Feb 2015**: (KBW) Change to load module based on MPL
    | **19 Mar 2015**: (KBW) Major changes to allow for more control
        over the directory structure containing the DRP data and for the
        DAP results; following changed in the drpfile and drpcomplete
        classes.
    | **27 Aug 2015**: (KBW) Sphinx documentation; prep for MPL-4;
        accommodate changes in drpcomplete; removed arg as an attribute
        because it was only used in
        :func:`mangadap.scripts.rundap._read_arg`; despite not doing
        anything before, commented out _write_module_commands for the
        time-being; removed file_root() and parameter_file() functions
        in favor of adding/using functions from
        :mod:`mangadap.config.defaults`; version number set to 1_1_0.
    | **06 Oct 2015**: (KBW) Added functionality to select the types of
        additional output.  Minor changes to allow for DR13 QA plots.
    | **22 Oct 2015**: (KBW) Added check that the selected MPL versions match
        the current environmental versions.  Includes edit to
        :class:`mangadap.survey.mangampl.MaNGAMPL` to include
        MANGACORE_VER.  Changed default number of nodes to 9.  Added
        definition of idl command for use on sciama cluster at the ICG,
        Portsmouth.
    | **04 Nov 2015**: (KBW) Added capability to turn off main DAP
        processing call so that :class:`rundap` can be used to, e.g.,
        just create the plots for existing data.
    | **17 Feb 2016**: (KBW) Added try/except block for importing
        pbs.queue; converted old drpfile to new
        :class:`mangadap.drpfits.DRPFits`, and drpcomplete to new
        :class:`mangadap.survey.drpcomplete.DRPComplete`.
    | **13 May 2016**: (KBW) Substantial changes to accommodate
        conversion of DAP from IDL to python
    | **19 May 2016**: (KBW) Added command-line options for including
        log file
    | **11 Oct 2016**: (KBW) Version checking moved to
        :class:`mangadap.survey.mangampl.MaNGAMPL`.
    | **23 Aug 2017**: (KBW) Include ppxffit QA plots.
    | **29 Aug 2017**: (KBW) Include DAPall post-processing.
    | **16 Jan 2019**: (KBW) Major changes: Only one operation mode
        (daily, all, redo removed; code always run in "redo" mode).
        Remove old "prior" options.  Setup done using DRPall by default,
        but option to use mangacore/platetargets.  'CUBE' vs. 'RSS' no
        longer an option, only analyzes CUBEs.
    | **23 Jan 2019**: (KBW) Add fit QA plots

----

.. include license and copyright
.. include:: ../include/copy.rst

----

.. include common links, assuming primary doc root is up one directory
.. include:: ../include/links.rst
"""
import sys
import time
import shutil
import os
import glob
import warnings
import argparse

try:
    import pbs.queue
except:
    warnings.warn('Could not import pbs.queue!  Any cluster submission will fail!', ImportWarning)

import numpy

from mangadap import __version__
from mangadap.config import defaults
from mangadap.survey.drpcomplete import DRPComplete
from mangadap.survey.mangampl import MaNGAMPL
from mangadap.util.exception_tools import print_frame
from mangadap.util.parser import arginp_to_list
from mangadap.util.fileio import create_symlink
from mangadap.par.analysisplan import AnalysisPlanSet
from mangadap.proc.spatiallybinnedspectra import SpatiallyBinnedSpectra
from mangadap.proc.stellarcontinuummodel import StellarContinuumModel
from mangadap.proc.emissionlinemodel import EmissionLineModel


def parse_args(options=None, return_parser=False):

    # Declare the ArgumentParser object
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    # Read the optional run-mode arguments
    parser.add_argument("--clobber",
                        help="if all selected, will run dap for all plates/ifudesigns/modes "
                        " regardless of state", action="store_true", default=False)
    parser.add_argument('-v', '--verbose', action='count',
                        help='Set verbosity level for manga_dap; can be omitted and set up '
                             'to -vv', default=0)
    parser.add_argument("--quiet", help="suppress screen output", action="store_true",
                        default=False)
    parser.add_argument("--print_version", help="print DAP version and stop",
                        action="store_true", default=False)
       
    # These arguments are used to override default behavior
    parser.add_argument("--loose", help="Only throw warnings if the versioning is "
                        "not identically as it should be for the designated MPL",
                        action="store_true", default=False)
    parser.add_argument("--mplver", type=str, help="select MPL version to analyze",
                        default=None)
    parser.add_argument("--redux_path", type=str, help="main DRP output path", default=None)
    parser.add_argument("--dapver", type=str,
                        help="optional output version, different from product version",
                        default=None)
    parser.add_argument("--analysis_path", type=str, help="main DAP output path", default=None)

    parser.add_argument("--plan_file", type=str, help="parameter file with the MaNGA DAP "
                        "execution plan to use instead of the default" , default=None)
    parser.add_argument("--platelist", type=str, help="set list of plates to reduce",
                        default=None)
    parser.add_argument("--ifudesignlist", type=str, help="set list of ifus to reduce",
                        default=None)
    parser.add_argument("--list_file", type=str,
                        help='A file with the list of plates and ifudesigns to analyze',
                        default=None)
    parser.add_argument("--combinatorics", help="force execution of all permutations of the "
                        "provided lists", action="store_true", default=False)

    parser.add_argument('--sres_ext', type=str, default=None,
                        help='Spectral resolution extension to use.  Default set by '
                             'MaNGADataCube class.')
    parser.add_argument('--sres_fill', type=str, default=None,
                        help='If present, use interpolation to fill any masked pixels in the '
                             'spectral resolution vectors. Default set by MaNGADataCube '
                             'class.')
    parser.add_argument('--covar_ext', type=str, default=None,
                        help='Use this extension to define the spatial correlation matrix.  '
                             'Default set by MaNGADataCube class.')

    parser.add_argument("--use_plttargets", action="store_true", default=False,
                        help="Use platetargets files instead of the DRPall file to "
                             "generate the DRP complete database")
    parser.add_argument("--plttargets", type=str, help="path to plateTargets file(s); "
                        "if provided will force update to drpcomplete fits file", default=None)
    parser.add_argument('--on_disk', action='store_true', default=False,
                        help='When using the DRPall file to collate the data for input to '
                             'the DAP, search for available DRP files on disk instead of '
                             'using the DRPall file content.')

    parser.add_argument("--log", help="Have the main DAP executable produce a log file",
                        action="store_true", default=False)
    parser.add_argument("--no_proc", help="Do NOT perform the main DAP processing steps",
                        action="store_true", default=False)
    parser.add_argument("--no_plots", help="Do NOT create QA plots", action="store_true",
                        default=False)
    parser.add_argument("--post", help="Create/Submit the post-processing scripts",
                        action="store_true", default=False)
    parser.add_argument("--post_plots", action="store_true", default=False,
                        help="Create/Submit the post-processing plotting scripts")

    parser.add_argument("--dapall", action="store_true", default=False,
                        help='Wait for any individual plate-ifu processes to finish and then'
                             'update the DAPall file')

    # Read arguments specific to the cluster submission behavior
    parser.add_argument("--label", type=str, help='label for cluster job', default=None)
    parser.add_argument("--nodes", type=int, help='number of nodes to use in cluster',
                        default=1)
    parser.add_argument("--cpus", type=int,
                        help='number of cpus to use per node.  Default is to use all available'
                             '; otherwise, set to minimum of provided number and number of '
                             'processors per node', default=None)
    parser.add_argument("--fast", dest='qos', type=str, help='qos state', default=None)
    parser.add_argument("--umask", type=str, help='umask bit for cluster job', default='0027')
    parser.add_argument("--walltime", type=str, help='walltime for cluster job',
                        default='240:00:00')
    parser.add_argument("--toughness", dest='hard', action='store_false', default=True,
                        help='turn off hard keyword for cluster submission')

    parser.add_argument('--create', action='store_true', default=False,
                        help='use the pbs package to create the cluster scripts')
    parser.add_argument('--submit', action='store_true', default=False,
                        help='submit the scripts to the cluster')
    parser.add_argument('--progress', action='store_true', default=False,
                        help='instead of closing the script, report the progress of the '
                             'analysis on the cluster; this is required if you want to submit '
                             'the DAPall script immediately after completing the individual '
                             'cube analysis')

    parser.add_argument("--queue", dest='queue', type=str, help='set the destination queue',
                        default=None)

    if return_parser:
        return parser

    return parser.parse_args() if options is None else parser.parse_args(options)


class rundap:
    r"""
    Class used for automated execution of the MaNGA DAP via submission
    of scripts to a cluster.

    [ More complete explanation of the class here. ]

    .. note::
    
        This is purely a survey-level utility for running the DAP at
        Utah, including submitting jobs to the cluster.  No user-level
        usage of the DAP should be considered here, apart from not
        hindering such usage of the primary python programs!

    .. todo::

        - Have :func:`_read_arg` return the boolean for :attr:`version`
          instead of keeping it as an attribute.
        - Ditch ``console`` and :func:`_read_arg` in favor of the
          more typical ``parser`` and ``main`` pairing.

    Args:
        clobber (:obj:`bool`, optional):
            Clobber and redo any existing analysis.
        console (:obj:`bool`, optional):
            The class has been executed from the terminal with
            command-line arguments that should be parsed.
        quiet (:obj:`bool`, optional):
            Suppress output
        print_version (:obj:`bool`, optional):
            Print the class version and return.
        strictver (:obj:`bool`, optional):
            Strictly check the version requirements of the dependencies
            defined in :mangadap.survey.mangampl.MaNGAMPL` for this
            version of the DAP.  Default is True, meaning the code will
            raise an exception if the expected version are not the same
            as the environment.
        mplver (:obj:`str`, optional):
            MPL version to analyze.  Used with
            :class:`mangadap.survey.mangampl.MaNGAMPL` to set the DAP
            dependencies.  Ideally, this controls the modules that are
            loaded before execution of the DAP.  The default version
            (selected if mplver=None) is defined by
            :class:`mangadap.survey.mangampl.MaNGAMPL`.
        redux_path (:obj:`str`, optional):
            The top-level path with the DRP files used to override the
            default defined by
            :func:`mangadap.config.defaults.drp_redux_path`.
        dapver (:obj:`str`, optional):
            The DAP version to use for the analysis, used to override
            the default defined by
            :func:`mangadap.config.defaults.dap_version`.  *Only
            used to define the file paths, not used to switch to the
            specified code version!*
        analysis_path (:obj:`str`, optional):
            The top-level path for the DAP output files, used to
            override the default defined by
            :func:`mangadap.config.defaults.dap_analysis_path`.
        plan_file (:obj:`str`, optional):
            Name of the plan file to use for *all* DRP data in this run,
            used to override the default defined by
            :func:`mangadap.par.analysisplan.AnalysisPlanSet.default`.
        platelist (:obj:`int`, :obj:`str`, optional):
            List of plates to analyze; default is to search the
            :attr:`redux_path` for any DRP files; see
            :class:`mangadap.survey.drpcomplete.DRPComplete`.
        ifudesignlist (:obj:`int`, :obj:`str`, optional):
            List of ifudesigns to analyze: default is to search the
            :attr:`redux_path` for any DRP files; see
            :class:`mangadap.survey.drpcomplete.DRPComplete`.
        combinatorics (:obj:`bool`, optional): 
            Use all unique combinations of the entered
            plate/ifudesign/mode lists to create the full list of DRP
            files to analyze.
        list_file (:obj:`str`, optional):
            File with a list of plates and ifudesigns to analyze. The
            file must have two columns, with one plate-ifu
            combination per file row.
        sres_ext (:obj:`str`, optional):
            The extension to use when constructing the
            spectral resolution vectors for the MaNGA datacubes. See
            :func:`mangadap.datacube.manga.MaNGADataCube.spectral_resolution`.
        sres_fill (:obj:`bool`, optional):
            Fill masked spectral-resolution data by simple linear
            interpolation.
        covar_ext (:obj:`str`, optional):
            Extension in the MaNGA DRP CUBE file to use as the single
            spatial correlation matrix for all wavelength channels.
        use_platetargets (:obj:`bool`, optional): 
            Flag to use the plateTargets files, instead of the DRPall
            file, as the source of the input data need by the DAP.
        platetargets (:obj:`str`, :obj:`list`, optional):
            List of platetargets files to search through to find any
            given plate-ifudesign combination.  Default is returned as
            the first element in
            :func:`mangadap.config.defaults.plate_target_files`.
        on_disk (:obj:`bool`, optional):
            When searching for available files to analyze, search the
            DRP directory path instead of using the data in the DRPall
            file.  This is implicitly true when using the plateTargets
            files to collate the data used for the DAP input instead of
            the DRPall file (i.e., when `use_platetargets=True`).
        log (:obj:`bool`, optional):
            Flag to create the DAP log files.
        dapproc (:obj:`bool`, optional):
            Flag to execute the main DAP processing.
        pltifu_plots (:obj:`bool`, optional):
            Create the QA plate-ifu specific plots.
        post_process (:obj:`bool`, optional):
            Prepare (and possibly perform) the post processing steps:
            create the DAPall file and its associated QA plots, and the
            fit QA plots for each plate.
        post_plots (:obj:`bool`, optional):
            Create the post-processing QA plots.
        report_progress (:obj:`bool`, optional):
            Report the progress of the analysis.
        vebose (:obj:`int`, optional):
            Verbosity level between 0 and 2.
        label (:obj:`str`, optional):
            Label to use in cluster queue.
        nodes (:obj:`int`, optional):
            Number of cluster nodes to use.
        qos (:obj:`str`, optional):
            Select a specific processor (e.g. sdss-fast when the manga
            user).  This option is ignored if :attr:`q` is set.
        umask (:obj:`str`, optional):
            umask to set for output.
        walltime (:obj:`str`, optional):
            Maximum wall time for cluster job.
        hard (:obj:`bool`, optional):
            Same as hard keyword in cluster submission; see
            :func:`pbs.queue.commit`.  This keyword is in the
            **opposite** sense of the "--toughness" command-line
            options; i.e. including --toughness on the command line sets
            hard=False.
        create (:obj:`bool`, optional):
            Use the pbs package to create the scripts.
        submit (:obj:`bool`, optional):
            Submit all jobs to the cluster.
        queue (:obj:`str`, optional):
            Name of the destination queue.  When submitting jobs at
            Utah, this should **not** be set (leaving it at the default
            of None).  When submitting jobs at Portsmouth, this can be
            used to select either sciama1.q, cluster.q (default), or
            sciama3.q.

    Attributes:
        clobber (bool): **Only used with :attr:`all`**.  Clobber and
            redo any existing analysis.
        quiet (bool): Suppress output
        print_version (bool): Print the version and then return
        strictver (bool): Strictly check the version requirements for
            this version of the dap.  If True, an exception is raised if
            the expected version are not the same as the environment.
        mpl (:class:`mangadap.survey.mangampl.MaNGAMPL`): MPL version;
            see above.
        redux_path (str): The top-level path with the DRP files.
        dapver (str): The DAP version.
        analysis_path (str): The top-level path for the DAP output
            files.
        plan_file (str): Name of the plan file to use for ALL DRP data
            in this run.  Will be None if using the default plan.
        plan (:class:`AnalysisPlanSet`): The set of analysis plans
        daptypes (:obj:`list`): This list of daptypes, one per plan
        platelist (list): List of plates to analyze.
        ifudesignlist (list): List of ifudesigns to analyze.
        combinatorics (bool): Use all unique combinations of the entered
            plate/ifudesign/mode lists to create the full list of DRP
            files to analyze.
        use_platetargets (bool): See above
        platetargets (list): List of platetargets files to search
            through to find any given plate-ifudesign combination.
        on_disk (bool): See above
        log (str): See above
        dapproc (bool): Flag to execute the main DAP processing
        pltifu_plots (bool): Create the QA plate-ifu specific plots
        post_process (:obj:`bool`): See above
        post_plots (:obj:`bool`): See above
        report_progress (:obj:`bool`): See abve
        label (str): Label to use in cluster queue.
        nodes (int): Number of cluster nodes to use.
        qos (str): Specific processor to use (Utah specific)
        umask (str): umask to set for output.
        walltime (str): Maximum wall time for cluster job.
        hard (bool): Same as hard keyword in cluster submission; see
            :func:`pbs.queue.commit`.
        create (bool): Use the pbs package to create the submission
            scripts.
        submit (bool): Submit all jobs to the cluster.
        q (str): Queue destination (Portsmouth specific)
        drpc (:class:`mangadap.survey.drpcomplete.DRPComplete`):
            Database of the available DRP files and the parameters
            necessary to write the DAP par files.

    Raises:
        FileNotFoundError:
            Raised if a plate-ifu list file is provided but it does not
            exist, or if the plan file does not exist.
        ValueError:
            Raised if all processing steps are turned off, if the
            selected MPL is undefined; if the fast node is selected with
            nodes > 1; if multiple run modes are selected; or if the
            plate list, ifudesign list, or mode list are not provided
            during a redo mode or if they are provided with any other
            mode.
    """
    def __init__(self, clobber=None, console=None, quiet=False, print_version=False,
                 strictver=True, mplver=None, redux_path=None, dapver=None, analysis_path=None,
                 plan_file=None, platelist=None, ifudesignlist=None, combinatorics=False,
                 list_file=None, sres_ext=None, sres_fill=None, covar_ext=None,
                 use_platetargets=False, platetargets=None, on_disk=False, log=False,
                 dapproc=True, pltifu_plots=True, post_process=False, post_plots=False,
                 report_progress=False, verbose=0, label='mangadap', nodes=1, cpus=None, qos=None,
                 umask='0027', walltime='240:00:00', hard=True, create=False, submit=False,
                 queue=None):

        # Save run-mode options
        self.clobber = clobber
        self.quiet = quiet
        self.print_version = print_version

        # Override environment
        self.strictver = strictver
        self.mpl = mplver
        self.redux_path = None if redux_path is None else os.path.abspath(redux_path)
        self.dapver = dapver
        self.analysis_path = None if analysis_path is None else os.path.abspath(analysis_path)

        # List of files to analyze
        self.plan_file = None if plan_file is None else os.path.abspath(plan_file)
        if list_file is not None and not os.path.isfile(list_file):
            raise FileNotFoundError('No file: {0}'.format(list_file))
        if list_file is None:
            self.platelist = arginp_to_list(platelist, evaluate=True)
            self.ifudesignlist = arginp_to_list(ifudesignlist, evaluate=True)
        else:
            if platelist is not None or ifudesignlist is not None:
                warnings.warn('Provided file with list of files supercedes other input.')
            self.platelist, self.ifudesignlist = self._read_file_list(list_file)
        self.combinatorics = combinatorics
        self.sres_ext = sres_ext
        self.sres_fill = sres_fill
        self.covar_ext = covar_ext

        # Select if plateTargets files should be used to generate the
        # DRPComplete database, and possibly provide them directly 
        # TODO: This is likely now outdated
        self.use_platetargets = use_platetargets
        self.platetargets = arginp_to_list(platetargets)
        self.on_disk = on_disk

        # Set the options for output
        self.log = log
        self.dapproc = dapproc
        self.pltifu_plots = pltifu_plots
        self.post_process = post_process
        self.post_plots = post_plots

        # Report the progess of the cluster, otherwise just finish
        # script
        self.report_progress = report_progress
        self.verbose = verbose

        # Cluster queue keywords
        self.label = label
        self.nodes = nodes
        self.cpus = cpus
        self.qos = qos
        self.umask = umask
        self.walltime = walltime
        self.hard = hard
        self.submit = submit

        self.q = queue

        # Read and parse command-line arguments
        if console:
            self._read_arg()

        # Only print the version of the DAP
        if self.print_version:
            print('DAP Version: {0}'.format(__version__))
            return

        # If submitting to the cluster, make sure that scripts will be
        # created!
        if not self.create and self.submit:
            self.create = True

        # If setting qos, the number of nodes must be 1
        if self.qos is not None and self.nodes > 1:
            raise ValueError('When selecting the fast node, must have nodes=1!')

        # Check processing steps
        nproc = 0
        processes = []
        if self.dapproc:
            nproc += 1
            processes += [ 'Main DAP blocks' ]
        if self.pltifu_plots:
            nproc += 1
            processes += [ 'PLATEIFU-specific QA plots' ]
        if self.post_process:
            nproc += 1
            processes += [ 'DAPall' ]
        if self.post_plots:
            nproc += 1
            processes += [ 'DAPall and PLATE-specific QA plots' ]
        if nproc < 1:
            raise ValueError('No processing steps requested!')

        if (self.dapproc or self.pltifu_plots) and self.submit and self.post_process \
                and not self.report_progress:
            warnings.warn('If both processing and post processing, need to report progress.  '
                          'Setting to report progress.')
            self.report_progress = True

        # Make sure the selected MPL version is available
        try:
            self.mpl = MaNGAMPL(version=self.mpl, strictver=self.strictver)
        except:
            e = sys.exc_info()
            print_frame(e[0])
            raise ValueError('MPL is undefined: {0}'.format(e[1]))
        self.dapver = defaults.dap_version() if self.dapver is None else self.dapver

        # Set the output paths
        self.redux_path = defaults.drp_redux_path(self.mpl.drpver) if self.redux_path is None \
                                    else os.path.abspath(self.redux_path)
        self.analysis_path = defaults.dap_analysis_path(self.mpl.drpver, self.dapver) \
                                    if self.analysis_path is None \
                                    else os.path.abspath(self.analysis_path)

        # Set the subdirectory from which the scripts are sourced and
        # the various log files are written.  Currently based on start
        # UTC; e.g., analysis_path/log/19May2016T09.17.42UTC
        self.calling_path = os.path.join(self.analysis_path, 'log',
                                         time.strftime('%d%b%YT%H.%M.%SUTC',time.gmtime()))

        # Create the calling path. The existence test is executed, even
        # though time stamp should mean that calling_path never exists.
        # TODO: allow argument to set calling_path directly?
        if not os.path.isdir(self.calling_path):
            os.makedirs(self.calling_path)

        # Allow for the default plan.
        if self.plan_file is None:
            warnings.warn('No analysis plan was provided.  DAP executions will adopt the default '
                          'plan.')
        else:
            if not os.path.isfile(self.plan_file):
                raise FileNotFoundError('No file: {0}'.format(self.plan_file))
            _plan_file = self.plan_file
            self.plan_file = os.path.join(self.calling_path, os.path.split(self.plan_file)[1])
            if not os.path.isfile(self.plan_file):
                shutil.copy2(_plan_file, self.plan_file)

        # Construct the DAPTYPEs and check that the plans yeild a fully
        # unique list
        # TODO: This should probably be check done in AnalysisPlanSet
        self.plan = AnalysisPlanSet.default() if self.plan_file is None \
                        else AnalysisPlanSet.from_par_file(self.plan_file)
        self.daptypes = []
        for p in self.plan:
            bin_method = SpatiallyBinnedSpectra.define_method(p['bin_key'])
            sc_method = StellarContinuumModel.define_method(p['continuum_key'])
            el_method = None if p['elfit_key'] is None \
                                else EmissionLineModel.define_method(p['elfit_key'])
            self.daptypes += [defaults.dap_method(bin_method['key'],
                                                  sc_method['fitpar']['template_library_key'],
                                                  'None' if el_method is None else
                                                        el_method['continuum_tpl_key'])]
        if len(numpy.unique(self.daptypes)) != self.plan.nplans:
            raise ValueError('{0} do not yield unique DAPTYPEs.'.format(
                                'Default plans' if self.plan_file is None 
                                    else 'Plans in {0} '.format(self.plan_file)))

        # Alert the user of the versions to be used
        if self.q is not None:
            print('Attempting to submit to queue: {0}'.format(self.q))
        print('Versions: DAP:{0}, {1}'.format(self.dapver, self.mpl.mplver))
        self.mpl.show()
        self.mpl.verify_versions(quiet=self.quiet)
        print('Paths:')
        print('      REDUX: {0}'.format(self.redux_path))
        print('   ANALYSIS: {0}'.format(self.analysis_path))
        print('    CALLING: {0}'.format(self.calling_path))
        if self.plan_file is not None:
            print('  PLAN FILE: {0}'.format(self.plan_file))
        print('Processing steps added to scripts:')
        print('        {0}'.format(processes))

        # Create and update the DRPComplete file if necessary
        self.drpc = DRPComplete(platetargets=self.platetargets, drpver=self.mpl.drpver,
                                redux_path=self.redux_path, dapver=self.dapver,
                                analysis_path=self.analysis_path)
        
        # Update the DRPComplete list; force an update if platetarget
        # files are provided
        self.drpc.update(platelist=self.platelist, ifudesignlist=self.ifudesignlist,
                         combinatorics=self.combinatorics, force=(self.platetargets is not None),
                         use_platetargets=self.use_platetargets, on_disk=self.on_disk)

        # Hereafter, self.platelist and self.ifuplatelist are the INPUT
        # values, whereas self.drpc.platelist and self.drpc.ifuplatelist
        # are the actual DRP reductions to analyze (including
        # combinatorics and file existence tests).  See, e.g.,
        # selected_drpfile_list().

        # Find the rows in the DRPComplete database with the DRP file
        # meta data
        drpc_rows = self._selected_drpfile_rows()

        print('Number of DRP files to process: {0}'.format(len(drpc_rows)))
        if len(drpc_rows) == 0:
            # If coded correctly, this function should never return here.
            return

        # Construct the post-processing scripts before
        # constructing/submitting datacube processing scripts
        if self.post_process or self.post_plots:
            post_queue = self._build_post_queue(drpc_rows)

        # Process the relevant datacubes
        if self.dapproc or self.pltifu_plots:
            # Create the script files for the set of drpfiles, and
            # submit the scripts to the queue if self.submit is true
            proc_queue = self._build_proc_queue(drpc_rows)

            # Submit the queue to the cluster
            if self.create:
                proc_queue.commit(hard=self.hard, submit=self.submit)

            # If not reporting progress, done
            if not self.report_progress:
                return

            # Wait until the queue is finished
            running = True
            while running and self.submit:
                time.sleep(60) 
                percent_complete = proc_queue.get_percent_complete()
                print('Percent complete ... {0:.1f}%'.format(percent_complete), end='\r')
                if percent_complete == 100:
                    running = False
            print('Percent complete ... {0:.1f}%'.format(percent_complete))

        if not self.post_process and not self.post_plots:
            # Don't perform the post-processing
            return

        # Submit the post-processing queue to the cluster
        if self.create:
            post_queue.commit(hard=self.hard, submit=self.submit)

    # ******************************************************************
    #  UTILITY FUNCTIONS
    # ******************************************************************
    def _read_file_list(self, list_file):
        db = numpy.genfromtxt(list_file, dtype=object)
        if db.shape[1] != 2:
            raise ValueError('Input file should contain 2 columns, plate and ifu.')
        return db[:,0].astype(int).tolist(), db[:,1].astype(int).tolist()

    def _read_arg(self):
        """Read and interpret the terminal command-line arguments.

        Function uses the argparse module, which provide a --help option
        to list all available arguments.  The arguments mimic those
        required by the class constructor.
        """
        arg = parse_args()

        ################################################################
        # Assign properties based on arguments; the suitability of these
        # arguments is checked in __init__()

        # Run-mode options
        # Will OVERWRITE existing input from __init__()
        if arg.clobber is not None:
            self.clobber = arg.clobber
        if arg.quiet is not None:
            self.quiet = arg.quiet
        if arg.print_version is not None:
            self.print_version = arg.print_version

        # Set the versions to use
        # Will OVERWRITE existing input from __init__()
        self.strictver = not arg.loose
        if arg.mplver is not None:
            self.mpl = arg.mplver
        if arg.redux_path is not None:
            self.redux_path = arg.redux_path
        if arg.dapver is not None:
            self.dapver = arg.dapver
        if arg.analysis_path is not None:
            self.analysis_path = arg.analysis_path

        if arg.plan_file is not None:
            self.plan_file = os.path.abspath(arg.plan_file)

        if arg.list_file is not None and not os.path.isfile(arg.list_file):
            raise FileNotFoundError('No file: {0}'.format(arg.list_file))

        if arg.list_file is None:
            if arg.platelist is not None:
                self.platelist = arginp_to_list(arg.platelist, evaluate=True)
            if arg.ifudesignlist is not None:
                self.ifudesignlist = arginp_to_list(arg.ifudesignlist, evaluate=True)
        else:
            if arg.platelist is not None or arg.ifudesignlist is not None:
                warnings.warn('Provided file with list of files supercedes other input.')
            self.platelist, self.ifudesignlist = self._read_file_list(arg.list_file)

        self.combinatorics = arg.combinatorics
        self.sres_ext = arg.sres_ext
        self.sres_fill = arg.sres_fill
        self.covar_ext = arg.covar_ext
   
        # Set the plateTargets and NSA catalog path
        if arg.plttargets is not None:
            self.platetargets = arginp_to_list(arg.plttargets)
        self.use_platetargets = arg.use_plttargets
        self.on_disk = arg.on_disk

        # Overwrite any existing output options
        self.log = arg.log
        self.dapproc = not arg.no_proc
        self.pltifu_plots = not arg.no_plots
        self.post_process = arg.post
        self.post_plots = arg.post_plots
        if self.verbose < arg.verbose:
            self.verbose = arg.verbose

        # Set to report the progress of the cluster.  This is needed if
        # waiting to submit the DAPall construction script
        self.report_progress = arg.progress

        # Set queue keywords
        if arg.label is not None:
            self.label = arg.label
        if arg.nodes is not None:
            self.nodes = arg.nodes
        if arg.cpus is not None:
            self.cpus = arg.cpus
        if arg.qos is not None:
            self.qos = arg.qos
            self.nodes = 1          # Force the number of nodes to be 1
        if arg.umask is not None:
            self.umask = arg.umask
        if arg.walltime is not None:
            self.walltime = arg.walltime
        if arg.hard is not None:
            self.hard = arg.hard
        if arg.create is not None:
            self.create = arg.create
        if arg.submit is not None:
            self.submit = arg.submit

        # Specify the destination queue
        if arg.queue is not None:
            self.q = arg.queue

    def _check_paths(self, plate, ifudesign):
        """
        Check if the output paths exists for a given plate and ifudesign.
        If not, create them.  The necessary output paths are:

            - The calling path that contains a copy of the plan, the
              script input and output files, and the touch files
            - The main common path defined by
              :func:`mangadap.config.defaults.dap_common_path`.
            - The plan-based subdirectories based on the plan file

        Args:
            plate (:obj:`int`):
                Plate number
            ifudesign (:obj:`int`):
                ifudesign number
        """
        # Generate the calling path for this plate/ifudesign
        path = os.path.join(self.calling_path, str(plate), str(ifudesign))
        if not os.path.isdir(path):
            os.makedirs(path)

        # Generate the main common path
        path = defaults.dap_common_path(plate=plate, ifudesign=ifudesign, drpver=self.mpl.drpver,
                                        dapver=self.dapver, analysis_path=self.analysis_path)
        if not os.path.isdir(path):
            os.makedirs(path)

        # Check the reference directories, which creates the full plan
        # path if necessary
        for daptype in self.daptypes:
            path = defaults.dap_method_path(daptype, plate=plate, ifudesign=ifudesign, ref=True,
                                            drpver=self.mpl.drpver, dapver=self.dapver,
                                            analysis_path=self.analysis_path)
            if not os.path.isdir(path):
                os.makedirs(path)


    def _create_queue(self):
        """
        Create a queue instance.
        
        This function is written to be functional for submission to the
        cluster queues both at Utah and using the sciama cluster at
        Portsmouth.  The available queues are:

        +-----------+--------+------+-----------+
        |     Queue |  Nodes |  PPN |   GB/core |
        +===========+========+======+===========+
        | sciama1.q |     80 |   12 |         2 |
        +-----------+--------+------+-----------+
        | sciama2.q |     96 |   16 |         4 |
        +-----------+--------+------+-----------+
        | sciama3.q |     48 |   20 |       3.2 |
        +-----------+--------+------+-----------+

        One should **not** submit to sciama3.q without permission from
        Will Percival.  The default queue is cluster.q (sciama2.q); so
        basically only use :attr:`q` to select sciama1.q.

        Submissions to Utah should have their queue destination set to
        None.  However, the fast node can be selected using :attr:`qos`.
        """
        # Instantiate
        queue = pbs.queue(verbose=not self.quiet)

        if self.q == 'ember':
            # Submitting to Utah ember cluster
            ppn = 12
            cpus = ppn if self.cpus is None else min(self.cpus, ppn)
            walltime = self.walltime if int(self.walltime.split(':')[0]) < 72 else '72:00:00'
            queue.create(label=self.label, nodes=self.nodes, qos=self.qos, umask=self.umask,
                         walltime=walltime, ppn=ppn, cpus=cpus, partition='ember', alloc='sdss')
        elif self.q is not None:
            # All other self.q values expected for Portsmouth cluster,
            # sciama.  In this case, the number of nodes is queue
            # dependent, and qos is not set
            if self.q == 'sciama1.q':
                ppn = 12
            elif self.q == 'sciama3.q':
                ppn = 20
            else:
                ppn = 16
            cpus = ppn if self.cpus is None else min(self.cpus, ppn)
            queue.create(label=self.label, nodes=self.nodes, umask=self.umask,
                         walltime=self.walltime, queue=self.q, ppn=ppn, cpus=cpus)
        else:
            # self.q can be None when submitting to both the Portsmouth
            # and Utah clusters.  In this case, the default queue
            # destination and ppn is correct.  qos is also set, but this
            # should only be used when submitting to Utah.
            ppn = 16
            cpus = ppn if self.cpus is None else min(self.cpus, ppn)
            queue.create(label=self.label, nodes=self.nodes, qos=self.qos, umask=self.umask,
                         walltime=self.walltime, ppn=ppn, cpus=cpus)

        return queue

    # Files:
    #       - mangadap-{plate}-{ifudesign}-LOG{mode} = script file = *
    #       - *.ready = script/directory is ready for queue submission
    #       - *.queued = script submitted to queue
    #       - *.started = script assigned a node/processor and has started running
    #       - *.done = completed execution of the full script
    # Other:
    #       - *.par = parameter file
    #       - *.out, *.err = stdout and stderr output from the script
    def _build_proc_queue(self, drpc_rows):
        """
        Create a queue instance, write the script files for all CPUs
        using :func:`prepare_for_analysis` and append the job to the
        queue, and then submit the jobs to the cluster queue.  If
        :attr:`submit` is False, this is essentially a wrapper for
        :func:`prepare_for_analysis`.

        This function is written to be functional for submission to the
        cluster queues both at Utah and using the sciama cluster at
        Portsmouth.  The available queues are:

        +-----------+--------+------+-----------+
        |     Queue |  Nodes |  PPN |   GB/core |
        +===========+========+======+===========+
        | sciama1.q |     80 |   12 |         2 |
        +-----------+--------+------+-----------+
        | sciama2.q |     96 |   16 |         4 |
        +-----------+--------+------+-----------+
        | sciama3.q |     48 |   20 |       3.2 |
        +-----------+--------+------+-----------+

        One should **not** submit to sciama3.q without permission from
        Will Percival.  The default queue is cluster.q (sciama2.q); so
        basically only use :attr:`q` to select sciama1.q.

        Submissions to Utah should have their queue destination set to
        None.  However, the fast node can be selected using :attr:`qos`.

        .. todo::
            - This algorithm is slow because it has to search through
              the drpcomplete fits file each time to find the
              plate/ifudesign combo.  Should make this more efficient.

        Args:
            drpc_rows (`numpy.ndarray`):
                A vector with the rows in the
                :class:`mangadap.survey.drpcomplete.DRPComplete`
                database with the DRP files to analyze.
        """
        # If submitting to the queue, create the queue object.  The
        # creation of the queue object is system dependent; however,
        # once created, the queue object is treated identically for both
        # the Utah and Portsmouth clusters.
        if self.create:
            queue = self._create_queue()
       
        # Create the script files, regardless of whether or not they are
        # submitted to the queue, appending the scripts to the queue if
        # they are to be submitted
        ndrp = len(drpc_rows)
        for i, index in enumerate(drpc_rows):
            # Write the compute script (write_compute_script also checks
            # the path exists!)
            scriptfile, stdoutfile, stderrfile \
                    = self.write_compute_script(index, dapproc=self.dapproc,
                                                plots=self.pltifu_plots, clobber=self.clobber)

            # Set the status to ready
            plate = self.drpc['PLATE'][index]
            ifudesign = self.drpc['IFUDESIGN'][index]
            self.set_pltifu_status(plate, ifudesign, status='ready')

            if self.create:
                queue.append('source {0}'.format(scriptfile), outfile=stdoutfile,
                                  errfile=stderrfile)
                self.set_pltifu_status(plate, ifudesign, status='queued')
            else:
                print('Preparing for analysis...{0:.1f}%'.format((i+1)*100/ndrp), end='\r')

        if not self.create:
            print('Preparing for analysis...{0:.1f}%'.format(100))

        return queue if self.create else None

    def _build_post_queue(self, drpc_rows):
        """
        Build the queue that performs the post-processing steps.
        """
        # Create the queue
        if self.create:
            queue = self._create_queue()

        if self.post_process:
            # Write the DAPall script and set its status
            dapall_scr, dapall_out, dapall_err \
                    = self.write_dapall_script(plots=self.post_plots, clobber=self.clobber)
            self.set_status(dapall_scr, status='ready')
            if self.create:
                # Add it to the queue
                queue.append('source {0}'.format(dapall_scr), outfile=dapall_out,
                             errfile=dapall_err)
                self.set_status(dapall_scr, status='queued')

        # Write the plate QA plot scripts
        if self.post_plots:
            # Find the unique plates
            plates = numpy.unique(self.drpc['PLATE'][drpc_rows])
            for plt in plates:
                # Write the plot script for this plate
                scriptfile, stdoutfile, stderrfile \
                        = self.write_plate_qa_script(plt, clobber=self.clobber)
                self.set_status(scriptfile, status='ready')
                
                if self.create:
                    queue.append('source {0}'.format(scriptfile), outfile=stdoutfile,
                                 errfile=stderrfile)
                    self.set_status(scriptfile, status='queued')

        return queue if self.create else None

    # ******************************************************************
    # Management
    # ******************************************************************
    def set_pltifu_status(self, plate, ifudesign, mode=None, status='queued'):
        """
        Generate a touch file that signifies the status for a given
        plate/ifudesign/mode.  The root of the touch file is set by
        :func:`mangadap.config.defaults.dap_file_root`.

        Args:
            plate (int): Plate number
            ifudesign (int): IFU design number
            mode (str): (**Optional**) DRP 3D mode, 'RSS' or 'CUBE'
            status (str): (**Optional**) That status signifier for the
                touch file.
        """
        # Get the name of the status file        
        root = defaults.dap_file_root(plate, ifudesign, mode=mode)
        self.set_status(os.path.join(self.calling_path, str(plate), str(ifudesign), root),
                        status)
        
    def set_status(self, root, status='queued'):
        """
        Generate a general touch file.

        Args:
            root (str): Full path and root name of the touch file
            status (str): (**Optional**) That status signifier for the
                touch file.
        """
        # Touch the file by opening and closing it
        file = open('{0}.{1}'.format(root, status), 'w')
        file.close()

    def _selected_drpfile_rows(self):
        """
        Find the rows in the
        :class:`mangadap.survey.drpcomplete.DRPComplete` database with
        the selected plates and IFUs to analyze.

        Returns:
            `numpy.ndarray`: A vector with the rows in the DRPComplete
            object with the relevant DRP file data.
        """
        # Find the rows in the DRPComplete database with the specified
        # plate and IFUS

        # All available PLATEIFUs
        drpc_pltifu = numpy.array(['{0}-{1}'.format(p,i) 
                                    for p,i in zip(self.drpc['PLATE'], self.drpc['IFUDESIGN'])])
        # The PLATEIFUs to analyze
        this_pltifu = numpy.array(['{0}-{1}'.format(p,i) 
                                    for p,i in zip(self.drpc.platelist, self.drpc.ifudesignlist)])

        # The rows with the PLATEIFUs to analyze
        rows = numpy.array([numpy.where(drpc_pltifu == pi)[0][0] if pi in drpc_pltifu else -1 \
                                for pi in this_pltifu])
        if numpy.any(rows < 0):
            raise ValueError('Should be able to find all plates and ifus.')
        return rows

    def write_compute_script(self, index, dapproc=True, plots=True, clobber=False,
                              relative_symlink=True):
        """
        Write the MaNGA DAP script file that is sent to a single CPU to
        analyze a single DRP file with a given plate and ifudesign.

        Args:
            index (:obj:`int`, optional):
                Row in the
                :class:`mangadap.survey.drpcomplete.DRPComplete`
                database with the DRP file to analysis.
            dapproc (:obj:`bool`, optional):
                Flag to execute the main DAP processing.
            plots (:obj:`bool`, optional):
                Create the QA plots.
            clobber (:obj:`bool`, optional):
                Flag to clobber any existing files.
            relative_symlink (:obj:`bool`, optional):
                Set the symlink to the par file to be a relative path
                instead of the absolute path.

        Returns:
            :obj:`str`: Three strings with the name of the written
            script file, the file for the output sent to STDOUT, and the
            file for the output sent to STDERR.
        """
        # Get the plate and IFU from the DRPComplete database; alway use
        # the CUBE file
        plate = self.drpc['PLATE'][index]
        ifudesign = self.drpc['IFUDESIGN'][index]
        mode = 'CUBE'

        # Check that the path exists, creating it if not
        # TODO: More efficient way to do this?
        self._check_paths(plate, ifudesign)

        # Create the parameter file
#        parfile = defaults.dap_par_file(plate, ifudesign, mode, drpver=self.mpl.drpver,
#                                        dapver=self.dapver, analysis_path=self.analysis_path)
        cfgfile = defaults.dap_config(plate, ifudesign, drpver=self.mpl.drpver, dapver=self.dapver,
                                      analysis_path=self.analysis_path)

        # Write the par file if it doesn't exist
#        if not os.path.isfile(parfile) or clobber:
        if not os.path.isfile(cfgfile) or clobber:
            # clobber defaults to True
#            self.drpc.write_par(parfile, mode, index=index)
            self.drpc.write_config(cfgfile, index=index, sres_ext=self.sres_ext,
                                   sres_fill=self.sres_fill, covar_ext=self.covar_ext)
            # and create symlinks to it
            for daptype in self.daptypes:
                # Generate the ref subdirectory for this plan
                path = defaults.dap_method_path(daptype, plate=plate, ifudesign=ifudesign,
                                                ref=True, drpver=self.mpl.drpver,
                                                dapver=self.dapver,
                                                analysis_path=self.analysis_path)
#                create_symlink(parfile, path, relative_symlink=relative_symlink, clobber=clobber,
#                               quiet=True)
                create_symlink(cfgfile, path, relative_symlink=relative_symlink, clobber=clobber,
                               quiet=True)

        # Set the root path for the scripts, inputs, outputs, and logs
        _calling_path = os.path.join(self.calling_path, str(plate), str(ifudesign))
        scr_file_root = defaults.dap_file_root(plate, ifudesign)

        # Set the names for the script, stdout, and stderr files
        scriptfile = os.path.join(_calling_path, scr_file_root)
        stdoutfile = '{0}.out'.format(scriptfile)
        stderrfile = '{0}.err'.format(scriptfile)

        ################################################################
        # Script file already exists, so just return
        if os.path.exists(scriptfile) and not clobber:
            return scriptfile, stdoutfile, stderrfile

        # Main script components are:
        #   - (dapproc) Run the DAP processing code
        #   - (plots) Run the plotting scripts to plot the output
        # ONE OF THESE 2 MUST BE TRUE

        ################################################################
        # Open the script file and write the date as a commented header
        # line
        file = open(scriptfile, 'w')
        file.write('# Auto-generated batch file\n')
        file.write('# {0}\n'.format(time.strftime("%a %d %b %Y %H:%M:%S",time.localtime())))
        file.write('\n')

        # Create the started touch file
        startfile = '{0}.started'.format(scriptfile)
        file.write('touch {0}\n'.format(startfile))
        file.write('\n')

        # Command that runs the DAP
        # TODO: Define a "default" plan?
        if dapproc:
            command = 'OMP_NUM_THREADS=1 manga_dap -c {0} -r {1} -a {2}'.format(
                            cfgfile, self.redux_path, self.analysis_path)
            if self.plan_file is not None:
                command += ' -p {0}'.format(self.plan_file)
            if self.log:
                command += ' --log {0}.log'.format(scriptfile)
            if self.verbose > 0:
                command += (' -'+'v'*self.verbose)
            file.write('{0}\n'.format(command))
            file.write('\n')

        # Plotting scripts
        if plots:
            command = 'OMP_NUM_THREADS=1 dap_ppxffit_qa {0} {1} --analysis_path {2}'.format(
                            plate, ifudesign, self.analysis_path)
            if self.plan_file is not None:
                command += ' --plan_file {0}'.format(self.plan_file)
            file.write('{0}\n'.format(command))
            file.write('\n')

            command = 'OMP_NUM_THREADS=1 spotcheck_dap_maps {0} {1} --analysis_path {2}'.format(
                            plate, ifudesign, self.analysis_path)
            if self.plan_file is not None:
                command += ' --plan_file {0}'.format(self.plan_file)
            file.write('{0}\n'.format(command))
            file.write('\n')

            command = 'OMP_NUM_THREADS=1 dap_fit_residuals {0} {1} --analysis_path {2}'.format(
                            plate, ifudesign, self.analysis_path)
            if self.plan_file is not None:
                command += ' --plan_file {0}'.format(self.plan_file)
            file.write('{0}\n'.format(command))
            file.write('\n')

        # Touch the done file
        donefile = '{0}.done'.format(scriptfile)
        file.write('touch {0}\n'.format(donefile))
        file.write('\n')

        # Done
        file.close()

        # Return the script file, file for stdout, and file for stderr
        return scriptfile, stdoutfile, stderrfile
    
    def write_dapall_script(self, plots=True, clobber=False):
        """
        Write the script used to construct the DAPall file and its
        associated quality assessment plots.

        Args:
            plots (:obj:`bool`, optional):
                Create the QA plots. Default is True.
            clobber (:obj:`bool`, optional):
                Flag to clobber any existing files.

        Returns:
            :obj:`str`: Three strings with the name of the written
            script file, the file for the output sent to STDOUT, and the
            file for the output sent to STDERR.
        """
        # Check that the path exists, creating it if not
        if not os.path.isdir(self.calling_path):
            os.makedirs(self.calling_path)
        
        # Set the names for the script, stdout, and stderr files
        scriptfile = os.path.join(self.calling_path, 'build_dapall')
        stdoutfile = '{0}.out'.format(scriptfile)
        stderrfile = '{0}.err'.format(scriptfile)

        # Script file already exists, so just return
        if os.path.exists(scriptfile) and not clobber:
            return scriptfile, stdoutfile, stderrfile

        # Open the script file and write the date as a commented header
        # line
        file = open(scriptfile, 'w')
        file.write('# Auto-generated batch file\n')
        file.write('# {0}\n'.format(time.strftime("%a %d %b %Y %H:%M:%S",time.localtime())))
        file.write('\n')

        # Create the started touch file
        startfile = '{0}.started'.format(scriptfile)
        file.write('touch {0}\n'.format(startfile))
        file.write('\n')

        # Command that constructs the DAPall file
        command = 'OMP_NUM_THREADS=1 '
        command += 'construct_dapall --drpver {0} -r {1} --dapver {2} -a {3}'.format(
                        self.mpl.drpver, self.redux_path, self.dapver, self.analysis_path)
        if self.plan_file is not None:
            command += ' --plan_file {0}'.format(self.plan_file)
        if self.verbose > 0:
            command += (' -'+'v'*self.verbose )
        file.write('{0}\n'.format(command))
        file.write('\n')

        # Add the plotting commands
        if plots:
            command = 'OMP_NUM_THREADS=1 dapall_qa --drpver {0} --redux_path {1}'.format(
                            self.mpl.drpver, self.redux_path) \
                        + ' --dapver {0} --analysis_path {1}'.format(self.dapver,
                                                                     self.analysis_path)
            if self.plan_file is not None:
                command += ' --plan_file {0}'.format(self.plan_file)
            file.write('{0}\n'.format(command))
            file.write('\n')

        # Touch the done file
        donefile = '{0}.done'.format(scriptfile)
        file.write('touch {0}\n'.format(donefile))
        file.write('\n')

        file.close()
        ################################################################

        # Return the script file, file for stdout, and file for stderr
        return scriptfile, stdoutfile, stderrfile
    
    def write_plate_qa_script(self, plate, clobber=False):
        """
        Write the script used to create the plate fit QA plot.

        Args:
            plate (:obj:`int`):
                The plate number for the plot.
            clobber (:obj:`bool`, optional):
                Flag to clobber any existing files.

        Returns:
            :obj:`str`: Three strings with the name of the written
            script file, the file for the output sent to STDOUT, and the
            file for the output sent to STDERR.
        """
        # Check that the path exists, creating it if not
        path = os.path.join(self.calling_path, str(plate))
        if not os.path.isdir(path):
            os.makedirs(path)
        
        # Set the names for the script, stdout, and stderr files
        scriptfile = os.path.join(self.calling_path, str(plate), '{0}_fitqa'.format(plate))
        stdoutfile = '{0}.out'.format(scriptfile)
        stderrfile = '{0}.err'.format(scriptfile)

        # Script file already exists, so just return
        if os.path.exists(scriptfile) and not clobber:
            return scriptfile, stdoutfile, stderrfile

        # Open the script file and write the date as a commented header
        # line
        file = open(scriptfile, 'w')
        file.write('# Auto-generated batch file\n')
        file.write('# {0}\n'.format(time.strftime("%a %d %b %Y %H:%M:%S",time.localtime())))
        file.write('\n')

        # Create the started touch file
        startfile = '{0}.started'.format(scriptfile)
        file.write('touch {0}\n'.format(startfile))
        file.write('\n')

        # Command that creates the plots
        command = 'OMP_NUM_THREADS=1 dap_plate_fit_qa {0} --analysis_path {1}'.format(
                        plate, self.analysis_path)
        if self.plan_file is not None:
            command += ' --plan_file {0}'.format(self.plan_file)
        file.write('{0}\n'.format(command))
        file.write('\n')

        # Touch the done file
        donefile = '{0}.done'.format(scriptfile)
        file.write('touch {0}\n'.format(donefile))
        file.write('\n')

        file.close()
        ################################################################

        # Return the script file, file for stdout, and file for stderr
        return scriptfile, stdoutfile, stderrfile
    


