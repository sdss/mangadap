# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-
"""
Defines the class used to automate the execution of the MaNGA DAP at
Utah.

*License*:
    Copyright (c) 2015, SDSS-IV/MaNGA Pipeline Group
        Licensed under BSD 3-clause license - see LICENSE.rst

*Source location*:
    $MANGADAP_DIR/python/mangadap/survey/rundap.py

*Imports and python version compliance*:
    ::


*Class usage examples*:
    Add some usage comments here!

.. todo::
    Add verbosity level

*Revision history*:
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
        :func:`mangadap.survey.rundap._read_arg`; despite not doing
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

.. _PEP 8: https://www.python.org/dev/peps/pep-0008
.. _PEP 257: https://www.python.org/dev/peps/pep-0257
.. _argparse.Namespace: https://docs.python.org/3/library/argparse.html#argparse.Namespace
.. _argparse.ArgumentParser: https://docs.python.org/3/library/argparse.html#argparse.ArgumentParser

"""

from __future__ import division
from __future__ import print_function
from __future__ import absolute_import
from __future__ import unicode_literals

import sys
if sys.version > '3':
    long = int

import subprocess
import time
import shutil
import os
import numpy
import glob
import warnings
from argparse import ArgumentParser

try:
    import pbs.queue
except:
    warnings.warn('Could not import pbs.queue!  Any cluster submission will fail!', ImportWarning)

from .drpcomplete import DRPComplete
from ..drpfits import DRPFits
from ..config.defaults import default_redux_path, default_drp_directory_path
from ..config.defaults import default_analysis_path, default_dap_common_path
from ..config.defaults import default_dap_method, default_dap_method_path, default_dap_file_root
from ..config.defaults import default_dap_par_file, default_dap_plan_file
from ..util.exception_tools import print_frame
from ..util.parser import arginp_to_list
from .mangampl import MaNGAMPL
from ..par.analysisplan import AnalysisPlanSet
from . import util

__author__ = 'Kyle Westfall'

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

    Args:
        daily (bool): (**Optional**) Execute daily update of analyses
        all (bool): (**Optional**) Analyze all existing DRP data.
        clobber (bool): (**Optional**) **Only used with :attr:`all`**.
            Clobber and redo any existing analysis.
        redo (bool): (**Optional**) Redo an existing set of analyses;
            must be accompanied by a list of plates, ifus, and modes to
            analyze.  **clobber** is implicitly true.
        console (bool): (**Optional**) The class has been executed from
            the terminal with command-line arguments that should be
            parsed.
        quiet (bool): (**Optional**) Suppress output
        version (bool): (**Optional**) Print the class version and
            return
        mplver (str): (**Optional**) MPL version to analyze.  Used with
            :class:`mangadap.survey.mangampl.MaNGAMPL` to set the
            idlutils version and the DRP version.  Ideally, this
            controls the modules that are loaded before execution of the
            DAP.  The default version (selected if mplver=None) is
            defined by :class:`mangadap.survey.mangampl.MaNGAMPL`.

            .. warning::

                This functionality is **NOT** actually implemented.
                Currently, the DAP will run with whatever set of modules
                are currently loaded when the jobs are submitted to the
                cluster.
                
        redux_path (str): (**Optional**) The top-level path with the DRP
            files used to override the default defined by
            :func:`mangadap.config.defaults.default_redux_path`.
        dapver (str): (**Optional**) The DAP version to use for the
            analysis, used to override the default defined by
            :func:`mangadap.config.defaults.default_dap_version`.
        analysis_path (str): (**Optional**) The top-level path for the
            DAP output files, used to override the default defined by
            :func:`mangadap.config.defaults.default_analysis_path`.
        plan_file (str): (**Optional**) Name of the plan file to use for
            *all* DRP data in this run, used to override the default
            defined by
            :func:`mangadap.config.defaults.default_dap_plan_file`.
        platelist (int, str): (**Optional**) List of plates to analyze;
            default is to search the :attr:`redux_path` for any DRP
            files; see :class:`mangadap.survey.drpcomplete.DRPComplete`.
        ifudesignlist (int, str): (**Optional**) List of ifudesigns to
            analyze: default is to search the :attr:`redux_path` for any
            DRP files; see
            :class:`mangadap.survey.drpcomplete.DRPComplete`.
        modelist (str): (**Optional**) DRP 3D modes ('RSS' or 'CUBE') to
            analyze; default is to analyze both.
        combinatorics (bool): (**Optional**) Use all unique combinations
            of the entered plate/ifudesign/mode lists to create the full
            list of DRP files to analyze.  See
            :func:`mangadap.drpfits.drpfits_list`.
        prior_mode (str): (**Optional**) Used to construct a plate-ifu
            specific prior fits file name to replace an existing value
            in the plan file.  This sets the mode ('RSS' or 'CUBE') of
            the fits file.
        prior_bin (str): (**Optional**) Used to construct a plate-ifu
            specific prior fits file name to replace an existing value
            in the plan file.  This sets the binning type of the fits
            file.
        prior_iter (int): (**Optional**) Used to construct a plate-ifu
            specific prior fits file name to replace an existing value
            in the plan file.  This sets the iteration number of the
            fits file.
        prior_old (str): (**Optional**) The value of the existing prior
            in the plan file that should be replaced by the prior fits
            file constructed using the other "prior_*" attributes.
        platetargets (str, list): (**Optional**) List of platetargets
            files to search through to find any given plate-ifudesign
            combination.  Default is returned as the first element in
            :func:`mangadap.config.defaults.default_plate_target_files`.
        dapproc (bool): (**Optional**) Flag to execute the main DAP
            processing. Default is True.
        plots (bool): (**Optional**) Create the QA plots. Default is
            True.
        label (str): (**Optional**) Label to use in cluster queue.
            Default is mangadap and the run mode (daily, etc).
        nodes (int): (**Optional**) Number of cluster nodes to use.
            Default is 9.
        qos (str): (**Optional**) Select a specific processor (only None
            for daily run).  This option is ignored if :attr:`q` is set.
        umask (str): (**Optional**) umask to set for output. Default is
            0027.
        walltime (str): (**Optional**) Maximum wall time for cluster
            job.  Default is 10 days ('240:00:00').
        hard (bool): (**Optional**) Same as hard keyword in cluster
            submission; see :func:`pbs.queue.commit`.

            .. note::

                This keyword is in the **opposite** sense of the
                "--toughness" command-line options; i.e. including
                --toughness on the command line sets hard=False.

        create (bool): (**Optional**) Use the pbs package to create the
            scripts. Default is False.
        submit (bool): (**Optional**) Submit all jobs to the cluster.
            Default is False.
        queue (str): (**Optional**) Name of the destination queue.  When
            submitting jobs at Utah, this should **not** be set (leaving
            it at the default of None).  When submitting jobs at
            Portsmouth, this can be used to select either sciama1.q,
            cluster.q (default), or sciama3.q.

    Attributes:
        daily (bool): Execute daily update of analyses
        all (bool): Analyze all existing DRP data.
        clobber (bool): **Only used with :attr:`all`**.  Clobber and
            redo any existing analysis.
        redo (bool): Redo an existing set of analyses; must be
            accompanied by a list of plates, ifus, and modes to analyze.
            In this case, :attr:`clobber` is implicitly true.
        quiet (bool): Suppress output
        version (bool): Print the class version and return.  Needed as
            an attribute so that it can be read and returned from the
            command line.
        mpl (:class:`mangadap.survey.mangampl.MaNGAMPL`): MPL version;
            see above.
        redux_path (str): The top-level path with the DRP files.
        dapver (str): The DAP version.
        analysis_path (str): The top-level path for the DAP output
            files.
        plan_file (str): Name of the plan file to use for ALL DRP data
            in this run.
        platelist (list): List of plates to analyze.
        ifudesignlist (list): List of ifudesigns to analyze.
        modelist (list): DRP 3D modes ('RSS' and/or 'CUBE') to analyze.
        combinatorics (bool): Use all unique combinations of the entered
            plate/ifudesign/mode lists to create the full list of DRP
            files to analyze.
        prior_mode (str): Used to construct a plate-ifu specific prior
            fits file name to replace an existing value in the plan
            file.  This sets the mode ('RSS' or 'CUBE') of the fits
            file.
        prior_bin (str): Used to construct a plate-ifu specific prior
            fits file name to replace an existing value in the plan
            file.  This sets the binning type of the fits file.
        prior_iter (int): Used to construct a plate-ifu specific prior
            fits file name to replace an existing value in the plan
            file.  This sets the iteration number of the fits file.
        prior_old (str): The value of the existing prior in the plan
            file that should be replaced by the prior fits file
            constructed using the other "prior_*" attributes.
        platetargets (list): List of platetargets files to search
            through to find any given plate-ifudesign combination.
        dapproc (bool): Flag to execute the main DAP processing
        plots (bool): Create the QA plots
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

    """
    def __init__(self,
                # Run mode options
                daily=None, all=None, clobber=None, redo=None,
                # STDIO options
                console=None, quiet=None, version=None,
                # Override default environmental variables
                mplver=None, redux_path=None, dapver=None, analysis_path=None, 
                # Definitions used to set files to process
                plan_file=None, platelist=None, ifudesignlist=None, modelist=None,
                combinatorics=False,
                prior_mode=None, prior_bin=None, prior_iter=None, prior_old=None,
                # Databases with input parameter information
                platetargets=None,
                # Flags for script contents
                log=False, dapproc=True, plots=True, verbose=0,
                # Cluster options
                label='mangadap', nodes=9, qos=None, umask='0027',walltime='240:00:00', hard=True,
                create=False, submit=False, queue=None):

        # Save run-mode options
        self.daily = daily
        self.all = all
        self.clobber = clobber
        self.redo = redo
        self.quiet = quiet
        self.version = version

        # Override environment
        self.mpl = mplver
        self.redux_path = redux_path

        self.dapver = util.product_version(simple=True, product='mangadap') if dapver is None \
                                                                            else dapver
        self.analysis_path = analysis_path

        # List of files to analyze
        self.plan_file = plan_file
        self.prior_mode = prior_mode
        self.prior_bin = prior_bin
        self.prior_iter = prior_iter
        self.prior_old = prior_old
        self.platelist = arginp_to_list(platelist, evaluate=True)
        self.ifudesignlist = arginp_to_list(ifudesignlist, evaluate=True)
        self.modelist = arginp_to_list(modelist)
        self.combinatorics = combinatorics

        # plateTargets file(s) to be used by the DRPComplete class,
        # which can be None
        self.platetargets = arginp_to_list(platetargets)

        # Set the options for output
        self.log = log
        self.dapproc = dapproc
        self.plots = plots
        self.verbose = verbose

        # Cluster queue keywords
        self.label = label
        self.nodes = nodes
        self.qos = qos
        self.umask = umask
        self.walltime = walltime
        self.hard = hard
        self.submit = submit

        self.q = queue

        # Read and parse command-line arguments
        if console:
            self._read_arg()

        # Make sure there's something to be done
        nproc = 0
        processes = []
        if self.dapproc:
            nproc += 1
            processes += [ 'Main DAP blocks' ]
        if self.plots:
            nproc += 1
            processes += [ 'QA plots' ]
        if nproc < 1:
            raise ValueError('No processing steps requested!')

        # Only print the version of the DAP
        if self.version:
            print('This is version 2.0')
            return

        # Make sure the selected MPL version is available
        try:
            self.mpl = MaNGAMPL(self.mpl)
        except:
            e = sys.exc_info()
            print_frame(e[0])
            raise ValueError('MPL is undefined: {0}'.format(e[1]))

        # Set the output paths
        self.redux_path = default_redux_path(self.mpl.drpver) if self.redux_path is None \
                                                              else str(self.redux_path)
        self.analysis_path = default_analysis_path(self.mpl.drpver, self.dapver) \
                             if self.analysis_path is None else str(self.analysis_path)

        # Set the subdirectory from which the scripts are sourced and
        # the various log files are written.  Currently based on start
        # UTC; e.g., analysis_path/log/19May2016T09.17.42UTC
        self.calling_path = os.path.join(self.analysis_path, 'log',
                                         time.strftime('%d%b%YT%H.%M.%SUTC',time.gmtime()))

        # Make sure there is an analysis plan
        if self.plan_file is None:
            self.plan_file = default_dap_plan_file(drpver=self.mpl.drpver, dapver=self.dapver,
                                                   analysis_path=self.analysis_path)
        if not os.path.isfile(self.plan_file):
            raise FileNotFoundError('No file: {0}'.format(self.plan_file))

        # Create the calling path and copy the plan file into the top
        # level of that directory; existence tests are run, even though
        # time stamp should mean that calling_path never exists
        # TODO: allow argument to set calling_path directly?
        if not os.path.isdir(self.calling_path):
            os.makedirs(self.calling_path)
        _plan_file = self.plan_file
        self.plan_file = os.path.join(self.calling_path, self.plan_file.split('/')[-1])
        if not os.path.isfile(self.plan_file):
            shutil.copy2(_plan_file, self.plan_file)

        # Alert the user of the versions to be used
        if self.q is not None:
            print('Attempting to submit to queue: {0}'.format(self.q))
        print('Versions: DAP:{0}, {1}'.format(self.dapver, self.mpl.mplver))
        self.mpl.show()
        print('Paths:')
        print('      REDUX: {0}'.format(self.redux_path))
        print('   ANALYSIS: {0}'.format(self.analysis_path))
        print('    CALLING: {0}'.format(self.calling_path))
        print('  PLAN FILE: {0}'.format(self.plan_file))
        print('Processing steps added to scripts:')
        print('        {0}'.format(processes))

        # Check the environment matches the selected MPL
        # These will throw KeyErrors if the appropriate environmental variables do not exist
        try:
            accessver_env = os.environ['SDSS_ACCESS_DIR'].split('/')[-1]
        except KeyError:
            accessver_env = None
        idlver_env = os.environ['IDLUTILS_DIR'].split('/')[-1]

        print('Current environment: ')
        print('    SDSS_ACCESS_VER: {0}'.format(accessver_env))
        print('       IDLUTILS_VER: {0}'.format(idlver_env))
        print('      MANGACORE_VER: {0}'.format(os.environ['MANGACORE_VER']))
        print('       MANGADRP_VER: {0}'.format(os.environ['MANGADRP_VER']))
        print('       MANGADAP_VER: {0}'.format(os.environ['MANGADAP_VER']))

#        print('mpl.accessver: {0}'.format(self.mpl.accessver))
#        print('mpl.idlver: {0}'.format(self.mpl.idlver))
#        print('mpl.corever: {0}'.format(self.mpl.corever))
#        print('mpl.drpver: {0}'.format(self.mpl.drpver))

        if self.mpl.accessver != accessver_env:
            print('mpl.accessver: {0}'.format(self.mpl.accessver))
            print('SDSS_ACCESS_VER: {0}'.format(accessver_env))
            raise EnvironmentError('MPL SDSS_ACCESS version does not match current environment;'
                                   ' see above')
        if self.mpl.idlver != idlver_env:
            print('mpl.idlver: {0}'.format(self.mpl.idlver))
            print('IDLUTILS_VER: {0}'.format(idlver_env))
            raise EnvironmentError('MPL IDLUTILS version does not match current environment;'
                                   ' see above')
        if self.mpl.corever != os.environ['MANGACORE_VER']:
            print('mpl.corever: {0}'.format(self.mpl.corever))
            print('MANGACORE_VER: {0}'.format(os.environ['MANGACORE_VER']))
            raise EnvironmentError('MPL CORE version does not match current environment; see above')
        if self.mpl.drpver != os.environ['MANGADRP_VER']:
            print('mpl.drpver: {0}'.format(self.mpl.drpver))
            print('MANGADRP_VER: {0}'.format(os.environ['MANGADRP_VER']))
            raise EnvironmentError('MPL DRP version does not match current environment; see above')

        # If submitting to the cluster, make sure that scripts will be
        # created!
        if not self.create and self.submit:
            self.create = True

        # If setting qos, the number of nodes must be 1
        if self.qos is not None and self.nodes > 1:
            raise ValueError('When selecting the fast node, must have nodes=1!')

        # Check that something is to be done
        nrun = 0
        if self.daily: nrun += 1
        if self.all: nrun += 1
        if self.redo: nrun += 1
        if nrun == 0 or nrun > 1:
            raise ValueError('Must request one (and only one) run mode!')

        # Check argument combination makes sense
        if self.clobber and not self.all:
            warnings.warn('Clobber keyword can only be set if all specified (-a,--all)!')
            self.clobber = None

        # If redo, platelist and ifudesignlist MUST be provided
        if self.redo and not self.platelist and not self.ifudesignlist and not self.modelist:
            raise ValueError('Must provide platelist, ifudesignlist, and modelist for redo!')

        # If not redo, platelist, ifudesignlist, and modelist should be
        # ignored
        # TODO: Does an error need to be thrown here, or just warn user.
        if not self.redo and (self.platelist or self.ifudesignlist or self.modelist):
            raise ValueError('platelist/ifudesignlist/modelist only allowed for redo!')

        # If a plan file is provided, make sure it exists
        if self.plan_file is not None and not os.path.exists(self.plan_file):
            raise ValueError('Provided plan file does not exist.')

        # If the prior is to be replaced, make sure all four arguments
        # are provided
        npri = 0
        if self.prior_mode is not None:
            npri += 1
        if self.prior_bin is not None:
            npri += 1
        if self.prior_iter is not None:
            npri += 1
        if self.prior_old is not None:
            npri += 1

        if npri > 0 and npri != 4:
            raise ValueError('To define prior, must provided mode, bin, iter, and old value!')
        # From now on can decide if prior exists by just testing one of
        # the four options

        if npri == 4 and self.plan_file is None:
            raise ValueError('To define prior, must provide a plan file to edit.')

        # If running all or daily, make sure lists are set to None such
        # that DRPComplete will search for all available DRP files and
        # make sure that the needed information is up-to-date
        if self.all or self.daily:
            self.platelist = None
            self.ifudesignlist = None
            self.modelist = None

        # Create and update the DRPComplete file if necessary
        self.drpc = DRPComplete(platetargets=self.platetargets, drpver=self.mpl.drpver,
                                redux_path=self.redux_path, dapver=self.dapver,
                                analysis_path=self.analysis_path)
        
        # Update the DRPComplete list; force an update if platetarget
        # files are provided
        self.drpc.update(platelist=self.platelist, ifudesignlist=self.ifudesignlist,
                         combinatorics=self.combinatorics, force=(self.platetargets is not None))

        # Hereafter, self.platelist and self.ifuplatelist are the INPUT
        # values, whereas self.drpc.platelist and self.drpc.ifuplatelist
        # are the actual DRP reductions to analyze (including
        # combinatorics and file existence tests).  These should be
        # combined with self.modelist to get the DRP files to analyze.
        # See, e.g., selected_drpfile_list().

        # Prep the verbosity of the queue if submitting
        if self.create:
            self.queue = pbs.queue(verbose=not self.quiet)

        # Run the selected mode
        if self.daily:
            self.run_daily()
        if self.all:
            self.run_all(clobber=self.clobber)
        if self.redo:
            self.run_redo()


    # ******************************************************************
    #  UTILITY FUNCTIONS
    # ******************************************************************
    def _read_arg(self):
        """Read and interpret the terminal command-line arguments.

        Function uses the argparse module, which provide a --help option
        to list all available arguments.  The arguments mimic those
        required by the class constructor.
        """
        # Declare the ArgumentParser object
        parser = ArgumentParser()

        # Set the run mode to be mutually exclusive
        mode = parser.add_mutually_exclusive_group(required=True)

        # Read the mode aruments
        mode.add_argument("--daily", help="runs dap for next mjd (as an after burner to drp)",
                          action="store_true")
        mode.add_argument("--all",
                          help="runs dap for all plates/ifudesigns/modes not currently running "
                          "or done", action="store_true")
        mode.add_argument("--redo",
                          help="runs dap for specified platelist/ifudesignlist/modelist "
                          "regardless of state", action="store_true")

        # Read the optional run-mode arguments
        parser.add_argument("--clobber",
                            help="if all selected, will run dap for all plates/ifudesigns/modes "
                            " regardless of state", action="store_true", default=False)
        parser.add_argument('-v', '--verbose', action='count',
                            help='Set verbosity level for manga_dap; can be omitted and set up '
                                 'to -vv', default=0)
        parser.add_argument("--quiet", help="suppress screen output", action="store_true",
                            default=False)
        parser.add_argument("--version", help="rundap version", action="store_true", default=False)
       
        # These arguments are used to override default behavior
        parser.add_argument("--mplver", type=str, help="select MPL version to analyze",
                            default=None)
        parser.add_argument("--redux_path", type=str, help="main DRP output path", default=None)
        parser.add_argument("--dapver", type=str,
                            help="optional output version, different from product version",
                            default=None)
        parser.add_argument("--analysis_path", type=str, help="main DAP output path", default=None)

        parser.add_argument("--plan_file", type=str, help="parameter file with the MaNGA DAP "
                            "execution plan to use instead of the default" , default=None)
        parser.add_argument("--prior_mode", type=str, help="mode type to use for prior",
                            default=None)
        parser.add_argument("--prior_bin", type=str, help="bin type to use for prior",
                            default=None)
        parser.add_argument("--prior_iter", type=str, help="iteration to use for prior",
                            default=None)
        parser.add_argument("--prior_old", type=str,
                            help="old value to replace in existing plan file", default=None)
        parser.add_argument("--platelist", type=str, help="set list of plates to reduce",
                            default=None)
        parser.add_argument("--ifudesignlist", type=str, help="set list of ifus to reduce",
                            default=None)
        parser.add_argument("--modelist", type=str, help="set list of DRP output modes to reduce"
                            " (CUBE or RSS)", default=None)
        parser.add_argument("--combinatorics", help="force execution of all permutations of the "
                            "provided lists", action="store_true", default=False)

        parser.add_argument("--plttargets", type=str, help="path to plateTargets file(s); "
                            "if provided will force update to drpcomplete fits file", default=None)

        parser.add_argument("--log", help="Have the main DAP executable produce a log file",
                            action="store_true", default=False)
        parser.add_argument("--no_proc", help="Do NOT perform the main DAP processing steps",
                            action="store_true", default=False)
        parser.add_argument("--no_plots", help="Do NOT create QA plots", action="store_true",
                            default=False)


        # Read arguments specific to the cluster submission behavior
        parser.add_argument("--label", type=str, help='label for cluster job', default='mangadap')
        parser.add_argument("--nodes", type=int, help='number of nodes to use in cluster',
                            default=18)
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

        parser.add_argument("--queue", dest='queue', type=str, help='set the destination queue',
                            default=None)

        
        # Finally parse the full set and set it to its own container
        arg = parser.parse_args()

        ################################################################
        # Assign properties based on arguments; the suitability of these
        # arguments is checked in __init__()

        # Run-mode
        # Will OVERWRITE existing input from __init__()
        if arg.daily is not None:
            self.daily = arg.daily
        if arg.all is not None:
            self.all = arg.all
        if arg.redo is not None:
            self.redo = arg.redo

        # Run-mode options
        # Will OVERWRITE existing input from __init__()
        if arg.clobber is not None:
            self.clobber = arg.clobber
        if arg.quiet is not None:
            self.quiet = arg.quiet
        if arg.version is not None:
            self.version = arg.version

        # Set the versions to use
        # Will OVERWRITE existing input from __init__()
        if arg.mplver is not None:
            self.mpl = arg.mplver
        if arg.redux_path is not None:
            self.redux_path = arg.redux_path
        if arg.dapver is not None:
            self.dapver = arg.dapver
        if arg.analysis_path is not None:
            self.analysis_path = arg.analysis_path

        if arg.plan_file is not None:
            self.plan_file = arg.plan_file

        if arg.prior_mode is not None:
            self.prior_mode = arg.prior_mode
        if arg.prior_bin is not None:
            self.prior_bin = arg.prior_bin
        if arg.prior_iter is not None:
            self.prior_iter = arg.prior_iter
        if arg.prior_old is not None:
            self.prior_old = arg.prior_old

        if arg.platelist is not None:
            self.platelist = arginp_to_list(arg.platelist, evaluate=True)
        if arg.ifudesignlist is not None:
            self.ifudesignlist = arginp_to_list(arg.ifudesignlist, evaluate=True)
        if arg.modelist is not None:
            self.modelist = arginp_to_list(arg.modelist)
        self.combinatorics = arg.combinatorics
   
        # Set the plateTargets and NSA catalog path
        if arg.plttargets is not None:
            self.platetargets = arginp_to_list(arg.plttargets)

        # Overwrite any existing output options
        # NOTE: This execution means that all 4 steps will default to
        # True if not specified on the command line.  This matches the
        # default set by __init__()
        self.log = arg.log
        self.dapproc = not arg.no_proc
        self.plots = not arg.no_plots
        if self.verbose < arg.verbose:
            self.verbose = arg.verbose

        # Set queue keywords
        if arg.umask is not None:
            self.umask = arg.umask
        if arg.nodes is not None:
            self.nodes = arg.nodes
        if arg.walltime is not None:
            self.walltime = arg.walltime
        if arg.qos is not None:
            self.qos = arg.qos
            self.nodes = 1          # Force the number of nodes to be 1
        if arg.umask is not None:
            self.umask = arg.umask
        if arg.hard is not None:
            self.hard = arg.hard
        if arg.create is not None:
            self.create = arg.create
        if arg.submit is not None:
            self.submit = arg.submit

        # Specify the destination queue
        if arg.queue is not None:
            self.q = arg.queue

#       print(self.q)
#       print('QOS: {0}'.format(self.qos))


    def _check_paths(self, plate, ifudesign):
        """
        Check if the output paths exists for a given plate and ifudesign.
        If not, create them.  The necessary output paths are:

            - The callind path that contains a copy of the plan, the
              script input and output files, and the touch files
            - The main common path defined by
              :func:`mangadap.config.defaults.default_dap_common_path`.
            - The plan-based subdirectories based on the plan file

        Args:
            plate (int): Plate number
            ifudesign (int): ifudesign number
        """
        # Generate the calling path for this plate/ifudesign
        path = os.path.join(self.calling_path, str(plate), str(ifudesign))
        if not os.path.isdir(path):
            os.makedirs(path)

        # Generate the main common path
        path = default_dap_common_path(plate=plate, ifudesign=ifudesign, drpver=self.mpl.drpver,
                                       dapver=self.dapver, analysis_path=self.analysis_path)
        if not os.path.isdir(path):
            os.makedirs(path)

        # Read the analysis plan
        plan = AnalysisPlanSet.from_par_file(self.plan_file)
        for i in range(plan.nplans):
            # Generate the ref subdirectory for this plan
            path = default_dap_method_path(default_dap_method(plan=plan.data[i]), plate=plate,
                                           ifudesign=ifudesign, ref=True, drpver=self.mpl.drpver,
                                           dapver=self.dapver, analysis_path=self.analysis_path)
            if not os.path.isdir(path):
                os.makedirs(path)


    # Files:
    #       - mangadap-{plate}-{ifudesign}-LOG{mode} = script file = *
    #       - *.ready = script/directory is ready for queue submission
    #       - *.queued = script submitted to queue
    #       - *.started = script assigned a node/processor and has started running
    #       - *.done = completed execution of the full script
    # Other:
    #       - *.par = parameter file
    #       - *.out, *.err = stdout and stderr output from the script
    def _fill_queue(self, drpfiles, clobber=False):
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
            drpfiles (list): The list of DRP file objects that define
                the DRP fits files to analyze.
            clobber (bool): (**Optional**) Overwrite any existing script
                files.
        """
        # If submitting to the queue, create the queue object.  The
        # creation of the queue object is system dependent; however,
        # once created, the queue object is treated identically for both
        # the Utah and Portsmouth clusters.
        if self.create:
            if self.q is not None:
                # Expect to select the queue only when submitting to the
                # Portsmouth cluster, sciama.  In this case, the number
                # of nodes is queue dependent, and qos is not set
                if self.q == 'sciama1.q':
                    ppn = 12
                elif self.q == 'sciama3.q':
                    ppn = 20
                else:
                    ppn = 16
                self.queue.create(label=self.label, nodes=self.nodes, umask=self.umask,
                                  walltime=self.walltime, queue=self.q, ppn=ppn)
            else:
                # self.q can be None when submitting to both the
                # Portsmouth and Utah clusters.  In this case, the
                # default queue destination and ppn is correct.  qos is
                # also set, but this should only be used when submitting
                # to Utah.
                self.queue.create(label=self.label, nodes=self.nodes, qos=self.qos,
                                  umask=self.umask, walltime=self.walltime)
       
        # Create the script files, regardless of whether or not they are
        # submitted to the queue, appending the scripts to the queue if
        # they are to be submitted
        ndrp = len(drpfiles)
        i = 0
        for drpf in drpfiles:
            scriptfile, stdoutfile, stderrfile = self.prepare_for_analysis(drpf, clobber=clobber)
            if self.create:
                self.queue.append('source {0}'.format(scriptfile), outfile=stdoutfile,
                                  errfile=stderrfile)
                self.set_status(drpf.plate, drpf.ifudesign, status='queued') #, mode=drpf.mode
            else:
                print('Preparing for analysis...{0:.1f}%'.format((i+1)*100/ndrp), end='\r')
                i += 1
        
        # Submit the queue to the cluster
        if self.create:
            self.queue.commit(hard=self.hard, submit=self.submit)
        else:
            print('Preparing for analysis...  DONE')


    # ******************************************************************
    #  MaNGA DAP RUN-MODE ROUTINES
    # ******************************************************************
    def run_daily(self):
        """
        Check for and analyze DRP files that have been produced through
        the 3D stage.  Ideally, this is automatically run by crontab.

        Raises:
            Exception: *Always* raised because daily mode is currently
                unavailable for the DAP.
        """

        # Until correctly implemented, raise an exception if this mode
        # is called.
        raise NotImplementedError('Daily mode (-d,--daily) not yet implemented.')

        # Daily analyses forced to use the sdss-fast node and only that
        # node
        self.nodes = 1
        self.qos='sdss-fast'
        self.label = '{0}_daily'.format(self.label)

        # Select the daily DRP files and write the script files (do not
        # clobber)
        drpfiles = self.select_daily()
        print('Number of DRP files to process: {0}'.format(len(drpfiles)))
        if len(drpfiles) == 0:
            return

        # Create the script files for the set of drpfiles, and 
        # submit the scripts to the queue if self.submit is true
        self._fill_queue(drpfiles)


    def run_all(self, clobber=False):
        """
        Analyze all existing DRP files.  This mode is called manually.
   
        Args:
            clobber (bool): (**Optional**) If True, all analyses will be
                run, regardless of their current state.  If false, only
                those analyses not currently running or done will be
                performed.

        Raises:
            ValueError: Raised if qos is requested; this mode is only
                allowed for the daily run.
        """
        self.label = '{0}_clobber'.format(self.label) if clobber else '{0}_all'.format(self.label)

        # In here, qos should never be anything but None; always set in daily
        if self.qos is not None:
            raise ValueError('This qos is reserved for single-node usage.')

        drpfiles = self.select_all(clobber=clobber)
        print('Number of DRP files to process: {0}'.format(len(drpfiles)))
        if len(drpfiles) == 0:
            return

        # Create the script files for the set of drpfiles, and 
        # submit the scripts to the queue if self.submit is true
        self._fill_queue(drpfiles, clobber=clobber)


    def run_redo(self):
        """
        Redo a selected set of analyses using the input plate,
        ifudesign, and mode lists.  This mode **automatically** clobbers
        any existing results and is called manually.

        Raises:
            ValueError: Raised if qos is requested; this mode is only
                allowed for the daily run.
        """
        self.label = '{0}_redo'.format(self.label)
#        # In here, qos should never be anything but None; always set in daily
#        if self.qos is not None:
#            raise ValueError('This qos is reserved for single-node usage.')

        drpfiles = self.select_redo()
        print('Number of DRP files to process: {0}'.format(len(drpfiles)))
        if len(drpfiles) == 0:
            return       # If coded correctly, this function should never return here.

        # Create the script files for the set of drpfiles, and 
        # submit the scripts to the queue if self.submit is true
        self._fill_queue(drpfiles, clobber=True)


    # ******************************************************************
    # Management
    # ******************************************************************
    def set_status(self, plate, ifudesign, mode=None, status='queued'):
        """
        Generate a touch file that signifies the status for a given
        plate/ifudesign/mode.  The root of the touch file is set by
        :func:`mangadap.config.defaults.default_dap_file_root`.

        Args:
            plate (int): Plate number
            ifudesign (int): IFU design number
            mode (str): (**Optional**) DRP 3D mode, 'RSS' or 'CUBE'
            status (str): (**Optional**) That status signifier for the
                touch file.
        """
        # Get the name of the status file        
        root = default_dap_file_root(plate, ifudesign, mode=mode)
        statfile = os.path.join(self.calling_path, str(plate), str(ifudesign),
                                '{0}.{1}'.format(root,status))

        # Touch the file by opening and closing it
        file = open(statfile,'w')
        file.close()


    def dap_complete(self, drpf):
        """
        Determine if the DAP has finished working on a given DRP file.

        Conditions that the DAP has completes its processing of the
        given DRP file:
            - the output path for the DAP file exists
            - within the path, the *.done touch file exits

        .. todo::

            - Add whether or not the DAP ended in error

        Args:
            drpf (:class:`mangadap.drpfits.DRPFits`) : DRP file object
                used to pass the plate, ifudesign, and mode values.

        Returns:
            bool: Flag that the DAP has finished (True) or not (False).
        """
        path = default_dap_common_path(plate=drpf.plate, drpver=self.mpl.drpver,
                                       dapver=self.dapver, analysis_path=self.analysis_path)
        if not os.path.isdir(path):
            return False

        root = default_dap_file_root(drpf.plate, drpf.ifudesign, drpf.mode)
#        print(root)
        donefile = '{0}.done'.format(os.path.join(path,root))
#        print(donefile)
        if os.path.exists(donefile):
            return True

        return False


    def select_daily(self):
        """
        Select the drp files created today.
       
        .. todo::

            - This selection is correctly implemented.  It should not
              just look for DRP files created in the same day as the
              execution of the script.  It should look for anything
              created based on newly created DRP files.  This could be
              very close to the combination of all=True, clobber=False.

        Raises:
            NotImplementedError: **Always** raised because function is not
                correctly implemented.
        """

        # TODO: This selection is incorrectly implemented!
        # "created_today" is not the correct selection criterion
        raise NotImplementedError('select_daily() routine not correctly implemented yet.')

        drplist = self.full_drpfile_list()
        n_drp = len(drplist)
#       print(n_drp)
#       for i in range(0,n_drp):
#           print(drplist[i].file_name())
#           print(drplist[i].created_today())

        # Select the files created today
        return [ drplist[i] for i in range(0,n_drp) if drplist[i].created_today() ]


    def select_all(self, clobber=False):
        """
        Select all the DRP files.  If clobber=False, omit files that
        have completed a run of DAP analysis.

        Args:
            clobber (bool): (**Optional**) Overwrite any existing
                analysis output.

        Returns:
            list: A list of :class:`mangadap.drpfits.DRPFits` objects
            defining the DRP files that should be analyzed.
        """

        # Get the full list of available DRP files that can be processed
        # by the DAP
        drplist = self.full_drpfile_list()

        n_drp = len(drplist)
#       print(n_drp)
#       for i in range(0,n_drp):
#           print(drplist[i].file_path())
#           print(self.dap_complete(drplist[i]))

        # If clobbering, just return the full list
        if clobber:
            return drplist

        # If clobbering is off, only return those that are not complete
        return [ drplist[i] for i in range(0,n_drp) if not self.dap_complete(drplist[i]) ]


    def select_redo(self):
        """
        Return a list of the DRP selected using the provided plates,
        ifudesigns, and modes.  By using this function, clobber is
        implicitly True.

        Returns:
            list: A list of :class:`mangadap.drpfits.DRPFits` objects
            defining the DRP files that should be analyzed.
        """
        return self.selected_drpfile_list()


    def full_drpfile_list(self):
        """
        Generate a list of all the DRP files based on the
        :class:`mangadap.survey.drpcomplete.DRPComplete` database.  DRP
        files constructed based on the
        :class:`mangadap.survey.drpcomplete.DRPComplete` database are
        expected to exist.

        Returns:
            list: A list of :class:`mangadap.drpfits.DRPFits` objects
            defining the available DRP files.
        """
        # Number of plates (CUBE only)
        n_plates = self.drpc.nobs
#        for i in range(n_plates):
#            print('{0:>5} {1:>5} {2:>10} {3:>10} {4:>10} {5:9.1f}'.format(self.drpc['PLATE'][i],
#                                                            self.drpc['IFUDESIGN'][i],
#                                                            self.drpc['MANGAID'][i],
#                                                            self.drpc['MANGA_TARGET1'][i],
#                                                            self.drpc['MANGA_TARGET3'][i],
#                                                            self.drpc['VEL'][i]))

        # Create the list of CUBE DRP files
        drplist = [ DRPFits(self.drpc['PLATE'][i], self.drpc['IFUDESIGN'][i], 'CUBE',
                            drpver=self.mpl.drpver, redux_path=self.redux_path) \
                            for i in range(n_plates) \
                                if self.drpc['MANGAID'][i] != 'NULL' \
                                   and (self.drpc['MANGA_TARGET1'][i] > 0 \
                                        or self.drpc['MANGA_TARGET3'][i] > 0) \
                                   and self.drpc['VEL'][i] > 0.0 ]

        # Add the list of RSS DRP files
        drplist = drplist + [ DRPFits(self.drpc['PLATE'][i], self.drpc['IFUDESIGN'][i],
                                      'RSS', drpver=self.mpl.drpver, redux_path=self.redux_path)
                  for i in range(n_plates)
                  if self.drpc['MANGAID'][i] != 'NULL' \
                     and (self.drpc['MANGA_TARGET1'][i] > 0 \
                          or self.drpc['MANGA_TARGET3'][i] > 0) \
                     and self.drpc['VEL'][i] > 0.0 \
                     and self.drpc['MODES'][i] == 2 ]
        return drplist


    def selected_drpfile_list(self):
        """
        Generate a list of all the DRP files based on the provided
        plates, ifudesigns, and modes, using the
        :class:`mangadap.survey.drpcomplete.DRPComplete` database.

        Returns:
            list: A list of :class:`mangadap.drpfits.DRPFits` objects
            defining the selected DRP files.
        """
        # Number of plates (CUBE only); self.drpc.platelist and
        # self.drpc.ifudesignlist should be the same size!
        n_plates = len(self.drpc.platelist)

        # Get the list of CUBE DRP files, if requested (implicitly or
        # otherwise)
        getcube = True
        if self.modelist is not None:
            try:
                index = self.modelist.index('CUBE')
            except ValueError: # as e:
                getcube = False

        drplist = []
        if getcube:
            for i in range(0,n_plates):
                j = self.drpc.entry_index(self.drpc.platelist[i], self.drpc.ifudesignlist[i])
                if self.drpc['MANGAID'][j] != 'NULL' \
                        and (self.drpc['MANGA_TARGET1'][j] > 0 \
                             or self.drpc['MANGA_TARGET3'][j] > 0) \
                        and self.drpc['VEL'][j] > 0.0:
                    drplist += [ DRPFits(self.drpc.platelist[i], self.drpc.ifudesignlist[i], 'CUBE',
                                 drpver=self.mpl.drpver, redux_path=self.redux_path) ]

        # Add the list of RSS DRP files, if requested (implicitly or
        # otherwise)
        if self.modelist is not None:
            try:
                index = self.modelist.index('RSS')
            except ValueError: #as e:
                return drplist                  # List complete

        for i in range(0,n_plates):
            j = self.drpc.entry_index(self.drpc.platelist[i], self.drpc.ifudesignlist[i])
            if self.drpc['MANGAID'][j] != 'NULL' \
                    and (self.drpc['MANGA_TARGET1'][j] > 0 \
                         or self.drpc['MANGA_TARGET3'][j] > 0) \
                    and self.drpc['VEL'][j] > 0.0 and self.drpc['MODES'][j] == 2:
                drplist += [ DRPFits(self.drpc.platelist[i], self.drpc.ifudesignlist[i], 'RSS',
                             drpver=self.mpl.drpver, redux_path=self.redux_path) ]

        return drplist


    def write_compute_script(self, plate, ifudesign, mode, dapproc=True, plots=True,
                             clobber=False, relative_symlink=True):
        """
        Write the MaNGA DAP script file that is sent to a single CPU to
        analyze a single DRP file with a given plate, ifudesign, and
        mode.

        Args:
            plate (int): Plate number
            ifudesign (int): IFU design number
            mode (int): DRP 3D mode ('RSS' or 'CUBE')
            dapproc (bool): (**Optional**) Flag to execute the main DAP
                processing. Default is True.
            plots (bool): (**Optional**) Create the QA plots. Default is
                True.
            clobber (bool): (**Optional**) Flag to clobber any existing
                files.
            relative_symlink (bool): (**Optional**) Set the symlink to
                the par file to be a relative path instead of the
                absolute path.

        Returns:
            str: Three strings with the name of the written script file,
            the file for the output sent to STDOUT, and the file for the
            output sent to STDERR.

        Raises:
            Exception: Raised if DAP module version is not correctly
                defined.
        """
        # Check that the path exists, creating it if not
        self._check_paths(plate, ifudesign)

        # Create the parameter file
        parfile = default_dap_par_file(plate, ifudesign, mode, drpver=self.mpl.drpver,
                                       dapver=self.dapver, analysis_path=self.analysis_path)
        # Write the par file if it doesn't exist
        if not os.path.isfile(parfile) or clobber:
            # clobber defaults to True
            self.drpc.write_par(parfile, mode, plate=plate, ifudesign=ifudesign)
            # and create symlinks to it
            plan = AnalysisPlanSet.from_par_file(self.plan_file)
            nplan = plan.nplans
            for i in range(nplan):
                # Generate the ref subdirectory for this plan
                path = default_dap_method_path(default_dap_method(plan=plan.data[i]), plate=plate,
                                               ifudesign=ifudesign, ref=True,
                                               drpver=self.mpl.drpver, dapver=self.dapver,
                                               analysis_path=self.analysis_path)
                olink_dest = os.path.join(path,parfile.split('/')[-1])
                olink_src = os.path.relpath(parfile, start=os.path.dirname(olink_dest)) \
                                if relative_symlink else parfile
                if os.path.isfile(olink_dest) or os.path.islink(olink_dest):
                    os.remove(olink_dest)
                os.symlink(olink_src, olink_dest)

#        # Generate the DRP input and DAP output paths, and the DAP
#        # source path
#        drppath = default_drp_directory_path(plate, drpver=self.mpl.drpver,
#                                             redux_path=self.redux_path)
#        output_path = default_dap_common_path(plate=plate, ifudesign=ifudesign,
#                                              drpver=self.mpl.drpver, dapver=self.dapver,
#                                              analysis_path=self.analysis_path)
#        dap_source = default_dap_source()

        # Set the root path for the scripts, inputs, outputs, and logs
        _calling_path = os.path.join(self.calling_path, str(plate), str(ifudesign))

        scr_file_root = default_dap_file_root(plate, ifudesign)#, mode)

        # TODO: Why is this check here?!?
        # Get module name and fault if no module version is available
        module_version = self.dapver
        if module_version is None:
            raise ValueError('No DAP module version!')
        
        # Set the names for the script, stdout, and stderr files
        scriptfile = os.path.join(_calling_path, scr_file_root)
        stdoutfile = '{0}.out'.format(scriptfile)
        stderrfile = '{0}.err'.format(scriptfile)

        ################################################################
        # Script file already exists, so just return
        if os.path.exists(scriptfile) and not clobber:
            return scriptfile, stdoutfile, stderrfile

        # Main script components are:
        #   - (dapproc) Run the DAP processing code:
        #       - copy/create the plan file
        #       - execute the manga_dap via IDL
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
        if dapproc:
            command = 'manga_dap {0} {1} -r {2} -a {3}'.format(parfile, self.plan_file,
                                                               self.redux_path, self.analysis_path)
            if self.log:
                command += (' --log {0}.log'.format(scriptfile))
            if self.verbose > 0:
                command += (' -'+'v'*self.verbose )
            file.write('{0}\n'.format(command))

#            # Create/Copy the plan file
#            if self.plan_file is None:
#                # Will create and use the default plan
#                file.write('echo \" manga_dap, par=\'{0}\', drppath=\'{1}\', dappath=\'{2}\', ' \
#                           '/nolog\" | {3} \n'.format(parfile, drppath, self.analysis_path,
#                           self.idl_cmnd))
#            else:
#                # Will use the provided plan file, but first copy it for
#                # documentation purposes
#                default_plan_file = default_dap_plan_file(self.mpl.drpver, self.dapver,
#                                                          self.analysis_path, None, plate,
#                                                          ifudesign, mode)
#                file.write('\cp -rf {0} {1}\n'.format(self.plan_file, default_plan_file))
#                file.write('\n')
#                # Change the prior if requested
#                if self.prior_mode is not None:
#                    prior_file = default_dap_file_name(plate, ifudesign, self.prior_mode,
#                                                       self.prior_bin, self.prior_iter)
#                    prior_file = os.path.join(output_path, prior_file)
#                    scr = os.path.join(dap_source, 'bin', 'edit_dap_plan')
#                    file.write('{0} {1} analysis_prior {2} {3}\n'.format(scr, default_plan_file,
#                                                                        self.prior_old, prior_file))
#                    file.write('\n')
#
#                file.write('echo \" resolve_all, resolve_procedure=\'manga_dap\' & manga_dap, ' \
#                           'par=\'{0}\', plan=\'{1}\', drppath=\'{2}\', dappath=\'{3}\', /nolog \"'
#                           ' | {4} \n'.format(parfile, default_plan_file, drppath,
#                                              self.analysis_path, self.idl_cmnd))
#
            file.write('\n')

        # Plotting scripts
        if plots:
            file.write('# No plotting scripts yet\n')
#            pylist_path = os.path.join(dap_source, 'python', 'mangadap', 'plot',
#                                       'make_qa_file_list.py')
#            file.write('python3 {0} {1} {2} {3}_files_to_plot.txt -overwrite \n'.format(pylist_path,
#                       output_path, mode, mode))
#            file.write('\n')
#
#            pyplot_path = os.path.join(dap_source, 'python', 'mangadap', 'plot', 'plotqa.py')
#            file.write('python3 {0} {1}/{2}_files_to_plot.txt dapqa_plottypes.ini \n'.format(
#                                    pyplot_path, output_path, mode))
            file.write('\n')

        # Touch the done file
        donefile = '{0}.done'.format(scriptfile)
        file.write('touch {0}\n'.format(donefile))

        file.write('\n')

        file.close()
        ################################################################

        # Return the script file, file for stdout, and file for stderr
        return scriptfile, stdoutfile, stderrfile
    

    def prepare_for_analysis(self, drpf, clobber=False):
        """
        Prepare the provided DRP file for analysis by setting up the
        output path, writing the input parameter file, writing the
        compute script, and setting the status to "ready", assuming the
        above steps are completed successfully.

        Args:
            drpf (:class:`mangadap.drpfits.DRPFits`): Object defining
                the DRP file to analyze.
            clobber (bool): (**Optional**) Flag to overwrite any
                existing files.
        
        Returns:
            str: Three strings with the name of the written script file,
            the file for the output sent to STDOUT, and the file for the
            output sent to STDERR, as provided by the call to
            :func:`write_compute_script`.
        """
        # TODO: Check that *.ready file exists?
        
        # Write the compute script (write_compute_script also checks the
        # path exists!)
#        sf, of, ef = self.write_compute_script(drpf.plate, drpf.ifudesign, drpf.mode, dapproc=self.dapproc, plots=self.plots,
#                                               clobber=clobber)
        try:
            sf, of, ef = self.write_compute_script(drpf.plate, drpf.ifudesign, drpf.mode,
                                                   dapproc=self.dapproc, plots=self.plots,
                                                   clobber=clobber)
        except:
            e = sys.exc_info()
            print_frame(e[0])
            warnings.warn('Problem writing compute script:: {0}: {1}.  Continuing...'.format(
                                                                                        e[1], e[2]))
            return None, None, None

        # Set the status to ready
        self.set_status(drpf.plate, drpf.ifudesign, status='ready') #, mode=drpf.mode

        # Return the list of script, stdout, and stderr files
        return sf, of, ef



