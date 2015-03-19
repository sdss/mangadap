from __future__ import division
from __future__ import print_function
from __future__ import absolute_import

import sys
if sys.version > '3':
    long = int

import subprocess
import time
import os.path
import pbs.queue
import numpy

from os import environ, makedirs
from argparse import ArgumentParser

# DAP imports
from mangadap.drpcomplete import drpcomplete
from mangadap.drpfile import drpfile
from mangadap.util.exception_tools import print_frame
from mangadap.util.parser import arginp_to_list
from mangadap.mangampl import mangampl
from mangadap.survey import util


__author__ = 'Kyle Westfall'

class rundap:
    """
    Automated procedures for MaNGA DAP processing at Utah.

    [ More complete explanation of the class here. ]

    IMPORTANT: This is purely a survey-level utility for running the DAP
    at Utah, including submitting jobs to the cluster.  No user-level
    usage of the DAP should be considered here, apart from not hindering
    such usage of the primary IDL/python programs!
        
    METHODS:
        run_daily       - process newly completed DRP reductions as
                          determined by the daily assessment

        run_all         - process all DRP reductions not currently
                          running or done unless clobber is specified to
                          force processing of all files

        run_redo        - redo the processing of selected DRP reductions


    UTILITY FUNCTIONS:
                
        _read_arg         - set command line arguments

    REVISION HISTORY:
        NOV-2014
            Architecture setup (Joel Brownstein
            <joelbrownstein@astro.utah.edu>)

        18 Nov 2014 (KBW):

            (Attempt to) conform style to PEP 8:

                https://www.python.org/dev/peps/pep-0008

            except that lines are 99 characters long (docstrings still
            limited to 72 characters).

            Also (attempt to) conform to docstrings style according to
            PEP 257:

                https://www.python.org/dev/peps/pep-0257

        01 Dec 2014: (KBW) Committed to SVN
        02 Dec 2014: (KBW) Changed drpver and idlutilsver to mangaver
        13 Jan 2015: (KBW) Changed it back, added platetargets
        05 Feb 2015: (KBW) Change to load module based on MPL
        19 Mar 2015: (KBW) Major changed to allow for more control over
                           the directory structure containing the DRP
                           data and for the DAP results; following
                           changed in the drpfile and drpcomplete
                           classes.
    """

    # TODO: outver is actually just the directory for output to use in
    # place of the current product version 
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
                # Databases with input parameter information
                platetargets=None, nsa_cat=None, nsa_catid=None,
                # Cluster options
                label='mangadap', nodes=18, qos=None, umask='0027',walltime='240:00:00', hard=True,
                submit=True):
        """
        Interprets the run-mode and cluster options.
    
        ARGUMENTS:
                daily bool 
                        run daily analyses
                all bool 
                        analyze all DRP data
                clobber bool 
                        when running all analyses, clobber existing
                        output
                redo bool 
                        redo existing analyses, clobber is implicitly
                        true

                console bool
                        executed from console, with command-line
                        arguments
                quiet bool
                        suppress output
                version bool
                        Only print version and return

                mplver string
                        MPL version to analyze.  Used with mangampl
                        class to set the idlutils version and the DRP
                        version.  THIS CONTROLS THE MODULES THAT ARE
                        LOADED WRT THE DRP AND ITS DEPENDENCIES.  The
                        default version (selected if mplver=None) is
                        defined by the mangampl class.
                redux_path string
                        Main reduction path with the DRP files.  See
                        drpfile class.
                dapver string
                        DAP version to use for analysis.  THIS CONTROLS
                        THE MODULES THAT ARE LOADED WRT THE DAP AND ITS
                        DEPENDENCIES (TBD if this can be fully
                        independent of the DRP/idlutils).
                analysis_path string
                        Main path for the DAP output, used to set
                        dappath in the call to manga_dap.  See dapfile
                        class.

                plan_file string
                        Name of the plan file to use for ALL DRP data in
                        this run
                platelist string
                        specified list of plates to analyze
                ifudesignlist string
                        specified list of ifudesigns to analyze
                modelist string
                        specified list of modes to analyze (CUBE or
                        RSS)
                combinatorics bool
                        use all unique combinations of the entered
                        plate/ifudesign/mode lists

                platetargets string
                        list of platetargets files to use for creating
                        the drpcomplete file. See drpcomplete.py.
                nsa_cat string
                        list of NSA catalogs to use for creating the
                        drpcomplete file.  See drpcomplete.py
                nsa_catid string
                        list of NSA catalog IDs.  This is expected to be
                        the first number in the MANGA ID, used to select
                        the data in the NSA catalog as specified by the
                        catalog index (the second number in the MANGA
                        ID)

                label string
                        label to use in cluster queue. Default is
                        mangadap and the run mode (daily, etc)
                nodes int
                        number of cluster nodes to use
                qos string
                        select processor (only None for daily run)
                umask string
                        umask to set for output
                walltime string
                        wall time for cluster job
                hard bool
                        turn OFF hard keyword for cluster submission
                submit bool
                        turn OFF submission of jobs to cluster
        """
        # Only print the version of the DAP
        if version:
            print('This is version 0.96.')
            return

        # Save run-mode options
        self.daily = daily
        self.all = all
        self.clobber=clobber
        self.redo=redo
        self.quiet = quiet

        # Override environment
        self.mpl = mplver
        self.redux_path = redux_path

        self.dapver = util.product_version(simple=True, product='mangadap') if dapver is None \
                                                                            else dapver
        self.analysis_path = dapfile._default_analysis_path() if analysis_path is None \
                                                        

        # List of files to analyze
        self.plan_file = plan_file
        self.platelist = arginp_to_list(platelist, evaluate=True)
        self.ifudesignlist = arginp_to_list(ifudesignlist, evaluate=True)
        self.modelist = arginp_to_list(modelist)
        self.combinatorics = combinatorics

        # plateTargets file(s) or NSA catalog to use for drpcomplete file
        # (drpcomplete can handle [platetargets = None, nsa_cat = None])
        self.platetargets = platetargets
        self.nsa_cat = nsa_cat
        self.nsa_catid = nsa_catid

        # Cluster queue keywords
        self.label = label
        self.nodes = nodes
        self.qos = qos
        self.umask = umask
        self.walltime = walltime
        self.hard = hard
        self.submit = submit

        # Read the analysis path, if it does not exist, exception is
        # raised.
        try:
            self.manga_spectro_analysis = environ['MANGA_SPECTRO_ANALYSIS']
        except:
            print_frame('Exception')
            raise Exception('Environmental variable MANGA_SPECTRO_ANALYSIS is undefined!')
       
        # Read and parse command-line arguments
        if console:
            self._read_arg()

#       # Make sure the selected MPL version is available
#       self.drpver = None
#       self._select_mpl()
        try:
            self.mpl = mangampl(self.mpl)
        except Exception as e:
            print_frame('Exception')
            raise Exception('Undefined MPL:'+e)

        # Alert the user of the versions to be used
#       print('Versions: DAP:{0}, IDLUTILS:{1}, DRP:{2}'.format(self.outver, self.idlutilsver, self.drpver))
#       print('Versions: DAP:{0}, {1}'.format(self.outver, self.mplver))
        print('Versions: DAP:{0}, {1}'.format(self.outver, self.mpl.mplver))

        # Check that something is to be done
        nrun = 0
        if self.daily: nrun += 1
        if self.all: nrun += 1
        if self.redo: nrun += 1
        if nrun == 0 or nrun > 1:
            raise Exception('Must request one (and only one) run mode!')

        # Check argument combination makes sense
        # TODO: Does an error need to be thrown here, or just warn user.
        # Use warnings class?
        if self.clobber and not self.all:
            print('Clobber keyword can only be set if all specified (-a,--all)!')
            self.clobber = None

        # If redo, platelist and ifudesignlist MUST be provided
        if self.redo and not self.platelist and not self.ifudesignlist and not self.modelist:
            raise Exception('Must provide platelist, ifudesignlist, and modelist for redo!')

        # If not redo, platelist, ifudesignlist, and modelist should be
        # ignored
        # TODO: Does an error need to be thrown here, or just warn user.
        # Use warnings class?
        if not self.redo and (self.platelist or self.ifudesignlist or self.modelist):
            raise Exception('platelist/ifudesignlist/modelist only allowed for redo!')

        # If a plan file is provided, make sure it exists
        if self.plan_file is not None and not os.path.exists(self.plan_file):
            raise Exception('Provided plan file does not exist.')

        # If running all or daily, make sure lists are set to None such
        # that drpcomplete will search for all available DRP files and
        # make sure that the needed information is up-to-date
        if self.all or self.daily:
            self.platelist = None
            self.ifudesignlist = None
            self.modelist = None

        # Create and update the drpcomplete file if necessary
#       self.drpc = drpcomplete(self.outver, self.drpver, platetargets=self.platetargets,
        self.drpc = drpcomplete(platetargets=self.platetargets, nsa_cat=self.nsa_cat,
                                dapver=self.outver, drpver=self.mpl.drpver)
        # Update the drpcomplete list; force an update if the
        # platetarget files or an NSA catalog file are provided
        if self.platetargets is not None or self.nsa_cat is not None:
            self.drpc.update(platelist=self.platelist, ifudesignlist=self.ifudesignlist,
                             combinatorics=self.combinatorics, force=True)
        else:
            self.drpc.update(platelist=self.platelist, ifudesignlist=self.ifudesignlist,
                             combinatorics=self.combinatorics)

        # Hereafter, self.platelist and self.ifuplatelist are the INPUT
        # values, whereas self.drpc.platelist and self.drpc.ifuplatelist
        # are the actual DRP reductions to analyze (including
        # combinatorics and file existence tests).  These should be
        # combined with self.modelist to get the DRP files to analyze.
        # See, e.g., selected_drpfile_list().

        # Prep the verbosity of the queue if submitting
        if self.submit:
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


    def _check_path(self, plate, ifudesign):
        """
        Check if the output path exists, creating it if it doesn't.
        """
        path = self.output_path(plate, ifudesign)
        if not os.path.isdir(path):
            makedirs(path)


    def _read_arg(self):
        """
        Interpret the command-line arguments to use during execution.
        """

        parser = ArgumentParser()
        mode = parser.add_mutually_exclusive_group(required=True)

        # Mode aruments
        mode.add_argument("-d", "--daily",
                          help="runs dap for next mjd (as an after burner to drp)",
                          action="store_true")
        mode.add_argument("-a", "--all",
                          help="runs dap for all plates/ifudesigns/modes not currently running "
                          "or done", action="store_true")
        mode.add_argument("-r", "--redo",
                          help="runs dap for specified platelist/ifudesignlist/modelist "
                          "regardless of state", action="store_true")

        # Run-mode optional arguments
        parser.add_argument("-c", "--clobber",
                            help="if all selected, will run dap for all plates/ifudesigns/modes "
                            " regardless of state", action="store_true", default=False)
        parser.add_argument("-q", "--quiet", help="suppress screen output", action="store_true",
                            default=False)
        parser.add_argument("-o", "--outver", type=str,
                            help="optional output version, different from product version",
                            default="trunk")

#       parser.add_argument("-i", "--idlutilsver", type=str, help="version of idlutils to use",
#                           default=None)
#       parser.add_argument("-v", "--drpver", type=str, help="version of mangadrp used to produce "
#                           "files to process", default=None)

        parser.add_argument("-v", "--mplver", type=str, help="MPL version to analyze", default=None)

        parser.add_argument("-p", "--platelist", type=str, help="set list of plates to reduce",
                            default=None)
        parser.add_argument("-b", "--ifudesignlist", type=str, help="set list of ifus to reduce",
                            default=None)
        parser.add_argument("-m", "--modelist", type=str, help="set list of DRP output modes to "
                            "reduce (CUBE or RSS)", default=None)
        parser.add_argument("-x", "--combinatorics", help="force execution of all permutations of "
                            "the provided lists", action="store_true", default=False)

        parser.add_argument("-g", "--plttargets", type=str, help="path to plateTargets file(s); "
                            "if provided will force update to drpcomplete file", default=None)
        parser.add_argument("-k", "--nsacat", type=str, help="path to NSA catalog to use; if "
                            "provided will force update to drpcomplete file", default=None)

        parser.add_argument("-e", "--plan_file", type=str, help="parameter file with the MaNGA "
                            "DAP execution plan to use instead of the default" , default=None)

        # Cluster arguments
        parser.add_argument('-l', '--label', type=str, help='label for cluster job',
                            default='mangadap')
        parser.add_argument('-n', '--nodes', type=int, help='number of nodes to use in cluster',
                            default=18)
        parser.add_argument('-f', '--fast', dest='qos', type=str, help='qos state',
                            default=None)
        parser.add_argument('-u', '--umask', type=str, help='umask bit for cluster job',
                            default='0027')
        parser.add_argument('-w', '--walltime', type=str, help='walltime for cluster job',
                            default='240:00:00')
        parser.add_argument('-t', '--toughness', dest='hard', action='store_false', default=True,
                            help='turn off hard keyword for cluster submission')
        parser.add_argument('-s', '--submit', help='turn off cluster submission',
                            action='store_false', default=True)
        
        # Parse the arguments and assign them to self
        self.arg = parser.parse_args()

        ################################################################
        # Assign properties based on arguments, if properties not set
        # previously; the suitability of these arguments is checked in
        # __init__()

        # Run-mode
        if self.daily is None:
            self.daily = self.arg.daily
        if self.all is None:
            self.all = self.arg.all
        if self.redo is None:
            self.redo = self.arg.redo

        # Run-mode options
        if self.clobber is None:
            self.clobber = self.arg.clobber
        if self.quiet is None:
            self.quiet = self.arg.quiet

        # Copy the platelist, ifudesignlist, and modelist, if provided.
        # Will OVERWRITE existing input from __init__()
        if self.arg.platelist is not None:
            self.platelist = arginp_to_list(self.arg.platelist, evaluate=True)
        if self.arg.ifudesignlist is not None:
            self.ifudesignlist = arginp_to_list(self.arg.ifudesignlist, evaluate=True)
        if self.arg.modelist is not None:
            self.modelist = arginp_to_list(self.arg.modelist)
        self.combinatorics = self.arg.combinatorics
   
        # Set the plateTargets and NSA catalog path
        if self.arg.plttargets is not None:
            self.platetargets = self.arg.plttargets
        if self.arg.nsacat is not None:
            self.nsa_cat = self.arg.nsacat

        # Save the plan file to use
        if self.arg.plan_file is not None:
            self.plan_file = self.arg.plan_file

        # Set the versions to use; versions will never be none because
        # of procedures in __init__().
        # Will OVERWRITE existing input from __init__()
        if self.arg.outver is not None:
            self.outver = self.arg.outver
#       if self.arg.idlutilsver is not None:
#           self.idlutilsver = self.arg.idlutilsver
#       if self.arg.drpver is not None:
#           self.drpver = self.arg.drpver
        if self.arg.mplver is not None:
            self.mpl = self.arg.mplver

        # Set queue keywords
        if self.arg.umask is not None:
            self.umask = self.arg.umask
        if self.arg.nodes is not None:
            self.nodes = self.arg.nodes
        if self.arg.walltime is not None:
            self.walltime = self.arg.walltime
        if self.arg.qos is not None:
            self.qos = self.arg.qos
        if self.arg.umask is not None:
            self.umask = self.arg.umask
        if self.arg.hard is not None:
            self.hard = self.arg.hard
        if self.arg.submit is not None:
            self.submit = self.arg.submit

#       print(self.submit)
#       self.submit = False


    # TODO: Files:
    #       - mangadap-{plate}-{ifudesign}-LOG{mode} = script file = *
    #       - *.ready = script/directory is ready for queue submission
    #       - *.queued = script submitted to queue
    #       - *.started = script assigned a node/processor and has started running
    #       - *.done = completed execution of the full script
    # TODO: Other:
    #       - *.par = parameter file
    #       - *.out, *.err = stdout and stderr output from the script
    def _fill_queue(self, drpfiles, clobber=False):

        # If submitting to the queue, create the queue object
        if self.submit:
            self.queue.create(label=self.label, nodes=self.nodes, qos=self.qos, umask=self.umask,
                              walltime=self.walltime)
       
        # Create the script files, regardless of whether or not they are
        # submitted to the queue, appending the scripts to the queue if
        # they are to be submitted
        for drpf in drpfiles:
            scriptfile, stdoutfile, stderrfile = self.prepare_for_analysis(drpf, clobber=clobber)
            if self.submit:
                self.queue.append('source {0}'.format(scriptfile), outfile=stdoutfile,
                                  errfile=stderrfile)
                self.set_status(drpf.plate, drpf.ifudesign, drpf.mode, stage='dap', status='queued')
        
        # Submit to queue to the cluster
        # hard = hard-coded script files that can be seen
        if self.submit:
            self.queue.commit(hard=self.hard, submit=self.submit)


    def _write_module_commands(self, file):
        """
        Write the module commands to the script file.  The file must
        already be open.

        """

        module = self.mpl.module_file()
#        file.write('module unload manga\n')
#        file.write('module load {0}\n'.format(module))

        # TODO: Is there a way to have the script catch errors in these
        # module calls (without proceeding to the remaining calls)?

        # TODO: Can I check that all the modules have been properly
        # loaded and then quit otherwise? (like the module_version
        # script above?)


    # ******************************************************************
    #  MaNGA DAP RUN-MODE ROUTINES
    # ******************************************************************


    def run_daily(self):
        """
        Run daily reductions for next mjd completed by the DRP.
        
        Automatically run by crontab.
        """

        # Until correctly implemented, raise an exception if this mode
        # is called.
        raise Exception('Daily mode (-d,--daily) not yet implemented.')

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
        Run all reductions.
    
        clobber: If True, all analyses will be run, regardless of their
        current state.  If false, only those analyses not currently
        running or done will be performed.

        This mode is manually called by humans.
        """

        self.label = '{0}_clobber'.format(self.label) if clobber else '{0}_all'.format(self.label)

        # In here, qos should never be anything but None; always set in daily
        if self.qos is not None:
            raise Exception('This qos is reserved for single-node usage.')

        drpfiles = self.select_all(clobber=clobber)
        print('Number of DRP files to process: {0}'.format(len(drpfiles)))
        if len(drpfiles) == 0:
            return

        # Create the script files for the set of drpfiles, and 
        # submit the scripts to the queue if self.submit is true
        self._fill_queue(drpfiles, clobber=clobber)

    def run_redo(self):
        """
        Run all reductions not currently running or done.
        
        Manually called by humans.
        """

        self.label = '{0}_redo'.format(self.label)
        # In here, qos should never be anything but None; always set in daily
        if self.qos is not None:
            raise Exception('This qos is reserved for single-node usage.')

        drpfiles = self.select_redo()
        print('Number of DRP files to process: {0}'.format(len(drpfiles)))
        if len(drpfiles) == 0:
            return       # If coded correctly, this function should never return here.

        # Create the script files for the set of drpfiles, and 
        # submit the scripts to the queue if self.submit is true
        self._fill_queue(drpfiles, clobber=True)



    # ******************************************************************
    #  Reduction Management
    # ******************************************************************

    
    def file_root(self, plate, ifudesign, mode, stage='dap'):
        """
        Generate the root name of the MaNGA DAP file for a given
        stage/plate/ifudesign.
        """

        return 'manga{0}-{1}-{2}-LOG{3}'.format(stage, plate, ifudesign, mode)


    def output_path(self, plate, ifudesign):
        """
        Set the path name for the ouptut of the MaNGA DAP for a given
        plate/ifudesign.
        """

        # TODO: Separate paths for RSS and CUBE?
        return os.path.join(self.manga_spectro_analysis, self.outver, str(plate), str(ifudesign))


    def set_status(self, plate, ifudesign, mode, stage='dap', status='queued'):
        """
        Generate a touch file that signifies the status for a given
        plate/ifudesign/mode.
        """
                
        # Touch status file
        root = self.file_root(plate, ifudesign, mode, stage)
        statfile = os.path.join(self.output_path(plate, ifudesign),'{0}.{1}'.format(root,status))
        file = open(statfile,'w')
        file.close()


    def dap_complete(self, drpf, stage='dap'):
        """
        Determine if the DAP was successfully completed on a given DRP
        file.

        Conditions that DAP is complete its processing of the given DRP file:
                - the output path for the DAP file exists
                - within the path, the *.done touch file exits
        """

        path = self.output_path(drpf.plate, drpf.ifudesign)
#        print(path)
        if not os.path.isdir(path):
            return False

        root = self.file_root(drpf.plate, drpf.ifudesign, drpf.mode, stage)
#        print(root)
        donefile = '{0}.done'.format(os.path.join(path,root))
#        print(donefile)
        if os.path.exists(donefile):
            return True

        return False


    def select_daily(self):
        """Select the drp files created today."""

        # TODO: This selection is incorrectly implemented!
        # "created_today" is not the correct selection criterion

        # raise an Exception until this method is correctly implemented
        raise Exception('select_daily() routine not correctly implemented yet.')

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
        Select all the drp files.  If clobber=False, omit files that
        have been successfully analyzed.
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

        # Get the selected list of DRP files based on the provided
        # plates/ifudesigns/and modes.  By selecting redo, clobber is
        # implicitly selected!
        return self.selected_drpfile_list()

#       if clobber:
#           return drplist

#       n_drp = len(drplist)
#       return [ drplist[i] for i in range(0,n_drp) if not self.dap_complete(drplist[i]) ]


    def full_drpfile_list(self):
        """
        Generate the list of drpfiles based on the drpcomplete database.
        
        DRP files constructed based on the drpcomplete database are
        expected to exist.
        """

#       # Number of plates (CUBE only)
#       n_plates = len(self.drpc.data['PLATE'])

#       # Create the list of CUBE DRP files
#       drplist = [ drpfile(self.drpc.data['PLATE'][i], self.drpc.data['IFUDESIGN'][i], 'CUBE',
#                           drpver=self.mpl.drpver) for i in range(0,n_plates) ]

        # Add the list of RSS DRP files
#       drplist = drplist + [ drpfile(self.drpc.data['PLATE'][i], self.drpc.data['IFUDESIGN'][i],
#                                     'RSS', drpver=self.mpl.drpver)
#                             for i in range(0,n_plates) if self.drpc.data['MODES'][i] == 2 ]

        # Edited to ignore drpfiles that are found but do not have an
        # entry in the platetargets files (mangaid = 'NULL')
        # Number of plates (CUBE only)
        n_plates = len(self.drpc.data['PLATE'])

        # Create the list of CUBE DRP files
        drplist = [ drpfile(self.drpc.data['PLATE'][i], self.drpc.data['IFUDESIGN'][i], 'CUBE',
                            drpver=self.mpl.drpver) for i in range(0,n_plates) 
                                                    if self.drpc.data['MANGAID'][i] != 'NULL' ]

        # Add the list of RSS DRP files
        drplist = drplist + [ drpfile(self.drpc.data['PLATE'][i], self.drpc.data['IFUDESIGN'][i],
                                      'RSS', drpver=self.mpl.drpver)
                  for i in range(0,n_plates)
                  if self.drpc.data['MANGAID'][i] != 'NULL' and self.drpc.data['MODES'][i] == 2 ]
        return drplist


    def selected_drpfile_list(self):
        """
        Generate the list of drpfiles, based on the provided plates,
        ifudesigns, and modes, using the drpcomplete database.
        
        See __init__ for the creation of the drpcomplete object
        (self.drpc).
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

        if getcube:
#           drplist = [ drpfile(self.drpc.platelist[i], self.drpc.ifudesignlist[i], 'CUBE', 
#                               drpver=self.mpl.drpver) for i in range(0,n_plates) ]
            drplist = [ drpfile(self.drpc.platelist[i], self.drpc.ifudesignlist[i], 'CUBE', 
                                drpver=self.mpl.drpver)
                      for i in range(0,n_plates)
                          if self.drpc.data['MANGAID'][self.drpc.entry_index(self.drpc.platelist[i],
                                                       self.drpc.ifudesignlist[i])] != 'NULL' ]
        else:
            drplist = list()

        # Add the list of RSS DRP files, if requested (implicitly or
        # otherwise)
        if self.modelist is not None:
            try:
                index = self.modelist.index('RSS')
            except ValueError: #as e:
                return drplist                  # List complete

#       drplist = drplist + [ drpfile(self.drpc.platelist[i], self.drpc.ifudesignlist[i], 'RSS',
#                                     drpver=self.mpl.drpver)
#                 for i in range(0,n_plates)
#                 if self.drpc.data['MODES'][self.drpc.entry_index(self.drpc.platelist[i],
#                                            self.drpc.ifudesignlist[i])] == 2 ]
        drplist = drplist + [ drpfile(self.drpc.platelist[i], self.drpc.ifudesignlist[i], 'RSS',
                                      drpver=self.mpl.drpver)
                  for i in range(0,n_plates)
                  if self.drpc.data['MANGAID'][self.drpc.entry_index(self.drpc.platelist[i],
                                               self.drpc.ifudesignlist[i])] != 'NULL' and
                     self.drpc.data['MODES'][self.drpc.entry_index(self.drpc.platelist[i],
                                             self.drpc.ifudesignlist[i])] == 2 ]

        return drplist


    def write_compute_script(self, plate, ifudesign, mode, stage='dap', clobber=False):
        """
        Write the MaNGA DAP script file for a given plate, ifudesign,
        and mode.
        """

        # Check the path exists, creating it if not
        self._check_path(plate, ifudesign)

        # TODO: Check that *.ready file exists?

        # Generate the path name and root name of the output files
        path = self.output_path(plate, ifudesign)
        root = self.file_root(plate, ifudesign, mode, stage)
            
        # Get module name
        # TODO: Set outver?
#        module_version = self.module_version()
        module_version = self.outver

        # Fault if no module version is available
        if module_version is None:
            raise Exception('No DAP module version!')
        
        # Set the names for the script, stdout, and stderr files
        scriptfile = os.path.join(path,root)
        stdoutfile = '{0}.out'.format(scriptfile)
        stderrfile = '{0}.err'.format(scriptfile)

        ################################################################
        # Script file already exists, so just return
        if os.path.exists(scriptfile) and not clobber:
            return scriptfile, stdoutfile, stderrfile


        ################################################################
        # Open the script file and write the date as a commented header
        # line
        file = open(scriptfile, 'w')
        file.write('# Auto-generated batch file {0}\n'.format(
                            time.strftime("%a, %d %b %Y %H:%M:%S +0000",time.localtime())))
        file.write('\n')

        # Append the module commands to the script
        self._write_module_commands(file)
        file.write('\n')

        # Create the started touch file
        startfile = '{0}.started'.format(scriptfile)
        file.write('touch {0}\n'.format(startfile))
        file.write('\n')

        # Set the output path manually
        dappath = os.path.join(self.manga_spectro_analysis, self.outver)

        # Command that runs the DAP
        parfile = self.parameter_file(plate, ifudesign, mode, stage)
        if self.plan_file is None:
            # Will create and use the default plan
            file.write('echo \" manga_dap, par=\'{0}\', dappath=\'{1}\', /nolog \" | idl \n'.format(parfile, dappath))
        else:
            # Will use the provided plan file, but first copy it for
            # documentation purposes
            default_plan_file = os.path.join(self.output_path(plate, ifudesign),
                                             'manga-{0}-{1}-LOG{2}-dapplan.par'.format(plate,
                                             ifudesign, mode))
            file.write('\cp -rf {0} {1}\n'.format(self.plan_file, default_plan_file))
            file.write('echo \" manga_dap, par=\'{0}\', plan=\'{1}\', dappath=\'{2}\', /nolog \"' \
                       ' | idl \n'.format(parfile, default_plan_file, dappath))
        file.write('\n')

        #file.write('setStatusDone -f "{0}" \n'.format(errfile))
        donefile = '{0}.done'.format(scriptfile)
        file.write('touch {0}\n'.format(donefile))
        file.write('\n')

        file.close()
        ################################################################

        # Return the script file, file for stdout, and file for stderr
        return scriptfile, stdoutfile, stderrfile
    

    def parameter_file(self, plate, ifudesign, mode, stage='dap'):
        """Get the name of the parameter file."""
        path = self.output_path(plate, ifudesign)
        root = self.file_root(plate, ifudesign, mode, stage)
        parfile = '{0}.par'.format(os.path.join(path,root))
        return parfile


    def prepare_for_analysis(self, drpfile, clobber=False):
        """
        Create the parameter file used by the DAP to process the
        provided DRP file, and create the script submitted to the queue
        to process the file.

            - drpfile is the file to analyze
            - clobber is used to overwrite any existing file(s)
        """
        
        # Check that the path exists, creating it if not
        self._check_path(drpfile.plate, drpfile.ifudesign)

        # Create the parameter file
        parfile = self.parameter_file(drpfile.plate, drpfile.ifudesign, drpfile.mode)
#       print(parfile)
        self.drpc.write_par(parfile, drpfile.mode, plate=drpfile.plate,
                            ifudesign=drpfile.ifudesign, clobber=clobber)

        # Write the compute script (write_compute_script also checks the
        # path exists!)
        try:
            sf, of, ef = self.write_compute_script(drpfile.plate, drpfile.ifudesign,
                                                   drpfile.mode, clobber=clobber)
        except Exception as e:
            print_frame('Exception')
            print("Exception: %s" % str(e))
            print("Skipping to next DRP file.")
            return None, None, None

        # Set the status to ready
        self.set_status(drpfile.plate, drpfile.ifudesign, drpfile.mode, stage='dap', status='ready')

        # Return the list of script, stdout, and stderr files
        return sf, of, ef

            


        

