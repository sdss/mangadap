from __future__ import division
from __future__ import print_function

import time
import os.path
#import pbs.queue

from os import environ, makedirs
from argparse import ArgumentParser
from drpcomplete import drpcomplete
from drpfile import drpfile, arginp_to_list

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
    """

    # TODO: outver is actually just the directory for output to use in
    # place of the current product version 
    def __init__(self,
                # Run mode options
                daily=None, all=None, clobber=None, redo=None, console=None, quiet=None,
                # Override defaults
                outver=None, mangaver=None, platelist=None, ifudesignlist=None, modelist=None,
                combinatorics=False, nsa_cat=None,
                # Cluster options
                label='mangadap', nodes=18, qos=None, umask='007',walltime='240:00:00', hard=True,
                submit=True):
        """
        Interprets the run-mode and cluster options.
    
        ARGUMENTS:
            daily - run daily analyses
            all - run all analyses
            clobber - when running all analyses, clobber existing output
            redo - redo existing analyses
            console - executed from console, with command-line arguments
            quiet - suppress output

            outver - output version
            mangaver - manga module version to use
            platelist - specified list of plates to analyze
            ifudesignlist - specified list of ifudesign to analyze
            modelist - specified list of modes to analysze (CUBE or RSS)
            combinatorics - use all unique combinations of the entered
                            plate/ifudesign/mode lists

            nsacat - specify the NSA catalog to use

            label - label to use in cluster queue
            nodes - number of cluster nodes to use
            qos - processor to use (?)
            umask - umask to set for output
            walltime - wall time for cluster job
            hard - turn OFF hard keyword for cluster submission
            submit - turn OFF submission of jobs to cluster
        """

        # Save run-mode options
        self.daily = daily
        self.all = all
        self.clobber=clobber
        self.redo=redo
        self.quiet = quiet

        # Determine current versions
#        self.idlutilsver = self.product_version(simple=True, product='idlutils')
#        self.corever = self.product_version(simple=True, product='mangacore')
#        self.drpver = self.product_version(simple=True, product='mangadrp')
        # TODO: Place holder!!!
        self.mangaver = self.product_version(simple=True, product='mangacore')
        self.dapver = self.product_version(simple=True, product='mangadap')

        # Use them or use input versions
        self.outver = self.dapver if outver is None else outver
        self.mangaver = self.mangaver if mangaver is None else mangaver

        # List of files to analyze
        self.platelist = arginp_to_list(platelist, evaluate=True)
        self.ifudesignlist = arginp_to_list(ifudesignlist, evaluate=True)
        self.modelist = arginp_to_list(modelist)
        self.combinatorics = combinatorics

        # NSA catalog to use for drpcomplete file (drpcomplete can
        # handle nsa_cat = None)
        self.nsa_cat = nsa_cat

        # Cluster queue keywords
        self.label = label
        self.nodes = nodes
        self.qos = qos
        self.umask = umask
        self.walltime = walltime
        self.hard = hard
        self.submit = submit

        # Read the analysis path
        try:
            self.manga_spectro_analysis = environ['MANGA_SPECTRO_ANALYSIS']
        except:
            raise Exception('Environmental variable MANGA_SPECTRO_ANALYSIS is undefined!')
#            print('MANGA_SPECTRO_ANALYSIS undefined!')
#            exit()
       
        # Read and parse command-line arguments
        if console:
            self._read_arg()

        # Check that something is to be done
        nrun = 0
        if self.daily: nrun += 1
        if self.all: nrun += 1
        if self.redo: nrun += 1
        if nrun == 0 or nrun > 1:
            raise Exception("Must request one (and only one) run mode!")

        # Check argument combination makes sense
        # TODO: Does an error need to be thrown here, or just warn user
        # that --clobber will be ignored?
        if self.clobber and not self.all:
            parser.error('Clobber keyword can only be set if all specified (-a,--all)!')

        # If redo, platelist and ifudesignlist MUST be provided
        # TODO: Is this right?
        if self.redo and (not self.platelist or not self.ifudesignlist or not self.modelist):
            parser.error('Must provide platelist, ifudesignlist, and modelist for redo!')

        # If not redo, platelist, ifudesignlist, and modelist should be ignored
        if not self.redo and (self.platelist or self.ifudesignlist or self.modelist):
            parser.error('platelist/ifudesignlist/modelist can only be set if redo specified!')

        # If running all or daily, make sure lists are set to None such
        # that drpcomplete will search for all available DRP files and
        # make sure that the needed information is up-to-date
        if self.all or self.daily:
            self.platelist = None
            self.ifudesignlist = None
            self.modelist = None

        # Create and update the drpcomplete file if necessary
        # TODO: Currently does not allow for force=True
        # TODO: Assumes mangaver = drpver!!
        self.drpc = drpcomplete(self.outver, self.mangaver, nsa_cat=self.nsa_cat)
        self.drpc.update(platelist=self.platelist, ifudesignlist=self.ifudesignlist,
                         combinatorics=self.combinatorics)

        # Hereafter, self.platelist and self.ifuplatelist are the INPUT
        # values, whereas self.drpc.platelist and self.drpc.ifuplatelist
        # are the actual DRP reductions to analyze (including
        # combinatorics and file existence tests).  These should be
        # combined with self.modelist to get the DRP files to analyze.
        # See, e.g., selected_drpfile_list().

        # TODO: Is this needed if submit=False?
#        self.queue = pbs.queue(verbose=not self.quiet)

        # Run the selected mode
        # TODO: Redo should automatically clobber, right?
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
        path = self.file_path(plate, ifudesign)
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
                          help="runs dap for all plates/ifudesigns not currently running or done",
                          action="store_true")
        mode.add_argument("-r", "--redo",
                          help="runs dap for all plates/ifudesigns in platelist/ifudesignlist "
                          "regardless of state", action="store_true")

        # Run-mode optional arguments
        parser.add_argument("-c", "--clobber",
                            help="if all selected, will run dap for all plates/ifudesigns in "
                            "platelist/ifudesignlist regardless of state", action="store_true",
                            default=False)
        parser.add_argument("-q", "--quiet", help="turn off verbosity", action="store_true",
                            default=False)
        parser.add_argument("-o", "--outver", type=str,
                            help="optional output version, different from product version",
                            default="trunk")
#        parser.add_argument("-i", "--idlutilsver", type=str, help="version of idlutils to use",
#                            default=None)
#        parser.add_argument("-v", "--drpver", type=str, help="version of mangadrp used to produce "
#                            "files to process", default=None)
        parser.add_argument("-v", "--mangaver", type=str, help="version of manga module to use for"
                            "processing (also sets mangadrp files to search for)", default=None)
        parser.add_argument("-p", "--platelist", type=str, help="set list of plates to reduce",
                            default=None)
        parser.add_argument("-b", "--ifudesignlist", type=str, help="set list of ifus to reduce",
                            default=None)
        parser.add_argument("-m", "--modelist", type=str, help="set list of DRP output modes to "
                            "reduce (CUBE or RSS)", default=None)
        parser.add_argument("-x", "--combinatorics", help="force execution of all permutations of "
                            "the provided lists", action="store_true", default=False)
        parser.add_argument("-k", "--nsacat", type=str, help="path to NSA catalog to use",
                            default=None)

        # Cluster arguments
        parser.add_argument('-l', '--label', type=str, help='label for cluster job',
                            default='mangadap')
        parser.add_argument('-n', '--nodes', type=int, help='number of nodes to use in cluster',
                            default=18)
        parser.add_argument('-f', '--fast', dest='qos', type=str, help='qos state',
                            default=None)
        parser.add_argument('-u', '--umask', type=str, help='umask bit for cluster job',
                            default='002')
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
        # previously

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
   
        # Set the NSA catalog path
        if self.arg.nsacat is not None:
            self.nsa_cat = self.arg.nsacat

        # Set the versions to use; self.outver and self.mangaver will
        # never be none here because of procedures in __init__().
        # Will OVERWRITE existing input from __init__()
        if self.arg.outver is not None:
            self.outver = self.arg.outver
        if self.arg.mangaver is not None:
            self.mangaver = self.arg.mangaver

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
    

    # ******************************************************************
    #  MaNGA DAP RUN-MODE ROUTINES
    # ******************************************************************


    def run_daily(self):
        """
        Run daily reductions for next mjd completed by the DAP.
        
        Automatically run by crontab.
        """

        # raise Exception('Daily mode (-d,--daily) not yet implemented.')
        self.nodes = 1
        self.qos='sdss-fast'
        self.label = '{0}_daily'.format(self.label)

        # Select the daily DRP files and write the script files (do not
        # clobber)
        drpfiles = self.select_daily()
        scriptfiles, stdoutfiles, stderrfiles = self.prepare_for_analysis(drpfiles)

        # If submit is true, submit the scripts to the queue?  What is
        # done otherwise?


    def run_all(self, clobber=False):
        """
        Run all reductions.
    
        clobber: If True, all analyses will be run, regardless of their
        current state.  If false, only those analyses not currently
        running or done will be performed.

        This mode is manually called by humans.
        """

        # raise Exception('All mode (-a,--all) not yet implemented.')

        # TODO: Should clobber override self.clobber?
        self.label = '{0}_clobber'.format(self.label) if clobber else '{0}_all'.format(self.label)

        drpfiles = self.select_all(clobber=clobber)
        scriptfiles, stdoutfiles, stderrfiles = self.prepare_for_analysis(drpfiles, clobber=clobber)

        # If submit is true and (*.done does not exist or clobber=True),
        # submit the scripts to the queue?  What is done otherwise?


    def run_redo(self, clobber=True):
        """
        Run all reductions not currently running or done.
        
        Manually called by humans.
        """

        #raise Exception('Redo mode (-r,--redo) not yet implemented.')

        self.label = '{0}_redo'.format(self.label)

        drpfiles = self.select_redo()
        scriptfiles, stdoutfiles, stderrfiles = self.prepare_for_analysis(drpfiles, clobber=clobber)

        # If submit is true, submit the scripts to the queue?  What is
        # done otherwise?
    
 
    # ******************************************************************
    #  Reduction Management
    # ******************************************************************


    # TODO: What is the use case for simple=False?  This is not really a
    # rundap-specific routine.  Should be in a lower-level class.
    def product_version(self, product='mangadap', simple=False):
        """
        Gets the version for the SDSS-III or SDSS-IV product.  The
        default product is mangadap.
        """

        # TODO: Where does subprocess get imported?  Always imported?

        # Expects to find an executable called {$product}_version that
        # reports the SDSS-III/SDSS-IV product version.  If simple=True,
        # only the first element of the reported version is set to
        # 'version'
        try:
            version = subprocess.check_output('%s_version' % product, shell=True)
            version = version.split(' ')[0].rstrip('\n') if simple else version.rstrip('\n')
        except:
            version = None

        return version


    # TODO: This is not really a rundap-specific routine.  Should be in
    # a lower-level class.
    def module_version(self, product='mangadap'):
        """
        Return the module version for the specified product.  The
        default product is mangadap.
        """
        
        try:
            modules = environ['LOADEDMODULES']
        except:
            modules = None
            return None
        # TODO: Re-raise the exception?
      
        # Parse the loaded version(s) of product
        versions = [module.split('/')[1] for module in modules.split(':')
                        if module.split('/')[0]==product]

        # If there is more than one version or no versions return None
        if len(versions) != 1:
            if len(versions) > 1:
                print('Multiple versions found for module {0}'.format(product))
            else:
                print('Module {0} is not loaded'.format(product))
            return None

        # Return the version
        return versions[0]


    def file_root(self, plate, ifudesign, mode, stage='dap'):
        """
        Generate the root name of the MaNGA DAP file for a given
        stage/plate/ifudesign.
        """

        return 'manga{0}-{1}-{2}-LOG{3}'.format(stage, plate, ifudesign, mode)


    def file_path(self, plate, ifudesign):
        """
        Set the path name for the ouptut of the MaNGA DAP for a given
        plate/ifudesign.
        """

        return os.path.join(self.manga_spectro_analysis, self.outver, str(plate), str(ifudesign))


    def set_status(self, plate, ifudesign, mode, stage='dap', status='queued'):
        """
        Generate a touch file that signifies the status for a given
        plate/ifudesign/mode.
        """
        
        # Check that the path exists, creating it if not
        self._check_path(plate, ifudesign)
                
        # Touch status file
        root = self.file_root(plate, ifudesign, mode, stage)
        statfile = os.path.join(self.file_path(plate, ifudesign),'{0}.{1}'.format(root,status))
        file = open(statfile,'w')
        file.close()


    def dap_complete(self, drpf, stage='dap'):
        """
        Determine if the DAP was successfully completed on a given DRP
        file.

        Condition that DAP is complete is that the path exists and that
        within the path, the *.done touch file exits.
        """

        path = self.file_path(drpf.plate, drpf.ifudesign)
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

        drplist = self.full_drpfile_list()
        n_drp = len(drplist)
        print(n_drp)
#        for i in range(0,n_drp):
#            print(drplist[i].file_name())
#            print(drplist[i].created_today())

        # Select the files created today
        return [ drplist[i] for i in range(0,n_drp) if drplist[i].created_today() ]


    def select_all(self, clobber=False):
        """
        Select all the drp files.  If clobber=False, omit files that
        have been successfully analyzed.
        """

        drplist = self.full_drpfile_list()
        n_drp = len(drplist)
        print(n_drp)
#        for i in range(0,n_drp):
#            print(drplist[i].file_name())
#            print(self.dap_complete(drplist[i]))

        if clobber:
            return drplist

        return [ drplist[i] for i in range(0,n_drp) if not self.dap_complete(drplist[i]) ]


    def select_redo(self, clobber=True):

        drplist = self.selected_drpfile_list()

        if clobber:
            return drplist

        n_drp = len(drplist)
        return [ drplist[i] for i in range(0,n_drp) if not self.dap_complete(drplist[i]) ]


    def full_drpfile_list(self):
        """
        Generate the list of drpfiles based on the drpcomplete database.
        
        DRP files constructed based on the drpcomplete database are
        expected to exist.
        """

        # Number of plates (CUBE only)
        n_plates = len(self.drpc.data['PLATE'])

        # Create the list of CUBE DRP files
        # TODO: Assumes mangaver = drpver!!
        drplist = [ drpfile(self.mangaver, self.drpc.data['PLATE'][i],
                            self.drpc.data['IFUDESIGN'][i], 'CUBE') for i in range(0,n_plates) ]

        # Add the list of RSS DRP files
        # TODO: Assumes mangaver = drpver!!
        drplist = drplist + [ drpfile(self.mangaver, self.drpc.data['PLATE'][i],
                                      self.drpc.data['IFUDESIGN'][i], 'RSS')
                              for i in range(0,n_plates) if self.drpc.data['MODES'][i] == 2 ]
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
            except ValueError, e:
                getcube = False

        if getcube:
            # TODO: Assumes mangaver = drpver!!
            drplist = [ drpfile(self.mangaver, self.drpc.platelist[i], self.drpc.ifudesignlist[i],
                                'CUBE') for i in range(0,n_plates) ]
        else:
            drplist = []

        # Add the list of RSS DRP files, if requested (implicitly or
        # otherwise)
        if self.modelist is not None:
            try:
                index = self.modelist.index('RSS')
            except ValueError, e:
                return drplist                  # List complete

        # TODO: Assumes mangaver = drpver!!
        drplist = drplist + [ drpfile(self.mangaver, self.drpc.platelist[i],
                                      self.drpc.ifudesignlist[i], 'RSS')
                for i in range(0,n_plates)
                    if self.drpc.data['MODES'][self.drpc.entry_index(self.drpc.platelist[i],
                                               self.drpc.ifudesignlist[i])] == 2 ]

        return drplist


    def write_compute_script(self, plate, ifudesign, mode, stage='dap', clobber=False):
        """
        Write the MaNGA DAP script file for a given plate, ifudesign,
        and mode.
        """

        # Check the path exists, creating it if not
        self._check_path(plate, ifudesign)

        # Generate the path and root name of the output files
        path = self.file_path(plate, ifudesign)
        root = self.file_root(plate, ifudesign, mode, stage)
            
        # Get module name
        # TODO: Set outver?
#        module_version = self.module_version()
        module_version = self.outver

        # Fault if no module version is available
        if module_version is None:
            raise Exception('No DAP module version!')
        
        # Open the script file        
        scriptfile = os.path.join(path,root)
        stdoutfile = '{0}.out'.format(scriptfile)
        stderrfile = '{0}.err'.format(scriptfile)

        # Script file already exists, so just return
        if os.path.exists(scriptfile) and not clobber:
            return scriptfile, stdoutfile, stderrfile

        file = open(scriptfile, 'w')

        # Write the date as a commented header line
        file.write('# Auto-generated batch file {0}\n'.format(
                            time.strftime("%a, %d %b %Y %H:%M:%S +0000",time.localtime())))
        file.write('\n')

        # Load the specified manga version (may just reload already loaded
        # version)
        if self.mangaver:
            file.write('module switch manga manga/{0}\n'.format(self.mangaver))

        # Load the version of the DAP that is current at the time when
        # rundap was called (may just reload already loaded version)
        # TODO: Load current version as opposed to outver?
        file.write('module switch mangadap mangadap/{0}\n'.format(module_version))

#        # Load the specified DRP version (may just reload already loaded
#        # version)
#        if self.drpver:
#            file.write('module unload mangadrp \n')
#            file.write('module switch mangadrp/{0} '.format(self.drpver))

        # Create the started and running touch files
        # TODO: Is there ever a case when *.started exists, but
        # *.running does not?  Why are they both necessary?
        startfile = '{0}.started'.format(scriptfile)
        runfile = '{0}.running'.format(scriptfile)
        file.write('\n')
        file.write('touch {0}'.format(startfile))
        file.write('\n')
        file.write('touch {0}'.format(runfile))
        file.write('\n')

        # Command that runs the DAP
        parfile = self.parameter_file(plate, ifudesign, mode, stage)
        file.write('echo \" manga_dap, par=\'{0}\', /nolog \" | idl \n'.format(parfile))
        file.write('\n')

        #file.write('setStatusDone -f "{0}" \n'.format(errfile))
        donefile = '{0}.done'.format(scriptfile)
        file.write('touch {0}'.format(donefile))
        file.write('\n')
        file.close()

        # Return the script file, file for stdout, and file for stderr
        return scriptfile, stdoutfile, stderrfile
    

    def parameter_file(self, plate, ifudesign, mode, stage='dap'):
        """Get the name of the parameter file."""
        path = self.file_path(plate, ifudesign)
        root = self.file_root(plate, ifudesign, mode, stage)
        parfile = '{0}.par'.format(os.path.join(path,root))
        return parfile


    def prepare_for_analysis(self, drpfiles, status=None, clobber=False):
        """
        Cycle through the provided DRP file objects, creating a script
        for each.
        """

        nn = len(drpfiles)
        scriptfiles = []
        stdoutfiles = []
        stderrfiles = []

        if nn == 0:
            return scriptfiles, stdoutfiles, stderrfiles

        for i in range(0,nn):

            # Set the status, if provided; be sure to check the path
            # otherwise
            if status is None:
                self._check_path(drpfiles[i].plate, drpfiles[i].ifudesign)
            else:
                self.set_status(drpfiles[i].plate, drpfiles[i].ifudesign, drpfiles[i].mode, status)

            # Create the parameter file
            parfile = self.parameter_file(drpfiles[i].plate, drpfiles[i].ifudesign,
                                          drpfiles[i].mode)
            print(parfile)
            self.drpc.write_par(parfile, drpfiles[i].mode, plate=drpfiles[i].plate,
                                ifudesign=drpfiles[i].ifudesign)

            # Write the compute script
            try:
                sf, of, ef = self.write_compute_script(drpfiles[i].plate, drpfiles[i].ifudesign,
                                                       drpfiles[i].mode, clobber=clobber)
            except Exception,e:
                print("Exception: %s" % str(e))
                if i < nn-1:
                    print("Skipping to next DRP file.")
            else:
                scriptfiles = scriptfiles + [sf]
                stdoutfiles = stdoutfiles + [of]
                stderrfiles = stderrfiles + [ef]

        # Return the list of script, stdout, and stderr files
        return scriptfiles, stdoutfiles, stderrfiles

            


        

