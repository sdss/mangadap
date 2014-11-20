from __future__ import division
from __future__ import print_function

import time
import os.path
import pbs.queue

from os import environ, makedirs
#from os.path import exists, join, isdir
#from pbs import queue

from argparse import ArgumentParser


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
        set_arg         - set command line arguments

        run_daily       - process newly completed DRP reductions as
                          determined by the daily assessment

        run_all         - process all DRP reductions not currently
                          running or done unless clobber is specified to
                          force processing of all files

        run_redo        - redo the processing of selected DRP reductions
                
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
                
    """

    # TODO: outver is actually just the directory for output to use in
    # place of the current product version 
    def __init__(self,
                # Run mode options
                daily=None, all=None, clobber=None, redo=None, console=None, quiet=None,
                # Override defaults
                outver=None, idlutilsver=None, drpver=None, ifudesignlist=None, platelist=None,
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
            idlutilsver - idlutils version to use
            drpver - DRP version to use
            ifudesignlist - specified list of ifudesign to analyze
            platelist - specified list of plates to analyze

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
        self.idlutilsver = self.product_version(simple=True, product='idlutils')
        self.corever = self.product_version(simple=True, product='mangacore')
        self.drpver = self.product_version(simple=True, product='mangadrp')
        self.dapver = self.product_version(simple=True, product='mangadap')

        # Use them or use input versions
        self.outver = self.dapver if outver is None else outver
        self.idlutilsver = self.idlutilsver if idlutilsver is None else idlutilsver
        self.drpver = self.drpver if drpver is None else drpver

        # List of files to analyze
        self.platelist = platelist
        self.ifudesignlist = ifudesignlist

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
       
        # Assign console keywords
        if console:
            self.set_arg()

        # TODO: Is this needed if submit=False?
        self.queue = pbs.queue(verbose=not self.quiet)

        # Run the selected mode
        # TODO: Redo should automatically clobber, right?
        if self.daily:
            self.run_daily()
        if self.all:
            self.run_all(clobber=self.clobber)
        if self.redo:
            self.run_redo()


    # ******************************************************************
    #  COMMAND-LINE UTILITIES
    # ******************************************************************
    
    def set_arg(self):
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
        parser.add_argument("-v", "--outver", type=str,
                            help="optional output version, different from product version",
                            default="trunk")
        parser.add_argument("-i", "--idlutilsver", type=str, help="version of idlutils to use",
                            default=None)
        parser.add_argument("-b", "--ifudesignlist", type=str, help="set list of ifus to reduce",
                            default=None)
        parser.add_argument("-p", "--platelist", type=str, help="set list of plates to reduce",
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
        
        # Check required arguments are present (now done using
        # required=True in mode)
#        if not self.arg.daily and not self.arg.all and not self.arg.redo:
#           parser.error('Must specify daily (-d,--daily), all (-a,--all), or redo (-r,--redo)!')

        # Check argument combination makes sense
        # TODO: Does an error need to be thrown here, or just warn user
        # that --clobber will be ignored?  What about usage with redo?
        if self.arg.clobber and not self.arg.all:
            parser.error('Clobber keyword can only be set if all specified (-a,--all)!')

        # If redo, platelist and ifudesignlist MUST be provided
        if self.arg.redo and (not self.arg.platelist or not self.arg.ifudesignlist):
            parser.error('Must provide platelist and ifudesignlist for redo (-r,--redo)!')

        # If not redo, platelist and ifudesignlist should be ignored
        if not self.arg.redo and (self.arg.platelist or self.arg.ifudesignlist):
            parser.error('platelist/ifudesignlist can only be set if redo specified (-r,--redo)!')

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
        if not self.outver:
            self.outver = self.arg.outver if self.arg.outver else self.dapver

        if self.ifudesignlist is None:
            self.ifudesignlist = eval(self.arg.ifudesignlist)
        if type(self.ifudesignlist == int:
            self.ifudesignlist = [self.ifudesignlist]

        if self.platelist is None:
            self.platelist = eval(self.arg.platelist)
        if type(self.platelist == int:
            self.platelist = [self.platelist]

        if self.idlutilsver is None:
            self.idlutilsver = self.arg.idlutilsver
            
        # Set queue keywords
        if self.umask is None:
            self.umask = self.arg.umask
        if self.nodes is None:
            self.nodes = self.arg.nodes
        if self.walltime is None:
            self.walltime = self.arg.walltime
        if self.qos is None:
            self.qos = self.arg.qos
        if self.umask is None:
            self.umask = self.arg.umask
        if self.hard is None:
            self.hard = self.arg.hard
        if self.submit is None:
            self.submit = self.arg.submit
        

    # ******************************************************************
    #  MaNGA DAP RUN-MODE ROUTINES
    # ******************************************************************

    def run_daily(self):
        """
        Run daily reductions for next mjd completed by the DAP.
        
        Automatically run by crontab.
        """

        raise Exception('Daily mode (-d,--daily) not yet implemented.')

        self.nodes = 1
        self.qos='sdss-fast'
        self.label = ''.join([self.label,'_daily'])

        # Run-mode daily body...

        pass


    def run_all(self, clobber=False):
        """
        Run all reductions.
    
        clobber: If True, all analyses will be run, regardless of their
        current state.  If false, only those analyses not currently
        running or done will be performed.

        This mode is manually called by humans.
        """

        # TODO: Should clobber override self.clobber?
        if clobber:
            self.label = ''.join([self.label,'_clobber'])
        else:
            self.label = ''.join([self.label,'_all'])

        # hack for one plate, ifudesign
        plate = {'plate':int(self.platelist[0]), 'ifudesign':int(self.ifudesignlist[0])}
        self.set_status(plate,status='manual')

        # TODO: Add try/except block
        self.generate_script(plate)


    def run_redo(self, clobber=True):
        """
        Run all reductions not currently running or done.
        
        Manually called by humans.
        """

        raise Exception('Redo mode (-r,--redo) not yet implemented.')

        self.label = ''.join([self.label,'_redo'])

        pass
    
 
    # ******************************************************************
    #  Reduction Management
    # ******************************************************************

    # TODO: What is the use case for simple=False?
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


    def module_version(self, product='mangadap'):
        """
        Return the module version for the specified product.  The
        default product is mangadap.
        """
        
        try:
            modules = environ['LOADEDMODULES']
        except:
            modules = None
        # TODO: Re-raise the exception?
      
        # No modules, so no version
        if modules is None:
            return None

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


    def set_file_root(self, plate, stage='dap'):
        """
        Generate the root name of the MaNGA DAP file for a given
        stage/plate/ifudesign.
        """

        return 'manga{0}-{1}-{2}'.format(stage, plate['plate'], plate['ifudesign'])


    def set_path(self, plate):
        """
        Set the path name for the ouptut of the MaNGA DAP for a given
        plate/ifudesign.
        """

        return os.path.join(self.manga_spectro_analysis, self.outver, str(plate['plate']), 
                            str(plate['ifudesign']))


    def set_status(self,plate,status='queued', stage='dap'):
        """
        Generate the touch files, setting the status for a given
        plate/ifudesign.
        """
        
        # Set path and check it exists
        path = set_path(plate)
        if not os.path.isdir(path):
            makedirs(path)
                
        # Touch status file
        root = set_file_root(plate,stage)
        statfile = os.path.join(path,'.'.join([root,status]))
        file = open(statfile,'w')
        file.close()


    def generate_script(self, plate, stage='dap'):
        """
        Generate the MaNGA DAP script files for a given plate and
        ifudesign.
        """

        # Generate the path, the script file name, and the spool files
        path = set_path(plate)
        root = set_file_root(plate,stage)
            
        scriptfile = os.path.join(path,root)
        stdoutfile = ''.join([scriptfile,'.out'])
        stderrfile = ''.join([scriptfile,'.err'])
        startfile = ''.join([scriptfile,'.started'])
        runfile = ''.join([scriptfile,'.running'])
        donefile = ''.join([scriptfile,'.done'])
        
        # Get module name
        # TODO: Set outver?
        module_version = self.module_version()

        if module_version is None:
            raise Exception('No DAP module version!')
        
        # Open the script file        
        file = open(scriptfile, 'w')

        # Write the date as a commented header line
        file.write(''.join(["# Auto-generated batch file ",
                            time.strftime("%a, %d %b %Y %H:%M:%S +0000",time.localtime()),"\n"]))
        file.write('\n')

        # TODO: Why is this necessary?  Just reloads already loaded
        # module
        file.write('module unload mangadap \n')
        file.write('module load mangadap/{0} \n'.format(module_version))

        # Load the specified idlutils (may just reload already loaded
        # version)
        if self.idlutilsver:
            file.write('module unload idlutils \n')
            file.write('module switch idlutils/{0} '.format(self.idlutilsver))

        # Load the specified DRP version (may just reload already loaded
        # version)
        if self.drpver:
            file.write('module unload mangadrp \n')
            file.write('module switch mangadrp/{0} '.format(self.drpver))

        # Create the started and running touch files
        # TODO: Is there ever a case when *.started exists, but
        # *.running does not?  Why are they both necessary?
        file.write('\n')
        file.write('touch {0}'.format(startfile))
        file.write('\n')
        file.write('touch {0}'.format(runfile))
        file.write('\n')

        # Command that runs the DAP
        # TODO: Replace hardcoded number!
        # .format(stage,plate['plate'],plate['ifudesign']))
        file.write("echo \"manga_dap, 38, /nolog \" | idl \n")

        file.write('\n')
        #file.write('setStatusDone -f "{0}" \n'.format(errfile))
        file.write('touch {0}'.format(donefile))
        file.close()
        
        return scriptfile, stdoutfile, stderrfile
        
        

                          
