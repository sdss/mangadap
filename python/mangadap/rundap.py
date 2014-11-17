from __future__ import division
from __future__ import print_function
from os.path import exists
from argparse import ArgumentParser
from pbs import queue

__author__ = 'Kyle Westfall'

'''
	NAME:
		rundap
		
	PURPOSE:
		Automated procedures for MaNGA DAP processing at Utah
	
	METHODS:
		set_arg				- set command line arguments
		run_daily			- process the daily drp stage
		run_all				- process all MJDs not currently running or done unless clobber is specified to force process all MJDs
		run_redo			- redoes the processing on some MJDs
		
	REVISION HISTORY:
		   NOV-2014 - Architecture setup (Joel Brownstein <joelbrownstein@astro.utah.edu>)
		
'''

class rundap:
    
    def __init__(self, daily=None, all=None, clobber=None, redo=None, console=None, quiet=None,
                outver=None, idlutilsver=None, ifudesignlist=None, platelist=None,
                label='mangadap', nodes=18, qos=None,umask='007',walltime='240:00:00', hard=True, submit=True):
        self.daily = daily
        self.all = all
        self.clobber=clobber
        self.redo=redo
        self.quiet = quiet
        self.dapver = self.product_version(simple=True, product='mangadap')
        self.outver = outver if outver else self.dapver
        self.drpver = self.product_version(simple=True, product='mangadrp')
        self.corever = self.product_version(simple=True, product='mangacore')
        self.idlutilsver = idlutilsver if idlutilsver else self.getversion(simple=True, product='idlutils')
        self.ifudesignlist = ifudesignlist
        self.platelist = platelist
        #queue keywords
        self.label = label
        self.nodes = nodes
        self.qos = qos
        self.umask = umask
        self.walltime = walltime
        self.hard = hard
        self.submit = submit
        #console keywords
        if console: self.set_arg()
        self.queue = queue(verbose=not self.quiet)
        #do keywords
        if self.daily:
            self.label += "_daily"
            self.nodes = 1
            self.qos='sdss-fast'
        if self.all:
            self.label += "_clobber" if self.clobber else self.label += "_all"
        if self.redo: self.label += "_redo"
        if self.daily: self.run_daily()
        elif self.all:
            if self.clobber: self.run_clobber()
            else: self.run_all()
        elif self.redo: self.run_redo()

    
    """ ********************************************
        COMMAND-LINE UTILITIES
        ********************************************"""
    
    def set_arg(self):
        parser = ArgumentParser()
        mode = parser.add_mutually_exclusive_group()
        mode.add_argument("-d", "--daily", help="runs dap for next mjd (as an after burner to drp)",action="store_true")
        mode.add_argument("-a", "--all", help="runs dap for all plates/ifudesigns not currently running or done",action="store_true")
        mode.add_argument("-c", "--redo", help="runs dap for all plates/ifudesigns in platelist/ifudesignlist regardless of state",action="store_true", default=False)
        parser.add_argument("-q", "--quiet", help="turn off verbosity", action="store_true", default=False)
        parser.add_argument("-v", "--outver", type=str, help="optional different output reduction version than product version", default="trunk")
        parser.add_argument("-u", "--idlutilsver", type=str, help="version of idlutils to use", default=None)
        parser.add_argument("-i", "--ifudesignlist", type=str, help="set list of ifus to reduce", default=None)
        parser.add_argument("-p", "--platelist", type=str, help="set list of plates to reduce", default=None)

        #queue
        parser.add_argument('-l', '--label', type=str, help='set label of cluster job', default='mangadap')
        parser.add_argument('-n', '--nodes', type=int,help='set number of nodes to use in cluster', default=18)
        parser.add_argument('-f', '--fast', dest='qos',type=str, help='set qos state', default=None)
        parser.add_argument('-u', '--umask', type=str, help='set umask bit for cluster job', default='002')
        parser.add_argument('-w', '--walltime', type=str, help='set walltime for cluster job', default='240:00:00')
        parser.add_argument('-t', '--toughness', dest='hard', action='store_false', default=True, help='turn off hard keyword for cluster submission')
        parser.add_argument('-s', '--submit', help='turn off cluster submission',action='store_false', default=True)
        
        self.arg = parser.parse_args()
        
        # One of these args required.  
        if not self.arg.daily and not self.arg.all and not self.arg.redo: parser.error('Must specifiy one option of Daily, All,or Redo!')
        if self.arg.clobber and not self.arg.all: parser.error('Clobber keyword can only be set with the all keyword (rundap --all)')
        
        if self.daily==None: self.daily = self.arg.daily
        if self.all==None: self.all = self.arg.all
        if self.clobber==None: self.clobber = self.arg.clobber
        if self.redo==None: self.redo = self.arg.redo
        if self.quiet==None: self.quiet = self.arg.quiet
        if not self.outver:
            self.outver = self.arg.outver if self.arg.outver else self.dapver
        if self.ifudesignlist == None: self.ifudesignlist = self.arg.ifudesignlist
        if self.platelist == None: self.platelist = self.arg.platelist
        if self.idlutilsver == None: self.idlutilsver = self.arg.idlutilsver
            
        # Set queue keywords
        if self.umask == None: self.umask = self.arg.umask
        if self.nodes == None: self.nodes = self.arg.nodes
        if self.walltime == None: self.walltime = self.arg.walltime
        if self.qos == None: self.qos = self.arg.qos
        if self.umask == None: self.umask = self.arg.umask
        if self.hard == None: self.hard = self.arg.hard
        if self.submit == None: self.submit = self.arg.submit
        

    """ ********************************************
        Manga DAP Reductions
        ********************************************"""

    def run_daily(self):
        ''' Run daily reductions for next mjd .  Automatically run by crontab'''
        pass

    def run_all(self, clobber=False):
        ''' Run all reductions not currently running or done (or regardless of whether running or done if self.clobber is specified).
            Manually called by humans. '''
        
        plate = {'plate':int(self.platelist[0]),'ifudesign':int(self.ifudesignlist[0])} #hack for one plate, ifudesign
        self.set_status(plate,status='manual')
        self.generate_script(plate)

    def run_redo(self, clobber=False):
        ''' Run all reductions not currently running or done.  Manually called by humans. '''
        pass
    
 
    """ ********************************************
        Reduction Management
        ********************************************"""

    def product_version(self, product='mangadap', simple=False):
        '''Gets the version for the SDSS-III or SDSS-IV product.  Default product is mangadap.'''
        
        try:
            version =  subprocess.check_output('%s_version' % product ,shell=True)
            if simple: version = version.split(' ')[0].rstrip('\n')
            else: version = version.rstrip('\n')
        except: version = None

        return version

    def module_version(self, product='mangadap'):
    ''' Return the module versions for the specified product.  Default product is mangadap.'''
        
        try: modules = os.environ['LOADEDMODULES']
        except: modules = None
        
        if modules:
            versions = [module.split('/')[1] for module in modules.split(':') if module.split('/')[0]==name] if modules else []
            if len(versions) == 1: version = versions[0]
            elif len(versions):
                print('Multiple versions found for module {0}'.format(name))
                version = None 
            else: 
                print('Module {0} is not loaded'.format(name))
                version = None
        else: version = None
        
        return version

    def set_status(self,plate,status='queued', stage='dap'):
        ''' Generate the touch status files for a given plate'''
        
        # Generate status path and name
        statfile = 'manga{0}-{1}-{2}.{3}'.format(stage, plate['plate'], plate['ifudesign'], status)
        path = os.path.join(os.getenv('MANGA_SPECTRO_ANALYSIS'), self.outver, str(plate['plate']), str(plate['ifudesign']))
        
        # Check for path existence
        if not os.path.isdir(path): os.makedirs(path)
                
        # Write status file
        fullfile = os.path.join(path,statfile)
        file = open(fullfile,'w')
        file.close()
                          
    def generate_script(self,plate, stage='dap'):
        ''' Generate the manga script files for a given plate'''
        
        # Generate script path and name
        statfile = 'manga{0}-{1}-{2}'.format(stage,plate['plate'],plate['ifudesign'])
        path = os.path.join(os.getenv('MANGA_SPECTRO_ANALYSIS'), self.outver, str(plate['plate']), str(plate['ifudesign']))
            
        fullfile = os.path.join(path,statfile)
        outfile = os.path.join(path,statfile+'.out')
        errfile = os.path.join(path,statfile+'.err')
        
        # Get module name
        module_version = self.module_version()
        
        # Write script file        
        file = open(fullfile,'w')
        file.write('# Auto-generated batch file '+time.strftime("%a, %d %b %Y %H:%M:%S +0000", time.localtime())+"\n")
        file.write('\n')
        file.write('module unload mangadap \n')
        file.write('module load mangadap/{0} \n'.format(module_version))
        if self.idlutilsver:
            file.write('module unload idlutils \n')
            file.write('module switch idlutils/{0} '.format(self.idlutilsver))
        file.write('\n')
        file.write('touch {0}.started'.format(fullfile))
        file.write('\n')
        file.write('touch {0}.running'.format(fullfile))
        file.write('\n')

        file.write("echo \"manga_dap, 38, /nolog \" | idl \n") # Eventually this will replace the hardcoded number: .format(stage,plate['plate'],plate['ifudesign']))

        file.write('\n')
        file.write('setStatusDone -f "{0}" \n'.format(errfile))
        file.close()
        
        return fullfile, outfile, errfile
        
        

                          
