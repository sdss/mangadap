from __future__ import division
from __future__ import print_function
from __future__ import absolute_import

import sys
if sys.version > '3':
    long = int

import os.path
from os import environ, makedirs, walk
import numpy
from scipy.spatial import KDTree
from astropy.io import fits
from astropy import constants

# DAP imports
from mangadap.util.yanny import yanny, read_yanny
from mangadap.drpfile import drpfile, drpfile_list, parse_drp_file_name
from mangadap.util.parser import arginp_to_list
from mangadap.util.exception_tools import print_frame

__author__ = 'Kyle Westfall'


def list_to_csl_string(flist):
    """
    Convert a list to a comma-separated string.
    """
    n = len(flist)
    out=str(flist[0])
    for i in range(1,n):
        out += ', '+str(flist[i])
    return out



# TODO: Let the paths be defined, using default otherwise

class drpcomplete:
    """
    drpcomplete object with creates, updates, and reads data from the
    drpcomplete file, uses as a reference file for the DAP.

    Provided files have priority over default ones.
    platetargets file has priority over NSA catalog.

    REVISION HISTORY:
        20 Nov 2014: (KBW) Started
        01 Dec 2014: (KBW) Committed to SVN
        12 Jan 2015: (KBW) Allow for use of plateTargets file using yanny.py
    """


    def __init__(self, platelist=None, ifudesignlist=None, platetargets=None, nsa_cat=None,
                 drpver=None, drppath=None, dapver=None, dappath=None, readonly=False):
        """
        Initializes the class and its members.
    
        ARGUMENTS:
            
            drpver - DRP version to use
            platelist - specified list of plates to analyze
            ifudesignlist - specified list of ifudesign to analyze

            platetargets - plate targets file to use for data collection
            nsa_cat - NSA catalog to use for data collection

            force - Force updates, even if the data already exists
        """

        # Input properties
        self.drpver = self._default_drp_version() if drpver is None else str(drpver)
        self.redux_path = self._default_redux_path() if drppath is None else str(drppath)

        self.dapver = self._default_dap_version() if dapver is None else str(dapver)
        self.analysis_path = self._default_analysis_path() if dappath is None else str(dappath)
        
        if os.path.exists(self.file_path()):
            self._read_data()
        else:
            self.data = None
            self.nrows = None

        if readonly:
            self.readonly=True
            self.platelist = None
            self.ifudesignlist = None
            self.platetargets = None
            self.nsa_cat=None
            return

        self.platelist = arginp_to_list(platelist, evaluate=True)
        self.ifudesignlist = arginp_to_list(ifudesignlist, evaluate=True)

        platetargets = arginp_to_list(platetargets)

        if platetargets is not None and all( os.path.exists(p) for p in platetargets ):
            self.platetargets = platetargets
            self.nsa_cat = None
        elif nsa_cat is not None and os.path.exists(nsa_cat):
            self.platetargets = None
            self.nsa_cat = nsa_cat
        else:
            self.platetargets = self._default_plate_targets_file()
            self.nsa_cat = self._default_nsa_catalog()


    # ******************************************************************
    #  Default values
    # ******************************************************************


    def _default_drp_version(self):
        """
        Return the DRP version defined by the environmental variable.
        """
        return environ['MANGADRP_VER']


    def _default_redux_path(self):
        """Return the directory path used by the DRP."""

        # Make sure the DRP version is set
        if self.drpver is None:
            self.drpver = self._default_drp_version()

        return os.path.join(environ['MANGA_SPECTRO_REDUX'], self.drpver)


    def _default_dap_version(self):
        """
        Return the DRP version defined by the environmental variable.
        """
        return environ['MANGADAP_VER']


    def _default_analysis_path(self):
        """Return the directory path used by the DAP."""

        # Make sure the DAP version is set
        if self.dapver is None:
            self.dapver = self._default_dap_version()

        return os.path.join(environ['MANGA_SPECTRO_ANALYSIS'], self.dapver)


    def _default_plate_targets_file(self):
        """
        Return the default plateTargets file used to get the NSA
        information
        """
        return [os.path.join(environ['MANGACORE_DIR'], 'platedesign', 'platetargets',
                                'plateTargets-12.par'),
                os.path.join(environ['MANGACORE_DIR'], 'platedesign', 'platetargets',
                            'plateTargets-1.par') ]


    def _default_nsa_catalog(self):
        """
        Return the default NSA catalog
        """
        if self.drpver == 'v1_0_0':
            return os.path.join(environ['MANGAWORK_DIR'], 'manga', 'target', 'temp',
                                '12-nsa_v1b_0_0_v2.fits.gz')

        return os.path.join(environ['MANGA_TARGET'], 'input', 'nsa_v1_0_0.fits')


    # ******************************************************************
    #  Utility functions
    # ******************************************************************


    def _drp_mangaid(self, drplist):
        """Grab the mangaids from the DRP fits files."""
        nn = len(drplist)
        mangaid = []
        print("Gathering MANGA-IDs for DRP files...")
        for i in range(0,nn):
            mangaid_, ra, dec = drplist[i].object_data()
            mangaid = mangaid + [mangaid_]
        print("... done")
        mangaid = numpy.array(mangaid)
        return mangaid


    def _drp_info(self, drplist):
        """
        Grab the mangaid, object RA, and object DEC from the DRP fits
        files.
        """
        nn = len(drplist)
        mangaid = []
        objra = numpy.empty(nn)
        objdec = numpy.empty(nn)
        print("Gathering DRP header data...")
        for i in range(0,nn):
            mangaid_, objra[i], objdec[i] = drplist[i].object_data()
            mangaid = mangaid + [mangaid_]
        print("... done")
        mangaid = numpy.array(mangaid)
        return mangaid, objra, objdec


    def _read_platetargets(self, verbose=False):
        """Read the platetargets file and the number of entries."""
        for p in self.platetargets:
            if not os.path.exists(p):
                raise Exception('Cannot open platetargets file: {0}'.format(p))

        nn = len(self.platetargets)
        par_data = yanny(self.platetargets[0])

        for i in range(1,nn):
            tmp_par = read_yanny(self.platetargets[1])
            par_data.append(tmp_par, add_to_file=False)

#       n_par = numpy.shape(par_data['PLTTRGT']['mangaid'])[0]
        par_data.convert_to_numpy_array()
        n_par = par_data.size('PLTTRGT')

        return par_data, n_par


    def _match_platetargets(self, def_veldisp):
        """
        Read the platetargets file to find the observed DRP objects and
        their ancillary information.
        """

        # Read the platetargets file
        print('Reading platetarget file(s)...')
        par_data, n_par = self._read_platetargets(verbose=True)
        print('DONE')

        # Initialize the output arrays (needed in case some DRP targets not found)
        n_drp = len(self.platelist)
        mangaid = numpy.empty(n_drp, dtype=par_data['PLTTRGT']['mangaid'].dtype)
        objra = numpy.empty(n_drp, dtype=par_data['PLTTRGT']['target_ra'].dtype)
        objdec = numpy.empty(n_drp, dtype=par_data['PLTTRGT']['target_dec'].dtype)
        nsaid = numpy.full(n_drp, -1, dtype=numpy.int32)
        vel = numpy.zeros(n_drp, dtype=par_data['PLTTRGT']['nsa_redshift'].dtype)
        veldisp = numpy.full(n_drp, def_veldisp, dtype=par_data['PLTTRGT']['nsa_vdisp'].dtype)
        ell = numpy.zeros(n_drp, dtype=par_data['PLTTRGT']['nsa_sersic_ba'].dtype)
        pa = numpy.zeros(n_drp, dtype=par_data['PLTTRGT']['nsa_sersic_phi'].dtype)
        Reff = numpy.zeros(n_drp, dtype=par_data['PLTTRGT']['nsa_sersic_th50'].dtype)

        print('Searching platetargets file for observed galaxies...')
        for i in range(0,n_drp):
            indx = numpy.where((par_data['PLTTRGT']['plateid'] == self.platelist[i]) &
                               (par_data['PLTTRGT']['ifudesign'] == self.ifudesignlist[i]))

            if len(indx[0]) == 0:
                raise Exception('Could not find plate={0}, ifudesign={1} in {2}.'.format(
                                self.platelist[i], self.ifudesignlist[i], self.platetargets))
            mangaid[i] = par_data['PLTTRGT']['mangaid'][indx][0]
            objra[i] = par_data['PLTTRGT']['target_ra'][indx][0]
            objdec[i] = par_data['PLTTRGT']['target_dec'][indx][0]

            nsaid[i] = numpy.int32(mangaid[i].decode("utf-8").split('-')[1])
            
            vel[i] = par_data['PLTTRGT']['nsa_redshift'][indx][0] * constants.c.to('km/s').value
            ell[i] = 1.0-par_data['PLTTRGT']['nsa_sersic_ba'][indx][0]
            pa[i] = par_data['PLTTRGT']['nsa_sersic_phi'][indx][0]
            Reff[i] = par_data['PLTTRGT']['nsa_sersic_th50'][indx][0]

            if par_data['PLTTRGT']['nsa_vdisp'][indx][0] > 0.:
                veldisp[i] = par_data['PLTTRGT']['nsa_vdisp'][indx][0]

        print('DONE')

        return mangaid, objra, objdec, nsaid, vel, veldisp, ell, pa, Reff


    def _read_nsa(self, verbose=False):
        """Read the NSA catalog and the number of rows."""
        if not os.path.exists(self.nsa_cat):
            raise Exception('Cannot open NSA catalog: {0}'.format(self.nsa_cat))
        hdulist = fits.open(self.nsa_cat)
        nsa_data = hdulist[1].data
        n_nsa = hdulist[1].header['NAXIS2']

        try:
            i = hdulist[1].columns.names.index('VELDISP')
        except ValueError:
            #print_frame('ValueError')
            has_veldisp = False
        else:
            has_veldisp = True

        if verbose and not has_veldisp:
            print("WARNING: NSA catalog does not have velocity dispersion.  Using default.")

        hdulist.close()

        return nsa_data, n_nsa, has_veldisp


    def _match_nsa(self, objra, objdec, def_veldisp, match_r):
        """
        Read the NSA catalog to find the observed DRP objects and their
        ancillary information.
        """

        # Read the NSA catalog
        print('Reading NSA...')
        nsa_data, n_nsa, nsa_has_veldisp = self._read_nsa(verbose=True)
        print('DONE')

        n = len((nsa_data['RA']).ravel())

        # Set the coordinates to a KDTree
        data = numpy.zeros((n, 2), dtype=numpy.float64)
        data[:,0] = (nsa_data['RA']).ravel()
        data[:,1] = (nsa_data['DEC']).ravel()
        print(numpy.shape(data))
        print(len(data))

        kdcoo = KDTree(data)
        #zip( (nsa_data['RA']).ravel(), (nsa_data['DEC']).ravel() ))
        #print(kdcoo.data.shape)

        # Initialize the output arrays (needed in case some DRP targets not found)
        n_drp = len(objra)
        nsaid = numpy.zeros(n_drp, dtype=numpy.int32)
        vel = numpy.zeros(n_drp, dtype=numpy.float64)
        veldisp = numpy.full(n_drp, def_veldisp, dtype=numpy.float64)
        ell = numpy.zeros(n_drp, dtype=numpy.float64)
        pa = numpy.zeros(n_drp, dtype=numpy.float64)
        Reff = numpy.zeros(n_drp, dtype=numpy.float64)

        # Match the DRP object coordinates to the NSA catalog
        d2r = numpy.pi/180.
        match_d = (match_r/3600.)**2.0

        # Search for the observed galaxies in the NSA using the KD tree
        print('Searching NSA for observed galaxies...')
        nsa_index = numpy.empty(n_drp, dtype=numpy.int32)
        for i in range(0,n_drp):
            # Get the 10 nearest neighbors using the KDTree
            if n_nsa > 10:
                #print('Finding neighbors')
                nneigh_d, nneigh_i = kdcoo.query( numpy.array([objra[i],objdec[i]]), k=10)
            else:
                nneigh_i = numpy.array([i for i in range(0,n_nsa)])

            # For the nearest neighbors, find the ones that are within
            # match_d (correcting for the cos(DEC) term)
            #print('Calculating distances')
            cosd=numpy.cos(objdec[i]*d2r)
            dd = (((kdcoo.data[nneigh_i,0]-objra[i])*cosd)**2.0
                   + (kdcoo.data[nneigh_i,1]-objdec[i])**2.0)
            if len(numpy.where(dd < match_d)) > 1:
                raise Exception('Multiple matches in the NSA catalog: ra={0},dec={1}'.format(
                                                                              objra[i], objdec[i]))

            di = numpy.argsort(dd)
            #print(dd, di)
            if dd[di[0]] < match_d:
                nsa_index[i] = nneigh_i[di[0]]
            else:
                print('Could not find match in NSA catalog for plate={0},ifudesign={1}!'.format(
                                    self.platelist[i], self.ifudesignlist[i]))
                nsa_index[i] = -1

        #print(nsa_index)
        print('DONE')

        noneg = numpy.where(nsa_index >= 0)
        #print(noneg)

        # Pull out the NSA data
        nsaid[noneg] = nsa_data['NSAID'][nsa_index[noneg]]
        vel[noneg] = nsa_data['Z'][nsa_index[noneg]]*constants.c.to('km/s').value
        ell[noneg] = 1.0-nsa_data['SERSIC_BA'][nsa_index[noneg]]
        pa[noneg] = nsa_data['SERSIC_PHI'][nsa_index[noneg]]
        Reff[noneg] = nsa_data['SERSIC_TH50'][nsa_index[noneg]]
        if nsa_has_veldisp:
            veldisp[noneg] = nsa_data['VELDISP'][nsa_index[noneg]]

        return nsaid, vel, veldisp, ell, pa, Reff


    def _all_data_exists(self, quiet=True):
        """
        Determine if the data for all the plates/ifudesigns is already
        present in the current drpcomplete file.
        """
        nn = len(self.platelist)
        for i in range(0,nn):
            try:
                index = self.entry_index(self.platelist[i], self.ifudesignlist[i])
            except Exception as e:
                if not quiet:
                    print_frame('Exception')
                    print(e)
                return False
        return True


    def _read_data(self):
        """Read the data and the number of rows."""
        hdulist = fits.open(self.file_path())
        self.data = hdulist[1].data
        self.nrows = hdulist[1].header['NAXIS2']
        hdulist.close()
        print('Read data: %d rows' % self.nrows)


    def _confirm_access(self, reread):
        """
        Check drpcomplete file is accessible.  Read the data if not yet
        read.
        """
        inp = self.file_path()
        if not os.path.exists(inp):
            raise Exception('Cannot open file: {0}'.format(inp))

        if self.data is None or reread:
            self._read_data()


    def _find_completed_reductions(self, mindesign=19, combinatorics=False):
        """
        Search the DRP path for reduced CUBE files.

        Function allows for one or both of self.platelist and
        self.ifudesignlist to be None:
            If both are None
                all available CUBE files are used to create the
                platelist and ifudesignlist.
            If one is None
                all available CUBE files within the provided list of the
                one that is not None are used to create the platelist
                and ifudesignlist.  E.g., if platelist=[7443] and
                ifudesignlist is None, all CUBE files with platelist =
                7443 are chosen.

            If both are not None and they have different lengths (or the
            same length and combinatorics=True)
                all available CUBE files within the constraints of both
                lists are selected. E.g., if platelist=[7443,7495] and
                ifudesignlist=[12704], the CUBE files with
                (plate,ifudesign)=[(7443,12704),(7459,12704) are chosen

            If both are not None and they have the same length and
            combinatorics=False
                all available CUBE files in the matched lists are
                chosen.  E.g. if platelist=[7443,7495] and
                ifudesignlist=[12703,12703], the CUBE files with
                (plate,ifudesign)=[(7443,12704),(7459,12704) are chosen
        """
        path = self.redux_path
        matchedlist = (self.platelist is not None and self.ifudesignlist is not None and
                       (len(self.platelist) == len(self.ifudesignlist) and not combinatorics))

        plates = []
        ifudesigns = []
        if matchedlist:
            n_plates=len(self.platelist)
            for i in range(0,n_plates):
                drpf = drpfile(self.drpver, self.platelist[i], self.ifudesignlist[i], 'CUBE')
                if os.path.exists(drpf.file_path()):
                    plates = plates + [self.platelist[i]]
                    ifudesigns = ifudesigns + [self.ifudesignlist[i]]
            return plates, ifudesigns

        for root, dir, files in walk(path):
            for file in files:
                if file.endswith('-LOGCUBE.fits.gz'):
                    p, b, m = parse_drp_file_name(file)

                    ip = 0
                    ib = 0

                    # Include only those files with fibers above a
                    # threshold ifu size
                    if mindesign is not None and b/100 < mindesign:
                        continue

                    if self.platelist is not None:
                        try:
                            ip = self.platelist.index(p)
                        except ValueError:
                            print_frame('ValueError')
                            ip = -1
                    if self.ifudesignlist is not None:
                        try:
                            ib = self.ifudesignlist.index(b)
                        except ValueError:
                            print_frame('ValueError')
                            ib = -1

                    if ip != -1 and ib != -1:
                        plates = plates + [p]
                        ifudesigns = ifudesigns + [b]

        return plates, ifudesigns


    def _find_modes(self, drplist):
        """
        Using the provided list of CUBE DRP files, find if RSS mode is
        also available.

        Currently only two modes are possible:
            modes[i] = 1 -> only 'CUBE' file are available
            modes[i] = 2 -> both 'CUBE' and 'RSS' files are available
        """
        n_drp = len(drplist)
        modes = numpy.empty(n_drp, dtype=numpy.uint8)
        for i in range(0,n_drp):
            drpf = drpfile(drplist[i].drpver, drplist[i].plate, drplist[i].ifudesign, 'RSS')
            modes[i] = 2 if os.path.exists(drpf.file_path()) else 1
        return modes


    def _write_parameter_struct(self, ostream):
        ostream.write('typedef struct {\n')
        ostream.write('    long plate;\n')
        ostream.write('    long ifudesign;\n')
        ostream.write('    char[4] mode;\n')
        ostream.write('    double vel;\n')
        ostream.write('    double vdisp;\n')
        ostream.write('    double ell;\n')
        ostream.write('    double pa;\n')
        ostream.write('    double reff;\n')
        ostream.write('} DAPPAR;\n')


    def _write_parameter_list(self, ostream, index, mode):
        ostream.write('DAPPAR {0:4d} {1:5d} {2:4s} {3:14.7e} {4:14.7e} {5:14.7e} {6:14.7e}'
                      ' {7:14.7e}\n'.format(self.data['PLATE'][index],
                                            self.data['IFUDESIGN'][index], mode,
                                            self.data['VEL'][index], self.data['VDISP'][index],
                                            self.data['ELL'][index], self.data['PA'][index],
                                            self.data['REFF'][index]))


    # ******************************************************************
    #  User functions
    # ******************************************************************

    def file_name(self):
        """Return the name of the DRP complete database."""
        return ('drpcomplete_{0}.fits'.format(self.drpver))


    def file_path(self):
        """Return the full pat to the DRP complete database."""
        return os.path.join(self.analysis_path, self.file_name())


    def update(self, platelist=None, ifudesignlist=None, combinatorics=False, force=False,
               alldrp=False, def_veldisp=100.0, match_r=10.0):
#              , null_veldisp=-999.0):
        """
        Update the DRP complete file.
        
        If platelist, ifudesignlist are provided, the existing (self)
        values are replaced.

        If platelist and ifudesignlist do not have the same length or
        combinatorics is True, the lists are expanded to include all
        combinations of their elements.

        If platelist and ifudesignlist are None or alldrp is True, all
        available plates/ifudesigns are collected from within the drpver
        directory structure.

        If self.file_path() does not exist, it is created.

        If self.file_path() does exist, the available plates and
        ifudesigns in the file are compared against the list (provided
        or collected) to update.  If the lists differ, self.file_path()
        is re-created from scratch.  If all the plates and ifudesigns
        are available, nothing is done, unless force=True.
        """

        if self.readonly:
            raise Exception('drpcomplete file was opened as read-only!')

        if platelist is not None:
            self.platelist = arginp_to_list(platelist, evaluate=True)
        if ifudesignlist is not None:
            self.ifudesignlist = arginp_to_list(ifudesignlist, evaluate=True)

        if alldrp:
            self.platelist = None
            self.ifudesignlist = None

        # TODO: _find_completed_reductions() only searches for CUBE files
        self.platelist, self.ifudesignlist = self._find_completed_reductions(
                                                                    combinatorics=combinatorics)

        # Setup the list of modes such that the combinatorics in
        # drpfile_list is correct.  After running
        # _find_completed_reductions(), the length of platelist and
        # ifudesignlist should be the same!
        n_plates = len(self.platelist)
        modelist = ['CUBE' for i in range(0,n_plates)]
        drplist = drpfile_list(self.drpver, self.platelist, self.ifudesignlist, modelist)

        # By running _all_data_exists() the self.data and self.nrows
        # attributes will be populated if they haven't been already!
        if self._all_data_exists() and not force:
            print('{0} up to date.'.format(self.file_path()))
            return
        else:
            print('Updating {0}.'.format(self.file_path()))

        # Past this point, the drpcomplete file will be ovewritten if it
        # already exists.  Only DRP files created using self.platelist
        # and self.ifudesignlist will be used to create the file, EVEN
        # IF OTHER FILES EXIST.

        modes = self._find_modes(drplist)

        nn = len(drplist)
        print('Number of DRP files: {0}'.format(nn))

        if self.platetargets is not None and all(os.path.exists(p) for p in self.platetargets):
            mangaid, objra, objdec, nsaid, vel, veldisp, ell, pa, Reff \
                    = self._match_platetargets(def_veldisp)

            self.write(self.dapver, self.drpver, self.platetargets, None,
                       [drplist[i].plate for i in range(0,nn)],
                       [drplist[i].ifudesign for i in range(0,nn)], modes, mangaid, nsaid, objra,
                       objdec, vel, veldisp, ell, pa, Reff)

            self._read_data()

        elif self.nsa_cat is not None and os.path.exists(self.nsa_cat):

            mangaid, objra, objdec = self._drp_info(drplist)
            nsaid, vel, veldisp, ell, pa, Reff = self._match_nsa(objra, objdec, def_veldisp,
                                                                 match_r)
            self.write(self.dapver, self.drpver, None, self.nsa_cat,
                       [drplist[i].plate for i in range(0,nn)],
                       [drplist[i].ifudesign for i in range(0,nn)], modes, mangaid, nsaid, objra,
                       objdec, vel, veldisp, ell, pa, Reff)
                       
            self._read_data()

        else:
            raise Exception('Cannot update because catalogs undefined or unavailable.');


    def write(self, dapver, drpver, platetargets, nsa_cat, platelist, ifudesignlist, modes,
              mangaid, nsaid, objra, objdec, vel, veldisp, ell, pa, Reff, clobber=True):
        """Write the drpcomplete file."""

        if self.readonly:
            raise Exception('drpcomplete file was opened as read-only!')

        out=self.file_path()
        if os.path.exists(out) and not clobber:
            raise Exception('DRP complete file already exists: {0}'.format(out))

        # Create the primary header
        hdr = fits.Header()
        hdr['VERSDAP'] = dapver
        hdr['VERSDRP'] = drpver
        if nsa_cat is not None:
            hdr['NSACAT'] = nsa_cat
        if platetargets is not None:
            hdr['PLTARG'] = list_to_csl_string(platetargets)
        hdr['AUTHOR'] = 'K.B. Westfall <kbwestfall@gmail.com>'

        hdu0 = fits.PrimaryHDU(header=hdr)

        # Create the Binary Table
        col1 = fits.Column(name='PLATE', format='1J', array=numpy.array(platelist))
        col2 = fits.Column(name='IFUDESIGN', format='1J', array=numpy.array(ifudesignlist))
        col3 = fits.Column(name='MODES', format='1B', array=modes)
        col4 = fits.Column(name='MANGAID', format='{0}A'.format(mangaid.dtype.itemsize),
                           array=mangaid)
        col5 = fits.Column(name='NSAID', format='1J', array=nsaid)
        col6 = fits.Column(name='OBJRA', format='1D', array=objra)
        col7 = fits.Column(name='OBJDEC', format='1D', array=objdec)
        col8 = fits.Column(name='VEL', format='1D', array=vel)
        col9 = fits.Column(name='VDISP', format='1D', array=veldisp)
        col10 = fits.Column(name='ELL', format='1D', array=ell)
        col11 = fits.Column(name='PA', format='1D', array=pa)
        col12 = fits.Column(name='REFF', format='1D', array=Reff)

        hdu1=fits.BinTableHDU.from_columns([col1, col2, col3, col4, col5, col6, col7, col8, col9,
                                           col10, col11, col12])
        hdulist = fits.HDUList([hdu0, hdu1])
        hdulist.writeto(out, clobber=True)


    def grab_data(self, plate=None, ifudesign=None, index=None, reread=False):

        if (plate is None or ifudesign is None) and index is None:
            raise ValueError('Must provide plate and ifudesign or row index!')

        if index is None:
            index = self.entry_index(plate, ifudesign, reread=reread)
        else:
            self._confirm_access(reread)
            if index >= self.nrows:
                raise ValueError('Selected row index does not exist')

        return [ self.data['PLATE'][index], self.data['IFUDESIGN'][index],
                 self.data['MODES'][index], self.data['MANGAID'][index], self.data['NSAID'][index],
                 self.data['OBJRA'][index], self.data['OBJDEC'][index], self.data['VEL'][index],
                 self.data['VDISP'][index], self.data['ELL'][index], self.data['PA'][index],
                 self.data['REFF'][index] ]


    def write_par(self, ofile, mode, plate=None, ifudesign=None, index=None, reread=False,
                  clobber=True):
        """
        Write the SDSS parameter (Yanny) file for use with the MaNGA DAP.
        """

        if os.path.exists(ofile) and not clobber:
            raise IOError('Parameter file already exists.  Set clobber=True to overwrite.')

        if (plate is None or ifudesign is None) and index is None:
            raise ValueError('Must provide plate and ifudesign or row index!')

        if mode != 'CUBE' and mode != 'RSS':
            raise ValueError('Mode must be either CUBE or RSS.')

        if index is None:
            index = self.entry_index(plate, ifudesign, reread=reread)
        else:
            self._confirm_access(reread)
            if index >= self.nrows:
                raise ValueError('Selected row index does not exist')

        if mode is 'RSS' and self.data['MODES'][index] != 2:
            raise ValueError('RSS mode not available for plate={0},ifudesign={1} !'.format(
                                    self.data['PLATE'][index], self.data['IFUDESIGN'][index]))

        # Write the SDSS parameter file
        ostream = open(ofile, 'w')
        ostream.write('\n')
        self._write_parameter_struct(ostream)
        ostream.write('\n')
        ostream.write('\n')
        self._write_parameter_list(ostream, index, mode)
        ostream.write('\n')
        ostream.close()


    def entry_index(self, plate, ifudesign, reread=False):
        """
        Find the index of the row with the information for the specified
        plate and ifudesign.
        """
        self._confirm_access(reread)

        for i in range(0,self.nrows):
            if self.data['PLATE'][i] == plate and self.data['IFUDESIGN'][i] == ifudesign:
                return i

        raise Exception('Could not find plate={0},ifudesign={1} in drpcomplete file!'.format(plate,
                        ifudesign))



