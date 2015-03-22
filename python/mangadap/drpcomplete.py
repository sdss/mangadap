from __future__ import division
from __future__ import print_function
from __future__ import absolute_import

import sys
if sys.version > '3':
    long = int

import os.path
from os import environ, makedirs, walk
import numpy
from astropy.io import fits
from astropy import constants

# DAP imports
from mangadap.util.yanny import yanny, read_yanny
from mangadap import drpfile
from mangadap import dapfile
from mangadap.util.parser import arginp_to_list, list_to_csl_string
from mangadap.util.exception_tools import print_frame

__author__ = 'Kyle Westfall'


# MOVED to util/parser.py
#def list_to_csl_string(flist):
#    """
#    Convert a list to a comma-separated string.
#    """
#    n = len(flist)
#    out=str(flist[0])
#    for i in range(1,n):
#        out += ', '+str(flist[i])
#    return out



class drpcomplete:
    """
    drpcomplete object with creates, updates, and reads data from the
    drpcomplete file, uses as a reference file for the DAP.

    Provided files have priority over default ones.
    platetargets file has priority over NSA catalog.

    REVISION HISTORY:
        20 Nov 2014: Started implementation by K. Westfall (KBW)
        01 Dec 2014: (KBW) Committed to SVN
        12 Jan 2015: (KBW) Allow for use of plateTargets file using yanny.py
        19 Mar 2015: (KBW) Adjustments for changes to drpfile and
                           drpfile_list.  Now always sets the names of
                           the NSA catalogs and the plateTargets files;
                           added the nsa_catid.  Changed drppath to
                           redux_path.  Added catid and catindx to
                           output file.  Major change to matching when
                           using NSA catalog(s).  Imports full drpfile
                           and dapfile, so now functions called as
                           drpfile.*, dapfile.*; uses default paths,
                           names, etc from there.
    """


    def __init__(self, platelist=None, ifudesignlist=None, platetargets=None, nsa_cat=None,
                 nsa_catid=None, drpver=None, redux_path=None, dapver=None, analysis_path=None,
                 readonly=False):
        """
        Initializes the class and its members.
    
        ARGUMENTS (all optional):
            
            platelist - specified list of plates to analyze
            ifudesignlist - specified list of ifudesign to analyze

            platetargets - plate targets file to use for data collection
            nsa_cat - NSA catalog to use for data collection
            nsa_catid - The ID numbers of the NSA catalog used to match
                        with the MaNGA ID

            drpver - DRP version:
                        - used to define the default DRP redux path
                        - used when declaring a drpfile instance
                        - used in the name of the drpcomplete file
                        - included as a header keyword in the output
                          file
                     If not provided, set to environ['MANGADRP_VER']

            redux_path  - The path to the top level directory containing
                          the DRP output; this is the same as the
                          redux_path in the drpfile class.  If not
                          provided, the default redux_path is:

                              os.path.join(environ['MANGA_SPECTRO_REDUX'],
                                self.drpver)

            dapver - DAP version:
                        - used to define the default DAP analysis path
                        - included as a header keyword in the output
                          drpcomplete file
                     If not provided, set to environ['MANGADAP_VER']

            analysis_path - The path to the top level directory for the
                            DAP output; this is DIFFERENT from the
                            directory_path in the dapfile class.  If not
                            provided, the default analysis_path is:

                              os.path.join(environ['MANGA_SPECTRO_ANALYSIS'],
                                self.dapver)
          
            readonly - Flag that the drpcomplete file is only opened for
                       reading, not for updating.
        """

        # Input properties
        self.drpver = drpfile.default_drp_version() if drpver is None else str(drpver)
        self.redux_path = drpfile.default_redux_path(self.drpver) if redux_path is None \
                                                                  else str(redux_path)

        self.dapver = dapfile.default_dap_version() if dapver is None else str(dapver)
        self.analysis_path = dapfile.default_analysis_path(self.dapver) if analysis_path is None \
                                                                        else str(analysis_path)
        
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
            self.nsa_cat = None
            self.nsa_catid = None
            return

        self.readonly = False

        self.platelist = arginp_to_list(platelist, evaluate=True)
        self.ifudesignlist = arginp_to_list(ifudesignlist, evaluate=True)

        self.platetargets = numpy.array( arginp_to_list(platetargets) ) \
                            if platetargets is not None else self._default_plate_target_files()

        if (nsa_cat is not None and nsa_catid is None) or \
           (nsa_cat is None and nsa_catid is not None):
            raise Exception('Must provide both nsa_cat and nsa_catid')

        if nsa_cat is not None:
            self.nsa_cat = numpy.array( arginp_to_list(nsa_cat) )
            self.nsa_catid = numpy.array( arginp_to_list(nsa_catid) )
        else:
            self.nsa_cat, self.nsa_catid = self._default_nsa_catalogs()


    # ******************************************************************
    #  Default values
    # ******************************************************************

# Use the drpfile versions
#   def _default_drp_version(self):
#       """
#       Return the DRP version defined by the environmental variable.
#       """
#       return environ['MANGADRP_VER']
#
#
#   def _default_redux_path(self):
#       """Return the directory path used by the DRP."""
#
#       # Make sure the DRP version is set
#       if self.drpver is None:
#           self.drpver = self._default_drp_version()
#
#       return os.path.join(environ['MANGA_SPECTRO_REDUX'], self.drpver)


# Use the dapfile versions
#   def _default_dap_version(self):
#       """
#       Return the DRP version defined by the environmental variable.
#       """
#       return environ['MANGADAP_VER']
#
#
#   def _default_analysis_path(self):
#       """Return the directory path used by the DAP."""
#
#       # Make sure the DAP version is set
#       if self.dapver is None:
#           self.dapver = self._default_dap_version()
#
#       return os.path.join(environ['MANGA_SPECTRO_ANALYSIS'], self.dapver)


    def _default_plate_target_files(self):
        """
        Return the default plateTargets file used to get the NSA
        information
        """
        return numpy.array([ os.path.join(environ['MANGACORE_DIR'], 'platedesign', 'platetargets',
                                          'plateTargets-12.par'),
                             os.path.join(environ['MANGACORE_DIR'], 'platedesign', 'platetargets',
                                          'plateTargets-1.par') ])


    def _default_nsa_catalogs(self):
        """
        Return the default NSA catalogs and their respective catalog id
        numbers.

        REVISION HISTORY:
                19 Mar 2015: (KBW) Now returns two lists, one with the
                                   catalog name and the other with the
                                   catalog number.  Should match the
                                   plateTarget files.
        """

        catalog = numpy.array([ os.path.join(environ['MANGAWORK_DIR'], 'manga', 'target', 'temp',
                                             '12-nsa_v1b_0_0_v2.fits.gz'), 
                                os.path.join(environ['MANGA_TARGET'], 'input', 'nsa_v1_0_0.fits') ])
        catalogid = numpy.array([ 12, 1 ])
        return catalog, catalogid


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
        objra = numpy.empty(nn, dtype=numpy.float64)
        objdec = numpy.empty(nn, dtype=numpy.float64)
        print("Gathering DRP header data...")
        for i in range(0,nn):
            mangaid_, objra[i], objdec[i] = drplist[i].object_data()
            mangaid = mangaid + [mangaid_]
        print("... done")
        mangaid = numpy.array(mangaid)
        return mangaid, objra, objdec


    def _read_platetargets(self, verbose=False):
        """
        Read all platetargets files and return a single set of data.
        """

        # Check the files exist
        if not all([ os.path.exists(p) for p in self.platetargets ]):
            raise Exception('Cannot open one or more plateTargets files!')

        nn = len(self.platetargets)             # Number of files
        par_data = yanny(self.platetargets[0])  # Read the first file
        print(self.platetargets[0])

        # Append the rest
        for i in range(1,nn):
            tmp_par = read_yanny(self.platetargets[i])
            print(self.platetargets[i])
            par_data.append(tmp_par, add_to_file=False)

        par_data.convert_to_numpy_array()       # Convert to a numpy array
        return par_data


    def _read_nsa_idonly(self, verbose=False):
        """
        Read the NSA catalogs.

        Returned quantites are:
            nsa_id - a list of numpy.ndarrays with the NSA IDs
            nsa_ndata - a numpy.array with the number of elements in
                        each nsa_id element
        
        """
        # Check the files exist
        if not all([ os.path.exists(n) for n in self.nsa_cat ]):
            raise Exception('Cannot open one or more NSA catalogs!')

        nsa_id = []
        nsa_ndata = []
        for n in self.nsa_cat:
            print(n)
            hdu = fits.open(n)

            nsa_ndata = nsa_ndata + [ hdu[1].header['NAXIS2'] ]
            nsa_id = nsa_id + [ numpy.array(hdu[1].data['NSAID']) ]
            hdu.close()

        return nsa_id, numpy.array(nsa_ndata)


    def _read_nsa(self, verbose=False):
        """Read all the necessary data from the NSA catalogs."""
        # Check the files exist
        if not all([ os.path.exists(n) for n in self.nsa_cat ]):
            raise Exception('Cannot open one or more NSA catalogs!')

        nsa_id = []
        nsa_vel = []
        nsa_ell = []
        nsa_pa = []
        nsa_Reff = []
        nsa_ndata = []
        nsa_veldisp = []
        nsa_has_veldisp = []
        for n in self.nsa_cat:
            print(n)
            hdu = fits.open(n)

            nsa_ndata = nsa_ndata + [ hdu[1].header['NAXIS2'] ]
            nsa_id = nsa_id + [ numpy.array(hdu[1].data['NSAID']) ]
            nsa_vel = nsa_vel + [ numpy.array(hdu[1].data['Z']*constants.c.to('km/s').value) ]
            nsa_ell = nsa_ell + [ numpy.array(1.0-hdu[1].data['SERSIC_BA']) ]
            nsa_pa = nsa_pa + [ numpy.array(hdu[1].data['SERSIC_PHI']) ]
            nsa_Reff = nsa_Reff + [ numpy.array(hdu[1].data['SERSIC_TH50']) ]

            # Try to get the column number with the velocity dispersion
            # to test if the NSA catalog has it.
            try:
                i = hdu[1].columns.names.index('VELDISP')
            except ValueError:
                if verbose:
                    print('WARNING: {0} does not have vel. dispersion.  Using default.'.format(n))
                nsa_has_veldisp = nsa_has_veldisp + [ False ]
                nsa_veldisp = nsa_veldisp + [ numpy.empty(0) ]
            else:
                nsa_has_veldisp = nsa_has_veldisp + [ True ]
                nsa_veldisp = nsa_veldisp + [ numpy.array(hdu[1].data['VELDISP']) ]

            hdu.close()

        return nsa_id, nsa_vel, nsa_veldisp, nsa_ell, nsa_pa, nsa_Reff, numpy.array(nsa_ndata), \
               numpy.array(nsa_has_veldisp)


    def _match_platetargets(self, def_veldisp):
        """
        Read the platetargets file, and match the completed DRP
        reductions based on the plate and ifudesign numbers.  The NSA
        catalogs are also read so that the ID can be pulled from the
        appropriate catalog.
        """

        # Read the platetargets file
        print('Reading platetarget file(s)...')
        par_data = self._read_platetargets(verbose=True)
        print('DONE')

        # Read the NSA catalogs
        print('Reading NSA catalog(s)...')
        nsa_id, nsa_ndata = self._read_nsa_idonly(verbose=True)
        print('DONE')

        # Initialize the output arrays (needed in case some DRP targets not found)
        n_drp = len(self.platelist)
#       print(par_data['PLTTRGT']['mangaid'].dtype)
        mangaid = []
#       mangaid = numpy.empty(n_drp, dtype=par_data['PLTTRGT']['mangaid'].dtype)
        objra = numpy.zeros(n_drp, dtype=par_data['PLTTRGT']['target_ra'].dtype)
        objdec = numpy.zeros(n_drp, dtype=par_data['PLTTRGT']['target_dec'].dtype)
        catid = numpy.full(n_drp, -1, dtype=numpy.int32)
        catindx = numpy.full(n_drp, -1, dtype=numpy.int32)
        nsaid = numpy.full(n_drp, -1, dtype=numpy.int32)
        vel = numpy.zeros(n_drp, dtype=par_data['PLTTRGT']['nsa_redshift'].dtype)
        veldisp = numpy.full(n_drp, def_veldisp, dtype=par_data['PLTTRGT']['nsa_vdisp'].dtype)
        ell = numpy.zeros(n_drp, dtype=par_data['PLTTRGT']['nsa_sersic_ba'].dtype)
        pa = numpy.zeros(n_drp, dtype=par_data['PLTTRGT']['nsa_sersic_phi'].dtype)
        Reff = numpy.zeros(n_drp, dtype=par_data['PLTTRGT']['nsa_sersic_th50'].dtype)

#       for t in par_data['PLTTRGT']['plateid']:
#           if type(t) != numpy.int64:
#               print(type(t))

        print('Searching platetargets file for observed galaxies...')
        for i in range(0,n_drp):
            indx = numpy.where((par_data['PLTTRGT']['plateid'] == self.platelist[i]) &
                               (par_data['PLTTRGT']['ifudesign'] == self.ifudesignlist[i]))

            if len(indx[0]) == 0:
#               print(type(par_data['PLTTRGT']['plateid'][0]))
#               print(type(self.platelist[i]))
#               print(type(par_data['PLTTRGT']['ifudesign'][0]))
#               print(type(self.ifudesignlist[i]))
                print('WARNING: Could not find plate={0}, ifudesign={1} in {2}.'.format(
                                self.platelist[i], self.ifudesignlist[i], self.platetargets))
                mangaid = mangaid + ['NULL']
                continue
#               raise Exception('Could not find plate={0}, ifudesign={1} in {2}.'.format(
#                               self.platelist[i], self.ifudesignlist[i], self.platetargets))
            mangaid = mangaid + [par_data['PLTTRGT']['mangaid'][indx][0].decode("ascii")]
            objra[i] = par_data['PLTTRGT']['target_ra'][indx][0]
            objdec[i] = par_data['PLTTRGT']['target_dec'][indx][0]

            # From David Wake:
            #
            # "MANGAID consists of CATID-CATIND, where CATID identifies
            # a parent catalog, in the case of the main manga samples
            # CATID = 1 which refers to nsa_v1_0_0.fits, and CATIND is
            # the position within that catalog (zero indexed). So if you
            # strip out CATIND from MANGAID you have the position within
            # nsa_v1_0_0.fits without any matching required."
            catid[i] = numpy.int32(mangaid[i].split('-')[0])
            catindx[i] = numpy.int32(mangaid[i].split('-')[1])

            nsa_indx = numpy.where(self.nsa_catid == catid[i])
            if len(nsa_indx[0]) == 0:
                print('WARNING: No NSA catalog {0} available.  Setting NSAID=-1.'.format(catid[i]))
                nsaid[i] = -1
            elif catindx[i] >= nsa_ndata[nsa_indx]:
                print('WARNING: No index {0} in NSA catalog {1}.  Setting NSAID=-1.'.format(
                        catindx[i], catid[i]))
                nsaid[i] = -1
            else:
                nsaid[i] = numpy.int32(nsa_id[nsa_indx[0]][catindx[i]])

            vel[i] = par_data['PLTTRGT']['nsa_redshift'][indx][0] * constants.c.to('km/s').value
            ell[i] = 1.0-par_data['PLTTRGT']['nsa_sersic_ba'][indx][0]
            pa[i] = par_data['PLTTRGT']['nsa_sersic_phi'][indx][0]
            Reff[i] = par_data['PLTTRGT']['nsa_sersic_th50'][indx][0]

            if par_data['PLTTRGT']['nsa_vdisp'][indx][0] > 0.:
                veldisp[i] = par_data['PLTTRGT']['nsa_vdisp'][indx][0]

        print('DONE')

        return numpy.array(mangaid), objra, objdec, catid, catindx, nsaid, vel, veldisp, ell, \
               pa, Reff


    def _match_nsa(self, mangaid, def_veldisp):
        """
        Read the NSA catalogs and extract the necessary information for
        the DRP files.
        """

        # Check the number of MaNGA IDs provided.
        n_drp = len(self.platelist)
        if len(mangaid) != n_drp:
            raise Exception('Incorrect number of MaNGA IDs')

        # Read the NSA catalogs
        print('Reading NSA catalog(s)...')
        nsa_id, nsa_vel, nsa_veldisp, nsa_ell, nsa_pa, nsa_Reff, nsa_ndata, nsa_has_veldisp = \
                self._read_nsa(verbose=True)
        print('DONE')

        # Initialize the output arrays (needed in case some DRP targets not found)
        catid = numpy.full(n_drp, -1, dtype=numpy.int32)
        catindx = numpy.full(n_drp, -1, dtype=numpy.int32)
        nsaid = numpy.full(n_drp, -1, dtype=numpy.int32)
        vel = numpy.zeros(n_drp, dtype=numpy.float64)
        veldisp = numpy.full(n_drp, def_veldisp, dtype=numpy.float64)
        ell = numpy.zeros(n_drp, dtype=numpy.float64)
        pa = numpy.zeros(n_drp, dtype=numpy.float64)
        Reff = numpy.zeros(n_drp, dtype=numpy.float64)

        print('Using MaNGA ID to get information from NSA catalog(s)...')
        for i in range(0,n_drp):
           
            catid[i] = numpy.int32(mangaid[i].split('-')[0])
            catindx[i] = numpy.int32(mangaid[i].split('-')[1])

            indx = numpy.where(self.nsa_catid == catid[i])
            if len(indx[0]) == 0:
                print('WARNING: No NSA catalog {0} available.  Setting NSAID=-1.'.format(catid[i]))
                nsaid[i] = -1
            elif catindx[i] >= nsa_ndata[indx]:
                print('WARNING: No index {0} in NSA catalog {1}.  Setting NSAID=-1.'.format(
                        catindx[i], catid[i]))
                nsaid[i] = -1
            else:
                nsaid[i] = numpy.int32(nsa_id[indx[0]][catindx[i]])
                vel[i] = numpy.float64(nsa_vel[indx[0]][catindx[i]])
                ell[i] = numpy.float64(nsa_ell[indx[0]][catindx[i]])
                pa[i] = numpy.float64(nsa_pa[indx[0]][catindx[i]])
                Reff[i] = numpy.float64(nsa_Reff[indx[0]][catindx[i]])
                if nsa_has_veldisp[indx[0]]:
                    veldisp[i] = numpy.float64(nsa_veldisp[indx[0]][catindx[i]])

        return catid, catindx, nsaid, vel, veldisp, ell, pa, Reff


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
        hdu = fits.open(self.file_path())
        self.data = hdu[1].data
        self.nrows = hdu[1].header['NAXIS2']
        hdu.close()
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
                drpf = drpfile.drpfile(self.platelist[i], self.ifudesignlist[i], 'CUBE',
                                       drpver=self.drpver)
                if os.path.exists(drpf.file_path()):
                    plates = plates + [self.platelist[i]]
                    ifudesigns = ifudesigns + [self.ifudesignlist[i]]
            return numpy.array(plates), numpy.array(ifudesigns)

        for root, dir, files in walk(path):
            for file in files:
                if file.endswith('-LOGCUBE.fits.gz'):
                    p, b, m = drpfile.parse_drp_file_name(file)

#                   print('{0}, {1}, {2}'.format(p,b,m))

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
                            #print_frame('ValueError')
                            ip = -1
                    if self.ifudesignlist is not None:
                        try:
                            ib = self.ifudesignlist.index(b)
                        except ValueError:
                            #print_frame('ValueError')
                            ib = -1

                    if ip != -1 and ib != -1:
                        plates = plates + [p]
                        ifudesigns = ifudesigns + [b]

#       print(plates)
#       print(ifudesigns)
        return plates, ifudesigns #numpy.array(plates), numpy.array(ifudesigns)


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
            drpf = drpfile.drpfile(drplist[i].plate, drplist[i].ifudesign, 'RSS',
                                   drpver=drplist[i].drpver)
            modes[i] = 2 if os.path.exists(drpf.file_path()) else 1
        return modes


    def _write_parameter_struct(self, ostream):
        ostream.write('typedef struct {\n')
        ostream.write('    long plate;\n')
        ostream.write('    long ifudesign;\n')
        ostream.write('    char mode[4];\n')
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
               alldrp=False, def_veldisp=100.0, use_nsa=False):
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

        OPTIONAL:
                use_nsa bool
                        Force the matching to be done using the NSA
                        catalog(s).  Default is to use the platetargets
                        files.
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
        drplist = drpfile.drpfile_list(self.platelist, self.ifudesignlist, modelist,
                                       drpver=self.drpver)

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

        # Default to using the plateTargets files
        if use_nsa:
            mangaid, objra, objdec = self._drp_info(drplist)
            catid, catindx, nsaid, vel, veldisp, ell, pa, Reff = self._match_nsa(mangaid,
                                                                                 def_veldisp)
        else:
            mangaid, objra, objdec, catid, catindx, nsaid, vel, veldisp, ell, pa, Reff \
                    = self._match_platetargets(def_veldisp)

        self.write([drplist[i].plate for i in range(0,nn)],
                   [drplist[i].ifudesign for i in range(0,nn)], modes, mangaid, objra, objdec,
                   catid, catindx, nsaid, vel, veldisp, ell, pa, Reff)

        self._read_data()


    def write(self, platelist, ifudesignlist, modes, mangaid, objra, objdec, catid, catindx, nsaid,
              vel, veldisp, ell, pa, Reff, drpver=None, redux_path=None, dapver=None,
              analysis_path=None, platetargets=None, nsa_cat=None, nsa_catid=None, clobber=True):
        """
        Write the drpcomplete file.
        
        REVISION HISTORY:
                19 Mar 2015: (KBW) Changed order of input and allowed
                                   for use of self components as
                                   default.
        """

        if self.readonly:
            raise Exception('drpcomplete file was opened as read-only!')

        out=self.file_path()
        if os.path.exists(out) and not clobber:
            raise Exception('DRP complete file already exists: {0}'.format(out))

        # Create the primary header
        hdr = fits.Header()
        hdr['VERSDRP'] = self.drpver if drpver is None else drpver
        hdr['RDXPTH'] = self.redux_path if redux_path is None else redux_path
        hdr['VERSDAP'] = self.dapver if dapver is None else dapver
        hdr['SISPTH'] = self.analysis_path if analysis_path is None else analysis_path
        if nsa_cat is None:
            nsa_cat = self.nsa_cat
        if nsa_catid is None:
            nsa_catid = self.nsa_catid
        for i in range(0,nsa_cat.size):
            hdr['NSACAT{0}'.format(i+1)] = nsa_cat[i]
            hdr['NSACID{0}'.format(i+1)] = nsa_catid[i]
        if platetargets is None:
            platetargets = self.platetargets
        for i in range(0,platetargets.size):
            hdr['PLTARG{0}'.format(i+1)] = platetargets[i]
        hdr['AUTHOR'] = 'K.B. Westfall <kbwestfall@gmail.com>'

        hdu0 = fits.PrimaryHDU(header=hdr)

        # Create the Binary Table
        col1 = fits.Column(name='PLATE', format='1J', array=numpy.array(platelist))
        col2 = fits.Column(name='IFUDESIGN', format='1J', array=numpy.array(ifudesignlist))
        col3 = fits.Column(name='MODES', format='1B', array=modes)
        col4 = fits.Column(name='MANGAID', format='{0}A'.format(mangaid.dtype.itemsize),
                           array=mangaid)
        col5 = fits.Column(name='OBJRA', format='1D', array=objra)
        col6 = fits.Column(name='OBJDEC', format='1D', array=objdec)
        col7 = fits.Column(name='CATID', format='1J', array=catid)
        col8 = fits.Column(name='CATINDX', format='1J', array=catindx)
        col9 = fits.Column(name='NSAID', format='1J', array=nsaid)
        col10 = fits.Column(name='VEL', format='1D', array=vel)
        col11 = fits.Column(name='VDISP', format='1D', array=veldisp)
        col12 = fits.Column(name='ELL', format='1D', array=ell)
        col13 = fits.Column(name='PA', format='1D', array=pa)
        col14 = fits.Column(name='REFF', format='1D', array=Reff)

        hdu1=fits.BinTableHDU.from_columns([col1, col2, col3, col4, col5, col6, col7, col8, col9,
                                           col10, col11, col12, col13, col14])
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
                 self.data['MODES'][index], self.data['MANGAID'][index], self.data['OBJRA'][index],
                 self.data['OBJDEC'][index], self.data['CATID'][index], self.data['CATINDX'][index],
                 self.data['NSAID'][index], self.data['VEL'][index], self.data['VDISP'][index],
                 self.data['ELL'][index], self.data['PA'][index], self.data['REFF'][index] ]


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



