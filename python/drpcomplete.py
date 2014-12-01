from __future__ import division
from __future__ import print_function

import os.path
from os import environ, makedirs, walk
import numpy

from drpfile import drpfile, drpfile_list, parse_drp_file_name, strcsv_to_list
from astropy.io import fits
from astropy import constants
from scipy.spatial import KDTree

__author__ = 'Kyle Westfall'


# TODO: Let the paths be defined, using default otherwise

class drpcomplete:
    """
    drpcomplete object with creates, updates, and reads data from the
    drpcomplete file, uses as a reference file for the DAP.

    REVISION HISTORY:
        20 Nov 2014: (KBW) Started
    """


    def __init__(self, dapver, drpver, platelist=None, ifudesignlist=None, nsa_cat=None,
                 force=False):
        """
        Initializes the class and its members.
    
        ARGUMENTS:
            
            drpver - DRP version to use
            platelist - specified list of plates to analyze
            ifudesignlist - specified list of ifudesign to analyze
            nsa_cat - NSA catalog to use for data collection
            force - Force updates, even if the data already exists
        """

        # Input properties
        if (type(dapver) != str):
            raise Exception('DAP version must be a string')
        self.dapver = dapver

        if (type(drpver) != str):
            raise Exception('DRP version must be a string')
        self.drpver = drpver

        self.platelist = platelist
        if self.platelist is not None and type(self.platelist) == str:
            self.platelist = strcsv_to_list(self.platelist, evaluate=True)
        if self.platelist is not None and type(self.platelist) != list:    
            self.platelist = [self.platelist]

        self.ifudesignlist = ifudesignlist
        if self.ifudesignlist is not None and type(self.ifudesignlist) == str:
            self.ifudesignlist = strcsv_to_list(self.ifudesignlist, evaluate=True)
        if self.ifudesignlist is not None and type(self.ifudesignlist) != list:    
            self.ifudesignlist = [self.ifudesignlist]

        self.nsa_cat = self.default_nsa_catalog() if nsa_cat is None else nsa_cat
        self.force = force

        if os.path.exists(self.file_path()):
            self._read_data()
        else:
            self.data = None
            self.nrows = None


    # ******************************************************************
    #  Utility functions
    # ******************************************************************


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


    def _read_nsa(self, verbose=False):
        """Read the NSA catalog and the number of rows."""
        if not os.path.exists(self.nsa_cat):
            raise Exception(''.join('Cannot open NSA catalog: ',self.nsa_cat))
        hdulist = fits.open(self.nsa_cat)
        nsa_data = hdulist[1].data
        n_nsa = hdulist[1].header['NAXIS2']

        try:
            i = hdulist[1].columns.names.index('VELDISP')
        except ValueError, e:
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

        # Set the coordinates to a KDTree
        kdcoo = KDTree(zip( (nsa_data['RA']).ravel(), (nsa_data['DEC']).ravel() ))
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

        # Brute force approach
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
            dd = ((kdcoo.data[nneigh_i,0]-objra[i])*cosd)**2.0 + (kdcoo.data[nneigh_i,1]-objdec[i])**2.0
            if len(numpy.where(dd < match_d)) > 1:
                raise Exception(''.join(['Multiple matches in the NSA catalog: ra=',
                                str(objra[i]), ',dec=', str(objdec[i]), '!']))

            di = numpy.argsort(dd)
            #print(dd, di)
            if dd[di[0]] < match_d:
                nsa_index[i] = nneigh_i[di[0]]
            else:
                print(''.join(['Could not find match in NSA catalog for plate=',
                      str(self.platelist[i]), ',ifudesign=', str(self.ifudesignlist[i]), '!']))
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

    def _all_data_exists(self):
        """
        Determine if the data for all the plates/ifudesigns is already
        present in in the current drpcomplete file.
        """
        nn = len(self.platelist)
        for i in range(0,nn):
            try:
                index = self.entry_index(self.platelist[i], self.ifudesignlist[i])
            except Exception, e:
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


#    def _check_pd(self, plate, ifudesign):
#        """Check plate and ifudesign are integers."""
#        if type(plate) != int or type(ifudesign) != int:
#            raise Exception("plate and ifudesign must be integers to access data in drpcomplete")


    def _check_access(self, reread):
        """Check drpcomplete file is accessible."""
        inp = self.file_path()
        if not os.path.exists(inp):
            raise Exception(''.join(['Cannot open file: ', inp]))

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
                one this not None are used to create the platelist and
                ifudesignlist.  E.g., if platelist=[7443] and
                ifudesignlist is None, all CUBE files with platelist =
                7443 are chosen.
            If both are not None
        
        
        """
        path = self.redux_path()
        matchedlist = (self.platelist is not None and self.ifudesignlist is not None and
                       (len(self.platelist) == len(self.ifudesignlist) and not combinatorics))
        plates = []
        ifudesigns = []
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
                        except ValueError, e:
                            ip = -1
                    if self.ifudesignlist is not None:
                        try:
                            ib = self.ifudesignlist.index(b)
                        except ValueError, e:
                            ib = -1

                    # Only accept the indices for a matched list if they
                    # are the same
                    if matchedlist and ip != ib:
                        ip = -1
       
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


    def default_nsa_catalog(self):
        return os.path.join(environ['MANGA_TARGET'],'input','nsa_v1_0_0.fits')


    def file_name(self):
        """Return the name of the DRP complete database."""
        return ''.join(['drpcomplete_',self.drpver,'.fits'])


    def redux_path(self):
        """Return the path to search for completed DRP reductions."""
        return os.path.join(environ['MANGA_SPECTRO_REDUX'], self.drpver)


    def file_path(self):
        """Return the full pat to the DRP complete database."""
        return os.path.join(environ['MANGA_SPECTRO_ANALYSIS'], self.dapver, self.file_name())


    def update(self, platelist=None, ifudesignlist=None, combinatorics=False, force=False,
               def_veldisp=100.0, match_r=10.0):
        """
        Update the DRP complete file.
        
        If platelist, ifudesignlist, or force are provided, the existing
        (self) values are replaced.

        If platelist and ifudesignlist are None, all available
        plates/ifudesigns are collected from within the drpver directory
        structure.

        If self.file_path() does not exist, it is created.

        If self.file_path() does exist, the available plates and
        ifudesigns in the file are compared against the list (provided
        or collected) to update.  If the lists differ, self.file_path()
        is re-created from scratch.  If all the plates and ifudesigns
        are available, nothing is done, unless force=True.
        """

        if platelist is not None:
            self.platelist = platelist
        if self.platelist is not None and type(self.platelist) == str:
            self.platelist = strcsv_to_list(self.platelist, evaluate=True)
        if self.platelist is not None and type(self.platelist) != list:    
            self.platelist = [self.platelist]

        if ifudesignlist is not None:
            self.ifudesignlist = ifudesignlist
        if self.ifudesignlist is not None and type(self.ifudesignlist) == str:
            self.ifudeisnglist = strcsv_to_list(self.ifudesignlist, evaluate=True)
        if self.ifudesignlist is not None and type(self.ifudesignlist) != list:    
            self.ifudesignlist = [self.ifudesignlist]

        # TODO: _find_completed_reductions() only searches for CUBE files
        self.platelist, self.ifudesignlist = \
            self._find_completed_reductions(combinatorics=combinatorics)

        # TODO: Warn the user if one of their plates/ifudesigns have no
        # completed reductions?

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
            print(''.join([self.file_path(), ' up to date.']))
            return

        # Past this point, the drpcomplete file will be ovewritten if it
        # already exists.  Only DRP files created using self.platelist
        # and self.ifudesignlist will be used to create the file, EVEN
        # IF OTHER FILES EXIST.

        modes = self._find_modes(drplist)

        nn = len(drplist)
        print(''.join(['Number of DRP files: ',str(nn)]))

        mangaid, objra, objdec = self._drp_info(drplist)
        nsaid, vel, veldisp, ell, pa, Reff = self._match_nsa(objra, objdec, def_veldisp, match_r)
        self.write(self.dapver, self.drpver, self.nsa_cat, [drplist[i].plate for i in range(0,nn)],
                   [drplist[i].ifudesign for i in range(0,nn)], modes, mangaid, nsaid, objra,
                   objdec, vel, veldisp, ell, pa, Reff)
        self._read_data()


    def write(self, dapver, drpver, nsa_cat, platelist, ifudesignlist, modes, mangaid, nsaid, objra,
              objdec, vel, veldisp, ell, pa, Reff, clobber=True):
        """Write the drpcomplete file."""

        out=self.file_path()
        if os.path.exists(out) and not clobber:
            raise Exception(''.join(['DRP complete file already exists: ',out]))

        # TODO: Check array types
        # platelist = list
        # ifudesignlist = list
        # modes = numpy.ndarray
        # mangaid = numpy.ndarray
        # nsaid = numpy.ndarray
        # objra = numpy.ndarray
        # objdec = numpy.ndarray
        # vel = numpy.ndarray
        # veldisp = numpy.ndarray
        # ell = numpy.ndarray
        # pa = numpy.ndarray
        # Reff = numpy.ndarray

        # Create the primary header
        hdr = fits.Header()
        hdr['VERSDAP'] = dapver
        hdr['VERSDRP'] = drpver
        hdr['NSACAT'] = nsa_cat
        hdr['AUTHOR'] = 'K.B. Westfall <kbwestfall@gmail.com>'

        hdu0 = fits.PrimaryHDU(header=hdr)

        # Create the Binary Table
        col1 = fits.Column(name='PLATE', format='1J', array=numpy.array(platelist))
        col2 = fits.Column(name='IFUDESIGN', format='1J', array=numpy.array(ifudesignlist))
        col3 = fits.Column(name='MODES', format='1B', array=modes)
        col4 = fits.Column(name='MANGAID', format=''.join([str(mangaid.dtype.itemsize),'A']),
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
            self._check_access(reread)
            if index >= self.nrows:
                raise ValueError('Selected row index does not exist')

        return [ self.data['PLATE'][index], self.data['IFUDESIGN'][index],
                 self.data['MODES'][index], self.data['MANGAID'][index], self.data['NSAID'][index],
                 self.data['OBJRA'][index], self.data['OBJDEC'][index], self.data['VEL'][index],
                 self.data['VDISP'][index], self.data['ELL'][index], self.data['PA'][index],
                 self.data['REFF'][index] ]


    def write_par(self, ofile, mode, plate=None, ifudesign=None, index=None, reread=False,
                  force=True):

        if os.path.exists(ofile) and not force:
            raise IOError('Parameter file already exists.  Set force=True to overwrite.')

        if (plate is None or ifudesign is None) and index is None:
            raise ValueError('Must provide plate and ifudesign or row index!')

        if mode != 'CUBE' and mode != 'RSS':
            raise ValueError('Mode must be either CUBE or RSS.')

        if index is None:
            index = self.entry_index(plate, ifudesign, reread=reread)
        else:
            self._check_access(reread)
            if index >= self.nrows:
                raise ValueError('Selected row index does not exist')

        if mode is 'RSS' and self.data['MODES'][index] != 2:
            raise ValueError(''.join(['RSS mode not available for plate=',
                                     str(self.data['PLATE'][index]), ',ifudesign=',
                                     str(self.data['IFUDESIGN'][index]), '!']))

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
#        self._check_pd(plate, ifudesign)
        self._check_access(reread)

        for i in range(0,self.nrows):
            if self.data['PLATE'][i] == plate and self.data['IFUDESIGN'][i] == ifudesign:
                return i

        raise Exception('Could not find plate={0},ifudesign={1} in drpcomplete file!'.format(plate,
                        ifudesign))



