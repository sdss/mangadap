"""

Defines a class used to produce and provide an interface for data
required to create the input `SDSS parameter files`_ for the DAP.

The drpcomplete file/data is primarily requred for the survey-level
execution of the `DAP at Utah`_.  However, it does provide useful
information regarding the completed DRP files, and it can be used to
create DAP par files.

Until further documented here, the description of the `DAP par
file`_ is available in the SDSS-IV/MaNGA `MPL-3 Technical Reference Manual`_.

*Source location*:
    $MANGADAP_DIR/python/mangadap/drpcomplete.py

*Python2/3 compliance*::

    from __future__ import division
    from __future__ import print_function
    from __future__ import absolute_import
    
    import sys
    if sys.version > '3':
        long = int

*Imports*::

    import os.path
    from os import environ, makedirs, walk
    import numpy
    from astropy.io import fits
    from astropy import constants
    from mangadap.util.yanny import yanny, read_yanny
    from mangadap.util import defaults
    from mangadap import drpfile
    from mangadap import dapfile
    from mangadap.util.parser import arginp_to_list, list_to_csl_string
    from mangadap.util.exception_tools import print_frame

*Class usage examples*:

    .. todo::

        Add some usage comments here!

*Revision history*:
    | **20 Nov 2014**: Started implementation by K. Westfall (KBW)
    | **01 Dec 2014**: (KBW) Committed to SVN
    | **12 Jan 2015**: (KBW) Allow for use of plateTargets file using yanny.py
    | **19 Mar 2015**: (KBW) Adjustments for changes to drpfile and
        drpfile_list.  Now always sets the names of the NSA catalogs and
        the plateTargets files; added the nsa_catid.  Changed drppath to
        redux_path.  Added catid and catindx to output file.  Major
        change to matching when using NSA catalog(s).  Imports full
        :class:`mangadap.drpfile` and :class:`mangadap.dapfile` classes,
        require a change to the function calls. Changed parameter order
        in :func:`write` and allowed for use of attributes as default.
    | **20 May 2015**: (KBW) Documentation and Sphinx tests. Prep for
        conversion from NSA to TRG nomenclatature, but still need to
        make the full conversion.

.. todo::
    - Need accommodate more than just the CATID=1,12 (NSA) catalogs
    - Move default plateTargets and catalog file names to defaults.py

.. _DAP par file: https://trac.sdss.org/wiki/MANGA/TRM/TRM_MPL-3/dap/Summary/parFile
.. _MPL-3 Technical Reference Manual: https://trac.sdss.org/wiki/MANGA/TRM/TRM_MPL-3
.. _SDSS parameter files: http://www.sdss.org/dr12/software/par/
.. _DAP at Utah: https://trac.sdss.org/wiki/MANGA/TRM/TRM_MPL-3/dap/Summary/Utah

"""

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
from mangadap.util.yanny import yanny, read_yanny
from mangadap.util import defaults
from mangadap import drpfile
from mangadap import dapfile
from mangadap.util.parser import arginp_to_list, list_to_csl_string, parse_drp_file_name
from mangadap.util.exception_tools import print_frame

__author__ = 'Kyle Westfall'

class drpcomplete:
    """
    Class is used to determine which DRP files are ready for analysis
    with the DAP, and collect the data necessary for and write the DAP
    input parameter file.

    Provided files have priority over default ones.
    platetargets file has priority over NSA catalog.

    Args:
        platelist (str or list): (Optional) List of plates to search
            for.  Default is to search the full DRP path.
        ifudesignlist (str or list): (Optional) List of ifudesign to
            search for.  Default is to search the full DRP path.
        platetargets (str or list): (Optional) List of platetargets
            files to search through to find any given plate ifudesign
            combination. Default is::
            
                numpy.array([ os.path.join(environ['MANGACORE_DIR'],
                                           'platedesign',
                                           'platetargets',
                                           'plateTargets-12.par'),
                              os.path.join(environ['MANGACORE_DIR'],
                                           'platedesign',
                                           'platetargets',
                                           'plateTargets-1.par') ])

        trg_cat (str or list): (Optional) List of catalogs to search
            through to grab the target ID from that catalog. Default
            is::
            
                numpy.array([ os.path.join(environ['MANGA_TARGET'],
                                           'catalogues',
                                           '12-nsa_v1b_0_0_v2.fits'), 
                              os.path.join(environ['MANGA_TARGET'],
                                           'catalogues',
                                           '1-nsa_v1_0_0.fits') ])
            
        trg_catid (str or list): (Optional) List of target catalog ID
            numbers to use.  Default is::
            
                numpy.array([ 12, 1 ])

        drpver (str): (Optional) DRP version, which is:
                - used to define the default DRP redux path
                - used when declaring a drpfile instance
                - used in the name of the drpcomplete file
                - included as a header keyword in the output file

            Default is defined by :func:`mangadap.util.defaults.default_drp_version`
        redux_path (str): (Optional) The path to the top level directory
            containing the DRP output files; this is the same as the
            *redux_path* in the :class:`mangadap.drpfile` class.
            Default is defined by
            :func:`mangadap.util.defaults.default_redux_path`.

        dapver (str): (Optional) DAP version, which is:
                - used to define the default DAP analysis path
                - included as a header keyword in the output drpcomplete
                  file

            Default is defined by
            :func:`mangadap.util.defaults.default_dap_version`
        analysis_path (str): The path to the top level directory for the
            DAP output files; this is **different** from the
            directory_path in the :class:`mangadap.dapfile` class.  Default is
            defined by
            :func:`mangadap.util.defaults.default_analysis_path`
        readonly (bool) : Flag that the drpcomplete fits file is only
            opened for reading, not for updating.

    Raises:
        Exception: Raised if user supplies only one of trg_cat or
            trg_catid instead of both.

    Attributes:
        platelist (list) : List of plates to search for, see above.
        ifudesignlist (list) : List of IFU designs to search for, see
            above.
        platetargets (numpy.ndarray, dtype=str) : List of platetargets
            files to search through, see above
        trg_cat (numpy.ndarray, dtype=str) : List of target catalogs to
            search through, see above
        trg_catid (numpy.ndarray, dtype=str) : List of target catalog
            IDs, see above
        drpver (str) : DRP version, see above.
        redux_path (str) : Path to the top level of directory for the
            DRP output, see above.
        dapver (str) : DAP version, see above.
        analysis_path (str) : Path to the top level directory for the
            DAP output files, see above.
        readonly (bool) : Flag that the drpcomplete fits file is only
            opened for reading, not for updating.
        data (dict) : Data storage structure.  *Need to check exact data
            type*.
        nrows (int) : Number of rows in the data structure.

    """
    def __init__(self, platelist=None, ifudesignlist=None, platetargets=None, trg_cat=None,
                 trg_catid=None, drpver=None, redux_path=None, dapver=None, analysis_path=None,
                 readonly=False):
        # Input properties
        self.drpver = defaults.default_drp_version() if drpver is None else str(drpver)
        self.redux_path = defaults.default_redux_path(self.drpver) if redux_path is None \
                                                                  else str(redux_path)

        self.dapver = defaults.default_dap_version() if dapver is None else str(dapver)
        self.analysis_path = defaults.default_analysis_path(self.drpver, self.dapver) \
                             if analysis_path is None else str(analysis_path)
        
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
            self.trg_cat = None
            self.trg_catid = None
            return

        self.readonly = False

        self.platelist = arginp_to_list(platelist, evaluate=True)
        self.ifudesignlist = arginp_to_list(ifudesignlist, evaluate=True)

        self.platetargets = numpy.array( arginp_to_list(platetargets) ) \
                            if platetargets is not None else self._default_plate_target_files()

        if (trg_cat is not None and trg_catid is None) or \
           (trg_cat is None and trg_catid is not None):
            raise Exception('Must provide both trg_cat and trg_catid')

        if trg_cat is not None:
            self.trg_cat = numpy.array( arginp_to_list(trg_cat) )
            self.trg_catid = numpy.array( arginp_to_list(trg_catid) )
        else:
            self.trg_cat, self.trg_catid = self._default_trg_catalogs()


    # ******************************************************************
    #  Default values
    # ******************************************************************

    def _default_plate_target_files(self):
        """
        Return the default plateTargets file used to get the target
        information

        .. todo::
            Move this to :file:`mangadap/util/defaults.py`

        """
        return numpy.array([ os.path.join(environ['MANGACORE_DIR'], 'platedesign', 'platetargets',
                                          'plateTargets-12.par'),
                             os.path.join(environ['MANGACORE_DIR'], 'platedesign', 'platetargets',
                                          'plateTargets-1.par') ])


    def _default_trg_catalogs(self):
        """
        Return the default target catalogs and their respective catalog id
        numbers.

        .. todo::
            Move this to :file:`mangadap/util/defaults.py`

        """
        catalog = numpy.array([ os.path.join(environ['MANGA_TARGET'], 'catalogues',
                                             '12-nsa_v1b_0_0_v2.fits'), 
                                os.path.join(environ['MANGA_TARGET'], 'catalogues',
                                             '1-nsa_v1_0_0.fits') ])
        catalogid = numpy.array([ 12, 1 ])
        return catalog, catalogid


    # ******************************************************************
    #  Utility functions
    # ******************************************************************


    def _drp_mangaid(self, drplist):
        """
        Grab the MaNGA IDs from the DRP fits files.
        
        Args:
            drplist (drpfile) : List :class:`mangadap.drpfile` objects

        Returns:
            numpy.ndarray : Array with the MaNGA IDs from the headers of
                the DRP fits files.

        .. note::
            This takes far too long; astropy.io.fits must be reading the
            entire file instead of just the header.

        """
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
        Grab the MaNGA IDs and object coordinates from the DRP fits
        files.
        
        Args:
            drplist (list) : List of :class:`mangadap.drpfile` objects

        Returns:
            numpy.ndarray, numpy.ndarray, numpy.ndarray : Arrays with
                the MaNGA IDs from the headers of the DRP fits files,
                the object right ascension, and the object declination.

        .. note::
            This takes far too long; astropy.io.fits must be reading the
            entire file instead of just the header.

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
        Read all the platetargets files using
        :class:`mangadap.util.yanny.yanny` and return a single set of
        data.

        Args:
            verbose (bool): (Optional) Suppress terminal output

        Returns:
            numpy.ndarray : Single array with the platetargets data.

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


    def _read_trg_idonly(self, verbose=False):
        """
        Read the target catalogs.

        Args:
            verbose (bool): (Optional) Suppress terminal output

        Returns:
            list, numpy.ndarray : A list of numpy.ndarrays with the
                target IDs and a single numpy.array with the number of
                elements in each trg_id element
        
        """
        # Check the files exist
        if not all([ os.path.exists(n) for n in self.trg_cat ]):
            raise Exception('Cannot open one or more NSA catalogs!')

        trg_id = []
        trg_ndata = []
        for n in self.trg_cat:
            print(n)
            hdu = fits.open(n)

            trg_ndata = trg_ndata + [ hdu[1].header['NAXIS2'] ]
#           trg_id = trg_id + [ numpy.array(hdu[1].data['TARGID']) ]
            trg_id = trg_id + [ numpy.array(hdu[1].data['NSAID']) ]
            hdu.close()

        return trg_id, numpy.array(trg_ndata)


    def _read_trg(self, verbose=False):
        """
        Read all the data from the target catalogs necessary to create
        the drpcomplete data array.  This is only necessary if the
        platetargets files are not available.

        Args:
            verbose (bool): (Optional) Suppress terminal output
        
        Returns:
            list, list, list, list, list, list, numpy.ndarray,
                numpy.ndarray : Lists of numpy.ndarrays with the target
                IDs, velocities, velocity dispersions, ellipticities,
                position angles, effective radii, and single
                numpy.arrays with the number of elements in each list
                element and flags signifying if the target catalog has
                velocity dispersions
        
        """
        # Check the files exist
        if not all([ os.path.exists(n) for n in self.trg_cat ]):
            raise Exception('Cannot open one or more of the target catalogs!')

        trg_id = []
        trg_vel = []
        trg_ell = []
        trg_pa = []
        trg_Reff = []
        trg_ndata = []
        trg_veldisp = []
        trg_has_veldisp = []
        for n in self.trg_cat:
            print(n)
            hdu = fits.open(n)

            trg_ndata = trg_ndata + [ hdu[1].header['NAXIS2'] ]
#           trg_id = trg_id + [ numpy.array(hdu[1].data['TARGID']) ]
            trg_id = trg_id + [ numpy.array(hdu[1].data['NSAID']) ]
            trg_vel = trg_vel + [ numpy.array(hdu[1].data['Z']*constants.c.to('km/s').value) ]
            trg_ell = trg_ell + [ numpy.array(1.0-hdu[1].data['SERSIC_BA']) ]
            trg_pa = trg_pa + [ numpy.array(hdu[1].data['SERSIC_PHI']) ]
            trg_Reff = trg_Reff + [ numpy.array(hdu[1].data['SERSIC_TH50']) ]

            # Try to get the column number with the velocity dispersion
            # to test if the NSA catalog has it.
            try:
                i = hdu[1].columns.names.index('VELDISP')
            except ValueError:
                if verbose:
                    print('WARNING: {0} does not have vel. dispersion.  Using default.'.format(n))
                trg_has_veldisp = trg_has_veldisp + [ False ]
                trg_veldisp = trg_veldisp + [ numpy.empty(0) ]
            else:
                trg_has_veldisp = trg_has_veldisp + [ True ]
                trg_veldisp = trg_veldisp + [ numpy.array(hdu[1].data['VELDISP']) ]

            hdu.close()

        return trg_id, trg_vel, trg_veldisp, trg_ell, trg_pa, trg_Reff, numpy.array(trg_ndata), \
               numpy.array(trg_has_veldisp)


    def _match_platetargets(self, def_veldisp):
        """
        Read the platetargets file, and match the completed DRP
        reductions based on the plate and ifudesign numbers.  The NSA
        catalogs are also read so that the ID can be pulled from the
        appropriate catalog.

        Args:
            def_veldisp (float): Default value to use for the velocity
                dispersion

        Returns:
            numpy.ndarray : Arrays with the MaNGA ID, object right
                ascension, object declination, catalog ID, index of the
                entry in the catalog, target ID, velocity, velocity
                dispersion, ellipticity, position angle, and effective
                radius

        """

        # Read the platetargets file
        print('Reading platetarget file(s)...')
        par_data = self._read_platetargets(verbose=True)
        print('DONE')

        # Read the NSA catalogs
        print('Reading the target catalog(s)...')
        trg_id, trg_ndata = self._read_trg_idonly(verbose=True)
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
        trgid = numpy.full(n_drp, -1, dtype=numpy.int32)
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

            trg_indx = numpy.where(self.trg_catid == catid[i])
            if len(trg_indx[0]) == 0:
#               print('WARNING: No target catalog {0} available.  Setting TARGID=-1.'.format(catid[i]))
                print('WARNING: No target catalog {0} available.  Setting NSAID=-1.'.format(catid[i]))
                trgid[i] = -1
            elif catindx[i] >= trg_ndata[trg_indx]:
#               print('WARNING: No index {0} in target catalog {1}.  Setting TARGID=-1.'.format(
                print('WARNING: No index {0} in target catalog {1}.  Setting NSAID=-1.'.format(
                        catindx[i], catid[i]))
                trgid[i] = -1
            else:
                trgid[i] = numpy.int32(trg_id[trg_indx[0]][catindx[i]])

            vel[i] = par_data['PLTTRGT']['nsa_redshift'][indx][0] * constants.c.to('km/s').value
            ell[i] = 1.0-par_data['PLTTRGT']['nsa_sersic_ba'][indx][0]
            pa[i] = par_data['PLTTRGT']['nsa_sersic_phi'][indx][0]
            Reff[i] = par_data['PLTTRGT']['nsa_sersic_th50'][indx][0]

            if par_data['PLTTRGT']['nsa_vdisp'][indx][0] > 0.:
                veldisp[i] = par_data['PLTTRGT']['nsa_vdisp'][indx][0]

        print('DONE')

        return numpy.array(mangaid), objra, objdec, catid, catindx, trgid, vel, veldisp, ell, \
               pa, Reff


    def _match_trg(self, mangaid, def_veldisp):
        """
        Read the target catalogs and extract the necessary data.

        Args:
            mangaid (list) : List of MaNGA IDs pulled from the DRP file headers

            def_veldisp (float): Default value to use for the velocity
                dispersion

        Returns:
            numpy.ndarray : Arrays with the catalog ID, index of the
                entry in the catalog, target ID, velocity, velocity
                dispersion, ellipticity, position angle, and effective
                radius

        Raises:
            Exception: Raised if the number of MaNGA IDs does not match
                the number of plates in :attr:`platelist`, see above.

        """

        # Check the number of MaNGA IDs provided.
        n_drp = len(self.platelist)
        if len(mangaid) != n_drp:
            raise Exception('Incorrect number of MaNGA IDs')

        # Read the NSA catalogs
        print('Reading target catalog(s)...')
        trg_id, trg_vel, trg_veldisp, trg_ell, trg_pa, trg_Reff, trg_ndata, trg_has_veldisp = \
                self._read_trg(verbose=True)
        print('DONE')

        # Initialize the output arrays (needed in case some DRP targets not found)
        catid = numpy.full(n_drp, -1, dtype=numpy.int32)
        catindx = numpy.full(n_drp, -1, dtype=numpy.int32)
        trgid = numpy.full(n_drp, -1, dtype=numpy.int32)
        vel = numpy.zeros(n_drp, dtype=numpy.float64)
        veldisp = numpy.full(n_drp, def_veldisp, dtype=numpy.float64)
        ell = numpy.zeros(n_drp, dtype=numpy.float64)
        pa = numpy.zeros(n_drp, dtype=numpy.float64)
        Reff = numpy.zeros(n_drp, dtype=numpy.float64)

        print('Using MaNGA ID to get information from target catalog(s)...')
        for i in range(0,n_drp):
           
            catid[i] = numpy.int32(mangaid[i].split('-')[0])
            catindx[i] = numpy.int32(mangaid[i].split('-')[1])

            indx = numpy.where(self.trg_catid == catid[i])
            if len(indx[0]) == 0:
#               print('WARNING: No target catalog {0} available.  Setting TARGID=-1.'.format(catid[i]))
                print('WARNING: No target catalog {0} available.  Setting NSAID=-1.'.format(catid[i]))
                trgid[i] = -1
            elif catindx[i] >= trg_ndata[indx]:
#               print('WARNING: No index {0} in target catalog {1}.  Setting TARGID=-1.'.format(
                print('WARNING: No index {0} in NSA catalog {1}.  Setting NSAID=-1.'.format(
                        catindx[i], catid[i]))
                trgid[i] = -1
            else:
                trgid[i] = numpy.int32(trg_id[indx[0]][catindx[i]])
                vel[i] = numpy.float64(trg_vel[indx[0]][catindx[i]])
                ell[i] = numpy.float64(trg_ell[indx[0]][catindx[i]])
                pa[i] = numpy.float64(trg_pa[indx[0]][catindx[i]])
                Reff[i] = numpy.float64(trg_Reff[indx[0]][catindx[i]])
                if trg_has_veldisp[indx[0]]:
                    veldisp[i] = numpy.float64(trg_veldisp[indx[0]][catindx[i]])

        return catid, catindx, trgid, vel, veldisp, ell, pa, Reff


    def _all_data_exists(self, quiet=True):
        """
        Determines if the data for all the plates/ifudesigns is already
        present in the current drpcomplete file.

        Args:
            quiet (bool): (Optional) Suppress terminal output

        Returns:
            bool : Flag that all data has already been collated for the
                requiested plate/ifudesign list.

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
        """
        Read the data in the existing file at :func:`file_path`.

        """
        hdu = fits.open(self.file_path())
        self.data = hdu[1].data
        self.nrows = hdu[1].header['NAXIS2']
        hdu.close()
        print('Read data: %d rows' % self.nrows)


    def _confirm_access(self, reread):
        """
        Check the drpcomplete file at :func:`file_path` is accessible,
        and read the data if not yet read.

        Args:
            reread (bool): Force the file to be re-read

        """
        inp = self.file_path()
        if not os.path.exists(inp):
            raise Exception('Cannot open file: {0}'.format(inp))

        if self.data is None or reread:
            self._read_data()


    def _find_completed_reductions(self, mindesign=19, combinatorics=False):
        """
        Search the DRP path for reduced CUBE files.

        Function allows for one or both of :attr:`platelist` and
        :attr:`ifudesignlist` to be None.  The behavior is:
        
            - If both are None, all available CUBE files are used to
              create :attr:`platelist` and :attr:`ifudesignlist`.

            - If one is None, all available CUBE files within the
              provided list of the one that is not None are used to
              create :attr:`platelist` and :attr:`ifudesignlist`.  E.g.,
              if :attr:`platelist` =[7443] and :attr:`ifudesignlist` is
              None, all CUBE files with plate=7443 are chosen.

            - If both are not None and they have different lengths, or
              the same length and combinatorics is True, all available
              CUBE files within the constraints of both lists are
              selected. E.g., if :attr:`platelist` =[7443,7495] and
              :attr:`ifudesignlist` =[12704], the CUBE files with
              (plate,ifudesign)=[(7443,12704),(7459,12704)] are chosen.

            - If both are not None and they have the same length and
              combinatorics is False, all available CUBE files in the
              matched lists are chosen.  E.g. if :attr:`platelist`
              =[7443,7495] and :attr:`ifudesignlist` =[12703,12703], the
              CUBE files with
              (plate,ifudesign)=[(7443,12704),(7459,12704) are chosen

        Args:
            mindesign (int): (Optional) Minimum bundle design to
                consider.  The bundle design is determined by
                ifudesign/100, such that::

                    if mindesign is not None and b/100 < mindesign:
                        continue

                a DRP file is ignored.  Default only ignores the
                minibundles (design=7).

            combinatorics (bool): (Optional) Create :attr:`platelist`
                and :attr:`ifudesignlist` by determining all possible
                combinations of the input values. See above.

        Returns:
            list, list: Lists with the available plates and ifudesigns
                for which to collect data.

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
            return plates, ifudesigns
#            return numpy.array(plates), numpy.array(ifudesigns)

        for root, dir, files in walk(path):
            for file in files:
                if file.endswith('-LOGCUBE.fits.gz'):
                    p, b, m = parse_drp_file_name(file)

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
        Using the provided list of CUBE DRP files, find if the RSS mode
        is also available.

        Currently only two modes are possible:

            - modes[i] = 1 -> only 'CUBE' file are available

            - modes[i] = 2 -> both 'CUBE' and 'RSS' files are available

        Args:
            drplist (list) : List of :class:`mangadap.drpfile` objects

        Returns:
            numpy.ndarray : Array of modes for each input DRP file.

        """
        n_drp = len(drplist)
        modes = numpy.empty(n_drp, dtype=numpy.uint8)
        for i in range(0,n_drp):
            drpf = drpfile.drpfile(drplist[i].plate, drplist[i].ifudesign, 'RSS',
                                   drpver=drplist[i].drpver)
            modes[i] = 2 if os.path.exists(drpf.file_path()) else 1
        return modes


    def _write_parameter_struct(self, ostream):
        """
        Write the yanny file structure to a file.

        Args:
            ostream (_io.TextIOWrapper) : File stream for output
        
        """
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
        """
        Write a DAPPAR entry to the yanny file.

        Args:
            ostream (_io.TextIOWrapper) : File stream for output
            index (int) : Index of the entry in :attr:`data` to write
            mode (int) : Mode of the entry; see :func:`_find_modes`

        """
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
               alldrp=False, def_veldisp=100.0, use_trg=False):
        """
        Update the DRP complete file.
        
        If *platelist* and/or *ifudesignlist* are provided, the existing
        :attr:`platelist` and :attr:`ifudesignlist` attributes are
        replaced.

        If *platelist* and *ifudesignlist* do not have the same length
        or *combinatorics* is True, the lists are expanded to include
        all combinations of their elements.

        If *platelist* and *ifudesignlist* are None or *alldrp* is True,
        all available plates/ifudesigns are collected from within the
        DRP directory structure.

        If the result of :func:`file_path` does not exist, it is
        created.

        If the result of :func:`file_path` does exist, the available
        plates and ifudesigns in the file are compared against the list
        (provided or collected) to update.  If the lists differ, the
        drpcomplete file is re-created from scratch.  If all the plates
        and ifudesigns are available, nothing is done, unless
        force=True.

        Args:
            platelist (str or list): (Optional) List of plates to
                include in the drpcomplete file.
            ifudesignlist (str or list): (Optional) List of ifudesigns
                to include in the drpcomplete file.
            combinatorics (bool): (Optional) Determine all combinations
                of the entered plates and ifudesigns.
            force (bool): (Optional) Overwrite any existing drpcomplete
                file with a new one built from scratch.
            alldrp (bool): (Optional) Find the full list of available
                DRP files.
            def_veldisp (float): (Optional) Default value of the
                velocity dispersion if the sourced databases do not
                contain it.
            use_trg (bool): (Optional) Force the matching to be done
                using the target catalog(s).  Default is to use the
                platetargets files.

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
        if use_trg:
            mangaid, objra, objdec = self._drp_info(drplist)
            catid, catindx, trgid, vel, veldisp, ell, pa, Reff = self._match_trg(mangaid,
                                                                                 def_veldisp)
        else:
            mangaid, objra, objdec, catid, catindx, trgid, vel, veldisp, ell, pa, Reff \
                    = self._match_platetargets(def_veldisp)

        self.write([drplist[i].plate for i in range(0,nn)],
                   [drplist[i].ifudesign for i in range(0,nn)], modes, mangaid, objra, objdec,
                   catid, catindx, trgid, vel, veldisp, ell, pa, Reff)

        self._read_data()


    def write(self, platelist, ifudesignlist, modes, mangaid, objra, objdec, catid, catindx, trgid,
              vel, veldisp, ell, pa, Reff, drpver=None, redux_path=None, dapver=None,
              analysis_path=None, platetargets=None, trg_cat=None, trg_catid=None, clobber=True):
        """
        Write the drpcomplete fits binary table.

        Header keywords are:
            - VERSDRP: DRP version
            - RDXPTH: DRP reduction path
            - VERSDAP: DAP version
            - SISPTH: DAP analysis path
            - NSACAT *N*: Location of target catalog *N*
            - NSACID *N*: ID Number of target catalog *N*
            - PLTARG *N*: plateTargets file *N*
            - AUTHOR: Set to 'K.B. Westfall <kbwestfall@gmail.com>'

        Binary table columns:
            1. PLATE (1J): Plate number
            2. IFUDESIGN (1J): IFU design
            3. MODES (1B): Modes available 
                * MODES=1: only 'CUBE' file is available
                * MODES=2: both 'CUBE' and 'RSS' files are available

            4. MANGAID (*n* A): MaNGA ID (same as 'CATID-CATINDX')
            5. OBJRA (1D): Object right ascension
            6. OBJDEC (1D): Object declination
            7. CATID (1J): Catalog ID used for target selection
            8. CATINDX (1J): Index in catalog with target data
            9. TARGID (1J): Target ID in catalog
            10. VEL (1D): Velocity from catalog, if available
            11. VDISP (1D): Velocity dispersion from catalog, if
                available
            12. ELL (1D): Ellipticity (1-b/a) from catalog, if available
            13. PA (1D): Position angle from catalog, if available
            14. REFF (1D): Effective radius from catalog, if available

        Args:
            platelist (list) : List of plates
            ifudesignlist (list) : List of IFU designs
            modes (numpy.ndarray): Mode values, see above
            mangaid (numpy.ndarray): MaNGA IDs
            objra (numpy.ndarray): Object right ascensions
            objdec (numpy.ndarray): Object declinations
            catid (numpy.ndarray): Catalog ID used for target selection
            catindx (numpy.ndarray): Index in catalog with target data
            trgid (numpy.ndarray): Target ID in catalog
            vel (numpy.ndarray): Velocity from catalog, if available
            veldisp (numpy.ndarray): Velocity dispersion from catalog,
                if available
            ell (numpy.ndarray): Ellipticity (1-b/a) from catalog, if
                available
            pa (numpy.ndarray): Position angle from catalog, if
                available
            Reff (numpy.ndarray): Effective radius from catalog, if
                available
            drpver (str) : (Optional) DRP version, see above.
            redux_path (str) : (Optional) Path to the top level of
                directory for the DRP output, see above.
            dapver (str) : (Optional) DAP version, see above.
            analysis_path (str) : (Optional) Path to the top level
                directory for the DAP output files, see above.
            platetargets (numpy.ndarray) : (Optional) List of the
                plateTargets files
            trg_cat (numpy.ndarray): (Optional) List of catalogs to
                search through to grab the target ID from that catalog.
            trg_catid (numpy.ndarray): (Optional) List of target catalog
                ID numbers to use.
            clobber (bool): (Optional) Overwrite any existing file.

        Raises:
            Exception: Raised if the drpcomplete file was opened as read
                only or if the drpcomplete file exists and clobber=False.

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
        if trg_cat is None:
            trg_cat = self.trg_cat
        if trg_catid is None:
            trg_catid = self.trg_catid
        for i in range(0,trg_cat.size):
#           hdr['TARGCAT{0}'.format(i+1)] = trg_cat[i]
#           hdr['TARGCID{0}'.format(i+1)] = trg_catid[i]
            hdr['NSACAT{0}'.format(i+1)] = trg_cat[i]
            hdr['NSACID{0}'.format(i+1)] = trg_catid[i]
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
#       col9 = fits.Column(name='NSAID', format='1J', array=nsaid)
        col9 = fits.Column(name='TARGID', format='1J', array=trgid)
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
        """
        Return the data for a given plate-ifudesign, or index in the
        :attr:`data`.  Even though *plate*, *ifudesign*, and *index* are
        all optional, the arguments must contain either *plate* and
        *ifudesign* or *index*.

        Args:
            plate (int): (Optional) Plate number
            ifudesign (int): (Optional) IFU design
            index (int): (Optional) Index of the row in :attr:`data`
                with the data to return
            reread (bool): (Optional) Force the database to be re-read

        Returns:
            list : List of the 14 elements of :attr:`data`; see
                :func:`write`.

        Raises:
            ValueError: Raised if the row with the data is unknown
                because either *index* is not defined or one or both of
                *plate* and *ifudesign* is not defined.  Also raised if
                *index* does not exist in the data array.

        """

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
#                self.data['TARGID'][index], self.data['VEL'][index], self.data['VDISP'][index],
                 self.data['NSAID'][index], self.data['VEL'][index], self.data['VDISP'][index],
                 self.data['ELL'][index], self.data['PA'][index], self.data['REFF'][index] ]


    def write_par(self, ofile, mode, plate=None, ifudesign=None, index=None, reread=False,
                  clobber=True):
        """
        Write the SDSS parameter (Yanny) file for use with the MaNGA DAP.

        Args: 
            ofile (str): Output file name
            mode (str): Mode of the DRP file to analyze; must be either
                'RSS' or 'CUBE'
            plate (int): (Optional) Plate number
            ifudesign (int): (Optional) IFU design
            index (int): (Optional) Index of the row in :attr:`data`
                with the data to return
            reread (bool): (Optional) Force the database to be re-read
            clobber (bool): (Optional) Overwrite any existing parameter
                file

        Raises:
            IOError: Raised if the parameter file already exists and
                clobber is False.
            ValueError: Raised if

                - the row with the data is unknown because either
                  *index* is not defined or one or both of *plate* and
                  *ifudesign* is not defined.
                - *index* does not exist in the data array.
                - the 'RSS' mode is selected but unavailable


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

        Args:
            plate (int): Plate number
            ifudesign (int): IFU design
            reread (bool): (Optional) Force the database to be re-read

        Returns:
            int: Index of the row in :attr:`data` with the data for the
                given *plate* and *ifudesign*

        Raises:
            Exception: Raised if the given *plate* and *ifudesign* were
                not found.

        """
        self._confirm_access(reread)

        for i in range(0,self.nrows):
            if self.data['PLATE'][i] == plate and self.data['IFUDESIGN'][i] == ifudesign:
                return i

        raise Exception('Could not find plate={0},ifudesign={1} in drpcomplete file!'.format(plate,
                        ifudesign))



