# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-
"""
Defines a class used to produce and provide an interface for data
required to create :ref:`execution-config`.

The drpcomplete file/data is primarily required for the survey-level
:ref:`execution` of the DAP. However, it does provide useful information
regarding the completed DRP files.

----

.. include license and copyright
.. include:: ../include/copy.rst

----

.. include common links, assuming primary doc root is up one directory
.. include:: ../include/links.rst
"""

import os
import time
import warnings
import glob


from IPython import embed

import numpy

from astropy.io import fits
import astropy.constants

from pydl.pydlutils.yanny import yanny

from ..datacube import MaNGADataCube
from ..spectra import MaNGARSS
from ..config import defaults
from ..util.parser import arginp_to_list, list_to_csl_string, parse_drp_file_name

class DRPComplete:
    r"""
    Database with information needed by the DAP to analyze the completed
    DRP files.

    This class searches the defined paths for files resulting from the
    3D phase of the MaNGA DRP, and then collates the information
    necessary to run those files through the MaNGA DAP.  The necessary
    parameters are pulled from the provided platetargets files; see
    :func:`update`.

    Args:

        platelist (:obj:`str`, :obj:`list`, optional):
            List of plates to search for.  Default is to search the full
            DRP path.
        ifudesignlist (:obj:`str`, :obj:`list`, optional):
            List of ifudesign to search for.  Default is to search the
            full DRP path.
        drpall (:obj:`str`, optional):
            The full path to the DRPall fits file.  Default is set by
            :func:`mangadap.config.defaults.drpall_file`.
        platetargets (:obj:`str`, :obj:`list`, optional):
            List of platetargets files to search through to find any
            given plate ifudesign combination.  Default is returned as
            the first element in
            :func:`mangadap.config.defaults.plate_target_files`.
        catid (:obj:`str`, :obj:`list`, optional):
            List of target catalog ID numbers.  Default is returned as
            the second element in
            :func:`mangadap.config.defaults.plate_target_files`.
        drpver (:obj:`str`, optional):
            DRP version, which is:

                - used to define the default DRP redux path
                - used when declaring a drpfits instance
                - used in the name of the drpcomplete fits file
                - included as a header keyword in the output file

            Default is defined by
            :func:`mangadap.config.defaults.drp_version`.
        redux_path (:obj:`str`, optional): 
            The path to the top level directory containing the DRP
            output files; this is the same as the ``redux_path`` in
            the :class:`mangadap.util.drpfits.DRPFits` class. Default
            is defined by
            :func:`mangadap.config.defaults.drp_redux_path`.
        dapver (:obj:`str`, optional):
            DAP version, which is:

                - used to define the default DAP analysis path
                - included as a header keyword in the output drpcomplete
                  fits file

            Default is defined by
            :func:`mangadap.config.defaults.dap_version`.
        analysis_path (:obj:`str`, optional):
            The path to the top level directory for the DAP output
            files; this is **different** from the directory_path in the
            :class:`mangadap.dapfile` class.  Default is defined by
            :func:`mangadap.config.defaults.dap_analysis_path`
        directory_path (:obj:`str`, optional):
            Direct path to the output file produced using
            :func:`mangadap.config.defaults.dap_common_path`
        readonly (:obj:`bool`, optional):
            Flag that the drpcomplete fits file is only opened for
            reading, not for updating.  If True, many of the attributes
            will be set to None.

    Attributes:
        platelist (:obj:`list`):
            List of plates to search for, see above.
        ifudesignlist (:obj:`list`):
            List of IFU designs to search for, see above.
        drpall (:obj:`str`):
            The DRPall file, see above.
        platetargets (`numpy.ndarray`):
            List of platetargets files to search through, see above.
        catid (`numpy.ndarray`):
            List of target catalog IDs, see above
        drpver (:obj:`str`):
            DRP version, see above.
        redux_path (:obj:`str`):
            Path to the top level of directory for the DRP output, see
            above.
        dapver (:obj:`str`):
            DAP version, see above.
        analysis_path (:obj:`str`):
            Path to the top level directory for the DAP output files,
            see above.
        directory_path (:obj:`str`):
            Direct path to the output file produced using
            :func:`mangadap.config.defaults.dap_common_path`
        readonly (:obj:`bool`):
            Flag that the drpcomplete fits file is only opened for
            reading, not for updating.
        hdu (`astropy.io.fits.HDUList`_):
            Fits data with binary table data.
        nobs (:obj:`int`):
            Number of observations in the file

    """
    def __init__(self, platelist=None, ifudesignlist=None, drpall=None, platetargets=None,
                 catid=None, drpver=None, redux_path=None, dapver=None, analysis_path=None,
                 directory_path=None, readonly=False):

        # Input properties
        self.drpver = defaults.drp_version() if drpver is None else str(drpver)
        self.redux_path = defaults.drp_redux_path(self.drpver) \
                                if redux_path is None else str(redux_path)

        self.dapver = defaults.dap_version() if dapver is None else str(dapver)
        self.analysis_path = defaults.dap_analysis_path(self.drpver, self.dapver) \
                             if analysis_path is None else str(analysis_path)
        self.directory_path = defaults.dap_common_path(drpver=self.drpver, dapver=self.dapver,
                                                       analysis_path=self.analysis_path) \
                             if directory_path is None else str(directory_path)

        self.hdu = None
        self.nobs = None
        if os.path.exists(self.file_path()):
            print('READING')
            self._read()

        if readonly:
            self.readonly=True
            self.platelist = None
            self.ifudesignlist = None
            self.platetargets = None
            self.catid = None
            return

        self.readonly = False

        self.platelist = arginp_to_list(platelist, evaluate=True)
        self.ifudesignlist = arginp_to_list(ifudesignlist, evaluate=True)

        self.drpall = defaults.drpall_file(drpver=self.drpver,redux_path=self.redux_path) \
                            if drpall is None else drpall
        
        if (platetargets is not None and catid is None) or \
           (platetargets is None and catid is not None):
            raise ValueError('To use user-provided platetargets files, must provide both '
                             'platetargets and catid.')

        self.platetargets = None
        if platetargets is not None:
            self.platetargets = numpy.array( arginp_to_list(platetargets) )
            self.catid = numpy.array( arginp_to_list(catid) ).astype(int)
        else:
            try:
                self.platetargets, self.catid = defaults.plate_target_files()
            except:
                warnings.warn('Could not define platetargets files.  '
                              'Updates must use DRPall file.')

    def __getitem__(self, key):
        return self.hdu['DRPC'].data[key]

    def __len__(self):
        return self.nobs

    # ******************************************************************
    #  Utility functions
    # ******************************************************************

    def _read_platetargets(self):
        """
        Read all the platetargets files using `pydl.pydlutils.yanny`_
        and return a list of yanny structures.

        Returns:
            list: A list of yanny structures, one per platetargets file.

        Raises:
            FileNotFoundError: Raised if cannot open one or more
                platetargets files.

        .. todo::
            This should be made more efficient by collating the required
            plateTargets data into a single record array!
        """
        plttrg_data = []

        root = '/'+os.path.join(*(self.platetargets[0].split('/')[:-1]))
        print('plateTargets root directory: {0}'.format(root))

        for p in self.platetargets:
            pfile = p.split('/')[-1]
            print('Reading plateTargets file: {0}'.format(pfile), end='\r')
            # Check that the file exists
            if not os.path.exists(p):
                raise FileNotFoundError('Cannot open {0}!'.format(p))
            # Read and append the data
            plttrg_data.append( yanny(filename=p) )
        print('\nReading plateTargets file: DONE')

        return plttrg_data

    def _read_fix_data(self, fix_file):
        """
        Read fix data from the provided file.

        Args:
            fix_file (:obj:`str`):
                SDSS parameter file with the fix data.

        Returns:
            yanny: A yanny structure with the data fixes.  Returns
            None and raises a warning if the file does not exist.
        """
        if not os.path.isfile(fix_file):
            warnings.warn('No redshift fix file available.')
            return None
        return yanny(filename=fix_file)

    def _read_redshift_fix(self):
        """
        Wrapper for :func:`_read_fix_data` and
        :func:`~mangadap.config.defaults.redshift_fix_file` that
        returns the redshift fix data.
        """
        return self._read_fix_data(defaults.redshift_fix_file())

    def _read_photometry_fix(self):
        """
        Wrapper for :func:`_read_fix_data` and
        :func:`~mangadap.config.defaults.photometry_fix_file` that
        returns the photometry fix data.
        """
        return self._read_fix_data(defaults.photometry_fix_file())

    def _match_platetargets(self, quiet=True):
        """
        Read the platetargets files and match the data therein to the
        completed DRP files based on the plate and ifudesign.

        If a plate-ifudesign combination is not found in the
        plateTargets file, the other parameter values are set to -9999
        and the MaNGA ID is set to NULL.

        If the plate-ifudesign combination is found, the columns in the
        plateTargets files that the code expects to find are 'plateid',
        'ifudesign', 'mangaid', 'object_ra', and 'object_dec'; the
        values that will be replaced with 'NULL' (str) or -9999 (int or
        float) if they do not exist are 'nsa_version', 'nsa_nsaid',
        'manga_target1', 'manga_target3', 'z' or 'nsa_z',
        'nsa_elpetro_ba', 'nsa_elpetro_phi', 'nsa_elpetro_th50_R',
        'nsa_vdisp'.
        
        .. todo::

            - Instead of searching through all files for the correct
              plate-ifudesign, use the MaNGA ID to get the correct
              catalog?  Requires getting the MaNGA ID from somewhere,
              but it takes a long time to read from the fits files...
            - Is there some more concise method of performing the same
              thing as what is done below with many try/except blocks?

        Args:
            quiet (:obj:`bool`, optional):
                Suppress terminal output

        Returns:
            numpy.array: 14 arrays with: MaNGA ID, object right
            ascension, object declination, catalog ID, index of the
            entry in the catalog, catalog version, ID of object in the
            catalog, main MaNGA survey target bitmask, ancillary MaNGA
            survey target bitmask, velocity, velocity dispersion,
            ellipticity, position angle, and effective radius.
        """
        # Read the platetargets file
        plttrg_data = self._read_platetargets()
        ntrg = len(plttrg_data)
        print('Read {0} plateTargets file(s)'.format(ntrg))

        # Initialize the output arrays (needed in case some DRP targets not found)
        n_drp = len(self.platelist)
        mangaid = []
        # TODO: try mangaid = numpy.empty(n_drp, dtype=object) ?
        objra = numpy.full(n_drp, -9999.0, dtype=numpy.float64)
        objdec = numpy.full(n_drp, -9999.0, dtype=numpy.float64)
        catid = numpy.full(n_drp, -9999, dtype=numpy.int32)
        catindx = numpy.full(n_drp, -9999, dtype=numpy.int32)
        trg_version = []
        # TODO: try trg_version = numpy.empty(n_drp, dtype=object) ?
        trg_id = numpy.full(n_drp, -9999, dtype=numpy.int32)
        manga_trg1 = numpy.full(n_drp, -9999, dtype=numpy.int64)
        manga_trg3 = numpy.full(n_drp, -9999, dtype=numpy.int64)
        vel = numpy.full(n_drp, -9999.0, dtype=numpy.float64)
        veldisp = numpy.full(n_drp, -9999.0, dtype=numpy.float64)
        ell = numpy.full(n_drp, -9999.0, dtype=numpy.float64)
        pa = numpy.full(n_drp, -9999.0, dtype=numpy.float64)
        Reff = numpy.full(n_drp, -9999.0, dtype=numpy.float64)

        print('Searching platetargets file for observed galaxies...', end='\r')
        for i in range(n_drp):
            plttrg_j = 0
            mangaid_i = 'NULL'
            for j in range(ntrg):
                indx = numpy.where((plttrg_data[j]['PLTTRGT']['plateid'] == self.platelist[i]) &
                                (plttrg_data[j]['PLTTRGT']['ifudesign'] == self.ifudesignlist[i]))
                if len(indx[0]) > 1:
                    raise ValueError('Multiple instances of {0}-{1} (MaNGA ID={2}) in {3}'.format(
                                     self.platelist[i], self.ifudesignlist[i], mangaid_i,
                                     self.platetargets[j]))
                if len(indx[0]) == 1:
                    plttrg_j = j
                    mangaid_i = plttrg_data[j]['PLTTRGT']['mangaid'][indx][0].decode("ascii")
                    if not quiet:
                        print('Found {0}-{1} (MaNGA ID={2}) in {3} (CAT={4})'.format(
                                        self.platelist[i], self.ifudesignlist[i], mangaid_i,
                                        self.platetargets[j], self.catid[j]))
                    break

            if len(indx[0]) == 0:
                warnings.warn('Could not find {0}-{1} in any plateTargets file!'.format(
                                                        self.platelist[i], self.ifudesignlist[i]))
                mangaid = mangaid + ['NULL']
                continue
            mangaid.append(plttrg_data[plttrg_j]['PLTTRGT']['mangaid'][indx][0].decode("ascii"))
            objra[i] = plttrg_data[plttrg_j]['PLTTRGT']['object_ra'][indx][0]
            objdec[i] = plttrg_data[plttrg_j]['PLTTRGT']['object_dec'][indx][0]

            # From David Wake:
            #
            # "MANGAID consists of CATID-CATIND, where CATID identifies
            # a parent catalog, in the case of the main manga samples
            # CATID = 1 which refers to nsa_v1_0_0.fits, and CATIND is
            # the position within that catalog (zero indexed). So if you
            # strip out CATIND from MANGAID you have the position within
            # nsa_v1_0_0.fits without any matching required."
            catid[i], catindx[i] = map(lambda x: int(x), mangaid[i].split('-'))
            if catid[i] != self.catid[plttrg_j]:
                warnings.warn('{0}-{1} (MaNGA ID={2}) found in {3} (CAT={4})!'.format(
                                self.platelist[i], self.ifudesignlist[i], mangaid_i,
                                self.platetargets[j], self.catid[j]))

            try:
                trg_version.append(plttrg_data[plttrg_j]['PLTTRGT']['nsa_version'][indx][0].decode(
                                                                                        'ascii'))
            except:
                trg_version.append('NULL')

            try:
                trg_id[i] = plttrg_data[plttrg_j]['PLTTRGT']['nsa_nsaid'][indx][0]
            except:
                trg_id[i] = -9999
                
            try:
                manga_trg1[i] = plttrg_data[plttrg_j]['PLTTRGT']['manga_target1'][indx][0]
            except:
                manga_trg1[i] = -9999

            try:
                manga_trg3[i] = plttrg_data[plttrg_j]['PLTTRGT']['manga_target3'][indx][0]
            except:
                manga_trg3[i] = -9999

            # As of MPL-6, platetargets files now include the redshift
            # from the targeting catalog which is the combination of the
            # NSA data and the ancillary targets; the NSA only redshift
            # column is 'nsa_z'.  To be compatible with previous
            # versions, first try to grab the redshift using the keyword
            # 'z', then try with 'nsa_z', then just set it to -9999.
            try:
                vel[i] = plttrg_data[plttrg_j]['PLTTRGT']['z'][indx][0] \
                            * astropy.constants.c.to('km/s').value
            except:
                try:
                    vel[i] = plttrg_data[plttrg_j]['PLTTRGT']['nsa_z'][indx][0] \
                                * astropy.constants.c.to('km/s').value
                except:
                    vel[i] = -9999.0

            try:
                if plttrg_data[plttrg_j]['PLTTRGT']['nsa_elpetro_ba'][indx][0] < 0:
                    raise
                ell[i] = 1.0-plttrg_data[plttrg_j]['PLTTRGT']['nsa_elpetro_ba'][indx][0]
            except:
                ell[i] = -9999.0
                warnings.warn('Ellipticity for {0}-{1} is bogus!'.format(self.platelist[i],
                                                                         self.ifudesignlist[i]))

            try:
                # THERE'S A BUG IN THE MANGACORE/V1_2_0 PLATETARGETS FILES
                pa[i] = plttrg_data[plttrg_j]['PLTTRGT']['nsa_elpetro_phi'][indx][0]
                if pa[i] < 0:
                    raise
            except:
                pa[i] = -9999.0
                warnings.warn('PA for {0}-{1} is bogus!'.format(self.platelist[i],
                                                                self.ifudesignlist[i]))

            try:
                Reff[i] = plttrg_data[plttrg_j]['PLTTRGT']['nsa_elpetro_th50_r'][indx][0]
            except:
                Reff[i] = -9999.0

            try:
                veldisp[i] = plttrg_data[plttrg_j]['PLTTRGT']['nsa_vdisp'][indx][0]
            except:
                veldisp[i] = -9999.0

            # Correct for known nan issue
            if numpy.isnan(ell[i]):
                raise ValueError('nan encountered!')
            if numpy.isnan(pa[i]):
                raise ValueError('nan encountered!')

        print('Searching platetargets file for observed galaxies...DONE')

        return numpy.array(mangaid), objra, objdec, catid, catindx, numpy.array(trg_version), \
               trg_id, manga_trg1, manga_trg3, vel, veldisp, ell, pa, Reff

    def _match_drpall(self):
        """
        Find the data needed for the DRP complete databas in the DRPall
        file.

        If a plate-ifudesign combination is not found in the DRPall
        file, the other parameter values are set to -9999 and the MaNGA
        ID is set to NULL.

        If the plate-ifudesign combination is found, the columns in the
        plateTargets files that the code expects to find are 'plate',
        'ifudsgn', 'mangaid', 'objra', 'objdec'; the values that will be
        replaced with 'NULL' (str) or -9999 (int or float) if they do
        not exist are 'nsa_version', 'nsa_nsaid', 'mngtarg1',
        'mngtarg3', 'z', 'nsa_elpetro_ba', 'nsa_elpetro_phi',
        'nsa_elpetro_th50_r', 'nsa_vdisp'.
        
        Returns:
            numpy.array: 14 arrays with: MaNGA ID, object right
            ascension, object declination, catalog ID, index of the
            entry in the catalog, catalog version, ID of object in the
            catalog, main MaNGA survey target bitmask, ancillary MaNGA
            survey target bitmask, velocity, velocity dispersion,
            ellipticity, position angle, and effective radius.
        """
        # Open the drpall file
        hdu = fits.open(self.drpall)
        # Check the plateifus are unique
        pltifu = hdu[1].data['PLATEIFU']
        if len(numpy.unique(pltifu)) != len(pltifu):
            raise ValueError('The PLATEIFU column in the DRPall file is not unique!')

        # Find the rows in the DRPall file with the correct; assumes the 
        this_pltifu = numpy.array(['{0}-{1}'.format(p,i) 
                                    for p,i in zip(self.platelist,self.ifudesignlist)])
        rows = numpy.array([numpy.where(pltifu == pi)[0][0] if pi in pltifu else -1 \
                                for pi in this_pltifu])
        if numpy.any(rows < 0):
            raise ValueError('Could not find the following in the DRPall file: '
                             '{0}'.format(this_pltifu[rows < 0]))
        indx = rows > -1
        rows = rows[indx]

        # Initialize the data
        n_drp = len(self.platelist)
        mangaid = numpy.full(n_drp, 'NULL', dtype=object)
        objra = numpy.full(n_drp, -9999.0, dtype=float)
        objdec = numpy.full(n_drp, -9999.0, dtype=float)
        catid = numpy.full(n_drp, -9999, dtype=int)
        catindx = numpy.full(n_drp, -9999, dtype=int)
        trg_version = numpy.full(n_drp, 'NULL', dtype=object)
        trg_id = numpy.full(n_drp, -9999, dtype=int)
        manga_trg1 = numpy.full(n_drp, -9999, dtype=int)
        manga_trg3 = numpy.full(n_drp, -9999, dtype=int)
        vel = numpy.full(n_drp, -9999.0, dtype=float)
        veldisp = numpy.full(n_drp, -9999.0, dtype=float)
        ell = numpy.full(n_drp, -9999.0, dtype=float)
        pa = numpy.full(n_drp, -9999.0, dtype=float)
        Reff = numpy.full(n_drp, -9999.0, dtype=float)

        # Get the data
        mangaid[indx] = hdu[1].data['mangaid'][rows]
        objra[indx] = hdu[1].data['objra'][rows]
        objdec[indx] = hdu[1].data['objdec'][rows]
        for i in numpy.where(indx)[0]:
            catid[i], catindx[i] = map(lambda x: x.astype(int).tolist(),
                                        numpy.array(mangaid[i].split('-')))
        trg_version[indx] = hdu[1].data['nsa_version'][rows]
        trg_id[indx] = hdu[1].data['nsa_nsaid'][rows]
        manga_trg1[indx] = hdu[1].data['mngtarg1'][rows]
        manga_trg3[indx] = hdu[1].data['mngtarg3'][rows]
        vel[indx] = hdu[1].data['z'][rows] * astropy.constants.c.to('km/s').value
        # There is no vdisp
        ell[indx] = 1 - hdu[1].data['nsa_elpetro_ba'][rows]
        pa[indx] = hdu[1].data['nsa_elpetro_phi'][rows]
        Reff[indx] = hdu[1].data['nsa_elpetro_th50_r'][rows]

        # Done
        hdu.close()
        return mangaid.astype(str), objra, objdec, catid, catindx, trg_version.astype(str), \
                    trg_id, manga_trg1, manga_trg3, vel, veldisp, ell, pa, Reff

    def _all_data_exists(self):
        """
        Determines if the data for all the plates/ifudesigns selected in
        the current compilation is already present in the current
        drpcomplete fits file.

        Returns:
            :obj:`bool`: Flag that all data has already been collated
            for the requested plate/ifudesign list.
        """
        if not self._confirm_access():
            return False
        for p,i in zip(self.platelist, self.ifudesignlist):
            if numpy.sum((self['PLATE'] == p) & (self['IFUDESIGN'] == i)) == 0:
                return False
        return True

    def _read(self):
        """Read the data in the existing file at :func:`file_path`."""
        if self.hdu is not None:
            self.hdu.close()
            self.hdu = None
        self.hdu = fits.open(self.file_path())
        self.nobs = self.hdu['DRPC'].header['NAXIS2']
        print('Read data: {0} rows'.format(self.nobs))

    def _confirm_access(self, reread=False):
        """
        Check the drpcomplete fits file at :func:`file_path` is
        accessible, and read the data if not yet read.

        Args:
            reread (:obj:`bool`, optional):
                Force the file to be re-read
        
        Returns:
            :obj:`bool`: Flag that :attr:`hdu` is read and valid.
        """
        if not os.path.exists(self.file_path()):
            return False
        if self.hdu is None or reread:
            self._read()
        return True

    def _find_completed_reductions(self, mindesign=19, combinatorics=False, on_disk=False):
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
              =[7443,7495] and :attr:`ifudesignlist` =[12704,12703], the
              CUBE files with
              (plate,ifudesign)=[(7443,12704),(7459,12703) are chosen

        Args:
            mindesign (:obj:`int`, optional):
                Minimum bundle design to consider.  For example, to
                ignore all the 7-fiber bundles, set `mindesign=19`
                (default) to only select the bundles with 19 or more
                fibers.
            combinatorics (:obj:`bool`, optional):
                Create :attr:`platelist` and :attr:`ifudesignlist` by
                determining all possible combinations of the input
                values. See above.
            on_disk (:obj:`bool`, optional):
                When searching for available files to analyze, search
                the DRP directory path instead of using the data in the
                DRPall file.

        Returns:
            list: Two lists with the available plates and ifudesigns for
            which to collect data.
        """
        # Are the plate and IFU lists matched?
        matchedlist = (self.platelist is not None and self.ifudesignlist is not None \
                       and (len(self.platelist) == len(self.ifudesignlist) and not combinatorics))

        # Get the list of files
        if matchedlist:
            # Lists already matched, just construct the file names
            files = [os.path.join(*MaNGADataCube.default_paths(p, i, drpver=self.drpver,
                                                               redux_path=self.redux_path)) \
                        for p,i in zip(self.platelist, self.ifudesignlist)]
        elif on_disk:
            # Find the DRP LOGCUBE files on disk
            print('Searching for completed DRP CUBE files...', end='\r')
            files = glob.glob(os.path.join(self.redux_path, '*', 'stack', '*-LOGCUBE.fits.gz'))
            print('Searching for completed DRP CUBE files...DONE.')
        else:
            # Use the DRPall file
            drpall_hdu = fits.open(self.drpall)
            pltifu = drpall_hdu[1].data['PLATEIFU']
            if len(numpy.unique(pltifu)) != len(pltifu):
                raise ValueError('The PLATEIFU column in the DRPall file is not unique!')
            files = [os.path.join(*MaNGADataCube.default_paths(int(p), int(i), drpver=self.drpver,
                                                               redux_path=self.redux_path)) \
                        for p,i in zip(drpall_hdu[1].data['plate'], drpall_hdu[1].data['ifudsgn'])]

        # Only use those files that exist
        for i in range(len(files)):
            if not os.path.isfile(files[i]):
                warnings.warn('No such file: {0}'.format(files[i]))
                del files[i]

        # Get the list of plate and ifus
        pltifu = numpy.array(list(map(lambda x: os.path.split(x)[-1].split('-')[1:3],
                                      files))).astype(int)

        # Ignore any below the minimum ifusize
        pltifu = pltifu[numpy.invert(pltifu[:,1]//100 < mindesign),:]
        if matchedlist:
            return pltifu[:,0], pltifu[:,1]

        # Select those that match to the provided list
        indx = numpy.ones(len(pltifu), dtype=bool) if self.platelist is None \
                    else numpy.array([p in self.platelist for p in pltifu[:,0]])
        indx &= numpy.ones(len(pltifu), dtype=bool) if self.ifudesignlist is None \
                    else numpy.array([i in self.ifudesignlist for i in pltifu[:,1]])
        return pltifu[indx,0], pltifu[indx,1]
        
    def _find_modes(self):
        """
        Using the provided list of CUBE DRP files, find if the RSS mode
        is also available.

        .. warning::
            - This assumes everything in :attr:`platelist` and
              :attr:`ifudesignlist` has a 'CUBE' file.

        Currently only two modes are possible:

            - (1) only 'CUBE' file are available

            - (2) both 'CUBE' and 'RSS' files are available

        Returns:
            `numpy.ndarray`: Array of modes for each input DRP file.
        """
        print('Checking for RSS counterparts...', end='\r')
        has_rss = [os.path.isfile(os.path.join(*MaNGARSS.default_paths(p, i, drpver=self.drpver,
                                                                    redux_path=self.redux_path)))
                        for p,i in zip(self.platelist, self.ifudesignlist)]
        print('Checking for RSS counterparts...DONE.')
        modes = numpy.ones(len(has_rss), dtype=int)
        modes[has_rss] = 2
        return modes
        
    def _write_parameter_struct(self, ostream):
        """
        Write the structure of the DAP SDSS-style parameter set to a
        file.

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
        Write a DAPPAR entry to the SDSS-style parameter file.

        Args:
            ostream (_io.TextIOWrapper) : File stream for output
            index (int) : Index of the entry in :attr:`data` to write
            mode (int) : Mode of the entry; see :func:`_find_modes`

        """
        ostream.write('DAPPAR {0:4d} {1:5d} {2:4s} {3:14.7e} {4:14.7e} {5:14.7e} {6:14.7e}'
                      ' {7:14.7e}\n'.format(self.hdu['DRPC'].data['PLATE'][index],
                                            self.hdu['DRPC'].data['IFUDESIGN'][index], mode,
                                            self.hdu['DRPC'].data['VEL'][index],
                                            self.hdu['DRPC'].data['VDISP'][index],
                                            self.hdu['DRPC'].data['ELL'][index],
                                            self.hdu['DRPC'].data['PA'][index],
                                            self.hdu['DRPC'].data['REFF'][index]))

    # ******************************************************************
    #  User functions
    # ******************************************************************
    def file_name(self):
        """Return the name of the DRP complete database."""
        return ('drpcomplete_{0}.fits'.format(self.drpver))

    def file_path(self):
        """Return the full pat to the DRP complete database."""
        return os.path.join(self.directory_path, self.file_name())

    def update(self, platelist=None, ifudesignlist=None, combinatorics=False, force=False,
               alldrp=False, use_platetargets=False, on_disk=False, quiet=False):
        """
        Update the DRP complete file.
        
        If `platelist` and/or `ifudesignlist` are provided, the existing
        :attr:`platelist` and :attr:`ifudesignlist` attributes are
        replaced.

        If `platelist` and `ifudesignlist` do not have the same length
        or `combinatorics` is True, the lists are expanded to include
        all combinations of their elements.

        If `platelist` and `ifudesignlist` are None or `alldrp` is True,
        all available plates/ifudesigns are collected from within the
        DRP directory structure.

        If the result of :func:`file_path` does not exist, it is
        created.

        If the result of :func:`file_path` does exist, the available
        plates and ifudesigns in the file are compared against the list
        (provided or collected) to update.  If the lists differ, the
        drpcomplete fits file is re-created from scratch.  If all the
        plates and ifudesigns are available, nothing is done, unless
        `force=True`.

        Args:
            platelist (:obj:`str`, :obj:`list`, optional):
                List of plates to include in the drpcomplete fits file.
            ifudesignlist (:obj:`str`, :obj:`list`, optional):
                List of ifudesigns to include in the drpcomplete fits
                file.
            combinatorics (:obj:`bool`, optional):
                Determine all combinations of the entered plates and
                ifudesigns.
            force (:obj:`bool`, optional):
                Overwrite any existing drpcomplete fits file with a new
                one built from scratch.
            alldrp (:obj:`bool`, optional):
                Find the full list of available DRP files.
            use_platetargets (:obj:`bool`, optional):
                Generate the data using the platetargets files instead
                of the DRPall file.
            on_disk (:obj:`bool`, optional):
                When searching for available files to analyze, search
                the DRP directory path instead of using the data in the
                DRPall file.
            quiet (:obj:`bool`, optional):
                Suppress terminal output

        Raises:
            ValueError:
                Raised if drpcomplete fits file was opened in read-only
                mode.

        """
        if self.readonly:
            raise ValueError('drpcomplete fits file was opened as read-only!')

        if platelist is not None:
            self.platelist = arginp_to_list(platelist, evaluate=True)
        if ifudesignlist is not None:
            self.ifudesignlist = arginp_to_list(ifudesignlist, evaluate=True)

        if alldrp:
            self.platelist = None
            self.ifudesignlist = None

        # This *only* searches for CUBE files
        self.platelist, self.ifudesignlist \
                = self._find_completed_reductions(combinatorics=combinatorics, on_disk=on_disk)

        # Check if these plate-ifus are already in the table.
        if not force and self._all_data_exists():
            print('{0} up to date.'.format(self.file_path()))
            return
        else:
            print('Updating {0}.'.format(self.file_path()))

        # Past this point, the drpcomplete fits file will be ovewritten
        # if it already exists.  Only DRP files defined by 
        # self.platelist and self.ifudesignlist will be included; i.e.,
        # this does *not* append to the file but completely overwrites
        # it.

        # Check for the RSS files
        modes = self._find_modes()

        # Notify
        print('Number of DRP files for DRPComplete file: {0}'.format(len(self.platelist)))

        # Get the data
        mangaid, objra, objdec, catid, catindx, trg_version, trg_id, manga_trg1, manga_trg3, vel, \
                veldisp, ell, pa, Reff = (self._match_platetargets(quiet=quiet)
                                            if use_platetargets else self._match_drpall())

        # Apply any corrections to the redshifts using the redshift fix file
        redshift_fix_data = self._read_redshift_fix()
        if redshift_fix_data is not None:
            fix_pltifu = numpy.array(['{0}-{1}'.format(p,i) 
                                    for p,i in zip(redshift_fix_data['DAPZCORR']['plate'],
                                                   redshift_fix_data['DAPZCORR']['ifudesign'])])
            this_pltifu = numpy.array(['{0}-{1}'.format(p,i) 
                                       for p,i in zip(self.platelist,self.ifudesignlist)])
            rows = numpy.array([numpy.where(fix_pltifu == pi)[0][0] if pi in fix_pltifu else -1 \
                                for pi in this_pltifu])
            indx = rows > -1
            rows = rows[indx]
            vel[indx] = redshift_fix_data['DAPZCORR']['z'][rows] \
                            * astropy.constants.c.to('km/s').value

        photometry_fix_data = self._read_photometry_fix()
        if photometry_fix_data is not None:
            fix_pltifu = numpy.array(['{0}-{1}'.format(p,i) 
                                    for p,i in zip(photometry_fix_data['DAPPHOTCORR']['plate'],
                                                   photometry_fix_data['DAPPHOTCORR']['ifudesign'])])
            this_pltifu = numpy.array(['{0}-{1}'.format(p,i) 
                                       for p,i in zip(self.platelist,self.ifudesignlist)])
            rows = numpy.array([numpy.where(fix_pltifu == pi)[0][0] if pi in fix_pltifu else -1 \
                                for pi in this_pltifu])
            indx = rows > -1
            rows = rows[indx]
            ell[indx] = photometry_fix_data['DAPPHOTCORR']['ell'][rows]
            pa[indx] = photometry_fix_data['DAPPHOTCORR']['pa'][rows]
            Reff[indx] = photometry_fix_data['DAPPHOTCORR']['reff'][rows]

        # Write the data to disk
        self.write(self.platelist, self.ifudesignlist, modes, mangaid, objra, objdec, catid,
                   catindx, trg_version, trg_id, manga_trg1, manga_trg3, vel, veldisp, ell, pa,
                   Reff)

    def write(self, platelist, ifudesignlist, modes, mangaid, objra, objdec, catid, catindx,
              trg_version, trg_id, manga_trg1, manga_trg3, vel, veldisp, ell, pa, Reff, drpver=None,
              redux_path=None, dapver=None, analysis_path=None, clobber=True):
        r"""
        Write the drpcomplete fits binary table.

        Header keywords are:
            - VERSDRP: DRP version
            - RDXPTH: DRP reduction path
            - VERSDAP: DAP version
            - SISPTH: DAP analysis path
            - PLTTRG *N*: plateTargets file *N*
            - CATID *N*: ID Number of target catalog *N*
            - AUTHOR: Set to 'K.B. Westfall <kbwestfall@gmail.com>'

        Binary table columns:
            1. PLATE (1J): Plate number
            2. IFUDESIGN (1J): IFU design
            3. MODES (1B): Modes available 
                - MODES=1: only 'CUBE' file is available
                - MODES=2: both 'CUBE' and 'RSS' files are available

            4. MANGAID (*n* A): MaNGA ID (same as 'CATID-CATINDX')
            5. OBJRA (1D): Object right ascension
            6. OBJDEC (1D): Object declination
            7. CATID (1J): Catalog ID used for target selection
            8. CATINDX (1J): Index in catalog with target data
            9. TRG_VERSION (*n* A): Version of the catalog (e.g., NSA)
                used in targetting (def='NULL' if not available)
            10. TRG_ID (1J): Target ID in catalog (def=-9999 if not
                available)
            11. MANGA_TARGET1: Main MaNGA survey target bitmask
            12. MANGA_TARGET3: Ancillary MaNGA survey target bitmask
            13. VEL (1D): Systemic velocity from catalog (def=-9999.0 if
                not available)
            14. VDISP (1D): Velocity dispersion from catalog
                (def=-9999.0 if not available)
            15. ELL (1D): Ellipticity (1-b/a) from catalog (def=-9999.0
                if not available); Elliptical Petrosian value from NSA
            16. PA (1D): Position angle from catalog (def=-9999.0 if not
                available); Elliptical Petrosian value from NSA
            17. REFF (1D): Effective radius from catalog (def=-9999.0 if
                not available); Elliptical Petrosian value from NSA

        Args:
            platelist (list) : List of plates
            ifudesignlist (list) : List of IFU designs
            modes (numpy.array): Mode values, see above
            mangaid (numpy.array): MaNGA IDs
            objra (numpy.array): Object right ascensions
            objdec (numpy.array): Object declinations
            catid (numpy.array): Catalog ID used for target selection
            catindx (numpy.array): Index in catalog with target data
            trg_version (numpy.array): Version of the catalog (e.g.,
                NSA) used in targetting
            trg_id (numpy.array): Target ID in catalog
            manga_trg1 (numpy.array): Main MaNGA survey target bitmask
            manga_trg3 (numpy.array): Ancillary MaNGA survey target
                bitmask
            vel (numpy.array): Velocity from catalog, if available
            veldisp (numpy.array): Velocity dispersion from catalog,
                if available
            ell (numpy.array): Ellipticity (1-b/a) from catalog, if
                available
            pa (numpy.array): Position angle from catalog, if
                available
            Reff (numpy.array): Effective radius from catalog, if
                available
            drpver (str) : (**Optional**) DRP version, see above.
            redux_path (str) : (**Optional**) Path to the top level of
                directory for the DRP output, see above.
            dapver (str) : (**Optional**) DAP version, see above.
            analysis_path (str) : (**Optional**) Path to the top level
                directory for the DAP output files, see above.
            clobber (bool): (**Optional**) Overwrite any existing file.

        Raises:
            AttributeError: Raised if drpcomplete fits file was opened
                in read-only mode.
            FileExistsError: Raised if the drpcomplete file exists and
                clobber=False.
        """
        if self.readonly:
            raise AttributeError('drpcomplete fits file was opened as read-only!')

        out=self.file_path()
        if os.path.isfile(out) and not clobber:
            raise FileExistsError('DRP complete file already exists: {0}'.format(out))

        if not os.path.isdir(self.directory_path):
            os.makedirs(self.directory_path)

        # Create the primary header
        nplttrg = 0 if self.platetargets is None else len(self.platetargets)

        hdr = fits.Header()
        hdr['VERSDRP'] = self.drpver if drpver is None else drpver
        hdr['RDXPTH'] = self.redux_path if redux_path is None else redux_path
        hdr['VERSDAP'] = self.dapver if dapver is None else dapver
        hdr['SISPTH'] = self.analysis_path if analysis_path is None else analysis_path
        for i in range(nplttrg):
            hdr['PLTTRG{0}'.format(i+1)] = self.platetargets[i]
            hdr['CATID{0}'.format(i+1)] = self.catid[i]
        hdr['AUTHOR'] = 'K.B. Westfall <kbwestfall@gmail.com>'

        # Create the Binary Table
        cols = []
        cols.append(fits.Column(name='PLATE', format='1J', array=numpy.array(platelist)))
        cols.append(fits.Column(name='IFUDESIGN', format='1J', array=numpy.array(ifudesignlist)))
        cols.append(fits.Column(name='MODES', format='1B', array=modes))
        cols.append(fits.Column(name='MANGAID', format='{0}A'.format(mangaid.dtype.itemsize),
                           array=mangaid))
        cols.append(fits.Column(name='OBJRA', format='1D', array=objra))
        cols.append(fits.Column(name='OBJDEC', format='1D', array=objdec))
        cols.append(fits.Column(name='CATID', format='1J', array=catid))
        cols.append(fits.Column(name='CATINDX', format='1J', array=catindx))
        cols.append(fits.Column(name='TRG_VERSION',
                                format='{0}A'.format(trg_version.dtype.itemsize),
                                array=trg_version))
        cols.append(fits.Column(name='TRG_ID', format='1J', array=trg_id))
        cols.append(fits.Column(name='MANGA_TARGET1', format='1J', array=manga_trg1))
        cols.append(fits.Column(name='MANGA_TARGET3', format='1J', array=manga_trg3))
        cols.append(fits.Column(name='VEL', format='1D', array=vel))
        cols.append(fits.Column(name='VDISP', format='1D', array=veldisp))
        cols.append(fits.Column(name='ELL', format='1D', array=ell))
        cols.append(fits.Column(name='PA', format='1D', array=pa))
        cols.append(fits.Column(name='REFF', format='1D', array=Reff))

        # Create the HDUList and write it to the fits file
        self.hdu = fits.HDUList([ fits.PrimaryHDU(header=hdr),
                                  fits.BinTableHDU.from_columns(cols, name='DRPC') ])
        print('Writing to disk: {0}'.format(out))
        self.hdu.writeto(out, overwrite=clobber) #clobber=clobber)
        self.nobs = self.hdu['DRPC'].header['NAXIS2']

#    def grab_data(self, plate=None, ifudesign=None, index=None, reread=False):
#        """
#        Return the data for a given plate-ifudesign, or index in the
#        :attr:`data`.  Even though *plate*, *ifudesign*, and *index* are
#        all optional, the arguments must contain either *plate* and
#        *ifudesign* or *index*.
#
#        Args:
#            plate (int): (**Optional**) Plate number
#            ifudesign (int): (**Optional**) IFU design
#            index (int): (**Optional**) Index of the row in :attr:`data`
#                with the data to return
#            reread (bool): (**Optional**) Force the database to be re-read
#
#        Returns:
#            list: List of the 14 elements of :attr:`data`; see
#            :func:`write`.
#
#        Raises:
#            ValueError: Raised if the row with the data is unknown
#                because either *index* is not defined or one or both of
#                *plate* and *ifudesign* is not defined.  Also raised if
#                *index* does not exist in the data array.
#        """
#        if (plate is None or ifudesign is None) and index is None:
#            raise ValueError('Must provide plate and ifudesign or row index!')
#
#        if index is None:
#            index = self.entry_index(plate, ifudesign, reread=reread)
#        else:
#            self._confirm_access(reread=reread)
#            if index >= self.nobs:
#                raise ValueError('Selected row index does not exist')
#
#        # TODO: Can't I just select index, i.e.,
#        # self.hdu[1].data[index]?
#
#        return [ self.hdu['DRPC'].data['PLATE'][index], self.hdu['DRPC'].data['IFUDESIGN'][index],
#                 self.hdu['DRPC'].data['MODES'][index], self.hdu['DRPC'].data['MANGAID'][index],
#                 self.hdu['DRPC'].data['OBJRA'][index], self.hdu['DRPC'].data['OBJDEC'][index],
#                 self.hdu['DRPC'].data['CATID'][index], self.hdu['DRPC'].data['CATINDX'][index],
#                 self.hdu['DRPC'].data['TRG_VERSION'][index], self.hdu['DRPC'].data['TRG_ID'][indx],
#                 self.hdu['DRPC'].data['NSA_V100_ID'][index], self.hdu['DRPC'].data['VEL'][index],
#                 self.hdu['DRPC'].data['VDISP'][index], self.hdu['DRPC'].data['ELL'][index],
#                 self.hdu['DRPC'].data['PA'][index], self.hdu['DRPC'].data['REFF'][index] ]

    def write_par(self, ofile, mode, plate=None, ifudesign=None, index=None, reread=False,
                  clobber=True):
        """
        Write the SDSS-style parameter (Yanny) file for use with the
        MaNGA DAP.

        Args: 
            ofile (str): Output file name
            mode (str): Mode of the DRP file to analyze; must be either
                'RSS' or 'CUBE'
            plate (int): (**Optional**) Plate number
            ifudesign (int): (**Optional**) IFU design
            index (int): (**Optional**) Index of the row in :attr:`data`
                with the data to return
            reread (bool): (**Optional**) Force the database to be re-read
            clobber (bool): (**Optional**) Overwrite any existing parameter
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
        if not self._confirm_access(reread=reread):
            raise IOError('Could not access database!')

        if os.path.exists(ofile) and not clobber:
            raise IOError('Parameter file already exists.  Set clobber=True to overwrite.')

        if (plate is None or ifudesign is None) and index is None:
            raise ValueError('Must provide plate and ifudesign or row index!')

        if mode != 'CUBE' and mode != 'RSS':
            raise ValueError('Mode must be either CUBE or RSS.')

        if index is None:
            index = self.entry_index(plate, ifudesign, reread=reread)
        else:
            if index >= self.nobs:
                raise ValueError('Selected row index does not exist')

        if mode == 'RSS' and self.hdu['DRPC'].data['MODES'][index] != 2:
            raise ValueError('RSS mode not available for plate={0},ifudesign={1} !'.format(
                                                        self.hdu['DRPC'].data['PLATE'][index],
                                                        self.hdu['DRPC'].data['IFUDESIGN'][index]))

        # Write the SDSS parameter file
        ostream = open(ofile, 'w')
        ostream.write('\n')
        self._write_parameter_struct(ostream)
        ostream.write('\n')
        ostream.write('\n')
        self._write_parameter_list(ostream, index, mode)
        ostream.write('\n')
        ostream.close()

    def write_config(self, ofile, plate=None, ifudesign=None, index=None, sres_ext=None,
                     sres_fill=None, covar_ext=None, reread=False, overwrite=True):
        """
        Write a config file with the data used to instantiate a
        :class:`mangadap.datacube.manga.MaNGADataCube` datacube for
        analysis.

        Args: 
            ofile (:obj:`str`):
                Output file name.
            plate (:obj:`int`, optional):
                Plate number.
            ifudesign (:obj:`int`, optional):
                IFU design.
            index (:obj:`int`, optional):
                Index of the row in :attr:`data` with the data to
                return.
            sres_ext (:obj:`str`, optional):
                The extension to use when constructing the spectral
                resolution vectors for the MaNGA datacubes. See
                :func:`mangadap.datacube.manga.MaNGADataCube.spectral_resolution`.
            sres_fill (:obj:`bool`, optional):
                Fill masked spectral-resolution data by simple linear
                interpolation.
            covar_ext (:obj:`str`, optional):
                Extension in the MaNGA DRP CUBE file to use as the
                single spatial correlation matrix for all wavelength
                channels.
            reread (:obj:`bool`, optional):
                Force the database to be re-read
            overwrite (:obj:`bool`, optional):
                Overwrite any existing parameter file

        Raises:
            IOError:
                Raised if the parameter file already exists and
                clobber is False.
            ValueError:
                Raised if
                    - the row with the data is unknown because either
                      ``index`` is not defined or one or both of
                      ``plate`` and ``ifudesign`` is not defined.
                    - ``index`` does not exist in the data array.

        """
        if not self._confirm_access(reread=reread):
            raise IOError('Could not access database!')

        if (plate is None or ifudesign is None) and index is None:
            raise ValueError('Must provide plate and ifudesign or row index!')

        if index is None:
            index = self.entry_index(plate, ifudesign, reread=reread)
        elif index >= self.nobs:
            raise ValueError('Selected row index does not exist')

        MaNGADataCube.write_config(ofile, self['PLATE'][index], self['IFUDESIGN'][index], log=True,
                                   z=self['VEL'][index]/astropy.constants.c.to('km/s').value,
                                   vdisp=self['VDISP'][index], ell=self['ELL'][index],
                                   pa=self['PA'][index], reff=self['REFF'][index],
                                   sres_ext=sres_ext, sres_fill=sres_fill, covar_ext=covar_ext,
                                   drpver=self.drpver, redux_path=self.redux_path,
                                   overwrite=overwrite)

    def entry_index(self, plate, ifudesign, reread=False):
        """
        Find the index of the row with the parameter data for the
        specified plate and ifudesign.

        .. warning::
            - This is very inefficient if you're looking for multiple
              entries...

        Args:
            plate (:obj:`int`):
                Plate number
            ifudesign (:obj:`int`):
                IFU design
            reread (:obj:`bool`, optional):
                Force the database to be re-read

        Returns:
            int: Index of the row in :attr:`data` with the data for the
            given *plate* and *ifudesign*

        Raises:
            ValueError:
                Raised if the given `plate` and `ifudesign` were not
                found.
        """
        if not self._confirm_access(reread=reread):
            raise IOError('Could not access database!')
        indx = (self['PLATE'] == plate) & (self['IFUDESIGN'] == ifudesign)
        if numpy.sum(indx) == 0:
            raise ValueError('Could not find plate={0},ifudesign={1} in drpcomplete file!'.format(
                                plate, ifudesign))
        return numpy.where(indx)[0][0]

    def can_analyze(self, row=None):
        """
        Determine if the DAP can analyze a plate-ifu entry in the
        database.

        The selection is currently:

            - MaNGAID != 'NULL'
            - MANGA_TARGET1 > 0 or MANGA_TARGET3 > 0
            - VEL > -500

        Args:
            row (:obj:`int`, optional):
                The specific row to test.  By default, return a boolean
                vector for all the database rows.

        Returns:
            Either a single boolean or boolean `numpy.ndarray`_
            flagging that DAP can (True) or cannot (False) analyze
            the data associated with the database entry (or entries).
        """
        if row is None:
            return (self['MANGAID'] != 'NULL') \
                    & ((self['MANGA_TARGET1'] > 0) | (self['MANGA_TARGET3'] > 0)) \
                    & (self['VEL'] > -500.0)
        return self['MANGAID'][row] != 'NULL' \
                and (self['MANGA_TARGET1'][row] > 0 or self['MANGA_TARGET3'][row] > 0) \
                and self['VEL'][row] > -500.0

