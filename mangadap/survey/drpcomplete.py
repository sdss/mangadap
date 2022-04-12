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
from pathlib import Path
import os
import warnings

from IPython import embed

import numpy

from astropy.io import fits
import astropy.constants

from pydl.pydlutils.yanny import yanny

from ..datacube import MaNGADataCube
from ..spectra import MaNGARSS
from ..config import manga
from ..util.parser import arginp_to_list

class DRPComplete:
    r"""
    Database with information needed by the DAP to analyze the completed DRP
    files.

    This class searches the defined paths for files resulting from the 3D phase
    of the MaNGA DRP, and then collates the information necessary to run those
    files through the MaNGA DAP.  The necessary metadata are pulled from the
    DRPall file; see :func:`update`.

    Args:
        platelist (:obj:`str`, :obj:`list`, optional):
            List of plates to search for.  Default is to search the full DRP
            path.
        ifudesignlist (:obj:`str`, :obj:`list`, optional):
            List of IFUs to search for.  Default is to search the full DRP path.
        drpall (:obj:`str`, optional):
            The full path to the DRPall fits file.  Default is set by
            :func:`mangadap.config.manga.drpall_file`.
        drpver (:obj:`str`, optional):
            DRP version, which is:

                - used to define the default DRP redux path
                - used in the name of the drpcomplete fits file
                - included as a header keyword in the output file

            Default is defined by :func:`mangadap.config.manga.drp_version`.
        redux_path (:obj:`str`, optional): 
            The path to the top-level directory containing the DRP output files.
            Default is defined by :func:`mangadap.config.manga.drp_redux_path`.
        dapver (:obj:`str`, optional):
            DAP version, which is:

                - used to define the default DAP analysis path
                - included as a header keyword in the output drpcomplete
                  fits file

            Default is defined by :func:`mangadap.config.manga.dap_version`.
        analysis_path (:obj:`str`, optional):
            The path to the top level directory for the DAP output files.
            Default is defined by
            :func:`mangadap.config.manga.dap_analysis_path`
        directory_path (:obj:`str`, optional):
            Direct path for the output file.  Default is a directory called
            ``common`` within the ``analysis_path`` directory.
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
        drpver (:obj:`str`):
            DRP version, see above.
        redux_path (`Path`_):
            Path to the top level of directory for the DRP output, see
            above.
        dapver (:obj:`str`):
            DAP version, see above.
        analysis_path (`Path`_):
            Path to the top level directory for the DAP output files,
            see above.
        directory_path (`Path`_):
            Direct path for the output file.
        readonly (:obj:`bool`):
            Flag that the drpcomplete fits file is only opened for
            reading, not for updating.
        hdu (`astropy.io.fits.HDUList`_):
            Fits data with binary table data.
        nobs (:obj:`int`):
            Number of observations in the file

    """
    def __init__(self, platelist=None, ifudesignlist=None, drpall=None, drpver=None,
                 redux_path=None, dapver=None, analysis_path=None, directory_path=None,
                 readonly=False):

        # Input properties
        self.drpver = manga.drp_version() if drpver is None else str(drpver)
        self.redux_path = manga.drp_redux_path(self.drpver) \
                                if redux_path is None else Path(redux_path).resolve()

        self.dapver = manga.dap_version() if dapver is None else str(dapver)
        self.analysis_path = manga.dap_analysis_path(self.drpver, self.dapver) \
                             if analysis_path is None else Path(analysis_path).resolve()
        self.directory_path = self.analysis_path / 'common' if directory_path is None \
                                    else Path(directory_path).resolve()

        self.hdu = None
        self.nobs = None
        if self.file_path.exists():
            print('READING')
            self._read()

        if readonly:
            self.readonly=True
            self.platelist = None
            self.ifudesignlist = None
            return

        self.readonly = False

        self.platelist = arginp_to_list(platelist, evaluate=True)
        self.ifudesignlist = arginp_to_list(ifudesignlist, evaluate=True)

        self.drpall = manga.drpall_file(drpver=self.drpver, redux_path=self.redux_path) \
                            if drpall is None else drpall
        
    def __getitem__(self, key):
        return self.hdu['DRPC'].data[key]

    def __len__(self):
        return self.nobs

    # ******************************************************************
    #  Utility functions
    # ******************************************************************
    def _read_fix_data(self, fix_file):
        """
        Read "fix" data from the provided file.

        Fix data are generally corrections to the metadata provides by the NSA
        for specific MaNGA observations.

        Args:
            fix_file (:obj:`str`):
                SDSS parameter file with the fix data.

        Returns:
            yanny: A yanny structure with the data fixes.  Returns
            None and raises a warning if the file does not exist.
        """
        if not fix_file.exists():
            warnings.warn(f'DAP "fix" file does not exist: {fix_file}')
            return None
        return yanny(filename=str(fix_file))

    def _read_redshift_fix(self):
        """
        Wrapper for :func:`_read_fix_data` and
        :func:`~mangadap.config.defaults.redshift_fix_file` that
        returns the redshift fix data.
        """
        return self._read_fix_data(manga.redshift_fix_file())

    def _read_photometry_fix(self):
        """
        Wrapper for :func:`_read_fix_data` and
        :func:`~mangadap.config.defaults.photometry_fix_file` that
        returns the photometry fix data.
        """
        return self._read_fix_data(manga.photometry_fix_file())

    def _match_drpall(self):
        """
        Find the data needed for the DRP complete database in the DRPall
        file.

        If a plate-ifudesign combination is not found in the DRPall
        file, the other parameter values are set to -9999 and the MaNGA
        ID is set to NULL.

        If the plate-ifudesign combination is found, the columns in the
        plateTargets files that the code expects to find are ``'plate'``,
        ``'ifudsgn'``, ``'mangaid'``, ``'objra'``, ``'objdec'``; the values that
        will be replaced with ``'NULL'`` (str) or -9999 (int or float) if they
        do not exist are ``'nsa_version'``, ``'nsa_nsaid'``, ``'mngtarg1'``,
        ``'mngtarg3'``, ``'z'``, ``'nsa_elpetro_ba'``, ``'nsa_elpetro_phi'``,
        ``'nsa_elpetro_th50_r'``, ``'nsa_vdisp'``.
        
        Returns:
            :obj:`tuple`: 14 `numpy.ndarray`_ arrays with the MaNGA ID, object
            right ascension, object declination, catalog ID, index of the entry
            in the catalog, catalog version, ID of object in the catalog, main
            MaNGA survey target bitmask, ancillary MaNGA survey target bitmask,
            velocity, velocity dispersion, ellipticity, position angle, and
            effective radius.
        """
        if not self.drpall.exists():
            raise ValueError(f'DRPall file does not exist! {self.drpall}')
        # Open the drpall file
        hdu = fits.open(str(self.drpall))
        # Check the plateifus are unique
        pltifu = hdu[1].data['PLATEIFU']
        if len(numpy.unique(pltifu)) != len(pltifu):
            raise ValueError('The PLATEIFU column in the DRPall file is not unique!')

        # Find the rows in the DRPall file with the correct; assumes the 
        this_pltifu = numpy.array([f'{p}-{i}' for p,i in zip(self.platelist,self.ifudesignlist)])
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
        self.hdu = fits.open(str(self.file_path))
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
        if not self.file_path.exists():
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
                to only select the bundles with 19 or more fibers.
            combinatorics (:obj:`bool`, optional):
                Create :attr:`platelist` and :attr:`ifudesignlist` by
                determining all possible combinations of the input
                values. See above.
            on_disk (:obj:`bool`, optional):
                When searching for available files to analyze, search
                the DRP directory path instead of using the data in the
                DRPall file.

        Returns:
            :obj:`tuple`: Two lists with the available plates and ifudesigns for
            which to collect data.
        """
        # Are the plate and IFU lists matched?
        matchedlist = (self.platelist is not None and self.ifudesignlist is not None \
                       and (len(self.platelist) == len(self.ifudesignlist) and not combinatorics))

        # Get the list of files
        if matchedlist:
            # Lists already matched, just construct the file names
            files = [manga.MaNGAConfig(int(p), int(i), drpver=self.drpver,
                                       redux_path=self.redux_path).file_path
                        for p,i in zip(self.platelist, self.ifudesignlist)]
        elif on_disk:
            # Find the DRP LOGCUBE files on disk
            print('Searching for completed DRP CUBE files...', end='\r')
            files = sorted(list(self.redux_path.glob('*/stack/*-LOGCUBE.fits.gz')))
            print('Searching for completed DRP CUBE files...DONE.')
        else:
            # Use the DRPall file
            drpall_hdu = fits.open(str(self.drpall))
            pltifu = drpall_hdu[1].data['PLATEIFU']
            if len(numpy.unique(pltifu)) != len(pltifu):
                raise ValueError('The PLATEIFU column in the DRPall file is not unique!')

            files = [manga.MaNGAConfig(int(p), int(i), drpver=self.drpver,
                                       redux_path=self.redux_path).file_path
                        for p,i in zip(drpall_hdu[1].data['plate'], drpall_hdu[1].data['ifudsgn'])]

        files = [f for f in files if f.exists()]
        if len(files) == 0:
            raise ValueError('Could not find any files to analyze.')

        # Get the list of plate and ifus
        pltifu = numpy.array([manga.parse_plate_ifu_from_file(f.name) for f in files]).astype(int)

        # Ignore any below the minimum ifusize
        pltifu = pltifu[numpy.logical_not(pltifu[:,1]//100 < mindesign),:]
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
        has_rss = [manga.MaNGAConfig(int(p), int(i), mode='RSS', drpver=self.drpver,
                                       redux_path=self.redux_path).file_path.exists()
                        for p,i in zip(self.platelist, self.ifudesignlist)]
        print('Checking for RSS counterparts...DONE.')
        modes = numpy.ones(len(has_rss), dtype=int)
        modes[has_rss] = 2
        return modes
        
    # ******************************************************************
    #  User functions
    # ******************************************************************
    @property
    def file_name(self):
        """Return the name of the DRP complete database."""
        return f'drpcomplete_{self.drpver}.fits'

    @property
    def file_path(self):
        """Return the full pat to the DRP complete database."""
        return self.directory_path / self.file_name

    def update(self, platelist=None, ifudesignlist=None, combinatorics=False, force=False,
               alldrp=False, on_disk=False, quiet=False):
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
            print(f'{self.file_path} up to date.')
            return
        else:
            print(f'Updating {self.file_path}.')

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
                veldisp, ell, pa, Reff = self._match_drpall()

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
              redux_path=None, dapver=None, analysis_path=None, overwrite=True):
        r"""
        Write the drpcomplete fits binary table.

        Header keywords are:

            - VERSDRP: DRP version
            - RDXPTH: DRP reduction path
            - VERSDAP: DAP version
            - SISPTH: DAP analysis path

        Binary table columns:

            1. PLATE: Plate number
            2. IFUDESIGN: IFU design
            3. MODES: Modes available 

                - MODES=1: only 'CUBE' file is available
                - MODES=2: both 'CUBE' and 'RSS' files are available

            4. MANGAID: MaNGA ID (same as 'CATID-CATINDX')
            5. OBJRA: Object right ascension
            6. OBJDEC: Object declination
            7. CATID: Catalog ID used for target selection
            8. CATINDX: Index in catalog with target data
            9. TRG_VERSION: Version of the catalog (e.g., NSA)
                used in targetting (def='NULL' if not available)
            10. TRG_ID: Target ID in catalog (def=-9999 if not available)
            11. MANGA_TARGET1: Main MaNGA survey target bitmask
            12. MANGA_TARGET3: Ancillary MaNGA survey target bitmask
            13. VEL: Systemic velocity from catalog (def=-9999.0 if
                not available)
            14. VDISP: Velocity dispersion from catalog
                (def=-9999.0 if not available)
            15. ELL: Ellipticity (1-b/a) from catalog (def=-9999.0
                if not available); Elliptical Petrosian value from NSA
            16. PA: Position angle from catalog (def=-9999.0 if not
                available); Elliptical Petrosian value from NSA
            17. REFF: Effective radius from catalog (def=-9999.0 if
                not available); Elliptical Petrosian value from NSA

        Args:
            platelist (:obj:`list`):
                List of plates
            ifudesignlist (:obj:`list`):
                List of IFU designs
            modes (`numpy.ndarray`_):
                Mode values, see above
            mangaid (`numpy.ndarray`_):
                MaNGA IDs
            objra (numpy.ndarray`_):
                Object right ascensions
            objdec (numpy.ndarray`_):
                Object declinations
            catid (numpy.ndarray`_):
                Catalog ID used for target selection
            catindx (numpy.ndarray`_):
                Index in catalog with target data
            trg_version (numpy.ndarray`_):
                Version of the catalog (e.g., NSA) used in targetting
            trg_id (numpy.ndarray`_):
                Target ID in catalog
            manga_trg1 (numpy.ndarray`_):
                Main MaNGA survey target bitmask
            manga_trg3 (numpy.ndarray`_):
                Ancillary MaNGA survey target bitmask
            vel (numpy.ndarray`_):
                Velocity from catalog, if available
            veldisp (numpy.ndarray`_):
                Velocity dispersion from catalog, if available
            ell (numpy.ndarray`_):
                Ellipticity (1-b/a) from catalog, if available
            pa (numpy.ndarray`_):
                Position angle from catalog, if available
            Reff (numpy.ndarray`_):
                Effective radius from catalog, if available
            drpver (:obj:`str`, optional):
                DRP version, see above.
            redux_path (`Path`_, optional):
                Path to the top level of directory for the DRP output, see
                above.
            dapver (:obj:`str`, optional):
                DAP version, see above.
            analysis_path (`Path`_, optional):
                Path to the top level directory for the DAP output files, see
                above.
            overwrite (:obj:`bool`, optional):
                Overwrite any existing file.

        Raises:
            AttributeError:
                Raised if drpcomplete fits file was opened in read-only mode.
            FileExistsError:
                Raised if the drpcomplete file exists and overwrite=False.
        """
        if self.readonly:
            raise AttributeError('drpcomplete fits file was opened as read-only!')

        out = self.file_path
        if out.exists() and not overwrite:
            raise FileExistsError(f'DRP complete file already exists: {out}')

        if not self.directory_path.exists():
            self.directory_path.mkdir(parents=True)

        hdr = fits.Header()
        hdr['VERSDRP'] = self.drpver if drpver is None else drpver
        hdr['RDXPTH'] = str(self.redux_path) if redux_path is None else str(redux_path)
        hdr['VERSDAP'] = self.dapver if dapver is None else dapver
        hdr['SISPTH'] = str(self.analysis_path) if analysis_path is None else str(analysis_path)
        hdr['AUTHOR'] = 'K.B. Westfall <kbwestfall@gmail.com>'

        # TODO: Convert this into a DataTable?
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
        self.hdu = fits.HDUList([fits.PrimaryHDU(header=hdr),
                                 fits.BinTableHDU.from_columns(cols, name='DRPC')])
        print(f'Writing to disk: {out}')
        self.hdu.writeto(str(out), overwrite=overwrite)
        self.nobs = self.hdu['DRPC'].header['NAXIS2']

    def write_config(self, ofile, plate=None, ifudesign=None, index=None, sres_ext=None,
                     sres_fill=None, covar_ext=None, reread=False, overwrite=True):
        """
        Write a config file with the data used to instantiate a
        :class:`mangadap.datacube.manga.MaNGADataCube` datacube for
        analysis.

        Args: 
            ofile (:obj:`str`, `Path`_):
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
                overwrite is False.
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

        manga.MaNGAConfig.write_config(ofile, self['PLATE'][index], self['IFUDESIGN'][index],
                                       log=True,
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
            :obj:`int`: Index of the row in :attr:`data` with the data for the
            given ``plate`` and ``ifudesign``.

        Raises:
            ValueError:
                Raised if the given ``plate`` and ``ifudesign`` were not found.
        """
        if not self._confirm_access(reread=reread):
            raise IOError('Could not access database!')
        indx = (self['PLATE'] == plate) & (self['IFUDESIGN'] == ifudesign)
        if numpy.sum(indx) == 0:
            raise ValueError(f'plate={plate},ifudesign={ifudesign} not in drpcomplete file!')
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

