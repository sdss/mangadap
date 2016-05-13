# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-
"""
Defines a class used to produce and provide an interface for data
required to create the input `SDSS parameter files`_ for the DAP.

The drpcomplete file/data is primarily required for the survey-level
execution of the `DAP at Utah`_.  However, it does provide useful
information regarding the completed DRP files, and it can be used to
create DAP par files.

Preferrably, however, one should get this information from the DRPall,
or forthcoming DAPall, file for a given MPL.

Until further documented here, the description of the `DAP par
file`_ is available in the SDSS-IV/MaNGA `Technical Reference Manual`_.

*License*:
    Copyright (c) 2015, SDSS-IV/MaNGA Pipeline Group
        Licensed under BSD 3-clause license - see LICENSE.rst

*Source location*:
    $MANGADAP_DIR/python/mangadap/survey/drpcomplete.py

*Imports and python version compliance*:
    ::

        from __future__ import division
        from __future__ import print_function
        from __future__ import absolute_import
        from __future__ import unicode_literals

        import sys
        if sys.version > '3':
            long = int

        import os
        import numpy
        import warnings
        from astropy.io import fits
        import astropy.constants
        from pydl.pydlutils.yanny import yanny

        from ..config import defaults
        from . import drpfits
        from ..util.parser import arginp_to_list, list_to_csl_string, parse_drp_file_name
        from ..util.exception_tools import print_frame

*Class usage examples*:
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
        :class:`mangadap.drpfits` and :class:`mangadap.dapfile` classes,
        require a change to the function calls. Changed parameter order
        in :func:`write` and allowed for use of attributes as default.
    | **20 May 2015**: (KBW) Documentation and Sphinx tests. Prep for
        conversion from NSA to TRG nomenclatature, but still need to
        make the full conversion.
    | **21 Aug 2015**: (KBW) Major revisions: File moved from main
        module to survey submodule; now only reads plateTargets files,
        which have the NSA IDs in them; changed the number of arrays
        returned by :func:`DRPComplete._match_platetargets` and the
        default values; removed def_veldisp and use_trg from
        :func:`DRPComplete.update`; added target-catalog version, and
        NSAID v1_0_0 to the output fits file; removed options from
        :func:`DRPComplete.write` related to the target catalogs and
        platetargets files (these should be defined by the object);
    | **27 Aug 2015**: (KBW) Include MANGA_TARGET1 and MANGA_TARGET3 in
        data structure, pulled from plateTargets files; NSAID v1_0_0
        removed; changed from returning Sersic parameters to elliptical
        Petrosian parameters.
    | **06 Oct 2015**: (KBW) Changed to reading 'object_ra' and
        'object_dec', instead of target counter parts due to changes in
        MaNGA core
    | **17 Feb 2016**: (KBW) Converted the name of the class to
        DRPComplete
    | **29 Feb 2016**: (KBW) Import drpfits, not drpfile
    | **13 May 2016**: Switch to using `pydl.pydlutils.yanny`_ instead
        of internal yanny reader.  Incorporated changes to plateTargets
        column names defined for DR13.
    
.. _DAP par file: https://trac.sdss.org/wiki/MANGA/TRM/TRM_ActiveDev/dap/Summary/parFile
.. _Technical Reference Manual: https://trac.sdss.org/wiki/MANGA/TRM/TRM_ActiveDev
.. _SDSS parameter files: http://www.sdss.org/dr12/software/par/
.. _DAP at Utah: https://trac.sdss.org/wiki/MANGA/TRM/TRM_ActiveDev/dap/Summary/Utah
.. _pydl.pydlutils.yanny: http://pydl.readthedocs.io/en/stable/api/pydl.pydlutils.yanny.yanny.html

"""

from __future__ import division
from __future__ import print_function
from __future__ import absolute_import
from __future__ import unicode_literals

import sys
if sys.version > '3':
    long = int

import os
import numpy
import warnings
from astropy.io import fits
import astropy.constants
from pydl.pydlutils.yanny import yanny

from ..config import defaults
from .. import drpfits
from ..util.parser import arginp_to_list, list_to_csl_string, parse_drp_file_name
from ..util.exception_tools import print_frame

__author__ = 'Kyle Westfall'

class DRPComplete:
    """Find DRP files ready for analysis and write parameter files.

    This class searches the defined paths for files resulting from the
    3D phase of the MaNGA DRP, and then collates the information
    necessary to run those files through the MaNGA DAP.  The necessary
    parameters are pulled from the provided platetargets files; see
    :func:`update`.

    Args:
        platelist (str or list): (Optional) List of plates to search
            for.  Default is to search the full DRP path.
        ifudesignlist (str or list): (Optional) List of ifudesign to
            search for.  Default is to search the full DRP path.
        platetargets (str or list): (Optional) List of platetargets
            files to search through to find any given plate ifudesign
            combination.  Default is returned as the first element in
            :func:`mangadap.config.defaults.default_plate_target_files`.
        catid (str or list): (Optional) List of target catalog ID
            numbers.  Default is returned as the second element in
            :func:`mangadap.config.defaults.default_plate_target_files`.
        drpver (str): (Optional) DRP version, which is:
                - used to define the default DRP redux path
                - used when declaring a drpfits instance
                - used in the name of the drpcomplete fits file
                - included as a header keyword in the output file
            Default is defined by
            :func:`mangadap.config.defaults.default_drp_version`

        redux_path (str): (Optional) The path to the top level directory
            containing the DRP output files; this is the same as the
            *redux_path* in the :class:`mangadap.drpfits.DRPFits` class.
            Default is defined by
            :func:`mangadap.config.defaults.default_redux_path`.
        dapver (str): (Optional) DAP version, which is:
                - used to define the default DAP analysis path
                - included as a header keyword in the output drpcomplete
                  fits file

            Default is defined by
            :func:`mangadap.config.defaults.default_dap_version`
        analysis_path (str): The path to the top level directory for the
            DAP output files; this is **different** from the
            directory_path in the :class:`mangadap.dapfile` class.  Default is
            defined by
            :func:`mangadap.config.defaults.default_analysis_path`
        readonly (bool) : Flag that the drpcomplete fits file is only
            opened for reading, not for updating.

    Raises:
        Exception: Raised if user supplies only one of platetargets or
            catid instead of both.

    Attributes:
        platelist (list) : List of plates to search for, see above.
        ifudesignlist (list) : List of IFU designs to search for, see
            above.
        platetargets (numpy.ndarray) : List of platetargets files to
            search through, see above
        catid (numpy.ndarray) : List of target catalog IDs, see above
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
    def __init__(self, platelist=None, ifudesignlist=None, platetargets=None, catid=None,
                 drpver=None, redux_path=None, dapver=None, analysis_path=None, readonly=False):
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
            self.catid = None
            return

        self.readonly = False

        self.platelist = arginp_to_list(platelist, evaluate=True)
        self.ifudesignlist = arginp_to_list(ifudesignlist, evaluate=True)

        if (platetargets is not None and catid is None) or \
           (platetargets is None and catid is not None):
            raise ValueError('To use user-provided platetargets files, must provide both '
                             'platetargets and catid.')

        if platetargets is not None:
            self.platetargets = numpy.array( arginp_to_list(platetargets) )
            self.catid = numpy.array( arginp_to_list(catid) ).astype(numpy.int)
        else:
            self.platetargets, self.catid = defaults.default_plate_target_files()


    # ******************************************************************
    #  Utility functions
    # ******************************************************************
    def _drp_mangaid(self, drplist):
        """
        Grab the MaNGA IDs from the DRP fits files.
        
        Args:
            drplist (list) : List of :class:`mangadap.drpfits.DRPFits`
                objects

        Returns:
            numpy.array: Array with the MaNGA IDs from the headers of
            the DRP fits files.

        .. note::
            This takes far too long; either astropy.io.fits is reading
            the entire file instead of just the header, or the slowdown
            is because the DRP files are compressed.

        """
        nn = len(drplist)
        mangaid = []
        print('Gathering MANGA-IDs for DRP files...', end='\r')
        for i in range(0,nn):
            mangaid_, ra, dec = drplist[i].object_data()
            mangaid = mangaid + [mangaid_]
        print('Gathering MANGA-IDs for DRP files...DONE')
        mangaid = numpy.array(mangaid)
        return mangaid


    def _drp_info(self, drplist):
        """
        Grab the MaNGA IDs and object coordinates from the DRP fits
        files.
        
        Args:
            drplist (list) : List of :class:`mangadap.drpfits.DRPFits`
                objects

        Returns:
            numpy.array: Three arrays with, respectively, the MaNGA IDs
            from the headers of the DRP fits files, the object right
            ascension, and the object declination.

        .. note::
            This takes far too long; either astropy.io.fits is reading
            the entire file instead of just the header, or the slowdown
            is because the DRP files are compressed.

        """
        nn = len(drplist)
        mangaid = []
        objra = numpy.empty(nn, dtype=numpy.float64)
        objdec = numpy.empty(nn, dtype=numpy.float64)
        print('Gathering DRP header data...', end='\r')
        for i in range(0,nn):
            mangaid_, objra[i], objdec[i] = drplist[i].object_data()
            mangaid = mangaid + [mangaid_]
        print('Gathering DRP header data...DONE')
        mangaid = numpy.array(mangaid)
        return mangaid, objra, objdec


    def _read_platetargets(self, quiet=True):
        """
        Read all the platetargets files using `pydl.pydlutils.yanny`_
        and return a list of yanny structures.

        Args:
            quiet (bool): (Optional) Suppress terminal output (NOT USED)

        Returns:
            list: A list of yanny structures, one per platetargets file.

        Raises:
            FileNotFoundError: Raised if cannot open one or more
                platetargets files.

        .. todo::
            This should be made more efficient by collating the required
            plateTargets data into a single record array!
        """


#        from matplotlib import pyplot
#        print(self.platetargets[0])
#        plttrg_data = yanny(filename=self.platetargets[0])
#        pyplot.scatter(plttrg_data['PLTTRGT']['nsa_elpetro_phi'],
#                       plttrg_data['PLTTRGT']['nsa_sersic_phi'], marker='.', color='k', s=30)
#        pyplot.show()
#        pyplot.scatter(plttrg_data['PLTTRGT']['nsa_elpetro_ba'],
#                       plttrg_data['PLTTRGT']['nsa_sersic_ba'], marker='.', color='k', s=30)
#        pyplot.show()
#        exit()

        plttrg_data = []

        for p in self.platetargets:
            print('Reading plateTargets file: {0}'.format(p), end='\r')
            # Check that the file exists
            if not os.path.exists(p):
                raise FileNotFoundError('Cannot open {0}!'.format(p))
            # Read and append the data
            plttrg_data.append( yanny(filename=p) )
        print('\nReading plateTargets file: DONE')

        return plttrg_data


    def _match_platetargets(self, quiet=True):
        """
        Read the platetargets files and match the data therein to the
        completed DRP files based on the plate and ifudesign.

        If a plate-ifudesign combination is not found in the
        plateTargets file, the other parameter values are set to -9999
        and the MaNGA ID is set to NULL.

        If the plate-ifudesign combination is found, the columns in the
        plateTargets files that the code expects to find are 'plate',
        'ifudesign', 'mangaid', 'object_ra', and 'object_dec'; the
        values that will be replaced with 'NULL' (str) or -9999 (int or
        float) if they do not exist are 'nsa_version', 'nsa_id',
        'manga_target1', 'manga_target3', 'nsa_redshift',
        'nsa_sersic_ba', 'nsa_sersic_phi', 'nsa_sersic_th50',
        'nsa_vdisp'.
        
        .. todo::

            - Instead of searching through all files for the correct
              plate-ifudesign, use the MaNGA ID to get the correct
              catalog?  Requires getting the MaNGA ID from somewhere,
              but it takes a long time to read from the fits files...
            - Is there some more concise method of performing the same
              thing as what is done below with many try/except blocks?

        Args:
            quiet (bool): (Optional) Suppress terminal output

        Returns:
            numpy.array: 14 arrays with: MaNGA ID, object right
            ascension, object declination, catalog ID, index of the
            entry in the catalog, catalog version, ID of object in the
            catalog, main MaNGA survey target bitmask, ancillary MaNGA
            survey target bitmask, velocity, velocity dispersion,
            ellipticity, position angle, and effective radius.
        """
        # Read the platetargets file
        plttrg_data = self._read_platetargets(quiet=quiet)
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
                    raise Exception('Multiple instances of {0}-{1} (MaNGA ID={2}) in {3}'.format(
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
            catid[i] = numpy.int32(mangaid[i].split('-')[0])
            if catid[i] != self.catid[plttrg_j]:
                warnings.warn('{0}-{1} (MaNGA ID={2}) found in {3} (CAT={4})!'.format(
                                self.platelist[i], self.ifudesignlist[i], mangaid_i,
                                self.platetargets[j], self.catid[j]))
            catindx[i] = numpy.int32(mangaid[i].split('-')[1])

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


    def _all_data_exists(self, quiet=True):
        """
        Determines if the data for all the plates/ifudesigns selected in
        the current compilation is already present in the current
        drpcomplete fits file.

        Args:
            quiet (bool): (Optional) Suppress output

        Returns:
            bool: Flag that all data has already been collated for the
            requested plate/ifudesign list.
        """
        nn = len(self.platelist)
        print('Checking for plate-ifudesign entries in existing file ...', end='\r')
        for i in range(0,nn):
            try:
                index = self.entry_index(self.platelist[i], self.ifudesignlist[i])
            except Exception as e:
                if not quiet:
                    print_frame('Exception')
                    print(e)
                print('Checking for plate-ifudesign entries in existing file ... DONE')
                return False
        print('Checking for plate-ifudesign entries in existing file ... DONE')
        return True


    def _read_data(self):
        """Read the data in the existing file at :func:`file_path`."""
        hdu = fits.open(self.file_path())
        self.data = hdu[1].data
        self.nrows = hdu[1].header['NAXIS2']
        hdu.close()
        print('Read data: %d rows' % self.nrows)


    def _confirm_access(self, reread):
        """
        Check the drpcomplete fits file at :func:`file_path` is
        accessible, and read the data if not yet read.

        Args:
            reread (bool): Force the file to be re-read

        Raises:
            FileNotFoundError: Raised if the drpcomplete fits file does
                not exist.
        """
        inp = self.file_path()
        if not os.path.exists(inp):
            raise FileNotFoundError('Cannot open file: {0}'.format(inp))
        if self.data is None or reread:
            self._read_data()


    def _find_completed_reductions(self, mindesign=19, combinatorics=False):
        """Search the DRP path for reduced CUBE files.

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
            list: Two lists with the available plates and ifudesigns for
            which to collect data.
        """

        path = self.redux_path
        matchedlist = (self.platelist is not None and self.ifudesignlist is not None \
                       and (len(self.platelist) == len(self.ifudesignlist) and not combinatorics))

        plates = []
        ifudesigns = []
        # Lists already matched, just see if the pairs exist
        if matchedlist:
            n_plates=len(self.platelist)
            for i in range(0,n_plates):
                drpf = drpfits.DRPFits(self.platelist[i], self.ifudesignlist[i], 'CUBE',
                                       drpver=self.drpver)
                if os.path.exists(drpf.file_path()):
                    plates.append(self.platelist[i])
                    ifudesigns.append(self.ifudesignlist[i])
            return plates, ifudesigns

#        print(path)
        print('Searching for completed DRP CUBE files...', end='\r')
        # Otherwise generate the lists
        for root, dir, files in os.walk(path):
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
#                        print('{0}-{1}'.format(p, b))

        print('Searching for completed DRP CUBE files...DONE.')

#        print(plates)
#        print(ifudesigns)
        return plates, ifudesigns


    def _find_modes(self, drplist):
        """
        Using the provided list of CUBE DRP files, find if the RSS mode
        is also available.

        Currently only two modes are possible:

            - modes[i] = 1 -> only 'CUBE' file are available

            - modes[i] = 2 -> both 'CUBE' and 'RSS' files are available

        .. todo::
            - Allow for there to be *only* 'RSS' files.

        Args:
            drplist (list) : List of :class:`mangadap.drpfits.DRPFits`
                objects

        Returns:
            numpy.array: Array of modes for each input DRP file.
        """
        n_drp = len(drplist)
        modes = numpy.empty(n_drp, dtype=numpy.uint8)
        print('Checking for RSS counterparts...', end='\r')
        for i in range(0,n_drp):
            drpf = drpfits.DRPFits(drplist[i].plate, drplist[i].ifudesign, 'RSS',
                                   drpver=drplist[i].drpver)
            modes[i] = 2 if os.path.exists(drpf.file_path()) else 1
        print('Checking for RSS counterparts...DONE.')
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
               alldrp=False):
        """Update the DRP complete file.
        
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
        drpcomplete fits file is re-created from scratch.  If all the
        plates and ifudesigns are available, nothing is done, unless
        force=True.

        Args:
            platelist (str or list): (Optional) List of plates to
                include in the drpcomplete fits file.
            ifudesignlist (str or list): (Optional) List of ifudesigns
                to include in the drpcomplete fits file.
            combinatorics (bool): (Optional) Determine all combinations
                of the entered plates and ifudesigns.
            force (bool): (Optional) Overwrite any existing drpcomplete
                fits file with a new one built from scratch.
            alldrp (bool): (Optional) Find the full list of available
                DRP files.

        Raises:
            AttributeError: Raised if drpcomplete fits file was opened
                in read-only mode.
        """
        if self.readonly:
            raise AttributeError('drpcomplete fits file was opened as read-only!')

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
        # drpfits_list is correct.  After running
        # _find_completed_reductions(), the length of platelist and
        # ifudesignlist should be the same!
        n_plates = len(self.platelist)
        modelist = ['CUBE' for i in range(0,n_plates)]
        drplist = drpfits.drpfits_list(self.platelist, self.ifudesignlist, modelist,
                                       drpver=self.drpver)

        # By running _all_data_exists() the self.data and self.nrows
        # attributes will be populated if they haven't been already!
        if not force and self._all_data_exists():
            print('{0} up to date.'.format(self.file_path()))
            return
        else:
            print('Updating {0}.'.format(self.file_path()))

        # Past this point, the drpcomplete fits file will be ovewritten
        # if it already exists.  Only DRP files created using
        # self.platelist and self.ifudesignlist will be used to create
        # the file, EVEN IF OTHER FILES EXIST.

        modes = self._find_modes(drplist)

        nn = len(drplist)
        print('Number of DRP files: {0}'.format(nn))

        # Use the  plateTargets files to gather the needed data
        mangaid, objra, objdec, catid, catindx, trg_version, trg_id, manga_trg1, manga_trg3, vel, \
            veldisp, ell, pa, Reff = self._match_platetargets()

        # Write the data to disk
        self.write([drplist[i].plate for i in range(0,nn)],
                   [drplist[i].ifudesign for i in range(0,nn)], modes, mangaid, objra, objdec,
                   catid, catindx, trg_version, trg_id, manga_trg1, manga_trg3, vel, veldisp, ell,
                   pa, Reff)

        # Read the data in from disk and save it within the object
        self._read_data()


    def write(self, platelist, ifudesignlist, modes, mangaid, objra, objdec, catid, catindx,
              trg_version, trg_id, manga_trg1, manga_trg3, vel, veldisp, ell, pa, Reff, drpver=None,
              redux_path=None, dapver=None, analysis_path=None, clobber=True):
        """
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
                * MODES=1: only 'CUBE' file is available
                * MODES=2: both 'CUBE' and 'RSS' files are available

            4. MANGAID (*n* A): MaNGA ID (same as 'CATID-CATINDX')
            5. OBJRA (1D): Object right ascension
            6. OBJDEC (1D): Object declination
            7. CATID (1J): Catalog ID used for target selection
            8. CATINDX (1J): Index in catalog with target data
            9. TRG_VERSION (*n*A): Version of the catalog (e.g., NSA)
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
            drpver (str) : (Optional) DRP version, see above.
            redux_path (str) : (Optional) Path to the top level of
                directory for the DRP output, see above.
            dapver (str) : (Optional) DAP version, see above.
            analysis_path (str) : (Optional) Path to the top level
                directory for the DAP output files, see above.
            clobber (bool): (Optional) Overwrite any existing file.

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

        if not os.path.isdir(self.analysis_path):
            os.makedirs(self.analysis_path)

        # Create the primary header
        nplttrg = len(self.platetargets)

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
        hdulist = fits.HDUList([ fits.PrimaryHDU(header=hdr), fits.BinTableHDU.from_columns(cols) ])

        # TODO: Make directory if it doesn't exist?
        print('Writing to disk: {0}'.format(out))
        hdulist.writeto(out, clobber=clobber)


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
            list: List of the 14 elements of :attr:`data`; see
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
                 self.data['TRG_VERSION'][index], self.data['TRG_ID'][indx],
                 self.data['NSA_V100_ID'][index], self.data['VEL'][index],
                 self.data['VDISP'][index], self.data['ELL'][index], self.data['PA'][index],
                 self.data['REFF'][index] ]


    def write_par(self, ofile, mode, plate=None, ifudesign=None, index=None, reread=False,
                  clobber=True):
        """
        Write the SDSS-style parameter (Yanny) file for use with the
        MaNGA DAP.

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
        Find the index of the row with the parameter data for the
        specified plate and ifudesign.

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

        for i in range(self.nrows):
            if self.data['PLATE'][i] == plate and self.data['IFUDESIGN'][i] == ifudesign:
                return i

        raise Exception('Could not find plate={0},ifudesign={1} in drpcomplete file!'.format(plate,
                        ifudesign))



