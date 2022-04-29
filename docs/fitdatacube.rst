
.. include:: include/links.rst

.. _fitdatacube:

How to Fit a non-MaNGA Datacube
===============================

After the initial ingestion of the MaNGA datacube, the main
:ref:`execution-mangadap` is effectively independent of the specific data
format.  Although all algorithms were fine-tuned using MaNGA datacubes, the
pipeline can be executed on data from other instruments.  This document details
how this is done and gives an example using data from Keck's KCWI.

Overview
--------

The critical components that the DAP needs to perform its analysis modules are:

    #. wavelengths are in vacuum,
    #. spectra that are geometrically binned in wavelength,
    #. spectral errors, pixel masks, and an estimate of the spectral resolution,
    #. an estimate of the redshift that is accurate to within :math:`\pm 2000` km/s,
    #. conventions for how to identify the datacube being analyzed,
    #. a plan for the analyses to be performed, and
    #. a naming convention for the DAP output.

The first four components are handled by the datacube container class,
:class:`~mangadap.datacube.datacube.DataCube`.  Therefore, to run the DAP on a
custom datacube, you need to build a custom
:class:`~mangadap.datacube.datacube.DataCube` subclass.

The last two components are handled by the parameter class
:class:`~mangadap.config.analysisplan.AnalysisPlanSet`.  Depending on the
desired output path structure, it is likely that DAP analysis of a custom
datacube can use this class directly, instead of having to create a new
subclass.

.. _datacube_subclass:

Creating a new DataCube subclass
--------------------------------

The primary function of the :class:`~mangadap.datacube.datacube.DataCube` is to
provide an interface for the DAP to read the data to be analyzed and provide
access to it in a way that the DAP understands.  All datacubes analyzed by the
DAP must subclass from this object.

To work directly with the ``manga_dap`` command-line script, the subclass must
instantiate the datacube object either directly from the datacube file (as its
``__init__`` method) or via a ``classmethod`` that reads a "configuration" file,
where the class method is called ``from_config``.  In the example below, we
focus on the former; however, the limitation is that the metadata required to
run the DAP (the redshift estimate) must be read from the file or hard-coded.
For reference, MaNGA datacubes analyzed at the survey level use configuration
files; see :func:`~mangadap.datacube.manga.MaNGADataCube.from_config`.

In its simplest implementation, a custom datacube subclass would look like this:

.. code-block:: python

    from pathlib import Path

    from astropy.io import fits
    from astropy.wcs import WCS

    from mangadap.datacube.datacube import DataCube

    class CustomDataCube(DataCube):
        r"""
        Container class for a custom datacube.

        Args:
            ifile (:obj:`str`):
                The name of the file to read.
        """

        instrument = 'custom'
        """
        Set the name of the instrument.  This is used to construct the
        output file names.
        """

        def __init__(self, ifile):

            _ifile = Path(ifile).resolve()
            if not _ifile.exists():
                raise FileNotFoundError(f'File does not exist: {_ifile}')

            # Set the paths
            self.directory_path = _ifile.parent
            self.file_name = _ifile.name

            # Collect the metadata into a dictionary
            meta = {}

            # Open the file and initialize the DataCube base class
            with fits.open(str(_ifile)) as hdu:
                print('Reading datacube ...', end='\r')
                prihdr = hdu[0].header
                wcs = WCS(header=prihdr, fix=True)
                ...
            print('Reading datacube ... DONE')

            # Instantiate the base class
            super().__init__(flux, ivar=ivar, mask=mask, sres=sres, wave=wave, meta=meta,
                             prihdr=prihdr, wcs=wcs, name=_ifile.name.split('_')[0]) 


In terms of the requirements above:

    #. Geometrically binned spectra: The binning is not checked by the DAP
       immediately when the code executed (yet).  Instead, the code will start
       and likely crash as soon as it tries to construct the template library to
       be fit.  See :ref:`resampling` and the KCWI example below for how to use
       DAP code to resample a datacube.  **NOTE**: The DAP expects datacube
       arrays to be oriented with the two spatial axes first and the spectral
       data along the last axis; i.e., the flux array should be ordered
       :math:`(x,y,$\lambda$)`.  This can be done internal to the instantiation
       by the base class; see the :class:`~mangadap.datacube.datacube.DataCube`
       ``axes`` argument.

    #. Ancillary data: The script provides spectral errors in the ``ivar``
       array, a pixel mask in the ``mask`` array, and an estimate of the
       spectral resolution in the ``sres`` array (this can be a single number if
       the resolution is constant as a function of wavelength and position).  It
       also provides the primary header, which is simply copied to output files,
       and the WCS read from the fits file (assuming it can be constructed from
       the primary header).

    #. Redshift estimate: The example above *does not* explicitly define the
       redshift.  This is assumed to be done within the ``...`` block.  The
       redshift must be included by defining ``meta['z']``.  The code will crash
       if this is not defined as soon as the DAP attempts to analyze the
       datacube.  All other metadata that can be used by the DAP is optional;
       see :func:`~mangadap.datacube.datacube.DataCube.populate_metadata`.
    
    #. Naming conventions: The code establishes the naming convention for output
       files using the defined ``instrument`` and ``name`` for each datacube.
       In this example, the ``instrument`` is defined as a class attribute, and
       the ``name`` is an instance attribute defined in the instantiation of the
       base class based on the input file name.  The actual value of these
       strings is free form, but at least one of them must be defined.  The name
       of the output files is then ``{instrument}-{name}``, if both are defined,
       and just the ``instrument`` or ``name`` string if only one is defined.

.. _fitdatacube-kcwi:

Example: KCWI
+++++++++++++

As an example, here's a fully implemented
:class:`~mangadap.datacube.datacube.DataCube` subclass for use with KCWI data,
where the data was reduced by the the Keck KCWI DRP (see `here
<https://kcwi-drp.readthedocs.io/en/latest/>`__).

.. code-block:: python

    from pathlib import Path

    import numpy

    from astropy.io import fits
    from astropy.wcs import WCS

    from mangadap.datacube.datacube import DataCube
    from mangadap.util.sampling import Resample, angstroms_per_pixel

    class KCWIDataCube(DataCube):
        r"""
        Container class for KCWI datacubes.

        Args:
            ifile (:obj:`str`):
                The name of the file to read.
        """

        instrument = 'kcwi'
        """
        Set the name of the instrument.  This is used to construct the
        output file names.
        """

        def __init__(self, ifile):

            _ifile = Path(ifile).resolve()
            if not _ifile.exists():
                raise FileNotFoundError(f'File does not exist: {_ifile}')

            # Set the paths
            self.directory_path = _ifile.parent
            self.file_name = _ifile.name

            # Collect the metadata into a dictionary
            meta = {}
            meta['z'] = 0.028       # Specific to the target observed
            sres = 1800             # Specific to this KCWI instrument setting

            # Open the file and initialize the DataCube base class
            with fits.open(str(_ifile)) as hdu:
                print('Reading KCWI datacube data ...', end='\r')
                prihdr = hdu[0].header
                wcs = WCS(header=prihdr, fix=True)
                flux = hdu[0].data
                err = hdu['UNCERT'].data
                mask = hdu['MASK'].data.astype(bool) | numpy.logical_not(err > 0.)
            print('Reading KCWI datacube data ... DONE')

            # Resample to a geometric sampling
            # - Get the wavelength vector
            spatial_shape = flux.shape[1:][::-1]
            nwave = flux.shape[0]
            coo = numpy.array([numpy.ones(nwave), numpy.ones(nwave), numpy.arange(nwave)+1]).T
            wave = wcs.all_pix2world(coo, 1)[:,2]*wcs.wcs.cunit[2].to('angstrom')
            # - Convert the fluxes to flux density
            dw = angstroms_per_pixel(wave, regular=False)
            flux /= dw[:,None,None]
            # - Set the geometric step to the mean value.  This means some
            # pixels will be oversampled and others will be averaged.
            dlogl = numpy.mean(numpy.diff(numpy.log10(wave)))
            # - Resample all the spectra.  Note that the Resample arguments
            # expect the input spectra to be provided in 2D arrays with the
            # last axis as the dispersion axis.
            r = Resample(flux.T.reshape(-1,nwave), e=err.T.reshape(-1,nwave),
                         mask=mask.T.reshape(-1,nwave), x=wave, inLog=False, newRange=wave[[0,-1]],
                         newLog=True, newdx=dlogl)
            # - Reshape and reformat the resampled data in prep for
            # instantiation
            ivar = r.oute.reshape(*spatial_shape,-1)
            mask = r.outf.reshape(*spatial_shape,-1) < 0.8
            ivar[mask] = 0.0
            gpm = numpy.logical_not(mask)
            ivar[gpm] = 1/ivar[gpm]**2
            _sres = numpy.full(ivar.shape, sres, dtype=float)
            flux = r.outy.reshape(*spatial_shape,-1)

            # Default name assumes file names like, e.g., '*_icubew.fits'
            super().__init__(flux, ivar=ivar, mask=mask, sres=_sres,
                             wave=r.outx, meta=meta, prihdr=prihdr, wcs=wcs,
                             name=_ifile.name.split('_')[0]) 

A few specifics of this implementation to note:

    - The redshift and spectral resolution for the data to reduce are
      *hard-coded*.  To be more general, one would rather have a ``from_config``
      classmethod for KCWI where the redshift can be file specific.

    - Much of the instantiation is concerned with resampling the data from a
      linear step in wavelength to a geometric step, which relies on the
      :class:`~mangadap.util.sampling.Resample` class.

    - Note that the ``wcs`` is still provided, unaltered by the change in the
      wavelength coordinates.  This is only possible because the wavelength
      vector (``wave``) is provided to the base-class instantiation directly,
      taking precedence over anything provided by the original WCS.  **Beware**,
      though, this means the input and output WCS may be inconsistent for the
      wavelength axis (GitHub issue raised).

Creating a new AnalysisPlan subclass
------------------------------------

For a detailed description of the components of the
:class:`~mangadap.config.analysisplan.AnalysisPlan` and its associated `toml`_
file, see :ref:`plan`.

The primary function of the :class:`~mangadap.config.analysisplan.AnalysisPlan`
class is to set how the DAP analyzes the data.  However, much of how this is
determined is by using abstracted methods and keywords, meaning that the
implementation of the :class:`~mangadap.config.analysisplan.AnalysisPlan` class
is largely independent of the data being analyzed.  The only exception to this
is the definition of the output paths.  Adopting a top-level output path
referred to as the ``analysis_path``, the default DAP output paths are:

    - ``{analysis_path}/common/``: The path where analysis products "common" to
      many analysis methods are placed.

    - ``{analysis_path}/{plan-key}/``: The path where analysis products specific
      to a single analysis method are placed.

    - ``{analysis_path}/{plan-key}/ref``: The path for DAP reference files

    - ``{analysis_path}/{plan-key}/qa``: The path for DAP quality assurance plots

These defaults are set by
:func:`~mangadap.config.analysisplan.AnalysisPlan.common_path` and
:func:`~mangadap.config.analysisplan.AnalysisPlan.method_path`.  The root of the
output file names is set by 
:func:`~mangadap.config.analysisplan.AnalysisPlan.dap_file_root`.  In most
cases, you will likely not want to change these.  However, the survey-level
execution of the DAP for the MaNGA Survey did change these to allow the files to
be organized by the plate and IFU number of the analyzed datacube.  See
:class:`~mangadap.config.manga.MaNGAAnalysisPlan` for an example of how to
create a subclass that can define specific output paths for your data (beyond
the ability of simply defining the top-level root directory).


Executing the DAP
-----------------

The main DAP scripts, ``manga_dap``, can accept the user-defined classes as
command-line arguments.  Say you have defined the class in a *local* file called
``my_datacube.py``, and the name of the implemented datacube class is
``MyDataCube``.  You can execute the DAP on a file (``input_file.fits``) with
the datacube to analyze as follows:

.. code-block:: console

    manga_dap -f input_file.fits --cube_module my_datacube MyDataCube --plan_module mangadap.config.analysisplan.AnalysisPlan -o dap

A few notes:

    - I've used the ``-f`` tag, which specifies the input file as the file that
      *contains* the datacube data.  To use a configuration file, you would use
      the ``-c`` option, and ``MyDataCube`` would have to have a ``from_config``
      class method.

    - In defining the cube module to use, note that I did not include the
      ``.py`` extension of the file name, and I included a space between
      ``my_datacube`` and ``MyDataCube``; these are the current way that you
      indicate that you're trying to import the relevant class from a local
      file.

    - Currently, the ``manga_dap`` defaults to using the MaNGA-specific classes.
      This means that you must force the code to use the
      :class:`~mangadap.config.analysisplan.AnalysisPlan` base class, if you
      haven't created a new subclass.  Here, I've defined the full object to
      import, connected by ``.``.

    - I did not include a plan file (using the ``-p`` option).  This means the
      code will use the "default" analysis plan; see
      :func:`~mangadap.config.analysisplan.AnalysisPlan.default`.  See
      `default_plan.toml
      <https://github.com/sdss/mangadap/blob/master/mangadap/config/default_plan.toml>`_

    - I set the top-level directory for the output to be ``./dap/``.

Problems?
---------

If you have problems, particularly those that you think may be a more general
problem, please `Submit an issue`_.



