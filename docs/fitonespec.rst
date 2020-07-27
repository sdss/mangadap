
.. _fitonespec:

How to Fit One Spectrum
=======================

The low-level functions in the DAP treat each spectrum independently.
As such, it's possible to use these low-level DAP functions to
effectively fit any spectrum you like, with the caveat that you'll
need to vet the results. Here, we guide you through the example
script ``$MANGADAP_DIR/examples/fit_one_spec.py``. It is not
sufficiently abstracted to be an executable (included in
``$MANGADAP_DIR/bin``); however, we may eventually provide such an
executable. For now, the script is meant to provide a worked example
for generating your own script.

Overall Workflow
----------------

    #. Read in desired spectrum and associated galaxy redshift.

    #. Generate a mask (i.e., a "SpectralPixelMask") of the 5577 Ang sky
       emission line and intrinsic galaxy emission lines.

    #. Read in the desired spectral template library (using
       :class:`mangadap.proc.templatelibrary.TemplateLibrary`).

    #. Use pPXF to fit a linear combination of these templates to the
       spectrum.

    #. Calculate moments of the emission lines after subtracting the
       stellar continuum (using "EmissionLineMoments").  (??? is this
       what is actually happening???)

    #. Generate a new mask including only the 5577 Ang sky emission line.

    #. Perform a simultaneous fit of the stellar continuum and emission
       lines while holding the stellar kinematics determined from the
       initial pPXF run fixed.  This is the model of the spectrum that
       is ultimately output by the DAP.

    #. Re-measure the emission line moments.

User-Supplied Inputs
--------------------

    #. Plate-ifu ("plt" and "ifu")

    #. Spaxel coordinates ("x" and "y")

    #. Stellar continuum template keyword (e.g., "sc_tpl_key = 'MILESHC'")

       This keyword specifies the template library that should be used
       in the initial run of pPXF (step 4 in the Overall Workflow).

       The desired template set must have an associated configuration
       file (``*.ini``) in
       ``$MANGADAP_DIR/python/mangadap/config/spectral_templates`` in
       which this keyword is specified.  The templates themselves must
       exist in a subdirectory of
       ``$MANGADAP_DIR/data/spectral_templates/``.

       We strongly recommend using the ``MILESHC`` keyword here, as we
       have found other template libraries are very likely to yield
       unreliable stellar kinematics.

    #. Emission plus stellar continuum template keyword (e.g.,
       "el_tpl_key = 'MILESHC'")

       This keyword specifies the template library that should be used
       in the second continuum+emission line fitting step (step 7 in the
       Overall Workflow).  This need not be the same as the stellar
       continuum template keyword described above.  However, as above,
       this library must have an associated ``*.ini`` configuration
       file, and must exist within the spectral_templates subdirectory.

    #. Emission line database to be used for calculating emission-line
       moments (e.g., ``elmom_key = 'ELBMILES'``)

       TBW

    #. Emission line database to be used for fitting emission lines
       (e.g., ``elfit_key = 'ELPMILES'``)

       TBW

    #. Desired spectral pixel scale of the templates relative to the
       data ("velscale_ratio")
        
       TBW

    #. Emission line database to be used for emission-line masking.
       This is specified in the lines

       .. code-block:: python

           # Mask the 5577 sky line and the emission lines
           sc_pixel_mask = SpectralPixelMask(artdb=ArtifactDB.from_key('BADSKY'),
                                             emldb=EmissionLineDB.from_key('ELPFULL'))

Examining the Results
---------------------

TBW

