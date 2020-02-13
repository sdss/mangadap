
=============  ==========  =================  ========  =====================================================================
Key            Type        Options            Default   Description                                                          
=============  ==========  =================  ========  =====================================================================
``plate``      int         ..                 ..        Plate number                                                         
``ifudesign``  int         ..                 ..        IFU designation                                                      
``mode``       str         ``CUBE``, ``RSS``  ``CUBE``  DRP 3D mode; see :func:`mangadap.drpfits.DRPFits.mode_options`.      
``vel``        int, float  ..                 ..        Systemic velocity (km/s)                                             
``vdisp``      int, float  ..                 100.0     Guess velocity dispersion (km/s)                                     
``ell``        int, float  ..                 0.0       Isophotal ellipticity (:math:`\varepsilon = 1-b/a`)                  
``pa``         int, float  ..                 0.0       Position angle (degrees) of the isophotal major axis from N through E
``reff``       int, float  ..                 1.0       Effective radius (arcsec)                                            
=============  ==========  =================  ========  =====================================================================

