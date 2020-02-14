Class Instantiation: :class:`mangadap.dapfits.DAPQualityBitMask`

=========  ===  ====================================================================================
Key        Bit  Description                                                                         
=========  ===  ====================================================================================
FORESTAR   0    There is a FORESTAR region within the data cube.                                    
BADZ       1    NSA redshift does not match derived redshift (*placeholder*)                        
LINELESS   2    No emission lines in data cube (*placeholder*)                                      
PPXFFAIL   3    pPXF failed to fit this object (*placeholder*)                                      
SINGLEBIN  4    Voronoi binning resulted in all spectra in a single bin                             
BADGEOM    5    Invalid input geometry; elliptical coordinates and effective radius are meaningless.
DRPCRIT    28   Critical failure in DRP                                                             
DAPCRIT    29   Critical failure in DAP                                                             
CRITICAL   30   Critical failure in DRP or DAP                                                      
=========  ===  ====================================================================================

