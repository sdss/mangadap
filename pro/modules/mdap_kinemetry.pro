;#######################################################################
;
; Copyright (C) 2005, Davor Krajnovic
;
; If you have found this software useful for your
; research, I would appreciate an acknowledgment to use of the
; `Kinemetry: a generalisation of photometry to the higher moments
; of the line-of-sight velocity distribution' by Krajnovic et al. (2006)'
;
; This software is provided as is without any warranty whatsoever.
; Permission to use, for non-commercial purposes is granted.
; Permission to modify for personal or internal use is granted,
; provided this copyright and disclaimer are included unchanged
; at the beginning of the file. All other rights are reserved.
;
;#######################################################################
;NAME:
;   KINEMETRY
;
;PURPOSE:
;   Perform harmonic expansion of 2D maps of observed kinematic
;   moments (velocity, velocity dispersion, h3, h4...) along the best
;   fitting ellipses (either fixed or free to change along the radii).
;
;EXPLANATION:
;   This program is a generalisation of the ellipse fitting method used
;   in photometry (e.g. Jedrzejewsky 1987) to the odd moments of the 
;   line-of-sight velocity distribution (LOSVD). The even moments are 
;   treated as in in photometry, while there is a modification for the 
;   odd moments. The method assumes that along the best fitting ellipse
;   the odd moment is well approximated by a simple cosine term, as in 
;   the tilted ring method by Begeman (1987). We use interpolation to 
;   sample the kinematic moment along the ellipses, as done in all 
;   ellipse-fitting implementations for photometry, in the few-pixels
;   regime.
;
;   For a given radius a best fit ellipse is described by flattening
;   'q' (= 1 - ellipticity) and its position angle 'PA'. The ellipse
;   parameters are found similarly to Jedrzejewsky (1987) photometry
;   approach, but in the case of odd moments by minimising 'a1', 'a3',
;   and 'b3' coefficients (or coeffs. next to sin(x), sin(3x), cos(3x))
;   of the Fourier expansion of a kinematic profile extracted along an 
;   ellipse. This is possible because a small error in q produces a 
;   non-zero b3 coefficient, and a small error in PA produces non-zero 
;   a1, a3, and b3 coefficients. Errors in the position of the centre 
;   produces nonzero even terms (a2,b2), which are ignored for the moment.
;
;   The determination of q and PA is done in two steps. First a
;   grid of q and PA is specified (not very dense but covering the
;   whole parameters space) and for each combination of the
;   ellipse parameters a kinematic profile is extracted. Using a
;   least-squares fit (singular value decomposition) the Fourier
;   coefficients are determined. The combination (q,PA), for which
;   a1^2+a3^2+b3^2 is the smallest, is used as the input for the
;   second step, which consists in a non-linear search for best
;   parameters performed by MPFIT.
;
;   After determination of the best fitting ellipse parameters
;   (described by q and PA for a given radius), a new Fourier expansion
;   is performed. The results are the Fourier coefficients and 
;   reconstructed kinematic moment maps.
;
;CALLING SEQUENCE:
;   KINEMETRY, xbin, ybin, moment, rad, pa, q, cf, [IMG=img, X0=x0, Y0=y0,$
;	  NTRM=NTRM, ERROR=error, SCALE=scale, NRAD=nrad, NAME=name, $
;	  PAQ=paq, NPA=npa, NQ=nq, RANGEQ=rangeq, RANGEPA=rangepa, ALL=all,$ 
;   EVEN=even, VSYS=VSYS, VELCIRC=velCirc, VELKIN=velkin, GASCIRC=gascirc, BMODEL=bmodel,$
;	  XC=xc, YC=yc, ER_CF=er_cf, ER_PA=er_pa, ER_q=er_q, ER_XC=er_xc, $
;	  ER_YC=er_yc, XELLIP=xellip, YELLIP=yellip, RING=ring, RADIUS=radius,$
;	  COVER=cover, VSYS=vsys, PLOT=plot, VERBOSE=verbose]
;
;INPUTS:
;   XBIN   - 1D array with X coordinates describing the map
;   YBIN   - 1D array with Y coordinates describing the map
;   MOMENT - 1D array with kin.moment (e.g. velocity) values 
;            at XBIN, YBIN positions
;
;OPTIONAL INPUTS KEYWORDS:   
;   /NTRM  - scalar specifying the number of terms for the harmonic analysis 
;	           of the profile extracted along the best fitting ellipse. Default
;	           value is 6 odd terms, which means the following terms will be 
;	           used in the expansion: a1, b1, a3, b3, a5 and a5, or the terms
;	           corresponding to sin(x),cos(x),sin(3x),cos(3x),sin(5x),cos(5x).
;   ERROR  - 1D array with errors to VELBIN values. If this
;            keyword is specified then the program will calculate
;	           formal (1sigma) errors of the coefficients which are
;	           returned in ER_PA, ER_Q (determined by MPFIT), ER_CF 
;	           (determined by SVD) variables. 
;            IF IMG keyword is set ERROR has to be a 2D array. If if is not
;	           supplied, it is creatred as 2D arrays with all values equal to 1. 
;   SCALE  - scalar specifying the pixel scale on the map. If not set, the
;            SAURON pixel scale (0.8 arcsec) is assumed.
;   IMG    - 2D array containing an image. This keyword was designed specifically
;	           for surface photometry (analysis of the zeroth moment of the LOSVD), 
;	           to increase the spead of calculations, but if kinematic map is large 
;	           (and has many pixels) it can be also passed through this keyword and 
;	           analysis will be quicker. To use kinemetry for photometry it is also 
;	           necessary to set keyword EVEN and it can be useful to use NTRM=10
; 	         in order to measure disky/boxy deviations (4th terms). 
;	           Images are very different from current kinematic maps. They are much larger
;            and usually not binned, so make regular 2D grids. This is the reason to
;	           treat the kinematic maps and images in a different way. If IMG is set,
;	           the treatment follows the ideas of Jedrzejewsky (1987): inner parts of
;	           (small radii, rad < 40 pixels) are interpolated while outer parts are 
;            binned in elliptical sectors (64 in angle to increase the signal and of 
;            thicknes equal to 10% of the given radius). When IMG is used center can 
;            be also fitted (currently only in photometry, so if keyword EVEN is not 
;            set - which usually means one is analysing a velocity map - center is not fitted). 
;            An estimate of the center is given through XC and YC keywords. ERROR should be a 2D array
;            of the same size as IMG. When IMG is used, VELCIRC and VELKIN keywords contain 
;            reconstructed images. It is assumed that image coordiantes are in pixels 
;            (not physical units). Keyword SCALE is automatically set to 1. If this is not 
;            the case, set SCALE to a required number. In order to be compatible with 
;            previous versions of kinemetry XBIN, YBIN and MOMENT still have to be passed
;            but they can be dummmy 1D variables, unless certain image areas are masked 
;            (as bad pixels). IF keyword BADPIX is used, XBIN and YBIN should be 1D arrays 
;            with the actual coordinates of the IMG. One can use the following set of lines 
;            to make the arrays: 
;                    s=size(img)  
;                    n=s[1]*s[2]
;                    yy=REPLICATE(1,s[1])#(indgen(s[2]))
;                    xx=(indgen(s[1]))#REPLICATE(1,s[2])
;                    xbin=REFORM(xx, n)
;                    ybin=REFORM(yy, n) 
;   X0     - an estimate of the X coordinate of the center (in pixels). If not given X0=0.
;	           For accurate determination of the center and other ellipse parameters at small
;	           radii it is important the ellipse includes the center of the galaxy.
;   Y0     - an estimate of the Y coordinate of the center (in pixels). If not given Y0=0.
;   FIXCEN - keyword, if set center will be fixed and no attempt will be made during the
;            harmonic decomposition to find new center. Center is fixed to X0 and Y0. 
;            This keyword is optional only for photometry (or even moments in general). 
;            For ODD moments, center is fixed always. 
;   NRAD   - scalar specifying the number of radii along which kinemetry should be run.
;            IF not specified, NRAD=100. Kinemetry will stop when the edge of the map
;	           is encountered and NRAD is not necessary achived. To force kinemetry to 
;	           do all radii, relax condition in keyword COVER.
;   NAME   - name of the object (used by VERBOSE keyword and for internal
;            plotting)
;   PAQ    - 2 element or 2*NRAD element vector specifying position angle (PA)
;	           and flattening (q) of the ellipses in that order (kept constant).
;	           It is possible to specify a set of PA and q values (that 
;	           correspond to given radii (see RADIUS keyword)), for which one
;	           wants to get the Fourier coefficients. In this case PAQ should
;	           be set as follows: PAQ=[PA1,Q1, PA2,Q2...., PAnrad,Qnrad]
;            It can be also used as an initial condition for determination of ellipses. 
;            In this case, it should be called together with /NOGRID keyword (currently  
;            implemented only for photometry). 
;            IF PAQ keyword is used to define ellipses along which harmonic
;            decomposition is made, then keyword NOGRID should not be used. In this case 
;            center is fixed (and should be defined via X0 and Y0 keywords if IMG keyword
;            is used).
;   NOGRID - keyword, if set it bypasses the direct minimisation via a grid in PA
;	           and Q values. It should be used together with PAQ parameters, when
;	           a good estimate of PA and Q are passed to the program, but not if PAQ
;            serve to pre-define ellipses for harmonic decomposition.
;	           It is desigend with photometry in mind, where the problem usually has
;	           only one well defined minimum (in PA,Q plane). It speeds up the
;	           calculation, but for the higher kinematic moments it is not as robust and 
;            it is advised not to be used (first use the grid to find the best fit 
;            PA and Q values.).
;   NPA    - scalar specifying the number of PA used to crudely estimate
;            the parameters of the best fit ellipse before entering
;            MPFIT. Default value is 21. To speed up the process and for
;	           quick tests it is useful to use a small number (e.g 5). 
;            Complicated maps may require a bigger number (e.g. 41).
;   NQ     - scalar specifying the number of q used to crudely estimate
;            the parameters of the best fit ellipse before entering
;            MPFIT. Default value is 21. To speed up the process and for
;	           quick tests it is useful to use a small number (e.g 5). 
;            Complicated maps may require a bigger number (e.g. 41).
;   RANGEQ - 2 element vector specifying the min and max value for 
;	           flattening Q. Default values are 0.2 and 1.0.
;   RANGEPA- 2 element vector specifying the min and max value for 
;	           position angle PA. Default values are -90 and 90 (degrees).
;   BADPIX - 1D array containing indices of pixels which should not be
;            used during harmonic fits. This keyword is used when data 
;            are passed via IMG. It is usefull for masking stars and
;            bad pixels. When used, XBIN and YBIN should be real 
;            coordinates of the pixels of IMG (see IMG for more details). 
;            The bad pixels are passed to the routine which defines/selects
;            the ellipse coordiantes (and values) to be fitted, and they are 
;            removed from the subsequent fits. All pixels of
;            the ellipse which are 2*da from the bad pixels are removed 
;            from the array, where da is the width of the ring. 
;   /ALL   - If this keyword is set then the harmonic analysis of the rings 
;	           will include both even and odd terms. If this keyword is set,
;	           and NTRM = n then the following terms are used in the expansion:
;	           a1, b2, a2, b2, a3, b3,...., an, bn (or coeffs nex to: sin(x),
;	           cos(x),sin(2x),cos(2x),sin(3x),cos(3x),...,sin(nx),cos(nx))
;   /EVEN  - set this keyword to do kinemetry on even kinematic moments. 
;	           In this case, kinemetry reduces to photometry and the best
;	           fitting ellipse is obtained by minimising a1, b1, a2, b2
;	           terms. When this keyword is set, keyword /ALL is automatically
;	           set and NTRM should be increased (e.g. NTRM=10 will use the 
;	           following terms in the expansion: a1, b2, a2, b2, a3, b3, a4, b4
;	           (or coeffs. next to sin(x),cos(x),sin(2x), cos(2x),sin(3x),
;	           cos(3x),sin(4x),cos(4x)))
;   /VSYS  - if this keyword is set the zeroth term (a0) is not extracted.
;	           (for odd moments).This might be useful for determinatio of
;	           rotation curves.One can first run kinemetry without setting 
;	           this keyword to find the systemic velocity (given as cf[*,0]). 
;	           Then subtract the systemic velocity form the velocity map and 
;	           re-run kinemetry with /vsys set. In this case the zeroth terms 
;	           will be zero. For completeness, it is also possible to input 
;	           VSYS, e.g. VSYS=10. The zeroth term will not be calculated, 
;	           but it will be set to 10 in output. Given that Fourier terms
;	           are orthogonal, it should not be necessary to set this keyword
;	           in general. 
;   RING   - scalar specifying desired radius of the first ring. Set this
;	           keyword to a value at which the extraction should begin. This 
;	           is useful in case of ring-like structures sometimes observed in 
;	           HI data.
;   RADIUS - 1D array with values specifying the lenght of the semi-major axis
;	           at which the data (kin.profile) should be extracted from the map 
;	           for the kinemetric analisys. The values should be in pixel space 
;            (not in physical units such as arcsec).  If this keyword is set, 
;	           the values are coopied into the output variable: RAD.
;   COVER  - Keyword controling the radius at which extraction of
;	           values from the map stops. Default value is 0.75, meaning
;	           that if less than 75% of the points along an ellipse are
;	           are not present, the program stops. The value of 75% is an 
;	           ad hoc, but a conservative value which makes sure that 
;	           kinemtric coeffs. are robust (at least the lower few orders).
;	           Sometimes it is necessary to relax this condition, especially
;	           for reconstruction of the maps.
;   /BMODEL- If this keyword is set, a model moment map is constructed. This keyword
;            should be set together with VELCIRC and VELKIN keywords, which will
;	           contain the reconstructed map, using the first dominant term and all 
;	           terms, respectively. IF IMG keyword is used, the outputs are 2D images,
; 	         otherwise BMODEL reconstructs the map at each input position XBIN,YBIN.
;	           If BMODEL is not set VELCIRC and VELKIN will contain reconstructed values
;	           at the positions of XELLIP and YELLIP.
;   /PLOT  - If this keyword is set, diagnostic plots are shown for each radii: 
;		       - the best ellipse (overploted on kin.moment map). If IMG
;		         keyword is set, the image is scaled to the size of the ellipse.
;		         The centering of the overplotted ellipse is good to 0.5 pixels
;		         so for small radii (r < a few pixels) it is possible that the 
;		         position of the center of the overplotted ellipse is not on the
;		         brightes pixel. 
;		       - PA-Q grid with the position of the minimum (green diamond) for 
;		         the parameters of the best fit ellipse determined by
;		         MPFIT, where the initial (input to MPFIT) values of PA
;		         and Q are presented by the grid of dots and the colours
;	           show the Chi2 square contours (linearly interpolated 
;		         between the PA,Q points),
;		       - fit to kin.profile (white are the DATA, red is the FIT, where
;		         FIT is given by a0+b1*cos(x) for odd, and a0 for even moments), 
;		       - residuals (DATA - FIT), and overplotted higher order
;		         terms (green: a1,a3 and b3, red: a1,a3,b3,a5 and b5; 
;		         for the /EVEN case - green: a1,b1,a2,b2, red:a1,b1,a2,b2,a4,b4)
;   /VERBOSE - set this keyword to print status of the fit on screen
;            including information on:
;               - Radius - number of the ring that was analysed
;               - RAD    - radius of the analysed ring (if SCALE is passed 
;                          it is given in the same units, otherwise in pixels)
;               - PA     - position angle of the best fitting ellipse
;               - Q      - flattening of the best fitting ellipse
;               - Xcen   - X coordinate of the centre of the ellipse
;               - Ycen   - Y coordinate of the centre of the ellipse
;               - # of ellipse elements - number of points to which the data
;                          points in the ring are sampled before derivation of
;                          the best fit parameters and harmonic analysis. It 
;                          varies between 20 and 64 (or 100 in non IMG model)
;                          depending on the ring size, giving a typical sampling 
;                          of 5.6 (3.6) degrees.
;
;
;OUTPUT:
;   RAD    - 1D array with radii at which kin.profiles were extracted
;   PA     - 1D array with position angle of the best fitting ellipses,
;            PA is first determined on an interval PA=[-90,90], where
;            PA=0 along positive x-axis. Above x-axis PA > 0 and below
;            x-axis Pa < 0. PA does not differentiate between receding
;            and approaching sides of (velocity) maps. This is
;            transformed to the usual East of North system, where the 
;	           East is the negative x-axis, and the North is the positive 
;	           y-axis. For odd kin.moments PA is measured from the North 
;	           to the receding (positive) side of the galaxy (which is 
;	           detected by checking the sign of the cos(theta) term. For
;	           the even terms it is left degenerate to 180 degrees rotation.
;   Q      - 1D array with flattening of the best fitting ellipses
;            (q=1-ellipticity), defined on interval q=[0.2,1]
;   CF     - 2D array containing coefficients of the Fourier expansion
;            for each radii cf=[Nradii, Ncoeff]. For example: 
;	           a0=cf[*,0], a1=cf[*,1], b1=cf[*,2]....
;
;OPTIONAL OUTPUT KEYWORDS:
;   VELKIN - 1D array of reconstructed kin.moment using NTRM harmonic
;            terms at positions XBIN,YBIN, obtained by linear interpolation
;	           from points given in XELLIP and YELLIP keywords (if BMODEL keyword
;	           is used, otherwise at positions XELLIP and YELLIP.
;   VELCIRC -1D array containinng 'circular velocity' or a0 + b1*cos(theta)
;            at positions XBIN, YBIN (velcirc = a0, in case of EVEN moments), 
;	           obtained by linear interpolation from points given in XELLIP 
;	           and YELLIP keywords (if BMODEL keyword is used, otherwise at 
;	           positions XELLIP and YELLIP.
;   GASCIRC -1D array containing circular velocity or Vcirc=cf[*,2]*cos(theta)
;            at positions XBIN and YBIN, obtained for fixed PA and q. 
;            PA and q are taken to be median values of the radial variation of PA and q. 
;            IF keyword PAQ is used than GASCIRC give the circular velocity (no systemic
;            velocity) for the median values of PA, Q values. Note that this is different
;            from VELCIRC (also VELKIN) which is obtained on the best
;            fitting ellipses and also includes Vsys (cf[*,0]) term. 
;            This keyowrd is useful for gas velocity maps, if one wants to obtain 
;            a quick disk model based on the circular velocity. 
;   ER_PA  - 1D array of 1sigma errors to the ellipse position angle
;   ER_Q   - 1D array of 1sigma errors to the ellipse axial ratio
;   ER_CF  - 2D array containing 1 sigma errors to the coefficients 
;	           of the Fourier expansion for each radii
;   XELLIP - 1D array with X coordintes of the best fitting ellipses
;   YELLIP - 1D array with Y coordintes of the best fitting ellipses
;    XC    - the X coordinate of the center (in pixels). If X0 not fit XC=0.
;    YC    - the Y coordinate of the center (in pixels). If Y0 not fit YC=0.
;
;RESTRICTIONS:
;   Speed and robustness of the program depends on the number of NQ
;   and NPA that define the initial grid. Small grid is fast, but
;   not precise in complicated cases. In case of IMG keyword, plotting
;   the results on the screen is a major contributor to the decrease
;   in speed.
;
;NOTE: determination of the centre and penalisation towards extraction 
;      on the circles is not yet included in this distribution. 
;
;REQUIRED ROUTINES:
;   MPFIT: by C.B. Markwardt from http://astrog.physics.wisc.edu/~craigm/idl/
;   RANGE: by M. Cappellari (included in kinemetry distribution)
;   SAURON_COLORMAP: by Michele Cappellari & Eric Emsellem, 2001
;			(included in kinemetry distribution)
;EXAMPLE:
;    1) 
;    Run kinemetry on a velocity map defined with coordinate arrays:
;    Xbin, Ybin, and a velocity array: Velbin. Desired outputs are: 
;    position angle and flattening of the ellipses, harmonic terms:
;    a0, a1, b1, a3, b3, a5, b5 and a reconstructed map of the circular
;    velocity: 
; 
;    KINEMETRY, Xbin, Ybin, Velbin, rad, pa, q, cf, NTRM=6, VELCIRC=velcirc,$
;	        /bmodel, /plot, /verbose 
;
;    2)
;    Run kinemetry on a velocity map starting at radius=5". Desired outputs 
;    are position angle and flattening of the ellipses, harmonic terms:
;    a0, a1, b1, a2, b2, a3, b3, a4, b4, a5, b5, a reconstructed map 
;    of the circular velocity and a map reconstructed with all terms: 
;
;    KINEMETRY, Xbin, Ybin, Velbin, rad, pa, q, cf, NTRM=10, /ALL, $
;    VELCIRC=velcirc, VELKIN=velkin, RING=5, /bmodel, /plot, /verbose 
;
;    3) 
;    Run kinemetry on a velocity dispersion map (given by array Sigma). 
;    Desired outputs are position angle and flattening of the ellipses,
;    and harmonic terms: a0, a1, b1, a2, b2, a3, b3, a4, b4:
;
;    KINEMETRY, Xbin, Ybin, Sigma, rad, pa, q, cf, NTRM=10, \EVEN, $
;               /plot, /verbose
;
;    4)
;    Run kinemetry on an image of a galaxy to perform surface photometry.
;    Image is given with 2D array IMG. X0 and Y0 are estimates of the center
;    of the galaxy. X,Y,FLUX are dummy 1D arrays. They are still required, but are
;    irrelevant if IMG keyword is used and can be anything. Note that plotting
;    slows down the program.
;
;    KINEMETRY, X, Y, FLUX, rad, pa, q, cf, IMG=IMG, X0=X0, Y0=Y0, $
;               NTRM=10, /EVEN, /verbose, /plot
;
;    5)
;    Run kinemetry on an image of a galaxy to perform surface photometry.
;    Image is given with 2D array IMG. X0 and Y0 are estimates of the center
;    of the galaxy. X,Y,FLUX are 1D dummy arrays. Determination of PA and Q 
;    is skipping the brute force via a grid (\NOGRID), but an initial estimate 
;    must be passed via PAQ keword. Outputs are for each ellipse of semi-major
;    axis length RAD (XELLIP, YELLIP): center coordinates (XC,YC), PA, Q, CF 
;    and associated errors. Models are made reconstruceted on each image pixel.
; 
;    KINEMETRY, X, Y, flux, rad, pa, q, cf, IMG=img, NTRM=10, /EVEN, $
;               X0=x0, Y0=y0, /nogrid, PAQ=[ang-90, 1-epsI], $
;               XC=xc, YC=yc, XELLIP=xellip, YELLIP=yellip,$
;               ER_CF=er_cf, ER_PA=er_pa, ER_q=er_q, ER_XC=er_xc, ER_YC=er_yc, $
;               VELCIRC=model, VELKIN=modelF, /bmodel, /verbose, /plot
;
;    6) 
;    Run kinemetry on an velocity map given as a 2D array (IMG_VEL). 
;    Firstly the minimisation is done on a grid of NPAxNQ values. 
;    Center is fixed at X0,Y0 and kept constant. Velocity maps are 
;    reconstructed at along the best fitting ellipses (XELLIP, YELLIP).
;    X,Y,V are dummy arrays. 
;
;    KINEMETRY, X, Y, V, rad, pa, q, cf, NPA=21, NQ=21, X0=x0, Y0=y0,$
;              IMG=IMG_VEL, VELCIRC=model, VELKIN=modelF, $
;              XELLIP=xellip, YELLIP=yellip, /verbose,  /plot, /bmodel
;
;
;REVISION HISTORY:
;   V1.0 - Written by Davor Krajnovic and Michele Cappellari (March, 2005)
;   V2.0 - First released version, Davor Krajnovic, Oxford, 7.12.2005. 
;   V2.1 - Add keyword RADIUS, DK, Oxford, 27.01.2006. 
;   v2.2 - Changed definition of PAQ keyword, DK, Oxford, 02.02.2006.
;   v2.3 - Corrected bug: keyword \even was not pass correctly to MPFIT. 
;	         Thanks to Roland Jesseit for prompting. Minor homogonisation 
;          of the code. DK, Oxford, 06.02.2006
;   V2.4 - Add keyword RANGEPA, DK, Oxfrod, 11.07.2006.
;   V2.5 - Add keyword COVER, DK, Zagreb, 02.11.2006.
;   v2.6 - Introduced standard definition of coeffs. for odd kinematic
;	         moments: they are measured from receding side of galaxy.
;	         Thanks to Anne-Marie Weijmans for asking. DK, La Palma, 24.04.2007.
;   v3.0 - Adaptation of kinemetry to efficiently work on images of galaxies
;          for doing surface photometry. New keywords included are IMG, X0,Y0,
;	         NRAD, BMODEL, XC,YC, ER_XC, ER_YC. kinem_extract_sector.pro was added, 
;	         while the main procedure and kinem_fitfunc_ellipse.pro were modified.
;	         Present distribution is fully compatible with the old one, with one 
;	         exception: VELCIRC and VELKIN are now returend as reconstructions at
;	         XELLIP and YELLIP positions unless keyword BMODEL is set as well. In 
;	         that case VELCIRC and VELKIN are reconstruceted as map. For the moment
;	         XBIN, YBIN and MOMENT still need to be passed at the same time as IMG, 
;	         if IMG is set they are not used at all and can be anything. This will 
;	         be change in future distributions. DK, Oxford, Queens College, 16.11.2007.
;   v3.1 - Clearing of a bug for using IMG for ODD moments. 2D images of both even (flux) 
;          and odd moments (e.g. velocity) can be analysied. Center is STILL not fitted 
;          for ODD moments. Thanks to Maxim Bois for prompting. DK, Hatfield, 03.09.2008.
;   v3.2 - Passing center values when no minimization is attempted and PAQ values are 
;          kept constant. DK, Oxford, 31.10.2008.
;   v3.3 - Introduced FIXCEN for EVEN moments. Sorted bug related to RADIUS and IMG
;          keywords. Thanks to Kristen Shapiro. DK, Hatfield, 03.11.2008.
;   v3.4 - Sort some minor bugs related to FIXCEN. Thanks to Kristen Shapiro. 
;          DK, Hatfield, 17.12.2008.
;   v3.5 - Introduced keyword BADPIX, Oxford, DK, 20.02.2009.
;   v3.6 - Introduced keyword GASCIRC and VSYS. Thanks to Witold Maciejewski. 
;          DK, Muechen, 08.12.2009.
;   v3.7 - Major testing. In IMG model ellipse sampling was changed to depend on 
;          the radius, but with the maximum of 64 points, affecting first 5 rings. 
;          For nonIMG mode the maximum wasleft to 100 as before. Cleaning of the
;          documentation. DK, Muenchen, 13.05.2010. 
;   v3.8 - Adjusted plotting of the map in /even case (but no IMG) for diagnostic
;          plots (when /plot set). DK, Garching, 20.05.2010.
;   v3.9 - Added SKY keyword to enable stopping at the sky level. DK, Munich, 26.10.2010.
;   v4.0 - Introduced GRIDDATA instead of a call to TRIGRID in the function kinem_trigrid_irregular, 
;          following a suggestion from Michele Cappellari. This speeded up the the evaulation
;          of the gird significantly (~6 times for a SAURON map). 
;        - Fixed fitting with a free centre for the EVEN case. 
;        - Re-structured calls to the initiatlization of parameters for MPFIT (parinfo)
;          for various input keywords. 
;          DK, Potsdam, 10.10.2013.
;#######################################################################
PRO sauron_colormap
;
; The official SAURON colormap: type "sauron_colormap" to load the colormap.
;
; Michele Cappellari & Eric Emsellem, Leiden, 10 July 2001
;
; Start with these 7 equally spaced coordinates, then add 4 additional points
; x = findgen(7)*255/6. + 1
; 1.0  43.5  86.0  128.5  171.0  213.5  256.0
;
x = [1.0, 43.5, 86.0, 86.0+20, 128.5-10, 128.5, 128.5+10, 171.0-20, 171.0, 213.5, 256.0]

red =   [0.01, 0.0, 0.4,  0.5, 0.3, 0.0, 0.7, 1.0, 1.0,  1.0, 0.9]
green = [0.01, 0.0, 0.85, 1.0, 1.0, 0.9, 1.0, 1.0, 0.85, 0.0, 0.9]
blue =  [0.01, 1.0, 1.0,  1.0, 0.7, 0.0, 0.0, 0.0, 0.0,  0.0, 0.9]

xnew = FINDGEN(256)+1

redV = INTERPOL(red, x, xnew)
greenV = INTERPOL(green, x, xnew)
blueV = INTERPOL(blue, x, xnew)

TVLCT, redV*255, greenV*255, blueV*255 ; load the SAURON colormap

END
;----------------------------------------------------------------------------
FUNCTION range, x1, x2, n
;
; RANGE(x1,x2) = x1,x1+1,...,x2. In this case x1, x2 should be integers.
; RANGE(x1,x2,n) = x1,x1+dx,...,x2 with N integer the result has length N.
; The result will have the type of x1, but if three parameters are used
; the result will be at least of type float.
;
; Michele cappellari, Leiden, 16 October 2001
;
compile_opt idl2
on_error, 2

t = SIZE(x1,/TYPE)
CASE n_params() of
    2: IF x1 LT x2 THEN $
            v = x1 + INDGEN(x2-x1+1, TYPE=t) $
        ELSE $
            v = x1 - INDGEN(x1-x2+1, TYPE=t)
    3: BEGIN
            v = x1 + (x2-x1)/(n-1.0)*INDGEN(n, TYPE=t)
    END
    ELSE: message, '2 or 3 parameters are needed'
ENDCASE

RETURN, v
END
;----------------------------------------------------------------------
PRO kinem_trigrid_irregular, xbin, ybin, moment, xnew, ynew, momNew, _EXTRA=extra
compile_opt idl2
;
; Given a set of irregular coordinates and corresponding values
; output coordinates, which are also irregular, it gives an interpolated values 
; at the output positions.
; The interpolation is done using TRIGRID and any option can be passed
; to this routine via the _EXTRA mechanism.
;
; V1.0: Michele Cappellari, Leiden, 17 December 2004
; V1.1: small cosmetic changes, Davor Krajnovic, Edinburgh, 24.11.2005.
; V2.0: introduced GRIDDATA instead of TRIGRID, following Michele's suggestion.
;----------
nbin = N_ELEMENTS(xBin)
IF nbin NE N_ELEMENTS(ybin) or nbin NE N_ELEMENTS(moment) THEN $
    message, 'XBIN, YBIN and VELBIN must have the same size'
nNew = N_ELEMENTS(xNew)
IF nNew NE n_elements(yNew) THEN message, 'XNEW and YNEW must have the same size'

TRIANGULATE, xbin, ybin, tr
momNew = GRIDDATA(xbin, ybin, moment, XOUT=xnew, YOUT=ynew, $
                   /LINEAR, TRIANGLES=tr, _EXTRA=extra)

END
;----------------------------------------------------------------------
PRO kinem_extract_sector, img, theta, a, da, q, pa, xc, yc, momell, $
	xell, yell, INTERP=interp, ERROR=error, ER_MOMELL=er_momell, $
	BADPIX=badpix, XBAR=xbar, YBAR=ybar, _EXTRA=extra
;
; Michele Cappellari, Oxford, 13 November 2007
; Expanded and incorporated in kinemetry, including interpolation
; and varios other keywords, DK, 15.11.2007

;
; ellipse along which the profile was extracted
;
x = a*COS(theta)
y = a*SIN(theta)*q
xEll = xc + x*COS(pa) - y*SIN(pa)
yEll = yc + x*SIN(pa) + y*COS(pa)

;
; definition of the sector apperture
;
thAper = [-0.5, -0.5, +0.5, +0.5]*(theta[1]-theta[0])
rAper  = [+0.5, -0.5, -0.5, +0.5]*da + a

IF KEYWORD_SET(interp) THEN BEGIN
	momell=INTERPOLATE(img, xell, yell, _EXTRA=extra)
	IF keyword_set(ERROR) THEN er_momell=INTERPOLATE(error, xell, yell, _EXTRA=extra)$
			                  ELSE er_momell=momell*0.+1
	;print, 'interpolating'		                  
ENDIF ELSE BEGIN
  ;print, 'ringing'
	s = size(img)
	momell = theta*0
	er_momell= momell +1.
;plot, 	xell, yell, /ynozero, psym=5;, /iso
;oplot, xbar, ybar, psym=3
	FOR j=0,n_elements(theta)-1 DO BEGIN
	   xaper = rAper * cos(theta[j]+thAper)
	   yaper = rAper*q * sin(theta[j]+thAper)    
	   xnew = xc + xaper*cos(pa) - yaper*sin(pa)
	   ynew = yc + xaper*sin(pa) + yaper*cos(pa)
	   w = POLYFILLV(xnew, ynew, s[1], s[2])
	   IF w[0] ne -1 THEN BEGIN
			 IF N_ELEMENTS(w) gt 5 THEN momell[j]=BIWEIGHT_MEAN(img[w]) ELSE momell[j] = mean(img[w])
			 IF keyword_set(ERROR) THEN er_momell[j] = TOTAL(SQRT(error[w]))/N_ELEMENTS(w)
	   ENDIF; ELSE print, j
;plots, xnew,ynew, psym=1
;oplot, xbar[w]+1, ybar[w]+1, psym=4
	ENDFOR
ENDELSE

;
; Check if there are bad pixels within the ring
;
IF keyword_set(badpix) THEN BEGIN
  xb=xbar[badpix]
  yb=ybar[badpix]
  bad=momell*0.
  FOR i=0, N_elements(momell)-1 DO BEGIN
      dist=SQRT( (xell[i]-xb)^2 + (yell[i]-yb)^2)
      w=where(dist le 2*da, nw)
      IF nw ne 0 THEN bad[i]=1
  ENDFOR  
  ;
  ; keep only good pixels  
  ;
  w=where(bad eq 0)
  momell=momell(w)       
  er_momell=er_momell[w]
  xell=xell[w]
  yell=yell[w]
  theta=theta[w]
ENDIF  ;(badpix)

  
END
;-----------------------------------------------------------------------------
FUNCTION kinem_SVD_Solve, a, b, ERRORS=errors
compile_opt idl2
;
; Solve a linear system using the method
; of the Singular Values Decomposition.
; IF error keyword set, a and b arrays are 
; considered to be properly divided by 
; Gaussian errors. Uncertainties of the returned
; parameters are calculated by computing the 
; covariance matrix and taking the SQRT of 
; the diagonal elements (following svdvar routine
; from Numerical Recepies, Press et al. 1992.). 
;
; This function returns coefficients of the linear 
; system. If /ERRORS is used then first N_ELEMENTS(w)
; values are the coefficients, while the second 
; N_ELEMENTS(w) values are the corresponding 
; uncertainties. 
;
; Written by Davor Krajnovic, Oxford, 25.05.2005, 
; (expanding on original routine by M.Cappellari)
; 

sa = SIZE(a)
sb = SIZE(b)
IF (sa[1] NE sb[1] or sa[0] NE 2 or sb[0] NE 1) THEN $
    message, 'SVD_Solve Error: incompatible dimensions'

SVDC, a, w, u, v, /COLUMN, /DOUBLE
small = WHERE(w LT MAX(w)*1d-13,m)
IF m NE 0 THEN w[small] = 0d
weights=SVSOL(u, w, v, b, /COLUMN, /DOUBLE)


IF keyword_set(errors) THEN BEGIN
	np = N_ELEMENTS(weights)
	wi=dblarr(np) & er = dblarr(np)
	cvm=dblarr(np,np)	; covariance matrix

	FOR i=0, np-1 DO BEGIN 
		FOR j=0, i DO BEGIN
			sum=0.
			FOR k=0, N_ELEMENTS(weights)-1 DO BEGIN
				IF w[k] ne 0. THEN wi[k] = 1./(w[k]*w[k])
				sum = sum + v[i,k]*v[j,k]*wi[k]
			ENDFOR
			cvm[i,j]=sum
			cvm[j,i]=sum
		ENDFOR
	ENDFOR

	FOR i=0, N_elements(weights)-1 DO er[i]=SQRT(cvm[i,i])

	param=[weights, er]
	RETURN, param
ENDIF

RETURN, weights

END
;----------------------------------------------------------------------
FUNCTION kinem_fit_trig_series, x, y, yfit, NTERMS=nterms, COEFF=coeff, $
			        ERR=err, ALL=all, EVEN=even, VSYS=vsys
compile_opt idl2
;
; v1.0 M. Cappellari, Leiden, 15 December 2004
; v1.1 D. Krajnovic, Paranal, 09.03.2004, adapted to be used for
;      reconstruction of kin.moment using arbitrary number of terms
; v1.2 expanded for error treatment, DK, Oxford, 26.05.2005
; v1.3 add vsys keyword and clean comments, DK, Edinburgh, 22.05.2005.
; v1.4 add even keyword, DK, Oxfrod, 06.02.2006.
;
;
; The harmonic series below can include 
;   1) odd J values:
;       1, [sin(x), cos(x)], [sin(3x), cos(3x)],...
;      The even trigonometric terms are zero in the point-symmetric case.
;      NTERMS=4 for all odd terms up to [sin(3x), cos(3x)] and
;      NTERMS=6 for all odd terms up to [sin(5x), cos(5x)]
;   2) even J values
;      1, [sin(x), cos(x)], [sin(2x), cos(2x)],...
;      This is used for a complete harmonic analysis, or to extract
;      EVEN terms. 
;      NTERMS=10 for all terms up to [sin(5x), cos(5x)]
;   3) odd or even, but without zero-th term (or systemic velocity)
;----------------------------------------------------------

arr = DBLARR(n_elements(x),nterms+1)
;
; if keyword vsys set, do not extract zero-th term
;
IF KEYWORD_SET(vsys) THEN arr[*,0] = 0d $
		     ELSE arr[*,0] = 1d


IF keyword_set(all) or keyword_set(even) THEN $ 
	FOR j=1,nterms,2 DO arr[*,[j,j+1]] = [sin(((j+1)/2)*x),cos(((j+1)/2)*x)] $
		    ELSE $
	FOR j=1,nterms,2 DO arr[*,[j,j+1]] = [sin(j*x),cos(j*x)]


;
; divide arrays with corresponding errors
;
IF keyword_set(err) THEN BEGIN
	y1 = y/err
	brr=arr*0.
	FOR j=0, Nterms DO brr[*,j] = arr[*,j]/err
	IF N_ELEMENTS(coeff) EQ 0 THEN coeff = KINEM_SVD_SOLVE(brr,y1, /ERRORS)
	koeff = coeff[0:nterms]
ENDIF ELSE BEGIN
	IF N_ELEMENTS(coeff) EQ 0 THEN coeff = KINEM_SVD_SOLVE(arr,y)
	koeff = coeff
ENDELSE
yfit = arr # koeff

RETURN, coeff
END
;----------------------------------------------------------------------
FUNCTION kinem_fitfunc_ellipse, p, $
	COEFF=coeff, NTERMS=nterms, XBAR=xbar, YBAR=ybar, MOMENT=moment, $
	RAD=r, XELL=xEll, YELL=yEll, MOMELL=momEll, THETA=theta, $
	MOMFIT=momFit, MPFIT=mpf, ELEM=elem, W=w, ERROR=error, $
	ER_MOMELL=er_momEll, ALL=all, EVEN=even, VSYS=vsys,$
	IMG=img, X0=x0, Y0=y0, BADPIX=badpix
compile_opt idl2
;
; Note that xell, yell and theta are passed further fully, while other
; products (momell, er_momell, momFit) are passed only for eccentricities
; (theta) that are on the map (selected by w): e.g. momell=momell[w] 
;

ang = p[0]/!RADEG

IF keyword_set(img) THEN BEGIN


  mi=(360./(180./10))
  theta = range(0.0,2.0*!pi, mi>10*r<64) 
  sp=SIZE(p)

	IF r lt 200 THEN BEGIN 
	  IF keyword_set(EVEN) AND sp[1] eq 4 THEN KINEM_EXTRACT_SECTOR, img, theta, r, r*0.1, p[1], ang, p[2], p[3], momell, xell, yell,$
 				                          /interp, MISSING=12345678, ERROR=error, er_momell=er_momell, BADPIX=badpix, XBAR=xbar, YBAR=ybar ELSE $
                              KINEM_EXTRACT_SECTOR, img, theta, r, r*0.1, p[1], ang, x0, y0, momell, xell, yell,$
                                  /interp, MISSING=12345678, ERROR=error, er_momell=er_momell, BADPIX=badpix, XBAR=xbar, YBAR=ybar
		w = WHERE(momEll NE 12345678, elem)
		IF elem eq 0 THEN RETURN, 1e30
		;theta=theta[w]
		momEll=momEll[w]
		IF KEYWORD_SET(error) THEN er_momEll=er_momEll[w] ELSE er_momELL=momEll*0.+1.
	ENDIF ELSE BEGIN
    IF keyword_set(EVEN) AND sp[1] eq 4 THEN KINEM_EXTRACT_SECTOR, img, theta, r, r*0.1, p[1], ang, p[2], p[3], momell, xell, yell,$
                                   MISSING=12345678, ERROR=error, er_momell=er_momell, BADPIX=badpix, XBAR=xbar, YBAR=ybar ELSE $
                              KINEM_EXTRACT_SECTOR, img, theta, r, r*0.1, p[1], ang, x0, y0, momell, xell, yell,$
                                   MISSING=12345678, ERROR=error, er_momell=er_momell, BADPIX=badpix, XBAR=xbar, YBAR=ybar
    w = WHERE(momEll NE 12345678, elem)
		IF elem eq 0 THEN RETURN, 1e30
		;theta=theta[w]
		momEll=momEll[w]
		IF KEYWORD_SET(error) THEN er_momEll=er_momEll[w] ELSE er_momELL=momEll*0.+1.
	ENDELSE

ENDIF ELSE BEGIN ;keyword_set(img) 
	
	;
	; construction of elliptical coordinates on which kin.moment is
	; interpolated; expansion of kin.profile in harmonic series;
	; used by both 'brute force' and MPFIT minimisation

	mi=(360./(180./10))
	theta = range(0.0,2.0*!pi, mi>10*r<100) 
	x = r*COS(theta)
	y = r*SIN(theta)*p[1]
	xEll = x*COS(ang) - y*SIN(ang)
	yEll = x*SIN(ang) + y*COS(ang)

	KINEM_TRIGRID_IRREGULAR, xBar, yBar, moment, xEll, yEll, momEll, MISSING=12345678
	w = WHERE(momEll NE 12345678, elem)
	IF elem eq 0 THEN RETURN, 1e30

	er_momEll=dblarr(N_ELEMENTS(w))
	xxe=xell[w] & yye=yell[w]
	;
	; Find the closest bin for the interpolated position 
	; (elliptical coordinates) and asign the error of 
	; the bin to that position 
	;
	IF KEYWORD_SET(error) THEN BEGIN
	   FOR i=0, N_ELEMENTS(w)-1 DO BEGIN
		     dist2 = (xbar-xxe[i])^2 + (ybar-yye[i])^2
		     wd = WHERE(dist2 eq MIN(dist2))	
		     er_momEll[i]=error[wd[0]]
	   ENDFOR
	ENDIF ELSE er_momELL=momEll*0.+1

	;theta=theta[w]
	momEll=momEll[w]
ENDELSE ;keyword_set(IMG)


;
; Calculation of coeffs. with options for ALL, EVEN or VSYS
;

coeff = KINEM_FIT_TRIG_SERIES(theta[w], momEll, momFit, NTERMS=nterms,$
				      ERR=er_momEll, ALL=all, EVEN=even, VSYS=vsys)


IF KEYWORD_SET(even) THEN BEGIN
	; 
	; Following eq.(1) of Jedrzejewski (1987), it tries to 
	; minimize a1,a2,a2,b2, which indicate incorrect PA and
	; flattening of the trial ellipse (terms are defined
	; differently than in 'odd' case	
	;
	IF KEYWORD_SET(mpf) THEN RETURN, coeff[[1,2,3,4]] $
			     ELSE return, TOTAL(coeff[[1,2,3,4]]^2)


ENDIF ELSE BEGIN
	;
	; Following eq.(1) of Jedrzejewski (1987), but for odd kinematic
	; moments, it tries to minimize a1,a3,b3, which indicate
	; incorrect PA and flattening of the trial ellipse
	;
	IF KEYWORD_SET(mpf) THEN RETURN, coeff[[1,3,4]] $
			    ELSE RETURN, TOTAL(coeff[[1,3,4]]^2)

ENDELSE

END
;----------------------------------------------------------------------
;        MAIN PROGRAM
;----------------------------------------------------------------------
PRO mdap_kinemetry, xbin, ybin, moment, rad, pa, q, cf, IMG=img, X0=x0, Y0=y0,$
	NTRM=NTRM, ERROR=error, SCALE=scale, NRAD=nrad, NAME=name, $
	PAQ=paq, NPA=npa, NQ=nq, RANGEQ=rangeq, RANGEPA=rangepa, ALL=all,$ 
  EVEN=even, VSYS=VSYS, VELCIRC=velCirc, VELKIN=velkin,  BMODEL=bmodel,$
	XC=xc, YC=yc, ER_CF=er_cf, ER_PA=er_pa, ER_q=er_q, ER_XC=er_xc, $
	ER_YC=er_yc, XELLIP=xellip, YELLIP=yellip, RING=ring, RADIUS=radius,$
	COVER=cover, PLOT=plot, VERBOSE=verbose, NOGRID=nogrid, FIXCEN=fixcen, $
	BADPIX=badpix, GASCIRC=gascirc, SKY=sky

compile_opt idl2

IF N_ELEMENTS(ntrm) EQ 0 THEN NTRM=6	    ; number of terms in F.expansion
IF N_ELEMENTS(cover) EQ 0 THEN cover = 0.65 ; required coverage of data

;
; make sure keywords are set correctly
;
ODD = 1.      ; per definition parity is set to odd and kinemetry works on velocity maps
IF KEYWORD_SET(even) THEN BEGIN                           ; even moments;
        even=1 
        odd=0                                             ; kinemetry doesn't work on velicity maps
ENDIF ELSE even = 0	          
IF KEYWORD_SET(all) THEN all=1 ELSE all = 0               ; all harmonic coeffs.;
IF KEYWORD_SET(X0) THEN X0=X0 ELSE X0=0.
IF KEYWORD_SET(Y0) THEN Y0=Y0 ELSE Y0=0.
IF KEYWORD_SET(img) THEN im=1 ELSE im=0
IF keyword_set(nogrid) THEN nogrid=1 ELSE nogrid=0
IF keyword_set(badpix) THEN badpix=badpix ELSE badpix=0

;
; check if using 2D image or 3x1D arrays, xbin,ybin,moment
;
IF KEYWORD_SET(img) THEN BEGIN
	simg=SIZE(img)
	IF KEYWORD_SET(error) THEN error=error ELSE error=img*0+1
	IF N_ELEMENTS(scale) EQ 0 THEN scale = 1  ; pixel size of input image
ENDIF ELSE BEGIN
;	img=0
	IF N_ELEMENTS(error) EQ 0 THEN error = 1. + 0*moment	  ; set errors
	IF N_ELEMENTS(scale) EQ 0 THEN scale = 0.8  ; pixel size of SAURON maps
ENDELSE

;
; change to pixel scale
;
xbar=xbin/scale
ybar=ybin/scale

;
; setting radii
;
IF KEYWORD_SET(RADIUS) THEN BEGIN
	rad = radius 
	nrad = N_ELEMENTS(rad)
ENDIF ELSE BEGIN
	IF keyword_set(NRAD) THEN nrad=nrad ELSE nrad = 100 ; Large value: it will be truncated anyway
	pix = FINDGEN(nrad)
	rad = pix + (1.1^pix)
ENDELSE

;
; The central pixel is left unchanged in the reconstruction.
; Shifting of the first radius in case of a central hole.
;
IF KEYWORD_SET(RING) THEN BEGIN
	rad = ring/scale + rad 
	xellip = 0.
	yellip = 0.
	vv = 0.
	vrec = 0.
ENDIF ELSE BEGIN
   IF KEYWORD_SET(IMG) THEN BEGIN
	xellip = X0
	yellip = Y0
	vv = img[round(x0), round(y0)]
	vrec = img[round(x0), round(y0)]
   ENDIF ELSE BEGIN
	mini=MIN(sqrt(xbar^2+ybar^2))
	ww = WHERE(sqrt(xbar^2+ybar^2) eq mini)
	xellip = xbin[ww]
	yellip = ybin[ww]
	vv = moment[ww]
	vrec = moment[ww]
   ENDELSE 
ENDELSE


;
; Initialised vectors of results
;
pa    = FLTARR(nrad)
q     = pa
cf    = FLTARR(nrad,ntrm+1)
xc    = pa
yc    = pa
er_cf = cf
er_pa = pa
er_q  = q
er_xc = pa
er_yc = pa

;IF keyword_set(FIXCEN) or not keyword_set(EVEN) THEN BEGIN
;  parinfo = REPLICATE({step:1.0,limits:[0d,0d],limited:[1,1]}, 2)
;  IF keyword_set(RANGEPA) THEN parinfo[0].limits = [RANGEPA[0],RANGEPA[1]] $
;                       ELSE parinfo[0].limits = [-95d,95d]  ; PA limits in degrees
;  IF keyword_set(RANGEQ) THEN parinfo[1].limits = [RANGEQ[0],RANGEQ[1]] $
;                       ELSE parinfo[1].limits = [0.2d,1d]   ; q limits
;  parinfo[0].step = 0.5d  ; Step in degrees (of the order of the expected accuracy)
;  parinfo[1].step = 0.01d ; q
;ENDIF ELSE BEGIN
;  parinfo = REPLICATE({step:1.0,limits:[0d,0d],limited:[1,1], fixed:0}, 4)
;  IF keyword_set(RANGEPA) THEN parinfo[0].limits = [RANGEPA[0],RANGEPA[1]] $
;                          ELSE parinfo[0].limits = [-95d,95d]  ; PA limits in degrees
;  IF keyword_set(RANGEQ) THEN parinfo[1].limits = [RANGEQ[0],RANGEQ[1]] $
;                         ELSE parinfo[1].limits = [0.01d,1d]   ; q limits
;  parinfo[0].step = 0.5d  ; Step in degrees (of the order of the expected accuracy)
;  parinfo[1].step = 0.01d ; q
;  parinfo[2].step = 0.1d
;  parinfo[3].step = 0.1d
;  parinfo[2].limits = [X0-100d,X0+100d] 
;  parinfo[3].limits = [Y0-100d,Y0+100d] 
;ENDELSE
;
; Initialises parameters for MPFIT
;
IF odd eq 1. THEN BEGIN
  parinfo = REPLICATE({step:1.0,limits:[0d,0d],limited:[1,1]}, 2)
  IF keyword_set(RANGEPA) THEN parinfo[0].limits = [RANGEPA[0],RANGEPA[1]] $
                       ELSE parinfo[0].limits = [-95d,95d]  ; PA limits in degrees
  IF keyword_set(RANGEQ) THEN parinfo[1].limits = [RANGEQ[0],RANGEQ[1]] $
                       ELSE parinfo[1].limits = [0.2d,1d]   ; q limits
  parinfo[0].step = 0.5d  ; Step in degrees (of the order of the expected accuracy)
  parinfo[1].step = 0.01d ; q
ENDIF 
;
; keep centre fixed in photometry
;
IF keyword_set(EVEN) and keyword_set(FIXCEN) THEN BEGIN
  parinfo = REPLICATE({step:1.0,limits:[0d,0d],limited:[1,1]}, 2)
  IF keyword_set(RANGEPA) THEN parinfo[0].limits = [RANGEPA[0],RANGEPA[1]] $
                       ELSE parinfo[0].limits = [-95d,95d]  ; PA limits in degrees
  IF keyword_set(RANGEQ) THEN parinfo[1].limits = [RANGEQ[0],RANGEQ[1]] $
                       ELSE parinfo[1].limits = [0.2d,1d]   ; q limits
  parinfo[0].step = 0.5d  ; Step in degrees (of the order of the expected accuracy)
  parinfo[1].step = 0.01d ; q
ENDIF 
;
; fit for centre in photometry
;
IF keyword_set(EVEN) and not keyword_set(FIXCEN) and keyword_set(img) THEN BEGIN
  parinfo = REPLICATE({step:1.0,limits:[0d,0d],limited:[1,1], fixed:0}, 4)
  IF keyword_set(RANGEPA) THEN parinfo[0].limits = [RANGEPA[0],RANGEPA[1]] $
                          ELSE parinfo[0].limits = [-95d,95d]  ; PA limits in degrees
  IF keyword_set(RANGEQ) THEN parinfo[1].limits = [RANGEQ[0],RANGEQ[1]] $
                         ELSE parinfo[1].limits = [0.01d,1d]   ; q limits
  parinfo[0].step = 0.5d  ; Step in degrees (of the order of the expected accuracy)
  parinfo[1].step = 0.01d ; q
  parinfo[2].step = 0.1d
  parinfo[3].step = 0.1d
  parinfo[2].limits = [X0-100d,X0+100d] 
  parinfo[3].limits = [Y0-100d,Y0+100d] 
ENDIF
;
; case for higher even moments, e.g. velocity dispersion
;
IF keyword_set(EVEN) and not keyword_set(FIXCEN) THEN BEGIN
  parinfo = REPLICATE({step:1.0,limits:[0d,0d],limited:[1,1]}, 2)
  IF keyword_set(RANGEPA) THEN parinfo[0].limits = [RANGEPA[0],RANGEPA[1]] $
                       ELSE parinfo[0].limits = [-95d,95d]  ; PA limits in degrees
  IF keyword_set(RANGEQ) THEN parinfo[1].limits = [RANGEQ[0],RANGEQ[1]] $
                       ELSE parinfo[1].limits = [0.2d,1d]   ; q limits
  parinfo[0].step = 0.5d  ; Step in degrees (of the order of the expected accuracy)
  parinfo[1].step = 0.01d ; q
ENDIF 




;
; Set grids for global minimization
;
IF N_ELEMENTS(npa) EQ 0 THEN npa = 21
IF N_ELEMENTS(nq) EQ 0 THEN nq = 21
IF N_ELEMENTS(rangeQ) EQ 0 THEN $
            IF KEYWORD_SET(IMG) THEN rangeQ=[0.1,1.0] ELSE rangeQ=[0.2,1.0] 
IF N_ELEMENTS(rangePA) EQ 0 THEN rangePA=[-90.,90.] 

;
; speed up the calcultions, without grid, but with initial PA and q values
; but put PAQ=0 
;
IF N_ELEMENTS(PAQ) GT 0 THEN BEGIN
	IF nogrid eq 1 THEN BEGIN
		pa_mpf=PAQ[0]
		q_mpf=PAQ[1]
		PAQ=0
	ENDIF 
ENDIF ELSE BEGIN
  PAQ=0
ENDELSE



pa_grid = RANGE(rangePA[0],rangePA[1],npa) # REPLICATE(1,nq)
q_grid = REPLICATE(1,npa) # RANGE(rangeQ[0],rangeQ[1],nq)
chi2_grid = FLTARR(nPa,nq)
k = 0
;
; loop over radii
;
FOR i=0, nrad-1 DO BEGIN
    ;
    ; check if PA and q are set constant for the whole map (PAQ)
    ;
    IF N_ELEMENTS(PAQ) EQ 1 THEN BEGIN
        ;
        ; Perform brute-force *global* minimization of chi2 on a regular grid. 
        ; This is needed as the problem can have multiple minima. The result
        ; are PA and Q which will be used as initial values for MPFIT.
        ; (nogrid eq 0 means that keyword nogrid was not used, hence
        ; minimization is done one the grid...)
       IF nogrid eq 0 THEN BEGIN
          FOR j=0,npa*nq-1 DO $
	            chi2_grid[j] = kinem_fitfunc_ellipse([pa_grid[j], q_grid[j], x0,y0],$
        	                   XBAR=xbar, YBAR=ybar, MOMENT=moment, NTERMS=4, RAD=rad[i], $
		                         EVEN=even, IMG=img, X0=x0, Y0=y0, BADPIX=badpix);ERROR=error, 
              tmp = MIN(chi2_grid,jmin)
	            pa_mpf=pa_grid[jmin]
	            q_mpf=q_grid[jmin]
       ENDIF ; (nogrid eq 0)
	      ;
        ; Perform least-squares minimization of the harmonic coefficients
        ; starting from the best values of the global minimization.
	      ; In case of the EVEN moment minimize a1,b1,a2,b2. In case of 
	      ; the ODD moment minimize a1,a3,b3.
        IF keyword_set(IMG) THEN BEGIN
           IF keyword_set(EVEN) THEN BEGIN
                IF keyword_set(FIXCEN) THEN BEGIN
                      par = [pa_mpf, q_mpf]
                      fa = {MOMENT:moment, XBAR:xbar, YBAR:ybar, NTERMS:4, RAD:rad[i], MPF:1, $
                            ERROR:error, EVEN:even, IMG:img, X0:x0, Y0:y0, BADPIX:badpix}
                      sol = MPFIT('kinem_fitfunc_ellipse', par, FUNCTARGS=fa, FTOL=1e-3, $
                            PARINFO=parinfo, NFEV=ncalls, BESTNORM=fmin, PERROR=err_sol, /QUIET)
                      PA_min = sol[0]
                      q_min = sol[1]
                      x0s = X0
                      y0s = Y0
                      er_PA_min = err_sol[0] 
                      er_q_min  = err_sol[1] 
                      er_x0s = 0
                      er_y0s = 0
                ENDIF ELSE BEGIN ; keyword(fixcen)
                      par = [pa_mpf, q_mpf, X0, Y0]
                      fa = {MOMENT:moment, XBAR:xbar, YBAR:ybar, NTERMS:4, RAD:rad[i], MPF:1, $
                           ERROR:error, EVEN:even, IMG:img, BADPIX:badpix}
                      sol = MPFIT('kinem_fitfunc_ellipse', par, FUNCTARGS=fa, FTOL=1e-3, $
                            PARINFO=parinfo, NFEV=ncalls, BESTNORM=fmin, PERROR=err_sol, /QUIET)
                      PA_min = sol[0]
                      q_min = sol[1]
                      x0s = sol[2]
                      y0s = sol[3]
                      er_PA_min = err_sol[0] 
                      er_q_min  = err_sol[1] 
                      er_x0s = err_sol[1]
                      er_y0s = err_sol[1]
                ENDELSE       ; keyword(fixcen) 
             ENDIF ELSE BEGIN ;keyword_set(EVEN)
                ; center fitting in odd case not implemented yet!!!!
 	              par = [pa_mpf, q_mpf]
                fa = {MOMENT:moment, XBAR:xbar, YBAR:ybar, NTERMS:4, RAD:rad[i], MPF:1, $
	                   ERROR:error, EVEN:even, IMG:img, X0:x0, Y0:y0, BADPIX:badpix}
                sol = MPFIT('kinem_fitfunc_ellipse', par, FUNCTARGS=fa, FTOL=1e-3, $
                     PARINFO=parinfo, NFEV=ncalls, BESTNORM=fmin, PERROR=err_sol, /QUIET)
                PA_min = sol[0]
                q_min = sol[1]
                x0s = x0;sol[2]
                y0s = y0;sol[3]
	        
                er_PA_min = err_sol[0]	
                er_q_min  = err_sol[1] 
                er_x0s = 0;err_sol[2]
                er_y0s = 0;err_sol[3]
           ENDELSE;keyword_set(EVEN)
	      ENDIF ELSE BEGIN ;keyword_set(IMG)
	      ;stop
	      ; when center fitting for kin data is implemented, this part can be removed
	      	      
          par = [pa_grid[jmin], q_grid[jmin]]
          fa = {MOMENT:moment, XBAR:xbar, YBAR:ybar, NTERMS:4, RAD:rad[i], MPF:1, $
  	            ERROR:error, EVEN:even}
          sol = MPFIT('kinem_fitfunc_ellipse', par, FUNCTARGS=fa, FTOL=1e-3, $
                PARINFO=parinfo, NFEV=ncalls, BESTNORM=fmin, PERROR=err_sol, /QUIET)
          PA_min = sol[0]
          q_min = sol[1]
	        er_PA_min = err_sol[0]	
	        er_q_min  = err_sol[1] 
	        X0s=0 & Y0s=0
	      ENDELSE ;keyword_set(IMG)

    ENDIF ELSE BEGIN ; PA and Q are set to fixed values: skip all minimization

	      ;
	      ; check if PAQ is an array of (nrad*2) values or has just 2 values
	      ;	      
	      IF N_ELEMENTS(PAQ) gt 2 THEN BEGIN
        	pa_min = paq[k]
        	q_min = paq[k+1]
		      k = k + 2
		      ;IF keyword_set(IMG) and keyword_set(EVEN) THEN BEGIN  ; keep center
                  X0s=X0 & Y0s=Y0
                  er_x0s=0 & er_y0s=0
		      ;ENDIF
	      ENDIF ELSE BEGIN
        	pa_min = paq[0]
        	q_min = paq[1]
          ;IF keyword_set(IMG) and keyword_set(EVEN) THEN BEGIN  ; keep center
                  X0s=X0 & Y0s=Y0
                  er_x0s=0 & er_y0s=0
          ;ENDIF
	      ENDELSE
	      er_PA_min = 0.
	      er_q_min  = 0. 

    ENDELSE ;(N_ELEMENTS(PAQ) EQ 1)

    ;
    ; Final harmonic expansion along the best fitting ellipse
    ; 	  
	  IF keyword_set(IMG) and keyword_set(EVEN) THEN BEGIN
        tmp = kinem_fitfunc_ellipse([pa_min, q_min, X0s, Y0s], $
              XBAR=xbar, YBAR=ybar, MOMENT=moment, NTERMS=ntrm, RAD=rad[i], COEFF=coeff, $
              XELL=xEll, YELL=yEll, MOMELL=momEll, MOMFIT=recon, THETA=theta, W=w, $
              ERROR=error, ER_MOMELL=er_momEll, ALL=all, EVEN=even, VSYS=vsys, IMG=img, BADPIX=badpix)
    ENDIF ELSE BEGIN
        tmp = kinem_fitfunc_ellipse([pa_min, q_min], $
              XBAR=xbar, YBAR=ybar, MOMENT=moment, NTERMS=ntrm, RAD=rad[i], COEFF=coeff, $
              XELL=xEll, YELL=yEll, MOMELL=momEll, MOMFIT=recon, THETA=theta, W=w, $
              ERROR=error, ER_MOMELL=er_momEll, ALL=all, EVEN=even, VSYS=vsys, IMG=img, X0=x0, Y0=y0)        
    ENDELSE
    ;	
    ; errors of harmonic coeffs
    ;
    er_coeff = coeff[ntrm+1 : 2*ntrm+1]
    coeff = coeff[0:ntrm]


    IF KEYWORD_SET(verbose) THEN BEGIN
           IF i eq 0 THEN print, '       Radius,      RAD,     PA,        Q,        Xcen,     Ycen,  # of ellipse elements'
	         print, format='(1i3, 1a11, 1f10.3, 4f10.3, 3i6)', i, '-th radius: ', $
		              rad[i]*scale, pa_min, q_min, x0s,y0s, N_elements(w)
    ENDIF
    ;
    ; Stops the fit when there are less than 3/4 of the pixels sampled
    ; along the best fitting ellipse. Use COVER keyword to relax this 
    ; condition (e.g. when also setting PAQ). In case of IMG keyword
    ; stop when lenght of the ellipse semi-major axis is 1.5 X size of 
    ; the image. Stop also if keyword sky is set and the intensity
    ; is 2xSKY.
    ;
    IF N_ELEMENTS(w) LE N_ELEMENTS(xEll)*cover THEN break
;    IF keyword_set(IMG) THEN IF rad[i+1] gt simg[1]/2. THEN break
    IF KEYWORD_SET(IMG) and not KEYWORD_SET(RADIUS) THEN IF rad[i+1] gt simg[1]/2. THEN break
    IF KEYWORD_SET(sky) THEN IF coeff[0] lt 0.5*sky THEN break;2*sky


    ;
    ; Assigning of outputs
    ;
    pa[i] = pa_min
    q[i] = q_min
    cf[i,*] = coeff

    er_cf[i,*] = er_coeff
    er_pa[i] = er_pa_min
    er_q[i] = er_q_min

    IF KEYWORD_SET(IMG) THEN BEGIN
	     xc[i]=x0s
	     yc[i]=y0s
	     er_xc[i]=er_x0s
    	 er_yc[i]=er_y0s
    ENDIF

    IF KEYWORD_SET(vsys) THEN BEGIN
			IF vsys ne 1 THEN cf[*,0] = Vsys ELSE cf[*,0] = 0.
			er_cf[*,0] = 0.
    ENDIF
    
    ;
    ; reconstruction of moments
    ;
    xellip = [xellip,xell[w]]
    yellip = [yellip,yell[w]]
    IF KEYWORD_SET(even) THEN vv = [vv, coeff[0] + 0.*coeff[2]*cos(theta[w])] $
                         ELSE vv = [vv, coeff[0] + coeff[2]*cos(theta[w])]
    vrec = [vrec, recon]

;plot, vv, /ynozero
;stop




    ;
    ; Optional plotting
    ;
    IF KEYWORD_SET(plot) THEN BEGIN
        !P.MULTI=[0,1,4]
        sauron_colormap

    IF KEYWORD_SET(EVEN) THEN BEGIN
		  IF KEYWORD_SET(IMG) THEN BEGIN
			   x1=ROUND(x0s-1.5*rad[i])
			   x2=ROUND(x0s+1.5*rad[i])
			   y1=ROUND(y0s-1.5*rad[i])
			   y2=ROUND(y0s+1.5*rad[i])
			   xcnew = (x0s+1.5*rad[i] - (x0s-1.5*rad[i]))/2.
			   ycnew = (y0s+1.5*rad[i] - (y0s-1.5*rad[i]))/2.
			   IF x1 ge 0 AND x2 le simg[1] and y1 ge 0 and y2 le simg[2] THEN BEGIN
				      img1=img[x1:x2,y1:y2]
	        		mn = MIN(SIGRANGE(img1),MAX=mx)
				      levels = mn + (mx-mn) * FINDGEN(10)/(10-1)
				      CONTOUR, img1>mn<mx, /ISO, LEVELS=levels, $
		    			       /XSTYLE, /YSTYLE, XTITLE='arcsec', YTITLE='arcsec'
		        	       PLOTS, xEll*scale- X0S + xcnew, yEll*scale- y0S + ycnew, COLOR=200;0;255;0
		        	       OPLOT, xcnew+[-Rad[i]*scale,Rad[i]*scale]*COS(PA[i]/!RADEG), ycnew+[-rad[i]*scale,rad[i]*scale]*SIN(pa[i]/!RADEG), COLOR=200;255
			   ENDIF ELSE BEGIN
	        		mn = MIN(SIGRANGE(img),MAX=mx)
				      levels = mn + (mx-mn) * FINDGEN(10)/(10-1)
				      CONTOUR, img>mn<mx, /ISO, LEVELS=levels, $
		    			       /XSTYLE, /YSTYLE, XTITLE='arcsec', YTITLE='arcsec'
		        	       PLOTS, xEll*scale, yEll*scale, COLOR=200;0;255;0
		        	       OPLOT, x0s+[-Rad[i]*scale,Rad[i]*scale]*COS(PA[i]/!RADEG), y0s+[-rad[i]*scale,rad[i]*scale]*SIN(pa[i]/!RADEG), COLOR=200;255
			   ENDELSE	 ;(x1 ge 0 AND x1 le simg[1] and y1 ge 0 and y2 le simg[2])	
		  ENDIF ELSE BEGIN;KEYWORD_SET(img)
	       mn = MIN(alog(moment),MAX=mx)
			   levels = mn + (mx-mn) * FINDGEN(64)/(64-1)
		 	   CONTOUR, alog(moment)>mn<mx, xbin, ybin, /FILL, /ISO, LEVELS=levels, $
		         	/XSTYLE, /YSTYLE, /IRREGULAR, XTITLE='arcsec', YTITLE='arcsec'       
         PLOTS, xEll*scale, yEll*scale, COLOR=200;0;255;0
         OPLOT, x0s+[-Rad[i]*scale,Rad[i]*scale]*COS(PA[i]/!RADEG), y0s+[-rad[i]*scale,rad[i]*scale]*SIN(pa[i]/!RADEG), COLOR=200;255
 		  ENDELSE ;KEYWORD_SET(img)
 ;stop		  
    ENDIF ELSE BEGIN ;KEYWORD_SET(EVEN) 
    
       IF KEYWORD_SET(IMG) THEN BEGIN
         x1=ROUND(x0s-1.5*rad[i])
         x2=ROUND(x0s+1.5*rad[i])
         y1=ROUND(y0s-1.5*rad[i])
         y2=ROUND(y0s+1.5*rad[i])
         xcnew = (x0s+1.5*rad[i] - (x0s-1.5*rad[i]))/2.
         ycnew = (y0s+1.5*rad[i] - (y0s-1.5*rad[i]))/2.
         ;IF x1 ge 0 AND x1 le simg[1] and y1 ge 0 and y2 le simg[2] THEN BEGIN
         ;     img1=img[x1:x2,y1:y2]
         ;     mn = MIN(SIGRANGE(img1),MAX=mx)
         ;     levels = mn + (mx-mn) * FINDGEN(10)/(10-1)
         ;     CONTOUR, img1>mn<mx, /ISO, LEVELS=levels, /FILL, $
         ;            /XSTYLE, /YSTYLE, XTITLE='arcsec', YTITLE='arcsec'
         ;            PLOTS, xEll*scale- X0S + xcnew, yEll*scale- y0S + ycnew, COLOR=200;0;255;0
         ;            OPLOT, xcnew+[-Rad[i]*scale,Rad[i]*scale]*COS(PA[i]/!RADEG), ycnew+[-rad[i]*scale,rad[i]*scale]*SIN(pa[i]/!RADEG), COLOR=200;255
         ;ENDIF ELSE BEGIN
              mn = MIN(SIGRANGE(img),MAX=mx)
              levels = mn + (mx-mn) * FINDGEN(10)/(10-1) 
              CONTOUR, img>mn<mx, /ISO, LEVELS=levels, /FILL, $
                     /XSTYLE, /YSTYLE, XTITLE='arcsec', YTITLE='arcsec'
                     PLOTS, xEll*scale, yEll*scale, COLOR=200;0;255;0
                     OPLOT, x0s+[-Rad[i]*scale,Rad[i]*scale]*COS(PA[i]/!RADEG), y0s+[-rad[i]*scale,rad[i]*scale]*SIN(pa[i]/!RADEG), COLOR=200;255
         ;ENDELSE   ;(x1 ge 0 AND x1 le simg[1] and y1 ge 0 and y2 le simg[2]) 
       ENDIF ELSE BEGIN;KEYWORD_SET(img)
	        mn = MIN(SIGRANGE(moment-median(moment)),MAX=mx)
          mx = mx < abs(mn)
		      levels = -mx + (mx-(-mx)) * FINDGEN(64)/(64-1)
		      CONTOUR, (moment-median(moment))>(-mx)<mx, xbin, ybin, /FILL, /ISO, LEVELS=levels, $
		               /XSTYLE, /YSTYLE, /IRREGULAR, XTITLE='arcsec', YTITLE='arcsec'
	        PLOTS, xEll*scale, yEll*scale, COLOR=200;0;255;0
	        OPLOT, x0s+[-Rad[i]*scale,Rad[i]*scale]*COS(PA[i]/!RADEG), y0s+[-rad[i]*scale,rad[i]*scale]*SIN(pa[i]/!RADEG), COLOR=200;255
	     ENDELSE;KEYWORD_SET(img)
	   ENDELSE;KEYWORD_SET(EVEN)

        LOADCT, 4, /silent
        IF keyword_set(nogrid) THEN plot, pa_grid, q_grid, PSYM=3, THICK=2, XTITLE='PA (deg)', YTITLE='q' $
        ELSE BEGIN
                  CONTOUR, ALOG(chi2_grid), pa_grid, q_grid, NLEVELS=15, XTITLE='PA (deg)', YTITLE='q', /FILL
                  PLOTS, pa_grid, q_grid, PSYM=3, THICK=2
        ENDELSE                            
        PLOTS, pa_min, q_min, PSYM=4, COLOR=125, THICK=3

        LOADCT, 12, /silent
        !Y.MARGIN=[0,0] ; show plots with shared X axis
        !Y.OMARGIN=[5,10] ; allow for space for the axis labels

	IF KEYWORD_SET(EVEN) THEN BEGIN
	  IF N_ELEMENTS(error) gt 1 THEN $
	        PLOTERROR, theta[w]*!RADEG, momEll, er_momEll, PSYM=-1, TITLE=rad[i]*scale, XTICKFORMAT='(A1)', YTITLE='V' $
	  ELSE $
        	PLOT, theta[w]*!RADEG, momEll, PSYM=-1, TITLE=rad[i], XTICKFORMAT='(A1)', YTITLE='V'
          OPLOT, theta*!RADEG, coeff[0]*(theta*0+1.), PSYM=-4, COLOR=200
	  mnp=min(SIGRANGE(momEll-coeff[0]), max=mxp)
          PLOT, theta[w]*!RADEG, SIGRANGE(momEll-coeff[0]), PSYM=-4, XTITLE='!7h!X!N', YTITLE='residuals', YRANGE=[mnp,mxp]
          OPLOT, theta*!RADEG, coeff[1]*SIN(theta) + coeff[2]*COS(theta) + coeff[3]*SIN(2*theta) + coeff[4]*COS(2*theta), PSYM=-1, COLOR=40
	  IF NTRM ge 8 THEN $
	          OPLOT, theta*!RADEG, coeff[1]*SIN(theta) + coeff[2]*COS(theta) + coeff[3]*SIN(2*theta) + coeff[4]*COS(2*theta) $
                       + coeff[7]*sin(4*theta) + coeff[8]*cos(4*theta), psym=-1, COLOR=200	

	ENDIF ELSE BEGIN

	  IF KEYWORD_SET(all) THEN BEGIN
		a1 = coeff[1] & b1 = coeff[2] & a3 = coeff[5] & b3 = coeff[6]
	  ENDIF ELSE BEGIN
		a1 = coeff[1] & b1 = coeff[2] & a3 = coeff[3] & b3 = coeff[4]
	  ENDELSE		

	  IF N_ELEMENTS(error) gt 1 THEN $
	        PLOTERROR, theta[w]*!RADEG, momEll, er_momEll, PSYM=-1, TITLE=rad[i]*scale, XTICKFORMAT='(A1)', YTITLE='V' $
	  ELSE $
        	PLOT, theta[w]*!RADEG, momEll, PSYM=-1, TITLE=rad[i], XTICKFORMAT='(A1)', YTITLE='V'
          OPLOT, theta*!RADEG, coeff[0]+b1*COS(theta), PSYM=-4, COLOR=200

          PLOT, theta[w]*!RADEG, SIGRANGE(momEll-coeff[0]-b1*COS(theta[w])), PSYM=-4, XTITLE='!7h!X!N', YTITLE='residuals', YRANGE=[-20,20]
          OPLOT, theta*!RADEG, a1*SIN(theta) + a3*SIN(3*theta) + b3*COS(3*theta), PSYM=-1, COLOR=40
	  
	  IF NTRM ge 6 THEN BEGIN
		IF KEYWORD_SET(all) and NTRM ge 9 THEN BEGIN
			a5 = coeff[9] & b5 = coeff[10]
		ENDIF ELSE BEGIN
	 		a5 = coeff[5] & b5 = coeff[6]
	  	ENDELSE		
          	OPLOT, theta*!RADEG, a1*SIN(theta) + a3*SIN(3*theta) + b3*COS(3*theta) $
                     + a5*sin(5*theta) + b5*cos(5*theta), psym=-1, COLOR=200	
	  ENDIF
        ENDELSE
        !Y.MARGIN=[4,2] ; back to default values
        !Y.OMARGIN=[0,0]
        !P.MULTI=0
    ENDIF
    ;stop
ENDFOR  ; end of loop over radii



;
; Final outputs (back to physical scale (arcsec))
;
xellip=xellip*scale
yellip=yellip*scale
rad   = rad[0:i-1]*scale
pa    = pa[0:i-1]
q     = q[0:i-1]
cf    = cf[0:i-1,*]
PAfix=(median(PA))/!radeg   ; used for reconstruction below, defined now before PA rotated N-E orientation

IF KEYWORD_SET(IMG) THEN BEGIN
	xc    = xc[0:i-1]
	xy    = yc[0:i-1]
ENDIF

IF N_ELEMENTS(error) gt 1 THEN BEGIN
	er_cf = er_cf[0:i-1,*]
	er_pa = er_pa[0:i-1]
	er_q  = er_q[0:i-1]
    IF KEYWORD_SET(IMG) THEN BEGIN
	xc    = xc[0:i-1]
	xy    = yc[0:i-1]
   ENDIF
ENDIF

;
; PA correction
;
IF N_ELEMENTS(PAQ) EQ 1 THEN BEGIN

	PA=270 + PA 		     	 ; PA measured East of North
	wb=WHERE(cf[*,2] lt 0., red) 	 ; PA measured from receding side of galaxy (positive vel)
	;
	; correction of coeffs. to standrd system: odd coeffs are measured from
	; the receding side of galaxy. This is done only for maps of odd kinematic 
	; moments and not when keyword EVEN is set
	;
	IF red ne 0 AND even eq 0 THEN BEGIN 		
		PA[wb]=pa[wb]-180
;		IF keyword_set(all) THEN BEGIN
;			FOR i=1, ntrm, 4 DO BEGIN
;				cf[*,i]   = -cf[*,i]
;				cf[*,i+1] = -cf[*,i+1]
;			ENDFOR
;		ENDIF ELSE cf = -cf
	ENDIF

ENDIF


;
; calculation of the circular velocity (not for the photometry mode): 
; (VELCIRC and VELKIN are calculated above before plotting and final outputs)
;   - fixed PA and q
;   - using only cf[*,2] terms (cosine terms)
;   - xellipF and yellipF are new ellipes along which gascirc is calculated
;     Each ellipse has 100 points - they are different from xellip,yellip
;
qfix=median(q)
xellipF = 0.
yellipF = 0.
vvF=0.
FOR i=0, N_elements(rad)-1 DO BEGIN
        theta = range(0.0,2.0*!pi,100) 
        xF = rad[i]*COS(theta) 
        yF = rad[i]*SIN(theta)*qfix
        xEllF = xF*COS(PAfix) - yF*SIN(PAfix) + X0s
        yEllF = xF*SIN(PAfix) + yF*COS(PAfix) + Y0s
        xellipF = [xellipF,xellF]
        yellipF = [yellipF,yellF]
        vvF = [vvF, cf[i,2]*cos(theta)]
ENDFOR 

;
; reconstruction of kinematic moment maps
;
IF KEYWORD_SET(bmodel) THEN BEGIN


   IF KEYWORD_SET(IMG) THEN BEGIN
	    TRIANGULATE, xellip, yellip, tr
	    velcirc= TRIGRID(xellip, yellip, vv, tr, xout=findgen(simg[1]), yout=findgen(simg[2]))
	    velKIn= TRIGRID(xellip, yellip, vrec, tr, xout=findgen(simg[1]), yout=findgen(simg[2]))
      TRIANGULATE, xellipF, yellipF, trF
      gascirc= TRIGRID(xellipF, yellipF, vvF, trF, xout=findgen(simg[1]), yout=findgen(simg[2]))

   ENDIF ELSE BEGIN 
	    KINEM_TRIGRID_IRREGULAR, xellip, yellip, vv, xbin, ybin, velcirc, MISSING=123456789
	    KINEM_TRIGRID_IRREGULAR, xellip, yellip, vrec, xbin, ybin, velkin, MISSING=12345789
	    KINEM_TRIGRID_IRREGULAR, xellipF, yellipF, vvF, xbin, ybin, GAScirc, MISSING=12345789
	    
	    
      ;
      ; Reconstruction cosmetic lines: un-commenting these lines will assign 
      ; original values to bins that are outside analysed region. *Not* recommended
      ; except for presentation purposes, since there is no science
      ; associated and may introduce confusion when interpreting results.
      ;
      ;wf = WHERE(velCirc eq -12345, elem)
      ;IF elem ne 0 THEN velCirc[wf] = velBin[wf] ; Leave unchanged values not fitted      
      ;IF elem ne 0 THEN velKin[wf] = velBin[wf]
      
   ENDELSE
       

ENDIF ELSE BEGIN ;KEYWORD_SET(bmodel)
	velcirc=vv	; reconstruction along ellipses
	velkin=vrec
	;gascirc=vvF ; disabled because xellipF, yellipF are different from xellip,yellip which are the actual outputs
  ;
  ; reconstruction of kinematic moment map from the best fitting ellipses
  ;
  ;KINEM_TRIGRID_IRREGULAR, xellip, yellip, vv, xbin, ybin, velCirc, MISSING=12345
  ;KINEM_TRIGRID_IRREGULAR, xellip, yellip, vrec, xbin, ybin, velKin
ENDELSE
;stop

END
;----------------------------------------------------------------------
