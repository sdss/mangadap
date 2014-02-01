;######################################################################
;
; Copyright (C) 2001-2006, Michele Cappellari
; E-mail: cappellari_at_astro.ox.ac.uk
;
; For details on the Voronoi binning method see:
;   Cappellari M., Copin Y., 2003, MNRAS, 342, 345
;
; Updated versions of the software are available from my web page
; http://www-astro.physics.ox.ac.uk/~mxc/idl/
;
; If you have found this software useful for your
; research, we would appreciate an acknowledgment to use of
; `the Voronoi 2D-binning method by Cappellari & Copin (2003)'.
;
; This software is provided as is without any warranty whatsoever.
; Permission to use, copy and distribute unmodified copies for
; non-commercial purposes is granted. Permission to modify for
; personal or internal use is granted, provided this copyright
; and disclaimer are included unchanged at the beginning of
; the file. All other rights are reserved.
;
;######################################################################
;+
; NAME:
;     VORONOI_2D_BINNING
;
; AUTHOR:
;       Michele Cappellari, University of Oxford
;       cappellari_at_astro.ox.ac.uk
;
; PURPOSE:
;       Perform adaptive spatial binning of Integral-Field Spectroscopic
;       (IFS) data to reach a chosen constant signal-to-noise ratio per bin.
;       This method is required for the proper analysis of IFS
;       observations, but can also be used for standard photometric
;       imagery or any other two-dimensional data.
;       This program precisely implements the algorithm described in
;       section 5.1 of the reference below.
;
; EXPLANATION:
;       Further information on VORONOI_2D_BINNING algorithm can be found in
;       Cappellari M., Copin Y., 2003, MNRAS, 342, 345
;
; CALLING SEQUENCE:
;       VORONOI_2D_BINNING, X, Y, Signal, Noise, targetSN, $
;               BinNumber, xBin, yBin, xBar, yBar, SN, NPixels, scale, $
;               /NO_CVT, /PLOT, /QUIET, /WVT
;
; INPUTS:
;            X: Vector containing the X coordinate of the pixels to bin.
;               Arbitrary units can be used (e.g. arcsec or pixels).
;               In what follows the term ?pixel? refers to a given
;               spatial element of the dataset (sometimes called ?spaxel? in
;               the IFS community): it can be an actual pixel of a CCD
;               image, or a spectrum position along the slit of a long-slit
;               spectrograph or in the field of view of an IFS
;               (e.g. a lenslet or a fiber).
;               It is assumed here that pixels are arranged in a regular
;               grid, so that the pixel size is a well defined quantity.
;               The pixel grid however can contain holes (some pixels can be
;               excluded from the binning) and can have an irregular boundary.
;               See the above reference for an example and details.
;            Y: Vector (same size as X) containing the Y coordinate
;               of the pixels to bin.
;       SIGNAL: Vector (same size as X) containing the signal
;               associated with each pixel, having coordinates (X,Y).
;               If the `pixels' are actually the apertures of an
;               integral-field spectrograph, then the signal can be
;               defined as the average flux in the spectral range under
;               study, for each aperture.
;               If pixels are the actual pixels of the CCD in a galaxy
;               image, the signal will be simply the counts in each pixel.
;        NOISE: Vector (same size as X) containing the corresponding
;               noise (1 sigma error) associated with each pixel.
;     TARGETSN: The desired signal-to-noise ratio in the final
;               2D-binned data. E.g. a S/N~50 per pixel may be a
;               reasonable value to extract stellar kinematics
;               information from galaxy spectra.
;
; KEYWORDS:
;      /NO_CVT: Set this keyword to skip the Centroidal Voronoi Tessellation
;               (CVT) step (vii) of the algorithm in Section 5.1 of
;               Cappellari & Copin (2003).
;               This may be useful if the noise is strongly non Poissonian,
;               the pixels are not optimally weighted, and the CVT step
;               appears to introduces significant gradients in the S/N.
;               A similar alternative consists of using the /WVT keyword below.
;        /PLOT: Set this keyword to produce a plot of the two-dimensional
;               bins and of the corresponding S/N at the end of the
;               computation.
;       /QUIET: by default the program shows the progress while accreting
;               pixels and then while iterating the CVT. Set this keyword
;               to avoid printing progess results.
;         /WVT: When this keyword is set, the routine bin2d_cvt_equal_mass is
;               modified as proposed by Diehl & Statler (2006, MNRAS, 368, 497).
;               In this case the final step of the algorithm, after the bin-accretion
;               stage, is not a modified Centroidal Voronoi Tessellation, but it uses
;               a Weighted Voronoi Tessellation.
;               This may be useful if the noise is strongly non Poissonian,
;               the pixels are not optimally weighted, and the CVT step
;               appears to introduces significant gradients in the S/N.
;               A similar alternative consists of using the /NO_CVT keyword above.
;               If you use the /WVT keyword you should also include a reference to
;               `the WVT modification proposed by Diehl & Statler (2006).'
;
; OUTPUTS:
;       BINNUMBER: Vector (same size as X) containing the bin number assigned
;               to each input pixel. The index goes from zero to Nbins-1.
;               This vector alone is enough to make *any* subsequent
;               computation on the binned data. Everything else is optional!
;         XBIN: Vector (size Nbins) of the X coordinates of the bin generators.
;               These generators uniquely define the Voronoi tessellation.
;         YBIN: Vector (size Nbins) of Y coordinates of the bin generators.
;         XBAR: Vector (size Nbins) of X coordinates of the bins luminosity
;               weighted centroids. Useful for plotting interpolated data.
;         YBAR: Vector (size Nbins) of Y coordinates of the bins luminosity
;               weighted centroids.
;           SN: Vector (size Nbins) with the final SN of each bin.
;      NPIXELS: Vector (size Nbins) with the number of pixels of each bin.
;        SCALE: Vector (size Nbins) with the scale length of the Weighted
;               Voronoi Tessellation, when the /WVT keyword is set.
;               In that case SCALE is *needed* together with the coordinates
;               XBIN and YBIN of the generators, to compute the tessellation
;               (but one can also simply use the BINNUMBER vector).
;
; PROCEDURES USED:
;       The following procedures are contained in the main VORONOI_2D_BINNING program.
;           WEIGHTED_CENTROID -- computes weighted centroid of one bin
;           BIN_ROUNDNESS     -- equation (5) of Cappellari & Copin (2003)
;           BIN_ACCRETION     -- steps (i)-(v) in section 5.1
;           REASSIGN_BAD_BINS -- steps (vi)-(vii) in section 5.1
;           CVT_EQUAL_MASS    -- the modified Lloyd algorithm in section 4.1
;           COMPUTE_USEFUL_BIN_QUANTITIES -- self explanatory
;           DISPLAY_PIXELS    -- plotting of colored pixels
;
; NOTE:
;       This program uses some features of IDL 5.4.
;       Please ask if you think support for IDL < 5.4 is useful.
;
; MODIFICATION HISTORY:
;       V1.0: First implementation. Michele Cappellari, Leiden, June 2001
;       V2.0: Major revisions. Stable version. MC, Leiden, 11 September 2001
;       V2.1: First released version. Written documentation.
;           MC, Vicenza, 13 February 2003
;       V2.2: Added computation of useful bin quantities in output. Deleted some
;           safety checks for zero size bins in CVT. Minor polishing of the code.
;           MC, Leiden, 11 March 2003
;       V2.3: Unified the three tests to stop the accretion of one bin.
;           This can improve some bins at the border. MC, Leiden, 9 April 2003
;       V2.31: Do *not* assume the first bin is made of one single pixel.
;           Added computation of S/N scatter and plotting of 1-pixel bins.
;           MC, Leiden, 13 April 2003
;       V2.4: Addedd basic error checking of input S/N. Reintroduced the
;           treatment for zero-size bins in CVT, which was deleted in V2.2.
;           Thanks to Robert Sharp and Kambiz Fathi for reporting problems.
;           MC, Leiden, 10 December 2003.
;       V2.41: Added /QUIET keyword and verbose output during the computation.
;           After suggestion by Richard McDermid. MC, Leiden, 14 December 2003
;       V2.42: Use LONARR instead of INTARR to define the CLASS vector,
;           to be able to deal with big images. Thanks to Tom Statler.
;           MC, Leiden, 4 August 2004
;       V2.43: Corrected bug introduced in version 2.31. It went undetected
;           for a long time because it could only happen in special conditions.
;           Now we recompute the index of the good bins after computing all
;           centroids of the reassigned bins in reassign_bad_bins. Many thanks
;           to Simona Ghizzardi for her clear analysis of the problem and
;           the solution. MC, Leiden, 29 November 2004
;       V2.44: Prevent division by zero for pixels with signal=0
;           and noise=sqrt(signal)=0, as can happen from X-ray data.
;           MC, Leiden, 30 November 2004
;       V2.45: Added BIN2D prefix to internal routines to avoid possible
;           naming conflicts. MC, Leiden, 3 December 2004
;       V2.46: Added /NO_CVT keyword to optionally skip the CVT step of
;           the algorithm. MC, Leiden, 27 August 2005
;       V2.47: Verify that SIGNAL and NOISE are non negative vectors.
;           MC, Leiden, 27 September 2005
;       V2.48: Use geometric centroid of a bin during the bin-accretion stage,
;           to allow the routine to deal with negative signal (e.g. in
;           background-subtracted X-ray images). Thanks to Steven Diehl for
;           pointing out the usefulness of dealing with negative signal.
;           MC, Leiden, 23 December 2005
;       V2.5: Added two new lines of code and the corresponding /WVT keyword
;           to implement the nice modification to the algorithm proposed by
;           Diehl & Statler (2006). MC, Leiden, 9 March 2006
;       V2.51: Updated documentation. MC, Oxford, 3 November 2006
;       V2.52: Print number of unbinned pixels. MC, Oxford, 28 March 2007
;       V2.53: Fixed program stop, introduced in V2.5, with /NO_CVT keyword.
;           MC, Oxford, 3 December 2007
; =================================================================================
;-      Nov 2013: ADAPTED FOR MANGA Data Analysis Pipeline by L. Coccato
;            - Recalculate targetSN to have at least 1 bin. L.Coccato Oct 2013
;            - Usage of empirical calibraion between estimated S/N
;              and real S/N (mdap_calibrate_sn.pro)
;
;----------------------------------------------------------------------------

pro bin2d_weighted_centroid, x, y, density, xBar, yBar
COMPILE_OPT IDL2, HIDDEN

; Computes weighted centroid of one bin.
; Equation (4) of Cappellari & Copin (2003)

mass = total(density)
xBar = total(x*density)/mass
yBar = total(y*density)/mass

end
;----------------------------------------------------------------------
function bin2d_roundness, x, y, pixelSize
COMPILE_OPT IDL2, HIDDEN

; Implements equation (5) of Cappellari & Copin (2003)

n = n_elements(x)
equivalentRadius = sqrt(n/!Pi)*pixelSize
xBar = mean(x) ; Geometric centroid here!
yBar = mean(y)
maxDistance = sqrt(max((x-xBar)^2 + (y-yBar)^2))
roundness = maxDistance/equivalentRadius - 1.0

return, roundness
end
;----------------------------------------------------------------------
PRO bin2d_accretion, x, y, signal, noise, targetSN, class, pixelSize, QUIET=quiet, SN_CALIBRATION=SN_CALIBRATION
COMPILE_OPT IDL2, HIDDEN

; Implements steps (i)-(v) in section 5.1 of Cappellari & Copin (2003)

n = n_elements(x)
class = lonarr(n) ; will contain the bin number of each given pixel
good = bytarr(n)  ; will contain 1 if the bin has been accepted as good

; For each point, find the distance to all other points and select the minimum.
; This is a robust but slow way of determining the pixel size of unbinned data.
;
dx = 1e30
for j=0,n-2 do dx = min((x[j]-x[j+1:*])^2 + (y[j]-y[j+1:*])^2) < dx
pixelSize = sqrt(dx)

SN = max(signal/noise, currentBin) ; Start from the pixel with highest S/N

; Rough estimate of the expected final bin number.
; This value is only used to have a feeling of the expected
; remaining computation time when binning very big dataset.
;
w = where(signal/noise lt targetSN, COMPLEMENT=w1, NCOMPLEMENT=nc)
maxnum = round(total((signal[w]/noise[w])^2)/targetSN^2) + nc


; The first bin will be assigned CLASS = 1
; With N pixels there will be at most N bins
;
for ind=1,n do begin

    if not keyword_set(quiet) then print, ind, maxnum, FORMAT='(%"Bin:  %d / %d")'

    class[currentBin] = ind ; Here currentBin is still made of one pixel
    xbar = x[currentBin]
    ybar = y[currentBin]    ; Centroid of one pixels

    while 1 do begin

        unBinned = where(class eq 0, m)
        if (m eq 0) then break ; Stops if all pixels are binned

        ; Find the unbinned pixel closest to the centroid of the current bin
        ;
        minDist = min((x[unBinned]-xBar)^2 + (y[unBinned]-yBar)^2, k)

        ; Find the distance from the closest pixel to the current bin
        ;
        minDist = min((x[currentBin]-x[unBinned[k]])^2 + (y[currentBin]-y[unBinned[k]])^2)

        ; Estimate the `roundness' of the POSSIBLE new bin
        ;
        nextBin = [currentBin,unBinned[k]]
        roundness = bin2d_roundness(x[nextBin], y[nextBin], pixelSize)

        ; Compute the S/N one would obtain by adding
        ; the CANDIDATE pixel to the current bin
        ;
        SNOld = SN
        ;--start original 
        ;SN = total(signal[nextBin])/sqrt(total(noise[nextBin]^2))
        ;--end original

        ;-- start manga
        SN_EST = total(signal[nextBin])/sqrt(total(noise[nextBin]^2))
        SN = SN_EST
        if n_elements(sn_calibration) ne 0 then SN = mdap_calibrate_sn(SN_EST,n_elements(nextBin),sn_calibration)
        
        ;-- end manga

        ; Test whether the CANDIDATE pixel is connected to the
        ; current bin, whether the POSSIBLE new bin is round enough
        ; and whether the resulting S/N would get closer to targetSN
        ;
        ;start original version
        ;if ((sqrt(minDist) gt 1.2*pixelsize) or (roundness gt 0.3) $
        ;     or (abs(SN-targetSN) gt abs(SNOld-targetSN))) then begin
        ;        if (SNOld gt 0.8*targetSN) then good[currentBin] = 1
        ;        break
        ; endif
        ;end original version

        ;start manga version
        if ((sqrt(minDist) gt 1.2*pixelsize) or (roundness gt 0.3) $
             or (abs(SN-targetSN) gt abs(SNOld-targetSN))) then begin
                if (SNOld gt 0.8*targetSN) then good[currentBin] = 1
                break
        endif
        ;end manga version

        ; If all the above tests are negative then accept the CANDIDATE pixel,
        ; add it to the current bin, and continue accreting pixels
        ;
        class[unBinned[k]] = ind
        currentBin = nextBin

        ; Update the centroid of the current bin
        ;
        xBar = mean(x[currentBin])
        yBar = mean(y[currentBin])

    endwhile

    ; Get the centroid of all the binned pixels
    ;
    unBinned = where(class eq 0, COMPLEMENT=binned, m)
    if (m eq 0) then break ; Stop if all pixels are binned
    xBar = mean(x[binned])
    yBar = mean(y[binned])

    ; Find the closest unbinned pixel to the centroid of all
    ; the binned pixels, and start a new bin from that pixel.
    ;
    minDist = min((x[unBinned]-xBar)^2 + (y[unBinned]-yBar)^2, k)
    currentBin = unBinned[k]    ; The bin is initially made of one pixel
    SN = signal[currentBin]/noise[currentBin]
    if n_elements(sn_calibration) ne 0 then sn = poly(sn^sn_calibration[0],sn_calibration[1:*])
endfor

class = class * good ; Set to zero all bins that did not reach the target S/N

END
;----------------------------------------------------------------------------
pro bin2d_reassign_bad_bins, x, y, signal, noise, targetSN, class, xnode, ynode
COMPILE_OPT IDL2, HIDDEN

; Implements steps (vi)-(vii) in section 5.1 of Cappellari & Copin (2003)

; Find the centroid of all succesful bins.
; CLASS = 0 are unbinned pixels which are excluded.
;
area = histogram(class, REVERSE_INDICES=r, MIN=1)
good = where(area gt 0, nnodes) ; Obtain the index of the good bins
xnode = fltarr(nnodes)
ynode = xnode
for j=0,nnodes-1 DO begin
    k = good[j]
    p = r[r[k]:r[k+1]-1] ; Find subscripts of pixels in bin k.
    xnode[j] = mean(x[p])
    ynode[j] = mean(y[p])
endfor

; Reassign pixels of bins with S/N < targetSN
; to the closest centroid of a good bin
;
bad = where(class eq 0, m)
for j=0,m-1 do begin
    tmp = min((x[bad[j]]-xnode)^2+(y[bad[j]]-ynode)^2, index)
    class[bad[j]] = good[index] + 1
endfor

; Recompute all centroids of the reassigned bins.
; These will be used as starting points for the CVT.
;
area = HISTOGRAM(class, REVERSE_INDICES=r)
good = where(area gt 0, nnodes) ; Re-obtain the index of the good bins
for j=0,nnodes-1 DO begin
    k = good[j]
    p = r[r[k]:r[k+1]-1] ; Find subscripts of pixels in bin k.
    xnode[j] = mean(x[p])
    ynode[j] = mean(y[p])
endfor

END
;----------------------------------------------------------------------------
PRO bin2d_cvt_equal_mass, x, y, signal, noise, xnode, ynode, scale, iter, $
    QUIET=quiet, WVT=wvt,sn_calibration=sn_calibration
COMPILE_OPT IDL2, HIDDEN

; Implements the modified Lloyd algorithm
; in section 4.1 of Cappellari & Copin (2003).
;
; NB: When the keyword /WVT is set this routine includes
; the modification proposed by Diehl & Statler (2006).

npixels = N_ELEMENTS(signal)
class = lonarr(npixels)  ; See beginning of section 4.1 of CC03
if keyword_set(wvt) then dens = 1.0 else dens = (signal/noise)^2
scale = 1.0 ; Start with the same scale length for all bins
sn = xnode

iter = 1
REPEAT BEGIN

    xnodeOld = xnode
    ynodeOld = ynode

    ; Computes (Weighted) Voronoi Tessellation of the pixels grid
    ;
    FOR j=0,npixels-1 DO BEGIN
        tmp = MIN(((x[j]-xnode)/scale)^2+((y[j]-ynode)/scale)^2, index)
        class[j] = index
    ENDFOR

    ; Computes centroids of the bins, weighted by dens^2.
    ; Exponent 2 on the density produces equal-mass Voronoi bins.
    ; The geometric centroids are computed if /WVT keyword is set.
    ;
    area = HISTOGRAM(class, REVERSE_INDICES=r)
    w = where(area gt 0, m) ; Check for zero-size Voronoi bins
    FOR j=0,m-1 DO BEGIN
        k = w[j]                 ; Only loop over nonzero bins
        index = r[r[k]:r[k+1]-1] ; Find subscripts of pixels in bin k.
        bin2d_weighted_centroid, x[index], y[index], dens[index]^2, xb, yb
        xnode[k] = xb
        ynode[k] = yb
        ;start original version
        ;sn[k] = total(signal[index])/sqrt(total(noise[index]^2))
        ;end original version

        ;start manga version
        sn_est = total(signal[index])/sqrt(total(noise[index]^2))
        sn[k] = sn_est
        if n_elements(sn_calibration) ne 0 then sn[k] =mdap_calibrate_sn(sn_est,n_elements(index),sn_calibration)
    ;end manga version


    ENDFOR

    if keyword_set(wvt) then scale = sqrt(area/sn) ; Eq. (4) of Diehl & Statler (2006)
    diff = TOTAL((xnode-xnodeOld)^2+(ynode-ynodeOld)^2)
    iter = iter + 1

    if not keyword_set(quiet) then $
        print, iter, diff, FORMAT='(%"Iter:  %d,  Diff:  %f")'

ENDREP UNTIL diff EQ 0

; Only return the generators of the nonzero Voronoi bins
;
xnode = xnode[w]
ynode = ynode[w]

END
;-----------------------------------------------------------------------
pro bin2d_compute_useful_bin_quantities, x, y, signal, noise, xnode, ynode, scale, $
    class, xbar, ybar, sn, area,SN_CALIBRATION=SN_CALIBRATION
COMPILE_OPT IDL2, HIDDEN

; Recomputes (Weighted) Voronoi Tessellation of the pixels grid to make sure
; that the class number corresponds to the proper Voronoi generator.
; This is done to take into account possible zero-size Voronoi bins
; in output from the previous CVT (or WVT).
;
npix = n_elements(x)
class = lonarr(npix) ; will contain the bin number of each given pixel
FOR j=0,npix-1 DO BEGIN
    tmp = MIN(((x[j]-xnode)/scale)^2+((y[j]-ynode)/scale)^2, index)
    class[j] = index
ENDFOR

; At the end of the computation evaluate the bin luminosity-weighted
; centroids (xbar,ybar) and the corresponding final S/N of each bin.
;
area = HISTOGRAM(class, REVERSE_INDICES=r)
nbins = n_elements(xnode)
xbar = fltarr(nbins)
ybar = xbar
sn = xbar
FOR j=0,nbins-1 DO BEGIN
    index = r[r[j]:r[j+1]-1] ; Find subscripts of pixels in bin j.
    bin2d_weighted_centroid, x[index], y[index], signal[index], xb, yb
    xbar[j] = xb
    ybar[j] = yb
    ;start original version
    ;sn[j] = total(signal[index])/sqrt(total(noise[index]^2))
    ;end original version

    ;start manga version
    sn_est = total(signal[index])/sqrt(total(noise[index]^2))
    sn[j] = sn_est
    if n_elements(sn_calibration) ne 0 then sn[j] =mdap_calibrate_sn(sn_est,n_elements(index),sn_calibration)
    ;if sn[j] le 50 then stop
    ;end manga version
ENDFOR

END
;-----------------------------------------------------------------------
pro bin2d_display_pixels, x, y, counts, pixelSize
COMPILE_OPT IDL2, HIDDEN

; Plots colored pixels with the same pixels size

PLOT, [MIN(x)-pixelSize,MAX(x)+pixelSize], [MIN(y)-pixelSize,MAX(y)+pixelSize], $
    /NODATA, /XSTYLE, /YSTYLE, XTITLE='arcsec', YTITLE='arcsec', /ISO
x1 = [-0.5, -0.5, +0.5, +0.5, -0.5] * pixelSize
y1 = [+0.5, -0.5, -0.5, +0.5, +0.5] * pixelSize
countMax = MAX(counts,MIN=countMin)
color = 255.0/(countMax-countMin)*(counts-countMin)
FOR j=0, N_ELEMENTS(counts)-1 DO POLYFILL, x[j]+x1, y[j]+y1, COLOR=color[j]

END
;----------------------------------------------------------------------
PRO mdap_voronoi_2d_binning, x, y, signal, noise, targetSN, $
    class, xNode, yNode, xBar, yBar, sn, area, scale, $
    NO_CVT=no_cvt, PLOT=plot, QUIET=quiet, WVT=wvt,SN_CALIBRATION=SN_CALIBRATION
COMPILE_OPT IDL2
ON_ERROR, 2
;
; This is the main program that has to be called from external programs.
; It simply calls in sequence the different steps of the algorithms
; and optionally plots the results at the end of the calculation.
;
; v0.1 adapted from original voronoi_2d_binning by Capellari & Copin
; 2003 for MANGA data

redo_with_less_sn_min:

npix = n_elements(x)
if n_elements(y) ne npix or n_elements(signal) ne npix $
    or n_elements(noise) ne npix then $
        message, 'Input vectors (x, y, signal, noise) must have the same size'
if n_elements(targetSN) ne 1 then message, 'targetSN must be a scalar'
if not array_equal(noise ge 0, 1) then message, 'NOISE cannot be negative'

; Perform basic tests to catch common input errors

if targetSN le 1 then goto, no_signal

if total(signal)/sqrt(total(noise^2)) lt targetSN then begin ;$
  ;  message, 'Not enough S/N in the whole set of pixels. ' $
  ;      + 'Many pixels may have noise but virtually no signal. ' $
  ;      + 'They should not be included in the set to bin, ' $
  ;      + 'or the pixels should be optimally weighted.' $
  ;      + 'See Cappellari & Copin (2003, Sec.2.1) and README file.'
  targetSN = targetSN*.8
 print, 'targetSN is now: ',targetSN

  goto,redo_with_less_sn_min
endif

;if min(signal/noise) gt targetSN then $
;    message, 'All pixels have enough S/N and binning is not needed'

; Prevent division by zero for pixels with signal=0 and
; noise=sqrt(signal)=0 as can happen with X-ray data
;
noise = noise > min(noise[where(noise gt 0)])*1e-9

if not keyword_set(quiet) then print, 'Bin-accretion...'
bin2d_accretion, x, y, signal, noise, targetSN, class, pixelSize, SN_CALIBRATION=SN_CALIBRATION,QUIET=quiet
if min(class) eq 0 and max(class) eq 0  then begin
   targetSN = targetSN*.8
   goto,redo_with_less_sn_min
endif
if not keyword_set(quiet) then print, strtrim(max(class),2), ' initial bins.'
if not keyword_set(quiet) then print, 'Reassign bad bins...'
bin2d_reassign_bad_bins, x, y, signal, noise, targetSN, class, xnode, ynode
if not keyword_set(quiet) then print, strtrim(n_elements(xnode),2), ' good bins.'
if not keyword_set(no_cvt) then begin
    if not keyword_set(quiet) then print, 'Modified Lloyd algorithm...'
    bin2d_cvt_equal_mass, x, y, signal, noise, xnode, ynode, scale, iter, QUIET=quiet, WVT=wvt, SN_CALIBRATION=SN_CALIBRATION
   if not keyword_set(quiet) then  print, STRTRIM(iter-1,2), ' iterations.'
endif else scale = 1.0
bin2d_compute_useful_bin_quantities, x, y, signal, noise, xnode, ynode, scale, class, xbar, ybar, sn, area, SN_CALIBRATION=SN_CALIBRATION
w1 = where(area eq 1, COMPLEMENT=w2, m)
if not keyword_set(quiet) then print, 'Unbinned pixels: ', m, ' / ', npix, FORMAT='(a,g0,a,g0)'
if not keyword_set(quiet) then print, 'Fractional S/N scatter (%):', stddev(sn[w2]-targetSN)/targetSN*100

if keyword_set(plot) then begin
    !p.multi=[0,1,2]
    rnd = randomu(seed,n_elements(xnode)) ; Randomize bin colors
    bin2d_display_pixels, x, y, rnd[class], pixelSize
    oplot, xnode, ynode, PSYM=1, SYMSIZE=0.5
    rad = sqrt(xbar^2+ybar^2) ; Use centroids, NOT generators
    plot, rad[w2], sn[w2], PSYM=6, SYMSIZE=0.5, $
        XTITLE="R (arcsec)", YTITLE="Bin S/N", /XSTYLE, $
        XRANGE=[min(rad),max(rad)], YRANGE=[0,max(sn)]
    if (m gt 0) then oplot, rad[w1], sn[w1], PSYM=1, SYMSIZE=0.5
    plots, [min(rad),max(rad)], [targetSN,targetSN]
    !p.multi=0
endif

no_signal:
END
;----------------------------------------------------------------------------
