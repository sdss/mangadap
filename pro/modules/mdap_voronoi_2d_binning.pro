;######################################################################
;
; Copyright (C) 2001-2014, Michele Cappellari
; E-mail: cappellari_at_astro.ox.ac.uk
;
; Updated versions of the software are available from my web page
; http://purl.org/cappellari/software
;
; If you have found this software useful for your
; research, we would appreciate an acknowledgment to use of
; `the Voronoi binning method by Cappellari & Copin (2003)'.
;
; This software is provided as is without any warranty whatsoever.
; Permission to use, for non-commercial purposes is granted.
; Permission to modify for personal or internal use is granted,
; provided this copyright and disclaimer are included unchanged
; at the beginning of the file. All other rights are reserved.
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
;	http://adsabs.harvard.edu/abs/2003MNRAS.342..345C
;
; CALLING SEQUENCE:
;       VORONOI_2D_BINNING, X, Y, Signal, Noise, targetSN, $
;               BinNumber, xBin, yBin, xBar, yBar, SN, NPixels, scale, $
;               /NO_CVT, /PLOT, /QUIET, /WVT, PIXSIZE=pixSize
;
;       The function BIN2D_SN_FUNC below returns the S/N of a bin and it can be
;       changed by the user if needed.
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
;      PIXSIZE: Optional pixel scale of the input data. 
;               This can be the size of a pixel of an image or the size 
;               of a spaxel or lenslet in an integral-field spectrograph.
;             - The value is computed automatically by the program, but 
;               this can take a long times when (X,Y) have many elements. 
;               In those cases the PIXSIZE keyword should be given.
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
;    BINNUMBER: Vector (same size as X) containing the bin number assigned
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
;           SN_FUNCTION       -- Example function to calculate S/N of a bin
;           WEIGHTED_CENTROID -- computes weighted centroid of one bin
;           BIN_ROUNDNESS     -- equation (5) of Cappellari & Copin (2003)
;           BIN_ACCRETION     -- steps (i)-(v) in section 5.1
;           REASSIGN_BAD_BINS -- steps (vi)-(vii) in section 5.1
;           CVT_EQUAL_MASS    -- the modified Lloyd algorithm in section 4.1
;           COMPUTE_USEFUL_BIN_QUANTITIES -- self explanatory
;           DISPLAY_PIXELS    -- plotting of colored pixels
;
; NOTE:
;       This program uses some features of IDL 6.1.
;       Please ask if you need to run it with IDL < 6.1.
;
; MODIFICATION HISTORY:
;       V1.0.0: First implementation. Michele Cappellari, Leiden, June 2001
;       V2.0.0: Major revisions. Stable version. MC, Leiden, 11 September 2001
;       V2.1.0: First released version. Written documentation.
;           MC, Vicenza, 13 February 2003
;       V2.2.0: Added computation of useful bin quantities in output. Deleted some
;           safety checks for zero size bins in CVT. Minor polishing of the code.
;           MC, Leiden, 11 March 2003
;       V2.3.0: Unified the three tests to stop the accretion of one bin.
;           This can improve some bins at the border. MC, Leiden, 9 April 2003
;       V2.3.1: Do *not* assume the first bin is made of one single pixel.
;           Added computation of S/N scatter and plotting of 1-pixel bins.
;           MC, Leiden, 13 April 2003
;       V2.4.0: Addedd basic error checking of input S/N. Reintroduced the
;           treatment for zero-size bins in CVT, which was deleted in V2.2.0.
;           Thanks to Robert Sharp and Kambiz Fathi for reporting problems.
;           MC, Leiden, 10 December 2003.
;       V2.4.1: Added /QUIET keyword and verbose output during the computation.
;           After suggestion by Richard McDermid. MC, Leiden, 14 December 2003
;       V2.4.2: Use LONARR instead of INTARR to define the CLASS vector,
;           to be able to deal with big images. Thanks to Tom Statler.
;           MC, Leiden, 4 August 2004
;       V2.4.3: Corrected bug introduced in V2.3.1. It went undetected for
;           a long time because it could only happen in special conditions.
;           Now we recompute the index of the good bins after computing all
;           centroids of the reassigned bins in reassign_bad_bins. Many thanks
;           to Simona Ghizzardi for her clear analysis of the problem and
;           the solution. MC, Leiden, 29 November 2004
;       V2.4.4: Prevent division by zero for pixels with signal=0
;           and noise=sqrt(signal)=0, as can happen from X-ray data.
;           MC, Leiden, 30 November 2004
;       V2.4.5: Added BIN2D prefix to internal routines to avoid possible
;           naming conflicts. MC, Leiden, 3 December 2004
;       V2.4.6: Added /NO_CVT keyword to optionally skip the CVT step of
;           the algorithm. MC, Leiden, 27 August 2005
;       V2.4.7: Verify that SIGNAL and NOISE are non negative vectors.
;           MC, Leiden, 27 September 2005
;       V2.4.8: Use geometric centroid of a bin during the bin-accretion stage,
;           to allow the routine to deal with negative signal (e.g. in
;           background-subtracted X-ray images). Thanks to Steven Diehl for
;           pointing out the usefulness of dealing with negative signal.
;           MC, Leiden, 23 December 2005
;       V2.5.0: Added two new lines of code and the corresponding /WVT keyword
;           to implement the nice modification to the algorithm proposed by
;           Diehl & Statler (2006). MC, Leiden, 9 March 2006
;       V2.5.1: Updated documentation. MC, Oxford, 3 November 2006
;       V2.5.2: Print number of unbinned pixels. MC, Oxford, 28 March 2007
;       V2.5.3: Fixed program stop, introduced in V2.5.0, with /NO_CVT keyword.
;           MC, Oxford, 3 December 2007
;       V2.5.4: Improved color shuffling for final plot. 
;           MC, Oxford, 30 November 2009
;       V2.5.5: Added PIXSIZE keyword. MC, Oxford, 28 April 2010
;       V2.5.6: Use IDL intrinsic function DISTANCE_MEASURE for 
;           automatic pixelSize, when PIXSIZE keyword is not given. 
;           MC, Oxford, 11 November 2011
;       V2.5.7: Included safety termination criterion of Lloyd algorithm
;           to prevent loops using /WVT. MC, Oxford, 24 March 2012
;       V2.5.8: Update Voronoi tessellation at the exit of bin2d_cvt_equal_mass.
;           This is only done when using /WVT, as DIFF may not be zero at the 
;           last iteration. MC, La Palma, 15 May 2012
;       V2.6.0: Included new SN_FUNCTION to illustrate the fact that the user can
;           define his own function to estimate the S/N of a bin if needed.
;           MC, London, 19 March 2014 
;-----
;       V3.0.0: Edited for MaNGA K. Westfall, 16 Mar 2015
;-
;----------------------------------------------------------------------------
;function bin2d_sn_function, signal, noise, index
;COMPILE_OPT IDL2, HIDDEN

; Generic function to calculate the S/N of a bin with spaxels "index".
; The Voronoi binning algorithm does not require this function to have a
; specific form and this generic one can be changed by the user if needed.
;
; The S/N returned by this function does not need to be an analytic function 
; of S and N. There is no need for this function to return the actual S/N. 
; Instead this function could return any quantity the user needs to optimize. 
;
; For example bin2d_sn_function could be a procedure which uses PPXF to 
; measure the velocity dispersion from the coadded spectrum of spaxels 
; "index" and returns the relative error in the dispersion. 
; Of course an analytic approximation of S/N speeds up the calculation.
;
;return, total(signal[index])/sqrt(total(noise[index]^2))
;end
;----------------------------------------------------------------------
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
PRO bin2d_accretion, x, y, signal, noise, targetSN, class, pixelSize, noise_calib=noise_calib, $
                     optimal_weighting=optimal_weighting, QUIET=quiet
COMPILE_OPT IDL2, HIDDEN

; Implements steps (i)-(v) in section 5.1 of Cappellari & Copin (2003)

n = n_elements(x)
class = lonarr(n) ; will contain the bin number of each given pixel
good = bytarr(n)  ; will contain 1 if the bin has been accepted as good

; For each point, find the distance to all other points and select the minimum.
; This is a robust but slow way of determining the pixel size of unbinned data.
;
if n_elements(pixelSize) eq 0 then $
    pixelSize = min(distance_measure( transpose([[x],[y]]) ))

SN = max(signal/noise, currentBin) ; Start from the pixel with highest S/N

; Rough estimate of the expected final bin number.
; This value is only used to give an idea of the expected
; remaining computation time when binning very big dataset.
;
w = where(signal/noise lt targetSN, COMPLEMENT=w1, NCOMPLEMENT=nc)
sn_total_eval = MDAP_CALCULATE_BIN_SN(signal[w], noise[w], noise_calib=noise_calib, $
                                      optimal_weighting=optimal_weighting)
; maxnum = round(total((signal[w]/noise[w])^2)/targetSN^2) + nc
maxnum = round(sn_total_eval^2/targetSN^2) + nc
;print, 'maxnum: ', maxnum
;stop

; The first bin will be assigned CLASS = 1
; With N pixels there will be at most N bins
;
for ind=1,n do begin

;    if not keyword_set(quiet) then print, ind, maxnum, FORMAT='(%"Bin:  %d / %d")'

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
;       SN = bin2d_sn_function(signal, noise, nextBin)
        SN = MDAP_CALCULATE_BIN_SN(signal[nextBin], noise[nextBin], noise_calib=noise_calib, $
                                   optimal_weighting=optimal_weighting)

        ; Test whether the CANDIDATE pixel is connected to the
        ; current bin, whether the POSSIBLE new bin is round enough
        ; and whether the resulting S/N would get closer to targetSN
        ;
        if ((sqrt(minDist) gt 1.2*pixelsize) or (roundness gt 0.3) $
             or (abs(SN-targetSN) gt abs(SNOld-targetSN))) then begin
                if (SNOld gt 0.9*targetSN) then good[currentBin] = 1
;               if (SNOld gt 0.8*targetSN) then good[currentBin] = 1
                break
             endif 

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

endfor

class = class * good ; Set to zero all bins that did not reach the target S/N

END
;----------------------------------------------------------------------------
pro bin2d_reassign_bad_bins, x, y, class, xnode, ynode
COMPILE_OPT IDL2, HIDDEN

; Implements steps (vi)-(vii) in section 5.1 of Cappellari & Copin (2003)

; Find the centroid of all succesful bins.
; CLASS = 0 are unbinned pixels which are excluded.
;
area = histogram(class, REVERSE_INDICES=r, MIN=1)
good = where(area gt 0, nnodes) ; Obtain the index of the good bins
xnode = dblarr(nnodes)
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
                          noise_calib=noise_calib, optimal_weighting=optimal_weighting, $
                          QUIET=quiet, WVT=wvt
COMPILE_OPT IDL2, HIDDEN

; Implements the modified Lloyd algorithm
; in section 4.1 of Cappellari & Copin (2003).
;
; NB: When the keyword /WVT is set this routine includes
; the modification proposed by Diehl & Statler (2006).

npixels = N_ELEMENTS(signal)
class = lonarr(npixels)  ; See beginning of section 4.1 of CC03
if keyword_set(wvt) then dens = signal*0+1.0 else dens = (signal/noise)^2
scale = xnode*0+1.0 ; Start with the same scale length for all bins
sn = xnode

FOR iter=1,n_elements(xnode) DO BEGIN

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
;       sn[k] = bin2d_sn_function(signal, noise, index)
        sn[k] = MDAP_CALCULATE_BIN_SN(signal[index], noise[index], noise_calib=noise_calib, $
                                      optimal_weighting=optimal_weighting)
    ENDFOR

    if keyword_set(wvt) then scale = sqrt(area/sn) ; Eq. (4) of Diehl & Statler (2006)
    diff = TOTAL((xnode-xnodeOld)^2+(ynode-ynodeOld)^2)

;    if not keyword_set(quiet) then $
;        print, iter, diff, FORMAT='(%"Iter:  %d,  Diff:  %f")'
        
    if diff eq 0 then break

ENDFOR

; If coordinates have changed, re-compute (Weighted) Voronoi Tessellation of the pixels grid
;
if diff gt 0 then begin
    FOR j=0,npixels-1 DO BEGIN
        tmp = MIN(((x[j]-xnode)/scale)^2+((y[j]-ynode)/scale)^2, index)
        class[j] = index
    ENDFOR
    area = HISTOGRAM(class, REVERSE_INDICES=r)
    w = where(area gt 0, m) ; Check for zero-size Voronoi bins
endif

; Only return the generators and scales of the nonzero Voronoi bins
;
xnode = xnode[w]
ynode = ynode[w]
scale = scale[w]

END
;-----------------------------------------------------------------------
pro bin2d_compute_useful_bin_quantities, x, y, signal, noise, xnode, ynode, scale, class, xbar, $
                                         ybar, sn, area, noise_calib=noise_calib, $
                                         optimal_weighting=optimal_weighting
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
xbar = dblarr(nbins)
ybar = xbar
sn = xbar
FOR j=0,nbins-1 DO BEGIN
    index = r[r[j]:r[j+1]-1] ; Find subscripts of pixels in bin j.
    bin2d_weighted_centroid, x[index], y[index], signal[index], xb, yb
    xbar[j] = xb
    ybar[j] = yb
;   sn[j] = bin2d_sn_function(signal, noise, index)
    sn[j] = MDAP_CALCULATE_BIN_SN(signal[index], noise[index], noise_calib=noise_calib, $
                                  optimal_weighting=optimal_weighting)
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
color = bytscl(counts)
FOR j=0, N_ELEMENTS(counts)-1 DO POLYFILL, x[j]+x1, y[j]+y1, COLOR=color[j]

END
;----------------------------------------------------------------------
; Added by KBW (3 Feb 2015); combines all pixels into a single bin
PRO BIN2D_ONE_BIN, $
                indx, x, y, signal, noise, xnode, ynode, scale, class, xbar, ybar, sn, area, $
                noise_calib=noise_calib, optimal_weighting=optimal_weighting
COMPILE_OPT IDL2, HIDDEN

        nbins = 1                                   ; Number of bins
        scale = 1.0                                 ; TODO: Is this right?

        xnode = make_array(nbins, /double, value=mean(x[indx]))   ; Set to mean X position
        ynode = make_array(nbins, /double, value=mean(y[indx]))   ; Set to mean Y position

        bin2d_weighted_centroid, x[indx], y[indx], signal[indx], xb, yb
        xbar = make_array(nbins, /double, value=xb) ; Luminosity-weighted center
        ybar = make_array(nbins, /double, value=yb)

        n = n_elements(x)                           ; Number of pixels
        class = make_array(n, /long, value=-1)      ; Initialize to not be in any bin
        class[indx] = 0                             ; All good cells set to bin=0
        area = make_array(nbins, /integer, value=n_elements(indx)) ; number of pixels in the bin
        sn = dblarr(nbins)                          ; Signal-to-noise
        sn[0] = MDAP_CALCULATE_BIN_SN(signal[indx], noise[indx], noise_calib=noise_calib, $
                                      optimal_weighting=optimal_weighting)
END
;----------------------------------------------------------------------
; Added by KBW (5 Mar 2015); produces bins with only individual spaxels
PRO BIN2D_NO_BINS, $
                x, y, signal, noise, xnode, ynode, scale, class, xbar, ybar, sn, area
COMPILE_OPT IDL2, HIDDEN

        nbins=n_elements(x)                         ; number of bins
        scale = 1.0                                 ; TODO: Is this right?

        xnode = x                                   ; Set to input X and Y positions
        ynode = y

        xbar = x                                    ; Set to input X and Y positions
        ybar = y

        class = indgen(nbins)                       ; All cells set to input bin index
        area = make_array(nbins, /integer, value=1) ; array with number of pixels in the bin
        sn = signal/noise                           ; Signal-to-noise
END
;----------------------------------------------------------------------
PRO mdap_voronoi_2d_binning, x, y, signal, noise, targetSN, $
    class, xNode, yNode, xBar, yBar, sn, area, scale, $
    noise_calib=noise_calib, optimal_weighting=optimal_weighting, $
    NO_CVT=no_cvt, PIXSIZE=pixelSize, PLOT=plot, QUIET=quiet, WVT=wvt
COMPILE_OPT IDL2
ON_ERROR, 2
;
; This is the main program that has to be called from external programs.
; It simply calls in sequence the different steps of the algorithms
; and optionally plots the results at the end of the calculation.

npix = n_elements(x)
if n_elements(y) ne npix or n_elements(signal) ne npix $
    or n_elements(noise) ne npix then $
        message, 'Input vectors (x, y, signal, noise) must have the same size'
if n_elements(targetSN) ne 1 then message, 'targetSN must be a scalar'
if not array_equal(noise ge 0, 1) then message, 'NOISE cannot be negative'

; Perform basic tests to catch common input errors
;
ston = signal/noise

indx = where(ston gt 0. and finite(ston), count)
if count eq 0 then $
    message, 'No bins with finite S/N and S/N > 0!'

sn_total_eval = MDAP_CALCULATE_BIN_SN(signal[indx], noise[indx], noise_calib=noise_calib, $
                                      optimal_weighting=optimal_weighting)
print, 'S/N of all data: ', sn_total_eval
if sn_total_eval lt targetSN then begin
;   message, 'Not enough S/N in the whole set of pixels. ' $
;       + 'Many pixels may have noise but virtually no signal. ' $
;       + 'They should not be included in the set to bin, ' $
;       + 'or the pixels should be optimally weighted.' $
;       + 'See Cappellari & Copin (2003, Sec.2.1) and README file.'

    print, 'WARNING: Combination of all spectra does not reach the requested S/N level.' $
           + '  Proceeding by combining all spectra.'

    BIN2D_ONE_BIN, indx, x, y, signal, noise, xnode, ynode, scale, class, xbar, ybar, sn, $
                   area, noise_calib=noise_calib, optimal_weighting=optimal_weighting
    return
endif

if min(ston) gt targetSN then begin
;   message, 'All pixels have enough S/N and binning is not needed'
    print, 'WARNING: All pixels have enough S/N, nothing to bin!'
    BIN2D_NO_BINS, x, y, signal, noise, xnode, ynode, scale, class, xbar, ybar, sn, area
    return
endif
    

; Prevent division by zero for pixels with signal=0 and
; noise=sqrt(signal)=0 as can happen with X-ray data
;
;noise = noise > min(noise[where(noise gt 0)])*1e-9
indx = where(noise gt 0, count)
if count ne 0 then begin
    noise = noise > min(noise[indx])*1e-9
endif else $
    noise = make_array(n_elements(noise), /double, value=1.)

if not keyword_set(quiet) then $
    print, 'Bin-accretion...'
bin2d_accretion, x, y, signal, noise, targetSN, class, pixelSize, noise_calib=noise_calib, $
                 optimal_weighting=optimal_weighting, QUIET=quiet
;print, 'class: ', class
indx = where(ston gt 0. and finite(ston), count)
;print, count, n_elements(class) 
if not keyword_set(quiet) then begin
    print, strtrim(max(class)+1,2), ' initial bins.'
endif

; All in a single bin
if max(class) eq 0 then begin
    indx = where(ston gt 0. and finite(ston), count)
    print, 'WARNING: All spectra assigned to a single bin!  Proceeding by combining all spectra.'
    BIN2D_ONE_BIN, indx, x, y, signal, noise, xnode, ynode, scale, class, xbar, ybar, sn, $
                   area, noise_calib=noise_calib, optimal_weighting=optimal_weighting
    return
endif

if not keyword_set(quiet) then begin
    print, 'Reassign bad bins...'
endif
bin2d_reassign_bad_bins, x, y, class, xnode, ynode
if not keyword_set(quiet) then $
    print, strtrim(n_elements(xnode),2), ' good bins.'
if not keyword_set(no_cvt) then begin
    if not keyword_set(quiet) then $
        print, 'Modified Lloyd algorithm...'
    bin2d_cvt_equal_mass, x, y, signal, noise, xnode, ynode, scale, iter, noise_calib=noise_calib, $
                          optimal_weighting=optimal_weighting, QUIET=quiet, WVT=wvt
    if not keyword_set(quiet) then $
        print, STRTRIM(iter-1,2), ' iterations.'
endif else scale = 1.0
bin2d_compute_useful_bin_quantities, x, y, signal, noise, xnode, ynode, scale, class, xbar, ybar, $
                                     sn, area, noise_calib=noise_calib, $
                                     optimal_weighting=optimal_weighting
if not keyword_set(quiet) then begin
    w1 = where(area eq 1, COMPLEMENT=w2, m)
    print, 'Unbinned pixels: ', m, ' / ', npix, FORMAT='(a,g0,a,g0)'
    print, 'Fractional S/N scatter (%):', stddev(sn[w2]-targetSN)/targetSN*100
endif

if keyword_set(plot) then begin
    !p.multi=[0,1,2]
    rnd = sort(randomu(seed,n_elements(xnode))) ; Randomize bin colors
    bin2d_display_pixels, x, y, rnd[class], pixelSize
    plots, xnode, ynode, PSYM=1, SYMSIZE=0.5
    rad = sqrt(xbar^2+ybar^2) ; Use centroids, NOT generators
    plot, rad[w2], sn[w2], PSYM=6, SYMSIZE=0.5, $
        XTITLE="R (arcsec)", YTITLE="Bin S/N", /XSTYLE, $
        XRANGE=[min(rad),max(rad)], YRANGE=[0,max(sn)]
    if (m gt 0) then oplot, rad[w1], sn[w1], PSYM=1, SYMSIZE=0.5
    plots, [min(rad),max(rad)], [targetSN,targetSN]
    !p.multi=0
endif

END
;----------------------------------------------------------------------------
