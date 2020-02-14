r"""
Eric EMSELLEM - ESO / CRAL
eric.emsellem@eso.org

The main goal of this script is to derive :math:`\lambda_R` and
:math:`V/\sigma` radial profiles, given a set of coordinates
:math:`(x,y)` and kinematic measurements (:math:`V`, :math:`\sigma`, but
also associated flux, :math:`F`)

The standard formulae for :math:`\lambda_R` comes from Emsellem et al.
2007, MNRAS 379, 401 and Emsellem et al. 2011, MNRAS 414, 888, while the
formulae for :math:`V/\sigma` done in the proper way with IFU come from
Binney 2005, 363, 937.

.. math::

    \lambda_R = \frac{\sum F\ |V|}{\sum F\ (V^2 + \sigma^2)^{1/2}}

and

.. math::

    V / \sigma = \left( \frac{\sum F\ V^2}{\sum F\ \sigma^2}
    \right)^{1/2}

The main difference between :math:`\lambda_R` and :math:`V/\sigma` is
that :math:`\lambda_R` is providing a radius-biased view of the
dynamical status and can be used as a proxy for apparent specific
stellar angular momentum, while :math:`V/\sigma` is flux-biased and more
prone to changes due to small central structures.  There is, however, a
nice relation between the two (see the above- mentioned papers) which
holds for regular kinematics.

Both quantities are derived within a given aperture (selected set of
spaxels) and are therefore "cumulative" in the sense that they should
include all spaxels within a certain "radius".  The advised selection is
to follow the moment ellipticity of the system: the derivation of
:math:`\lambda_R` and :math:`V/\sigma` thus require values for the
ellipticity (Eps; :math:`\epsilon`) and position angle (PA;
:math:`\phi_0`).

*Source location*:
    $MANGADAP_DIR/python/mangadap/contrib/LambdaR_2D_forMaNGA.py

*Python2/3 compliance*::

    from __future__ import division
    from __future__ import print_function
    from __future__ import absolute_import
    
    import sys
    if sys.version > '3':
        long = int
        xrange = range

*Imports*::

    import numpy as np
    from numpy import cos, sin, sum, sqrt
    import scipy
    from scipy import interpolate
    import copy

*Revision history*:

    | version 1.1.1 - June 15, 2015: K. Westfall - inclusion in MaNGA
        DAP, minor modifications to docstrings, python 2/3 compliance.
    | version 1.1.0 - March 30, 2015: transformed for MaNGA usage
    | version 1.0.3 - August 20, 2014: added SigmaR
    | version 1.0.2 - August 19, 2014
    | version 1.0.1 - April 23, 2014
    | version 1.0.0 - August 12, 2011 : creation after 2006 module

.. _scipy.interpolate.interp1d: http://docs.scipy.org/doc/scipy/reference/generated/scipy.interpolate.interp1d.html#scipy.interpolate.interp1d

"""

from __future__ import division
from __future__ import print_function
from __future__ import absolute_import

import sys
if sys.version > '3':
    long = int
    xrange = range

import numpy as np
from numpy import cos, sin, sum, sqrt
import scipy
from scipy import interpolate
import copy

def dist2(x1,y1,x2,y2, scale=1.0) :
    return ((x1-x2)**2 + (y1-y2)**2)/scale**2

def guess_regular_grid(xnodes, ynodes, pixelsize=None) :
    """
    Return a regular grid guessed on an irregular one (Voronoi).

    Args:
        xnodes (numpy.array): arrays of Voronoi bin x positions
        ynodes (numpy.array): arrays of Voronoi bin y positions

    Return:
        numpy.ndarray : (xunb, yunb) regular grid for x and y (unbinned)

    """
    ## First deriving a pixel size
    xn_rav, yn_rav = xnodes.ravel(), ynodes.ravel()
    if pixelsize is None :
        pixelsize = derive_pixelsize(xnodes, ynodes)
    minxn = np.int(np.min(xn_rav) / pixelsize) * pixelsize
    minyn = np.int(np.min(yn_rav) / pixelsize) * pixelsize
    xunb, yunb = np.meshgrid(np.arange(minxn, np.max(xn_rav)+pixelsize, pixelsize),
                           np.arange(minyn, np.max(yn_rav)+pixelsize, pixelsize))

    return xunb, yunb

def derive_unbinned_field(xnodes, ynodes, data, xunb=None, yunb=None, mask_array=None) :
    """
    Provide an array of the same shape as the input xunb, and yunb with
    the values derived from the Voronoi binned data

    Args:
        xnodes (numpy.array): arrays of Voronoi bin x positions
        ynodes (numpy.array): arrays of Voronoi bin y positions
        data (numpy.array): values for each node
        xunb (numpy.array): x coordinates of the unbinned data if not
            provided (default) they will be guessed from the nodes
        yunb (numpy.array): y coordinates of the unbinned data if not
            provided (default) they will be guessed from the nodes
        mask_array (numpy.ndarray): array with the same shape than xunb
            providing mask values for positions to ignore

    Returns:
        numpy.ndarray : (xunb, yunb, unbinned_data) arrays with the same
            shape as xunb,

    """
    if xunb is None :
        xunb, yunb = guess_regular_grid(xnodes, ynodes)

    x_rav, y_rav = xunb.ravel(), yunb.ravel()
    xnodes_rav, ynodes_rav = xnodes.ravel(), ynodes.ravel()
    data_rav = data.ravel()
    unbinned_data = np.zeros_like(x_rav, dtype=data_rav.dtype)
    for i in xrange(len(x_rav)) :
        indclosestBin = np.argmin(dist2(x_rav[i], y_rav[i], xnodes_rav, ynodes_rav))
        unbinned_data[i] = data_rav[indclosestBin]

    return xunb, yunb, unbinned_data.reshape(xunb.shape)

###################################################################
# Derive V / Sigma and LambdaR profiles from a set 
# of X, Y, I, V and S
# =================================================================
def Derive_LR_VS_Profiles(X, Y, F, V, S, Rprofile=None, R_EpsPA=None, Eps=None, PA=None,
                          maskDic=None, maskName=None, verbose=False, **kwargs) :
    r"""
    Compute the :math:`\lambda_R` and :math:`V/\sigma` parameters from a
    set of :math:`(x,y)` (positions) and :math:`V`, :math:`\sigma`
    (stellar kinematics).

    If an ellipticity profile is provided, it is used to select the
    points within growing apertures.  Otherwise it uses a default
    ellipse value, constant over the full field.

    Same with the position angle (in degrees), which can vary with
    radius or be constant.

    Rprofile and R_EpsPA are the effective (Geometric) radius.

    .. todo::

        - Documentation out of date.  E.g., not all kwargs listed.
        - Convert result to a true class
        - raise exceptions instead of writing errors and returning
    
    Args:
        X, Y (numpy.ndarray) : Cartesian coordinates. Units can be
            arbitrary as :math:`\lambda_R` and :math:`V/\sigma` are
            dimensionless (but should be linear, so not RA, DEC) *X*,
            and *Y* should be associated with square spaxels but do not
            need to fill in a full rectangular grid.  No default. 
        F, V, S (numpy.ndarray): flux, stellar velocity and velocity
            dispersion measured at *X*, *Y* positions.  Units are
            arbitrary as for *X*, *Y* (since the final result will be
            dimensionless) No default.
        Rprofile (numpy.ndarray): (Optional) Radius sampling imposed by
            the user.  These are the radii where the profiles will be
            derived. Default to None, means that Rsampling will be
            automatically set up.
        R_EpsPA, Eps, PA (numpy.ndarray): ellipticity and Position Angle
            profiles, given as radii (R_EpsPA), ellipticity values (Eps)
            and position angle (degrees). If given (default to None),
            this profile will be used to derive the apertures.
        maskDic (?): (Optional) ?
        maskName (?): (Optional) ?
        verbose (bool): (Optional) Print some additional information if
            relevant.
        Symmetrize (bool): (via kwargs: default is True) to specify
            whether or not the fields should be symmetrized.
        Re (float): (via kwargs: default is 1) effective radius (in
            units of input X, Y)
        Estimate_Vsys (bool): (via kwargs: default is True) If True,
            estimate the systemic velocity (Systemic_Velocity is then
            overwritten). If False, value of Systemic_Velocity will be
            used.
        Systemic_Velocity (float): (via kwargs: default is 0.) Systemic
            Velocity to subtract from V.
        Vaperture (float): (via kwargs: default is 3.) Radius (in units
            of input X, Y) within which V will be summed to estimate the
            systemic velocity.
        Maximum_Velocity (float): (via kwargs: default is 400.)
            Threshold to select out stellar velocities.
        Maximum_Dispersion (float): (via kwargs: default is 500.)
            Threshold to select out stellar velocity dispersion.
        Threshold_Cover (float): (via kwargs: default is 1.25) Fraction
            of pixels which should be covered within an ellipse. Default
            to 1.25, meaning that if more than values for more than 20%
            of the pixels within the ellipse are missing, we stop the
            calculation (as it means that there are not enough points).
        Threshold_Radius (float): (via kwargs: default is 1.15) Same as
            Threshold_Cover but this time in terms of the ratio between
            the actual effective radius, and the effective radius of the
            corresponding ellipse. 
        Min_Radius (float): (via kwargs: default is 1) Minimum radius to
            be considered in the profile.
        Max_Radius (float): (via kwargs: default is 50) Maximum radius
            to be considered in the profile.
        N_Radius (int): (via kwargs: default is 100) Number of points
            within the radius range.

    Returns:
        Result : A class that contains the results of the computation.
            May return None if errors occur.

    """

    ## Defining the Result Class ================================
    class Result :
       def __init__(self):
          self.type = "LambdaR Computation"
    ## ===========================================================

    ## Here is the result structure which will be filled in at the end
    result = Result()
    ## ===========================================================

    ## Default setup =============================================
    result.Symmetrize = kwargs.get('Symmetrize', True)
    result.ExtraCoverageCriterion = kwargs.get('ExtraCoverageCriterion', False)

    result.Re = np.float32(kwargs.get('Re', 1.))

    result.Estimate_Vsys = kwargs.get('Estimate_Vsys', True)
    result.Systemic_Velocity = np.int(kwargs.get('Systemic_Velocity', 0.))
    result.Vaperture = np.float32(kwargs.get('Vaperture', 3.))

    result.Maximum_Velocity = np.float32(kwargs.get('Maximum_Velocity', 500.))
    result.Maximum_Dispersion = np.float32(kwargs.get('Maximum_Dispersion', 500.))

    result.Threshold_Cover = np.float32(kwargs.get('Threshold_Cover', 1.25))
    result.Threshold_Radius = np.float32(kwargs.get('Threshold_Radius', 1.15))

    result.Min_Radius = np.float32(kwargs.get('Min_Radius', 1.))
    result.Max_Radius = np.float32(kwargs.get('Max_Radius', 50.))
    result.N_Radius = np.int(kwargs.get('N_Radius', 100))
    ## ===========================================================

    ## Checking array ============================================
    ref_ShapeX = X.shape
    Shapes = Y.shape, F.shape, V.shape, S.shape
    if not all(s == ref_ShapeX for s in Shapes) :
        print('ERROR: input arrays (X, Y, F, V, S) are not all of the same size')
        print('ERROR: please check the shapes/sizes of these data.')
        return

    if R_EpsPA != None :
        ref_ShapeEps = len(R_EpsPA)
        Shapes = len(Eps), len(PA)
        if not all(s == ref_ShapeEps for s in Shapes) :
            print('ERROR: input arrays (REps, Eps, PA) are not all of the same size')
            print('ERROR: please check the shapes/sizes of these data.')
            return
    ## ===========================================================

    ## Select the right spaxels
    ## All pixels with positive flux, all pixels with positive velocity and dispersion, 
    ## all pixels with velocities lower (amplitude) than the set up Maximum velocity
    ## all pixels with dispersion lower than the set up Maximum dispersion

    ## selFlux is the frst selection to just get an idea of the systemic V
    First_selFlux = (F > 0.) & (S >= 0.) & (np.abs(S) < result.Maximum_Dispersion)
    ## Finally getting the min and max of the velocity and the systemic value from an aperture median
    AperV, minV, maxV = find_centerodd(X[First_selFlux],Y[First_selFlux],V[First_selFlux], Radius=result.Vaperture)

    ## selFlux is the final selection for the spaxels
    selFlux = (F > 0.) & (S >= 0.) & (np.abs(V - AperV) < result.Maximum_Velocity) & (np.abs(S) < result.Maximum_Dispersion)
    ## We combine this with the Masks which may be defined for that galaxy
    mask = select_spaxels(maskDic, maskName, X, Y) * selFlux

    ## Now reducing the arrays to the selected ones
    ## These are the arrays we will use for the calculation
    sX = np.compress(mask, X, axis=0)
    sY = np.compress(mask, Y, axis=0)
    sF = np.compress(mask, F, axis=0)
    sV = np.compress(mask, V, axis=0)
    sS = np.compress(mask, S, axis=0)

    ## Getting the central values and normalising V with Vsys
    ## Finally getting the min and max of the velocity and the systemic value from an aperture median
    AperV, minV, maxV = find_centerodd(sX, sY, sV, Radius=result.Vaperture)
    ## Same for the dispersion (just min and max values)
    minS, maxS = find_centereven(sX,sY,sS)

    if result.Estimate_Vsys : 
        ## Use the estimate aperture value
        result.Systemic_Velocity = AperV
        if verbose:
            print('Central velocity measured is : {0}'.format(AperV))
    else : 
        if verbose:
            print('Systemic Velocity value use is : {0}'.format(result.Systemic_Velocity))

    # Return the symmetrised points if Symmetrize is true
    if result.Symmetrize :
        if verbose : 
            print('Symmetrizing the points')
        nsV = SymmetrizeField(sX, sY, sV - result.Systemic_Velocity, odd=1)
        nsS = SymmetrizeField(sX, sY, sS, odd=0)
    else :
        nsV = sV - result.Systemic_Velocity
        nsS = sS

    ## Define the Radii to be used for the profile
    ## If not provided by Rprofile, then build one using Min/Max_Radius and N_Radius
    if Rprofile == None :
       Rprofile = np.linspace(result.Min_Radius, result.Max_Radius, result.N_Radius)

    ## Creation of a big grid of step
    dist = sqrt((X - X[0])**2 + (Y - Y[0])**2)
    ## looking for the minimum radius (which is not zero) which should be the Size of the smallest Spaxel
    binstep = np.min(dist[dist>0])
    if verbose: 
        print('BINSTEP is {0}'.format(binstep))

    ## Computing the Map Grid Corners (after masking)
    X1, X2, Y1, Y2 = Get_CornersArray(sX, sY, binstep)

    ## And computing the large rectangular grid containing the full set of spaxels
    Nspaxels = np.max([X.min()/binstep, X.max()/binstep, Y.min()/binstep, Y.max()/binstep])
    XYsample = np.arange(-Nspaxels*binstep, (Nspaxels+1)*binstep, binstep)
    GX, GY = np.meshgrid(XYsample, XYsample)
    ## Now computing the Full grid corners
    GX1, GX2, GY1, GY2 = Get_CornersArray(GX, GY, binstep)

    ## Define the ellipticity profiles to be used ================
    if R_EpsPA == None :
        if (Eps == None)  :
            print('Error: no Ellipticity provided')
            return 
        elif (PA == None) :
            print('Error: no Position Angle provided')
            return 
        else :
            R_EpsPA = Rprofile
            Eps = np.zeros_like(Rprofile) + Eps
            PA = np.zeros_like(Rprofile) + PA

    ## Defining the function to interpolate Eps and PA
    profEps = interpolate.interp1d(R_EpsPA, Eps, kind="linear")
    profPA = interpolate.interp1d(R_EpsPA, PA, kind="linear")
    ## ===========================================================

    ## Contracting Rprofile with the input R_EpsPA ==================
    compactRprofile = Rprofile[(Rprofile > np.min(R_EpsPA)) & (Rprofile < np.max(R_EpsPA))]

    ## Start the Loop ============================================
    Npoints = len(compactRprofile)

    LambdaR = np.zeros(Npoints, dtype=np.float32)
    FlagRadius = np.ones(Npoints, dtype=np.int32)

    VS = np.zeros_like(LambdaR)
    sumF = np.zeros_like(LambdaR)
    FV2 = np.zeros_like(LambdaR)
    FS2 = np.zeros_like(LambdaR)
    FMU2 = np.zeros_like(LambdaR)
    RFV = np.zeros_like(LambdaR)
    RFMU = np.zeros_like(LambdaR)
    SemiMajorRadius = np.zeros_like(LambdaR)
    EffRadius = np.zeros_like(LambdaR)
    GEffRadius = np.zeros_like(LambdaR)

    ## Starting the actual loop on the profile radii
    for i in range(Npoints) :
       locRadius = compactRprofile[i]    ## This is the geometric radius
       Eps_Here = profEps(locRadius)     ## THe ellipticity at that radius
       if Eps_Here >= 1. :
          print('Error: we get an ellipticity of 1 at radius {0}'.format(locRadius))
          continue

       PA_Here = profPA(locRadius)       ## The PA at that radius
       SemiMajorRadius[i] = locRadius / sqrt(1. - Eps_Here)  ## Semi major axis radius

       ## Rotating the axis to have it aligned
       rsX, rsY = rotPA(sX,sY, PA_Here)
       rellipse = XY_toSemiMajor(rsX, rsY, Eps_Here) ## SEMI MAJOR AXIS

       ## Rotating the Grid in the same way
       rGX, rGY = rotPA(GX,GY, PA_Here)
       rGellipse = XY_toSemiMajor(rGX.ravel(), rGY.ravel(), Eps_Here) ## SEMI MAJOR AXIS

       ## Selection for the present Radius
       selection = (rellipse <= SemiMajorRadius[i]) 
       Gselection = (rGellipse <= SemiMajorRadius[i]) 

       ## Rotating the Corners
       rXa, rXb, rXc, rXd, rYa, rYb, rYc, rYd = Rotate_Corners(X1, X2, Y1, Y2, PA_Here)
       ## Deriving the corresponding Major-Axis radius
       rElla = XY_toSemiMajor(rXa, rYa, Eps_Here) 
       rEllb = XY_toSemiMajor(rXb, rYb, Eps_Here) 
       rEllc = XY_toSemiMajor(rXc, rYc, Eps_Here) 
       rElld = XY_toSemiMajor(rXd, rYd, Eps_Here) 
       ## Also just the radius
       ra = XY_toSemiMajor(rXa, rYa, 0.0) 
       rb = XY_toSemiMajor(rXb, rYb, 0.0) 
       rc = XY_toSemiMajor(rXc, rYc, 0.0) 
       rd = XY_toSemiMajor(rXd, rYd, 0.0) 

       ## Same with the Full Grid
       rGXa, rGXb, rGXc, rGXd, rGYa, rGYb, rGYc, rGYd = Rotate_Corners(GX1, GX2, GY1, GY2, PA_Here)
       rGElla = XY_toSemiMajor(rGXa, rGYa, Eps_Here) 
       rGEllb = XY_toSemiMajor(rGXb, rGYb, Eps_Here) 
       rGEllc = XY_toSemiMajor(rGXc, rGYc, Eps_Here) 
       rGElld = XY_toSemiMajor(rGXd, rGYd, Eps_Here) 

       Allr = np.concatenate((ra,rb,rc,rd))
       AllEllr = np.concatenate((rElla,rEllb,rEllc,rElld))
       AllGEllr = np.concatenate((rGElla,rGEllb,rGEllc,rGElld))
       MaxR = np.max(Allr)

       ## Selecting the coordinates within this region
       srsX = np.compress(selection, rsX, axis=0)
       srsY = np.compress(selection, rsY, axis=0)
       ssX = np.compress(selection, sX, axis=0) 
       ssY = np.compress(selection, sY, axis=0) 
       ssR = sqrt(ssX**2 + ssY**2)
       ssF = np.compress(selection, sF, axis=0)
       ssV = np.compress(selection, nsV, axis=0)
       ssS = np.compress(selection, nsS, axis=0)

       ## Calculate the full Area on the full Grid
       ssGX = np.compress(Gselection, GX.ravel(), axis=0)
       ## The effective radius is then just the SQRT of the AREA
       GEffRadius[i] = np.sqrt(len(ssGX) * binstep * binstep / np.pi)

       ## Calculate the effective radius for the selection on Input Data
       EffRadius[i] = sqrt(len(ssF) * binstep * binstep / np.pi)

       ## For V/S calculate I * V^2 and I * S^2 as in Binney 2005
       FV2[i] = sum(ssF * ssV**2, axis=0)
       FS2[i] = sum(ssF * ssS**2, axis=0)

       ## For Sigma_e
       FMU2[i] = FV2[i] + FS2[i]
       sumF[i] = sum(ssF)

       ## For LambdaR, derive R * I * ABS(V)
       ## and R * I * (V^2 +  S^2)
       RFV[i] = sum(ssF * np.abs(ssV) * ssR, axis=0)
       RFMU[i] = sum(ssF * sqrt(ssV**2 + ssS**2) * ssR, axis=0)

       #  Lp, Lm, momEps, momPA = Flux_Inertia(ssX, ssY, ssF)

       ## Check if the Effective radius passes the Threshold
       ## If not, Flag the point to 0
       if (GEffRadius[i] / EffRadius[i] > result.Threshold_Radius) :
          FlagRadius[i] = 0.

       ## Now check other criteria by looking at how many pixels are intersecting
       ## This criterion should in principle be ignored
       if result.ExtraCoverageCriterion :
          nPixEll = len(sX[check_cross_pixel(AllEllr,SemiMajorRadius[i])])
          nPixGEll = len(GX[check_cross_pixel(AllGEllr,SemiMajorRadius[i])])
          if nPixEll == 0 :
             FlagRadius[i] = 0
          elif (((nGPixEll*1.00)/(nPixEll*1.00)) > result.Threshold_Cover) :
             FlagRadius[i] = 0

    ## Finalising the calculation
    ## V / Sigma
    VS = sqrt(FV2 / FS2)
    ## LambdaR
    LambdaR = RFV / RFMU
    ## Aperture Sigma
    SigmaR = sqrt(FMU2 / sumF)

    ## Interpolation at Various Radii : Re/2, Re 
    selRadius = (FlagRadius == 1)
    EffRadiusN = EffRadius / result.Re
    MaxEffRadius = np.max(EffRadius[selRadius])
    MaxEffRadiusN = MaxEffRadius / result.Re

    ## Calculation at 0.5, 1 and 2 Re
    RadiusN_for_re = np.minimum(1.0,  MaxEffRadiusN)
    RadiusN_for_2re = np.minimum(2.0,  MaxEffRadiusN)
    RadiusN_for_re2 = np.minimum(0.5,  MaxEffRadiusN)

    VSre2 = interp_value_fromR(EffRadiusN, VS, RadiusN_for_re2)
    LambdaRre2 = interp_value_fromR(EffRadiusN, LambdaR, RadiusN_for_re2)
    Sigmare2 = interp_value_fromR(EffRadiusN, SigmaR, RadiusN_for_re2)
    VSre = interp_value_fromR(EffRadiusN, VS, RadiusN_for_re)
    LambdaRre = interp_value_fromR(EffRadiusN, LambdaR, RadiusN_for_re)
    Sigmare = interp_value_fromR(EffRadiusN, SigmaR, RadiusN_for_re)
    VS2re = interp_value_fromR(EffRadiusN, VS, RadiusN_for_2re)
    LambdaR2re = interp_value_fromR(EffRadiusN, LambdaR, RadiusN_for_2re)
    Sigma2re = interp_value_fromR(EffRadiusN, SigmaR, RadiusN_for_2re)

    ## Just adding all the relevant arrays into the structure
    result.LambdaR = LambdaR
    result.VS = VS
    result.GEffRadius = GEffRadius
    result.EffRadius = EffRadius
    result.EffRadiusN = EffRadiusN
    result.FlagRadius = FlagRadius
    result.MaxEffRadius = MaxEffRadius
    result.LambdaR_re = LambdaRre
    result.LambdaR_re2 = LambdaRre2
    result.LambdaR_2re = LambdaR2re
    result.RadiusN_for_re = RadiusN_for_re
    result.RadiusN_for_2re = RadiusN_for_2re
    result.RadiusN_for_re2 = RadiusN_for_re2
    result.Sigma_re = Sigmare
    result.Sigma_re2 = Sigmare2
    result.Sigma_2re = Sigma2re
    result.VS_re = VSre
    result.VS_re2 = VSre2
    result.VS_2re = VS2re
    result.binstep = binstep
    result.R_EpsPA = R_EpsPA
    result.Eps = Eps
    result.PA = PA
    result.SemiMajorRadius = SemiMajorRadius
    return result
# ===============================================================
###################################################################
def interp_value_fromR(R=None, val=None, radius=None) :
    """
    Interpolation of a profile at a certain radius.  Check if the radius
    is really within the given input range

    If radius out of range, return the closest (in radius) value. This
    is thus different from the option of bounds_error of
    `scipy.interpolate.interp1d`_ which is used here.

    Args:
        R (numpy.ndarray): (Optional) Input radii
        val (numpy.ndarray): (Optional) Input values. Should have the
            same size as R.
        radius (float): (Optional) Radius at which to interpolate

    Returns:
        float : Value interpolated from the input data.

    """
    if (radius < R[0]) :
        return val[0]

    elif (radius > R[-1]) :
        return val[-1]

    else :
        fval = interpolate.interp1d(R, val, kind="linear")
        return fval([radius])[0]
# =================================================================
####################################################################
# def Flux_Inertia(X, Y, Flux) :
#     """
#     Derive the moment of inertia from a flux map
#     Given X, Y and Flux values
# 
#     Return minor, major, ellipticity, and PA in degrees
#     """
#     momI = sum(Flux, axis=0)
#     if momI == 0. :
#         return 0., 0., 1., 0.
#     momIX = sum(Flux * X, axis=0) / momI
#     momIY = sum(Flux * Y, axis=0) / momI
#     a = sum(Flux * X * X, axis=0) / momI - momIX**2
#     b = sum(Flux * Y * Y, axis=0) / momI - momIY**2
#     c = sum(Flux * X * Y, axis=0) / momI - momIX * momIY
#     if c == 0 :
#         if a == 0. :
#             return b, a, 0., 0.
#         if b > a :
#             return b, a, 1. - sqrt(a / b), 0.
#         else :
#             if b == 0.:
#                 return a,b,0.,0.
#             else :
#                 return a, b, 1. - sqrt(b/a), 90.
# 
#     delta = (a - b)**2. + 4 * c * c 
#     Lp2 = ((a + b) + sqrt(delta)) / 2.
#     Lm2 = ((a + b) - sqrt(delta)) / 2.
#     Lp = sqrt(np.maximum(Lp2, 0.))
#     Lm = sqrt(np.maximum(Lm2, 0.))
#     eps = (Lp - Lm) / Lp
#     PA = np.rad2deg((np.arctan((b - Lp2) / c))
#     return Lp, Lm, eps, PA
# ==================================================================

####################################################################
# Function to see if the ellipse is basically 
#          in or out or intersects the spaxel
####################################################################
def check_cross_pixel(rpixels, rlimit) :
    """
    Check whether the ellipse is going THROUGH the pixel or not.  Return
    an array of True (intersecting) or False (Not intersecting).
    """
    signC = np.sign(rpixels - rlimit)
    # True means the pixel has the circle intersecting it
    # False means the ellipse is in or out
    return np.abs(signC.sum(axis=0)) < 4.0
# ==================================================================

####################################################################
# Define the central value of a map
# Case of an odd symmetry
def find_centerodd(X, Y, Z, Radius=3.0) :
    """
    Find central value for an odd sided field.

    Args:
        X (numpy.ndarray): x coordinates
        Y (numpy.ndarray): y coordinates
        Z (numpy.ndarray): values
        Radius (float): (Optional) Radius within which to derive the
            central value.

    Returns:
        float: The central value and some amplitude value (range about
            that value?).

    """
    # First select the points which are non-zero and compress
    ravZ = np.ravel(Z)
    sel = (ravZ != 0.)
    selz = np.compress(sel,ravZ)

    ## Select the central points for the central value
    selxy = np.ravel((np.abs(X) < Radius) & (np.abs(Y) < Radius)) & sel
    selzxy = np.compress(selxy, ravZ)
    cval = np.median(selzxy)

    sig = np.std(selz)
    sselz = np.compress(np.abs(selz - cval) < 3 * sig, selz - cval)
    ampl = np.max(np.abs(sselz)) / 1.1
    return cval, cval - ampl, cval + ampl
# ==================================================================

####################################################################
# Define the central value of a map
# Case of an even symmetry
def find_centereven(X, Y, Z, Radius=3.0, sel0=True) :
    """
    Find the central value for an even sided field.

    Args:
        X (numpy.ndarray): x coordinates
        Y (numpy.ndarray): y coordinates
        Z (numpy.ndarray): values
        Radius (float): (Optional) Radius within which to derive the
            central value.
        sel0 (bool): (Optional) ?

    Returns:
        float : The minimum and maximum value.

    """
    ravZ = np.ravel(Z)
    if sel0 : sel = (ravZ != 0.)
    else : sel = (np.abs(ravZ) >= 0.)
    sel = np.ravel((np.abs(X) < Radius) & (np.abs(Y) < Radius)) & sel
    selz = np.compress(sel, ravZ)
    maxz = np.max(selz)
    minz = np.min(selz)

    return minz, maxz
# ==================================================================

########## point symmetrization of the points #####################
# Use a simple symmetry since pixels are centred
# So reading the opposite pixel
###################################################################
def SymmetrizeField(X, Y, Z, odd=1, maskDic=None, maskName=None) :
    """
    Symmetrize a field. Only works when spaxels are located around a
    central spaxel.

    Args:
        X (numpy.ndarray): x coordinates
        Y (numpy.ndarray): y coordinates
        Z (numpy.ndarray): input values
        odd (bool): (Optional) True if it is odd w.r.t. centre
            (X=0,Y=0). False if even.
        maskDic, maskName: (Optional) Dictionary and name of masks

    Returns:
        numpy.ndarray : The symmetrized array (symmetrized Z values)

    """
    newz = copy.copy(Z)

    if odd == 1:
        weight = -1.
    else :
        weight = 1.

    ## Interpolation
    xp = np.zeros(1, np.float32)
    yp = np.zeros(1, np.float32)
    for i in range(len(Z)) :
        ## Look for the symmetric point (w.r.t centre)
        xp[0] = -X[i]
        yp[0] = -Y[i]
        ## Only one spaxel to test
        selmask = select_spaxels(maskDic, maskName, xp, yp) # Mask the bad regions
        ## If the symmetric spaxel is not good then I just use the value within the masked region
        ## Otherwise I average neighbouring spaxels
        if selmask[0] ==  True :
            nz = Z[(abs(X - xp[0]) + abs(Y - yp[0])) < 1.e-3]
            ## If found some spaxels, average the values
            if len(nz) > 0 :
                newz[i] = (nz[0] * weight + Z[i]) / 2.0

    return newz
# ===============================================================
###################################################################
def Get_CornersArray(X, Y, binstep) :
     """
     Return corners given a set of X and Y and a fixed binstep.
     """
     X1, X2 = X + binstep/2., X - binstep/2.
     Y1, Y2 = Y + binstep/2., Y - binstep/2.
     return X1, X2, Y1, Y2
# ===============================================================
####################################################################
## Rotating coordinates, this includes a 90 degrees rot to
## have the major axis along the Y axis transformed in X
## So rotangle is the POSITION ANGLE from Top, counter-clockwise
def rotPA(X,Y, PA) :
     """ Rotation of X and Y
         PA in degrees
     """ 
     PA_rad = np.deg2rad(PA)
     Xr = - X * sin(PA_rad) + Y * cos(PA_rad)
     Yr = - X * cos(PA_rad) - Y * sin(PA_rad)
     return Xr, Yr
# ==================================================================
###################################################################
def Rotate_Corners(X1, X2, Y1, Y2, PA) :
     """
     Rotate corners of a grid and return all corners
     Using a PA in degrees
     """
     rXa, rYa = rotPA(X1, Y1, PA)
     rXb, rYb = rotPA(X2, Y2, PA)
     rXc, rYc = rotPA(X1, Y2, PA)
     rXd, rYd = rotPA(X2, Y1, PA)
     return rXa, rXb, rXc, rXd, rYa, rYb, rYc, rYd
# ===============================================================
###################################################################
def XY_toSemiMajor(X, Y, Eps) :
     """
     Transform X, Y into Semi-Major axis radii using a given ellipticity

     Args:
        X (numpy.ndarray):
            Array of X coordinate
        Y (numpy.ndarray):
            Array of Y coordinate
        Eps (:obj:`float`):
            Ellipticity (1-b/a).

     """
     return sqrt(X**2 + Y**2 / (1. - Eps)**2)
# ===============================================================

####################################################################
# Function to select (Mask) good values from a map
# Using both Rectangle and Circular Zones from the Selection_Zone class
#
# Input is name of galaxy, and coordinates
####################################################################
def select_spaxels(maskDic, maskName, X, Y) :
   ## All spaxels are set to GOOD (True) first
   selgood = (X**2 >= 0)

   ## If no Mask is provided, we just return the full set of input X, Y
   if maskDic == None :
      return selgood

   ## We first check if the maskName is in the list of the defined Masks
   ## If the galaxy is not in the list, then the selection is all True
   if (maskDic.has_key(maskName)) :
      ## The mask is defined, so Get the list of Regions
      ## From the defined dictionary
      listRegions = maskDic[maskName]
      ## For each region, select the good spaxels
      for region in  listRegions :
         selgood = selgood & region.select(X, Y)

   return selgood
#=================================================================
####################################################################
# Parent class for the various types of Zones
####################################################################
class Selection_Zone :
    """
    Parent class for Rectangle_Zone and Circle_Zone
    """
    def __init__(self, params=None) :
        self.params = params

#=================================================================
####################################################################
class Rectangle_Zone(Selection_Zone) :
    type = "Rectangle"
    def select(self, xin, yin) :
        """ Define a selection within a rectangle
            It can be rotated by an angle theta (in degrees) """
        if self.params == None :
           return (xin**2 >=0)
        [x0,y0,length,width,theta] = self.params
        dx = xin - x0
        dy = yin - y0
        thetarad = theta * np.pi / 180.
        nx =   dx * cos(thetarad) + dy * sin(thetarad) 
        ny = - dx * sin(thetarad) + dy * cos(thetarad) 
        selgood = (np.abs(ny) > width/2.) | (np.abs(nx) > length/2.)
        return selgood
#=================================================================
####################################################################
class Circle_Zone(Selection_Zone) :
    type = "Circle"
    def select(self, xin, yin) :
        """ Define a selection within a circle """
        if self.params == None :
           return (xin**2 >=0)
        [x0,y0,radius] = self.params
        selgood = (sqrt((xin - x0)**2 + (yin - y0)**2) > radius)
        return selgood
#=================================================================

####################################################################
## Example of Selection Zone for a set of galaxies
## A3DMask = {'IC0598'      : [Circle_Zone([-13.6, +7.9, 2.6])],
##   'IC0676'    : [Circle_Zone([+10.5, +18.2, 6.2])]
## , 'NGC0448'   : [Circle_Zone([+19.0, -5.6, 2.2])]
## , 'NGC0509'   : [Circle_Zone([-19.1, -8.0, 3.0]), Circle_Zone([-0.9, +11.1, 3.0])]
## , 'NGC0661'   : [Circle_Zone([+11.9, -19.1, 4.2])]
## , 'NGC4486'   : [Rectangle_Zone([0.0,0.0,2.0]),Rectangle_Zone([12.0,4.8,12.0,4.0,23])]
## , 'NGC4486A'  : [Circle_Zone([+1.0, -2.3, 3.5]), Circle_Zone([+1.4, -9.8, 2.4]), Circle_Zone([-4.1, -13.3, 2.4]), Circle_Zone([+7.6, -14.1, 2.4])]
####################################################################
