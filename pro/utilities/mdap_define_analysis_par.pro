;+
; NAME:
;       MDAP_DEFINE_ANALYSIS_PAR
;
; PURPOSE:
;       Some of the analysis steps have parameters that can be set to
;       specify certain aspects of their execution.  These parameters
;       are collected into the structure defined here.  The currently
;       available parameters are:
;
;           AnalysisPar.moments
;               Number of *stellar* kinematics to fit.  The number of
;               emission-line kinematic moments is currently fixed at 2
;               (v, sigma).
;
;           AnalysisPar.degree
;               Degree of the *additive* polynomial used by pPXF when
;               fitting the stellar continuum.
;
;           AnalysisPar.mdegree
;               Degree of the *multiplicative* polynomial used by pPXF
;               and GANDALF when fitting the stellar continuum and full
;               spectrum, respectively.
;
;           AnalysisPar.reddening_order
;               The order of the reddening fit by GANDALF:  0 means no
;               reddening; 1 means fit a single dust-screen model; 2
;               means fit a dust-screen model plus a nebular only
;               component.  NOTE, if the reddening is fit, the value of
;               MDEGREE is ignored because fitting them both is
;               degenerate.
;
;           AnalysisPar.reddening
;               A two-element array used to set the initial guess for
;               the reddening coefficients.  The provided values must
;               match the input reddening_order, but the array must
;               always contain two elements.
;
; CALLING SEQUENCE:
;       result = MDAP_DEFINE_ANALYSIS_PAR()
;
; REVISION HISTORY:
;       22 Oct 2014: (KBW) Original implementation
;       10 Nov 2014: (KBW) Add the default values
;       05 Dec 2014: (KBW) Add some more comments
;-
;------------------------------------------------------------------------------
FUNCTION MDAP_DEFINE_ANALYSIS_PAR

        return, { AnalysisPar, moments:2, degree:-1, mdegree:-1, reddening:dblarr(2), $
                               reddening_order:0 }
END

