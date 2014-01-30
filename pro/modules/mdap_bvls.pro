;######################################################################
;
; Copyright (C) 1999-2001, Michele Cappellari
; E-mail: cappellari@strw.leidenuniv.nl
;
; Updated versions of the software are available from my web page
; http://www.strw.leidenuniv.nl/~mcappell/idl
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
;   BVLS
;
; AUTHOR:
;   Michele Cappellari, Leiden Observatory, The Netherlands
;   cappellari@strw.leidenuniv.nl
;
; PURPOSE:
;   Perform bound constrained linear least-squares minimization
;
; CATEGORY:
;   Least Squares
;
; CALLING SEQUENCE:
;   BVLS, A, B, BND, X, $
;       EPS=eps, /FASTNORM, IERR=ierr, INDEX=index, ITER=iter, $
;       ITMAX=itmax, NSETP=nsetp, RNORM=rnorm, W=w
;
; DESCRIPTION:
;
;   Given an M by N matrix, A, and an M-vector, B,  compute an
;   N-vector, X, that solves the least-squares problem A # X = B
;   subject to X(J) satisfying  BND(0,J) <= X(J) <= BND(1,J)
;
;   The values BND(0,J) = -(MACHAR()).XMAX and BND(1,J) = (MACHAR()).XMAX
;   are suggested choices to designate that there is no constraint in that
;   direction.
;
;   This algorithm is a generalization of  NNLS, that solves
;   the least-squares problem,  A # X = B,  subject to all X(J) >= 0.
;   The subroutine NNLS appeared in 'SOLVING LEAST SQUARES PROBLEMS,'
;   by Lawson and Hanson, Prentice-Hall, 1974.  Work on BVLS was started
;   by C. L. Lawson and R. J. Hanson at Jet Propulsion Laboratory,
;   1973 June 12.  Many modifications were subsequently made.
;   The Fortran 90 code was completed in April, 1995 by R. J. Hanson.
;   The BVLS package is an additional item for the reprinting of the book
;   by SIAM Publications and the distribution of the code package
;   using netlib and Internet or network facilities.
;
;   This IDL version was ported from the original Fortran 90 code
;   by Michele Cappellari, Leiden Observatory, The Netherlands
;   E-mail: cappellari@strw.leidenuniv.nl
;
; INPUT PARAMETERS:
;
;   A(:,:)     [INTENT(InOut)]
;       On entry A() contains the M by N matrix, A.
;       On return A() contains the product matrix, Q*A, where
;       Q is an M by M orthogonal matrix generated by this
;       subroutine.  The dimensions are M=size(A,1) and N=size(A,2).
;
;   B(:)     [INTENT(InOut)]
;       On entry B() contains the M-vector, B.
;       On return, B() contains Q*B.  The same Q multiplies A.
;
;   BND(:,:)  [INTENT(In)]
;       BND(1,J) is the lower bound for X(J).
;       BND(2,J) is the upper bound for X(J).
;       Require:  BND(1,J) <= BND(2,J).
;
; OUTPUT PARAMETER:
;
;   X(:)    [INTENT(Out)]
;       On entry X() need not be initialized.  On return,
;       X() will contain the solution N-vector.
;
; KEYWORD PARAMETERS:
;
;   RNORM    [INTENT(Out)]
;       The Euclidean norm of the residual vector, b - A*X.
;
;   NSETP    [INTENT(Out)]
;       Indicates the number of components of the solution
;       vector, X(), that are not at their constraint values.
;
;   W(:)     [INTENT(Out)]
;       An N-array.  On return, W() will contain the dual solution
;       vector.   Using Set definitions below:
;       W(J) = 0 for all j in Set P,
;       W(J) <= 0 for all j in Set Z, such that X(J) is at its
;       lower bound, and
;       W(J) >= 0 for all j in Set Z, such that X(J) is at its
;       upper bound.
;       If BND(1,J) = BND(2,J), so the variable X(J) is fixed,
;       then W(J) will have an arbitrary value.
;
;   INDEX(:)    [INTENT(Out)]
;       An INTEGER working array of size N.  On exit the contents
;       of this array define the sets P, Z, and F as follows:
;
;   INDEX(1)   through INDEX(NSETP) = Set P.
;   INDEX(IZ1) through INDEX(IZ2)   = Set Z.
;   INDEX(IZ2+1) through INDEX(N)   = Set F.
;   IZ1 = NSETP + 1 = NPP1
;       Any of these sets may be empty.  Set F is those components
;       that are constrained to a unique value by the given
;       constraints.   Sets P and Z are those that are allowed a non-
;       zero range of values.  Of these, set Z are those whose final
;       value is a constraint value, while set P are those whose
;       final value is not a constraint.  The value of IZ2 is not returned.
;...    It is computable as the number of bounds constraining a component
;...    of X uniquely.
;
;   IERR    [INTENT(Out)]
;   Indicates status on return.
;   = 0   Solution completed.
;       = 1   M <= 0 or N <= 0
;       = 2   B(:), X(:), BND(:,:), W(:), or INDEX(:) size or shape violation.
;       = 3   Input bounds are inconsistent.
;       = 4   Exceed maximum number of iterations.
;
;   EPS [real(kind(one))]
;       Determines the relative linear dependence of a column vector
;       for a variable moved from its initial value.  This is used in
;       one place with the default value EPS=(MACHAR()).EPS.  Other
;       values, larger or smaller may be needed for some problems.
;
;   ITMAX  [integer]
;       Set to 3*N.  Maximum number of iterations permitted.
;       This is usually larger than required.
;
;   ITER   [integer]
;       Iteration counter.
;
;   /FASTNORM
;       Perform Euclidean Norm computation without checking for over/underflows.
;       It can speed up the program considerably when M is large, but has to
;       be used with care since may lead to instabilities!
;
; MODIFICATION HISTORY:
;   V1.0: Written by Michele Cappellari, Padova, 2000
;   V1.1: Added /FASTNORM keyword, MC, Leiden, 20 September 2001
;   V1.2: Use MAKE_ARRAY to deal with float or double arrays,
;       MC, Leiden, 19 October 2001
;   V1.3: Added compilation options and converted to IDL V5.0,
;       MC, Leiden 20 May 2002
;   V1.4: Define optional parameters using optional keywords.
;       The new calling sequence is not directly compatible with
;       the previous versions. MC, Leiden, 20 March 2004
;-
;----------------------------------------------------------------------
FUNCTION NRM2, X, fastNorm
;
;   NRM2 returns the Euclidean norm of a vector so that
;
;   NRM2 := sqrt( x'*x )
;
COMPILE_OPT IDL2, HIDDEN

IF fastNorm THEN RETURN, SQRT(TOTAL(X^2)) ; brute force approach: use with care!

ZERO = 0.0
ONE = 1.0

N = N_ELEMENTS(X)
IF( N LT 1)THEN $
    NORM  = ZERO $
ELSE IF( N EQ 1 )THEN $
    NORM  = ABS( X[0] ) $
ELSE BEGIN
    SCALE = ZERO
    SSQ   = ONE
;
    FOR IX = 0L, N-1 DO BEGIN
        ABSXI = ABS( X[IX] )
        IF(ABSXI GT ZERO )THEN $
            IF( SCALE LT ABSXI )THEN BEGIN
                SSQ   = ONE + SSQ*( SCALE/ABSXI )^2
                SCALE = ABSXI
            ENDIF ELSE $
                SSQ   = SSQ + ( ABSXI/SCALE )^2
    ENDFOR
    NORM  = SCALE * SQRT( SSQ )
ENDELSE
;
RETURN, NORM
END
;----------------------------------------------------------------------
PRO TERMINATION
;
COMPILE_OPT IDL2, HIDDEN

COMMON bvls_global, FIND, HITBND, FREE1, FREE2, FREE, IERR, M, N, I, $
    IBOUND, II, IP, ITER, ITMAX, IZ, IZ1, IZ2, J, JJ, JZ, L, LBOUND, $
    NPP1, NSETP, INDEX, ZERO, ONE, TWO, A, B, S, X, W, Z, BND, ALPHA, $
    ASAVE, CC, EPS, RANGE, RNORM, NORM, SM, SS, T, UNORM, UP, ZTEST, FASTNORM
;
;    IF (IERR LE 0) THEN BEGIN
;
;   Compute the norm of the residual vector.
    SM = ZERO
    IF (NPP1 LE M) THEN $
        SM = NRM2(B[NPP1-1:M-1],fastNorm) $
    ELSE $
        W[0:N-1] = ZERO
    RNORM = SM
;    ENDIF; ( IERR...)

END ; ( TERMINATION )
;----------------------------------------------------------------------
PRO MOVE_COEF_I_FROM_SET_P_TO_SET_Z
;
COMPILE_OPT IDL2, HIDDEN

COMMON bvls_global

X[I-1] = BND[IBOUND-1,I-1]
IF (ABS(X[I-1]) GT ZERO AND JJ GT 0) THEN $
    B[0:JJ-1] = B[0:JJ-1]-A[0:JJ-1,I-1]*X[I-1]

;   The following loop can be null.
FOR J=JJ+1,NSETP DO BEGIN
    II = INDEX[J-1]
    INDEX[J-2] = II

    SCALE = TOTAL(ABS(A[J-2:J-1,II-1]))
    IF (SCALE GT ZERO) THEN BEGIN
        R = SCALE * SQRT(TOTAL( (A[J-2:J-1,II-1]/SCALE)^2 ))
        IF (ABS(A[J-2,II-1]) GT ABS(A[J-1,II-1])) THEN $
            ROE = A[J-2,II-1] $
        ELSE $
            ROE = A[J-1,II-1]
        IF (ROE LT ZERO) THEN R = -R
        CC = A[J-2,II-1]/R
        SS = A[J-1,II-1]/R
        A[J-2,II-1] = R
    ENDIF ELSE BEGIN
        CC = ONE
        SS = ZERO
    ENDELSE

    SM = A[J-2,II-1]
;
;   The plane rotation is applied to two rows of A and the right-hand
;   side.  One row is moved to the scratch array S and THEN the updates
;   are computed.  The intent is for array operations to be performed
;   and minimal extra data movement.  One extra rotation is applied
;   to column II in this approach.
    S = A[J-2,0:N-1]
    A[J-2,0:N-1] = CC*S+SS*A[J-1,0:N-1]
    A[J-1,0:N-1] = CC*A[J-1,0:N-1]-SS*S
    A[J-2,II-1] = SM
    A[J-1,II-1] = ZERO
    SM = B[J-2]
    B[J-2] = CC*SM+SS*B[J-1]
    B[J-1] = CC*B[J-1]-SS*SM
ENDFOR
;
NPP1 = NSETP
NSETP = NSETP-1
IZ1 = IZ1-1
INDEX[IZ1-1] = I

END ; ( MOVE COEF I FROM SET P TO SET Z )
;----------------------------------------------------------------------
PRO  SEE_IF_ALL_CONSTRAINED_COEFFS_ARE_FEASIBLE
;
COMPILE_OPT IDL2, HIDDEN

COMMON bvls_global
;
;   See if each coefficient in set P is strictly interior to its constraint region.
;   If so, set HITBND = false.
;   If not, set HITBND = true, and also set ALPHA, JJ, and IBOUND.
;   Then ALPHA will satisfy  0.  < ALPHA  <=  1.
;
ALPHA=TWO
FOR IP=1L,NSETP DO BEGIN
    L = INDEX[IP-1]
    IF  (Z[IP-1]  LE  BND[0,L-1]) THEN $
;   Z(IP) HITS LOWER BOUND
        LBOUND=1 $
    ELSE  IF  (Z[IP-1]  GE  BND[1,L-1]) THEN $
;   Z(IP) HITS UPPER BOUND
        LBOUND=2 $
    ELSE $
        LBOUND = 0
;
    IF  ( LBOUND   NE   0 ) THEN BEGIN
        T = (BND[LBOUND-1,L-1]-X[L-1])/(Z[IP-1]-X[L-1])
        IF  (ALPHA   GT  T) THEN BEGIN
            ALPHA = T
            JJ = IP
            IBOUND = LBOUND
        ENDIF; ( LBOUND )
    ENDIF; ( ALPHA   >  T )
ENDFOR
HITBND = ABS(ALPHA - TWO) GT ZERO

END ;( SEE IF ALL CONSTRAINED COEFFS ARE FEASIBLE )
;----------------------------------------------------------------------
PRO TEST_SET_P_AGAINST_CONSTRAINTS
;
COMPILE_OPT IDL2, HIDDEN

COMMON bvls_global

WHILE 1 DO BEGIN
;   The solution obtained by solving the current set P is in the array Z().
;
    ITER = ITER+1
    IF (ITER GT ITMAX) THEN BEGIN
        IERR = 4
        GOTO, fine_LOOPB
    ENDIF
;
    SEE_IF_ALL_CONSTRAINED_COEFFS_ARE_FEASIBLE
;
;   The above call sets HITBND.  If HITBND = true THEN it also sets
;   ALPHA, JJ, and IBOUND.
    IF (NOT HITBND) THEN GOTO, fine_LOOPB
;
;   Here ALPHA will be between 0 and 1 for interpolation
;   between the old X() and the new Z().
    FOR IP=1L,NSETP DO BEGIN
        L = INDEX[IP-1]
        X[L-1] = X[L-1]+ALPHA*(Z[IP-1]-X[L-1])
    ENDFOR
;
    I = INDEX[JJ-1]
;   Note:  The exit test is done at the end of the loop, so the loop
;   will always be executed at least once.
    WHILE 1 DO BEGIN
;
;   Modify A(*,*), B(*) and the index arrays to move coefficient I
;   from set P to set Z.
;
        MOVE_COEF_I_FROM_SET_P_TO_SET_Z
;
        IF (NSETP LE 0) THEN GOTO, fine_LOOPB
;
;   See if the remaining coefficients in set P are feasible.  They should
;   be because of the way ALPHA was determined.  If any are infeasible
;   it is due to round-off error.  Any that are infeasible or on a boundary
;   will be set to the boundary value and moved from set P to set Z.
;
        IBOUND = 0
        FOR JJ=1L,NSETP DO BEGIN
            I = INDEX[JJ-1]
            IF  (X[I-1] LE BND[0,I-1]) THEN BEGIN
                IBOUND = 1
                GOTO, fine_ciclo1
            ENDIF ELSE IF (X[I-1] GE BND[1,I-1]) THEN BEGIN
                IBOUND = 2
                GOTO, fine_ciclo1
            ENDIF
        ENDFOR
        fine_ciclo1:
        IF (IBOUND LE 0) THEN GOTO, fine_ciclo2
    ENDWHILE
    fine_ciclo2:
;
;   Copy B( ) into Z( ).  Then solve again and loop back.
    Z[0:M-1] = B[0:M-1]
;
    FOR I=NSETP,1,-1 DO BEGIN
        IF (I NE NSETP) THEN $
            Z[0:I-1] = Z[0:I-1]-A[0:I-1,II-1]*Z[I]
        II = INDEX[I-1]
        Z[I-1] = Z[I-1]/A[I-1,II-1]
    ENDFOR
ENDWHILE
fine_LOOPB:

;   The following loop can be null.
FOR IP=1L,NSETP DO BEGIN
    I = INDEX[IP-1]
    X[I-1] = Z[IP-1]
ENDFOR

END ; ( TEST SET P AGAINST CONSTRAINTS)
;----------------------------------------------------------------------
PRO MOVE_J_FROM_SET_Z_TO_SET_P
;
COMPILE_OPT IDL2, HIDDEN

COMMON bvls_global
;
;   The index  J=index(IZ)  has been selected to be moved from
;   set Z to set P.  Z() contains the old B() adjusted as though X(J) = 0.
;   A(*,J) contains the new Householder transformation vector.
B[0:M-1] = Z[0:M-1]
;
INDEX[IZ-1] = INDEX[IZ1-1]
INDEX[IZ1-1] = J
IZ1 = IZ1+1
NSETP = NPP1
NPP1 = NPP1+1
;   The following loop can be null or not required.
NORM = A[NSETP-1,J-1]
A[NSETP-1,J-1] = UP
IF(ABS(NORM) GT ZERO) THEN BEGIN
    FOR JZ=IZ1,IZ2 DO BEGIN
        JJ = INDEX[JZ-1]
        SM = TOTAL(A[NSETP-1:M-1,J-1]/NORM * A[NSETP-1:M-1,JJ-1])/UP
        A[NSETP-1:M-1,JJ-1] = A[NSETP-1:M-1,JJ-1]+SM*A[NSETP-1:M-1,J-1]
    ENDFOR
    A[NSETP-1,J-1] = NORM
ENDIF
;   The following loop can be null.
FOR L=NPP1,M DO A[L-1,J-1] = ZERO
;
W[J-1] = ZERO
;
;   Solve the triangular system.  Store this solution temporarily in Z().
FOR I=NSETP,1,-1 DO BEGIN
    IF (I NE NSETP) THEN Z[0:I-1] = Z[0:I-1]-A[0:I-1,II-1]*Z[I]
    II = INDEX[I-1]
    Z[I-1] = Z[I-1]/A[I-1,II-1]
ENDFOR

END ; ( MOVE J FROM SET Z TO SET P )
;----------------------------------------------------------------------
PRO TEST_COEF_J_FOR_DIAG_ELT_AND_DIRECTION_OF_CHANGE
;
COMPILE_OPT IDL2, HIDDEN

COMMON bvls_global
;
;   The sign of W(J) is OK for J to be moved to set P.
;   Begin the transformation and check new diagonal element to avoid
;   near linear dependence.
;
ASAVE = A[NPP1-1,J-1]
;
;   Construct a Householder transformation.

VNORM = NRM2(A[NPP1-1:M-1,J-1],fastNorm)
IF (A[NPP1-1,J-1] GT ZERO) THEN VNORM = -VNORM
UP = A[NPP1-1,J-1] - VNORM
A[NPP1-1,J-1] = VNORM

IF (NSETP LT 1) THEN UNORM = 0.0 ELSE UNORM = NRM2(A[0:NSETP-1,J-1],fastNorm)
IF (ABS(A[NPP1-1,J-1]) GT EPS*UNORM) THEN BEGIN
;
;   Column J is sufficiently independent.  Copy b into Z, update Z.
    Z[0:M-1] = B[0:M-1]
; Compute product of transormation and updated right-hand side.
    NORM = A[NPP1-1,J-1]
    A[NPP1-1,J-1] = UP
    IF (ABS(NORM) GT ZERO) THEN BEGIN
        SM = TOTAL(A[NPP1-1:M-1,J-1]/NORM * Z[NPP1-1:M-1])/UP
        Z[NPP1-1:M-1] = Z[NPP1-1:M-1]+SM*A[NPP1-1:M-1,J-1]
        A[NPP1-1,J-1] = NORM
    ENDIF

    IF (ABS(X[J-1]) GT ZERO) THEN $
        Z[0:NPP1-1] = Z[0:NPP1-1]+A[0:NPP1-1,J-1]*X[J-1]
;   Adjust Z() as though X(J) had been reset to zero.
    IF ( FREE ) THEN $
        FIND = 1 $
    ELSE BEGIN
;
;   Solve for ZTEST ( proposed new value for X(J) ).
;   Then set FIND to indicate if ZTEST has moved away from X(J) in
;   the expected direction indicated by the sign of W(J).
        ZTEST = Z[NPP1-1]/A[NPP1-1,J-1]
        FIND = ( W[J-1] LT ZERO AND ZTEST LT X[J-1] ) OR $
               ( W[J-1] GT ZERO AND ZTEST GT X[J-1] )
    ENDELSE
ENDIF
;
;   If J was not accepted to be moved from set Z to set P,
;   restore A(NNP1,J).  Failing these tests may mean the computed
;   sign of W(J) is suspect, so here we set W(J) = 0.  This will
;   not affect subsequent computation, but cleans up the W() array.
IF  ( NOT FIND ) THEN BEGIN
    A[NPP1-1,J-1] = ASAVE
    W[J-1] = ZERO
ENDIF; ( .not. FIND )

END ;TEST_COEF_J_FOR_DIAG_ELT_AND_DIRECTION_OF_CHANGE
;----------------------------------------------------------------------
PRO  SELECT_ANOTHER_COEFF_TO_SOLVE_FOR
;
COMPILE_OPT IDL2, HIDDEN

COMMON bvls_global
;
;   1. Search through set z for a new coefficient to solve for.
;   First select a candidate that is either an unconstrained
;   coefficient or ELSE a constrained coefficient that has room
;   to move in the direction consistent with the sign of its dual
;   vector component.  Components of the dual (negative gradient)
;   vector will be computed as needed.
;   2. For each candidate start the transformation to bring this
;   candidate into the triangle, and THEN do two tests:  Test size
;   of new diagonal value to avoid extreme ill-conditioning, and
;   the value of this new coefficient to be sure it moved in the
;   expected direction.
;   3. If some coefficient passes all these conditions, set FIND = true,
;   The index of the selected coefficient is J = INDEX(IZ).
;   4. If no coefficient is selected, set FIND = false.
;
FIND = 0
FOR IZ=IZ1,IZ2 DO BEGIN
    J = INDEX[IZ-1]
;
;   Set FREE1 = true if X(J) is not at the left end-point of its
;   constraint region.
;   Set FREE2 = true if X(J) is not at the right end-point of its
;   constraint region.
;   Set FREE = true if X(J) is not at either end-point of its
;   constraint region.
;
    FREE1 = X[J-1] GT BND[0,J-1]
    FREE2 = X[J-1] LT BND[1,J-1]
    FREE = FREE1 AND FREE2

    IF ( FREE ) THEN $
        TEST_COEF_J_FOR_DIAG_ELT_AND_DIRECTION_OF_CHANGE $
    ELSE BEGIN
;   Compute dual coefficient W(J).
        W[J-1] = TOTAL(A[NPP1-1:M-1,J-1] * B[NPP1-1:M-1])
;
;   Can X(J) move in the direction indicated by the sign of W(J)?
;
        IF ( W[J-1] LT ZERO ) THEN BEGIN
            IF ( FREE1 ) THEN $
                TEST_COEF_J_FOR_DIAG_ELT_AND_DIRECTION_OF_CHANGE
        ENDIF ELSE  IF ( W[J-1]  GT ZERO ) THEN BEGIN
            IF ( FREE2 ) THEN $
                TEST_COEF_J_FOR_DIAG_ELT_AND_DIRECTION_OF_CHANGE
        ENDIF
    ENDELSE
    IF ( FIND ) THEN RETURN
ENDFOR;  IZ

END ; ( SELECT ANOTHER COEF TO SOLVE FOR )
;----------------------------------------------------------------------
PRO INITIALIZE, eps1, itmax1
;
COMPILE_OPT IDL2, HIDDEN

COMMON bvls_global

IF (N LT 2 OR M LT 2 OR N_ELEMENTS(B) NE M $
  OR (SIZE(BND))[1] NE 2 OR (SIZE(BND))[2] NE N) THEN BEGIN
    IERR = 2
    MESSAGE, 'Wrong input arrays size'
ENDIF

IERR = 0
mch = MACHAR()
IF N_ELEMENTS(eps1) EQ 0 THEN EPS = mch.EPS ELSE EPS = eps1
IF N_ELEMENTS(itmax1) EQ 0 THEN ITMAX = 3L*N ELSE ITMAX = itmax1
ITER = 0L
;
IZ2 = N
IZ1 = 1L
NSETP = 0L
NPP1 = 1L
;
;   Begin:  Loop on IZ to initialize  X().
IZ = IZ1
WHILE 1 DO BEGIN
    IF (IZ GT IZ2) THEN GOTO, fine1
    J = INDEX[IZ-1]
    IF (BND[0,J-1] LE -mch.XMAX) THEN $
        IF (BND[1,J-1] GE mch.XMAX) THEN $
            X[J-1] = ZERO $
        ELSE $
            X[J-1] = ZERO < BND[1,J-1] $
    ELSE IF (BND[1,J-1] GE mch.XMAX) THEN $
        X[J-1] = ZERO > BND[0,J-1] $
    ELSE BEGIN
        RANGE = BND[1,J-1] - BND[0,J-1]
        IF (RANGE LE ZERO) THEN BEGIN
;
;   Here X(J) is constrained to a single value.
            INDEX[IZ-1] = INDEX[IZ2-1]
            INDEX[IZ2-1] = J
            IZ = IZ-1
            IZ2 = IZ2-1
            X[J-1] = BND[0,J-1]
            W[J-1] = ZERO
        ENDIF ELSE IF (RANGE GT ZERO) THEN $
;
;   The following statement sets X(J) to 0 if the constraint interval
;   includes 0, and otherwise sets X(J) to the endpoint of the
;   constraint interval that is closest to 0.
;
            X[J-1] = BND[0,J-1] > (BND[1,J-1] < ZERO) $
        ELSE BEGIN
            IERR = 3
            RETURN
        ENDELSE ; ( RANGE:.)
    ENDELSE
;
;   Change B() to reflect a nonzero starting value for X(J).
;
    IF (ABS(X[J-1]) GT ZERO) THEN $
        B[0:M-1] = B[0:M-1]-A[0:M-1,J-1]*X[J-1]
    IZ = IZ+1
ENDWHILE; ( IZ <= IZ2 )
fine1:

END ; ( INITIALIZE )
;----------------------------------------------------------------------
PRO mdap_bvls, A1, B1, BND1, X1, $
    EPS=eps1, FASTNORM=fastNorm1, IERR=ierr1, INDEX=index1, $
    ITER=iter1, ITMAX=itmax1, NSETP=nsetp1, RNORM=rnorm1, W=w1


COMPILE_OPT IDL2
ON_ERROR, 2

COMMON bvls_global
;
; Load needed input parameters into the COMMON block
;
A = TEMPORARY(A1)
B = TEMPORARY(B1)
BND = TEMPORARY(BND1)

siz = SIZE(A)
M = siz[1]
N = siz[2]
X = MAKE_ARRAY(N,TYPE=siz[3])
W = MAKE_ARRAY(N,TYPE=siz[3])
INDEX = LINDGEN(N)+1

S = MAKE_ARRAY(N,TYPE=siz[3])
Z = MAKE_ARRAY(M,TYPE=siz[3])
;
; Load some constants into the COMMON block
;
ZERO = 0.0
ONE = 1.0
TWO = 2.0
IF KEYWORD_SET(fastNorm1) THEN fastNorm = 1 ELSE fastNorm = 0

INITIALIZE, eps1, itmax1
;
;   The above call will set IERR.
;
WHILE 1 DO BEGIN
    ;
    ;   Quit on error flag, or if all coefficients are already in the
    ;   solution, .or. if M columns of A have been triangularized.
    IF (IERR NE 0 OR IZ1 GT IZ2 OR NSETP GE M) THEN GOTO, fine
    ;
    SELECT_ANOTHER_COEFF_TO_SOLVE_FOR
    ;
    ;   See if no index was found to be moved from set Z to set P.
    ;   Then go to termination.
    IF  ( NOT FIND ) THEN GOTO, fine
    ;
    MOVE_J_FROM_SET_Z_TO_SET_P
    ;
    TEST_SET_P_AGAINST_CONSTRAINTS
    ;
    ;   The above call may set IERR.
    ;   All coefficients in set P are strictly feasible.  Loop back.
ENDWHILE
fine:
;
TERMINATION

A1 = TEMPORARY(A)
B1 = TEMPORARY(B)
BND1 = TEMPORARY(BND)
X1 = TEMPORARY(X)
RNORM1 = RNORM
NSETP1 = NSETP
W1 = TEMPORARY(W)
INDEX1 = TEMPORARY(INDEX)
IERR1 = IERR
ITER1 = ITER

END
;----------------------------------------------------------------------