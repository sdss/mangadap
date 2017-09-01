################################################################################
#
# Copyright (C) 2017, Michele Cappellari
# E-mail: michele.cappellari_at_physics.ox.ac.uk
#
# Updated versions of the software are available from my web page
# http://purl.org/cappellari/software
#
# If you have found this software useful for your research,
# I would appreciate an acknowledgment to the use of the
# "CAPFIT program within the pPXF software distribution of Cappellari (2017),
#  which implements a Levenberg-Marquardt approach with the addition of
#  box constraints and fixed/tied variables".
#
# This software is provided as is without any warranty whatsoever.
# Permission to use, for non-commercial purposes is granted.
# Permission to modify for personal or internal use is granted,
# provided this copyright and disclaimer are included unchanged
# at the beginning of the file. All other rights are reserved.
#
################################################################################
#
# Implementation of the Levenberg-Marquardt (LM) algorithm with the addition
# of bound constraints and optionally fixed or tied variables.
#
# The meaning of the input and output parameters is fully consistent with
# the parameters with the same name in scipy.optimize.least_squares
# https://docs.scipy.org/doc/scipy/reference/generated/scipy.optimize.least_squares.html
#
# A general description of the unconstrained LM algorithm is given in
# - Chapter 5.2 of Fletcher R., 1987, Practical Methods of Optimization, 2nd ed., Wiley
#   http://doi.org/10.1002/9781118723203
# - Chapter 10.3 of Nocedal J. & Wright S.J., 2006, Numerical Optimization, 2nd ed., Springer
#   http://doi.org./10.1007/978-0-387-40065-5
#
# However, the present program differs from the standard algorithm as it
# includes bound constraints by solving a full quadratic programming problem at
# every iteration. This approach is robust and accurate, and is efficient when
# the function evaluation is at least as expensive to compute as the quadratic
# solution.
#
# The Jacobian scaling and convergence tests follow
# More', J.J., Garbow, B.S. & Hillstrom, K.E. 1980, User Guide for MINPACK-1,
# Argonne National Laboratory Report ANL-80-74 (http://cds.cern.ch/record/126569)
#
# V1.0.0: By Michele Cappellari, Oxford, 15 June 2017
#
################################################################################

from __future__ import print_function

import numpy as np
from scipy import optimize, linalg

################################################################################

def chi2(x):
    return x.dot(x)

################################################################################

def cov_err(jac):
    """
    Covariance and 1sigma formal errors calculation

    """
    U, s, Vh = linalg.svd(jac, full_matrices=False)
    w = s > np.finfo(float).eps*max(jac.shape)*s[0]
    cov = (Vh[w].T/s[w]**2).dot(Vh[w])
    perr = np.sqrt(np.diag(cov))

    return cov, perr

################################################################################

class capfit(object):

    def __init__(self, func, p1, abs_step=None, bounds=(-np.inf, np.inf),
                 diff_step=1e-4, fixed=None, ftol=1e-4, max_nfev=None, tied=None,
                 verbose=False, x_scale='jac', xtol=1e-4, args=(), kwargs={}):

        p1 = np.array(p1, dtype=float)  # Make copy to leave input unchanged
        bounds = np.asarray([np.resize(b, p1.size) for b in bounds])
        p1 = p1.clip(*bounds)   # Make initial guess feasible

        if fixed is None:
            fixed = np.full(p1.size, False)
        if tied is None:
            tied = np.full(p1.size, '')
        assert len(p1) == len(fixed) == len(tied), \
            "`x0`, `fixed` and `tied` must have the same size"

        self.nfev = 0
        self.njev = 0
        self.diff_step = diff_step
        self.abs_step = abs_step
        self.tied = np.asarray([a.strip() for a in tied])
        self.free = (np.asarray(fixed) == 0) & (self.tied == '')
        self.args = args
        self.kwargs = kwargs
        if max_nfev is None:
            self.max_nfev = 100*self.free.sum()
        self.ftol = ftol
        self.xtol = xtol
        self.verbose = verbose

        f1 = self.call(func, p1)
        J1 = self.fdjac(func, p1, f1)
        dd = linalg.norm(J1, axis=0)
        mx = np.max(dd)
        dd[dd < mx*np.finfo(float).eps] = 1  # As More'+80
        lam = 0.01*mx**2  # max(diag(J1.T @ J1))

        if verbose == 2:
            print("Start lam:", lam, "chi2:", chi2(f1), "\nStart p:", p1)

        while(1):

            if isinstance(x_scale, str) and x_scale == 'jac':
                dd = np.maximum(dd, linalg.norm(J1, axis=0))
            else:
                dd = np.ones_like(p1)/x_scale

            # Solve the LM system without explicitly creating J1.T @ J1
            D = np.diag(dd/np.max(dd))
            A = np.vstack([J1, np.sqrt(lam)*D])
            b = np.append(-f1, np.zeros_like(p1))
            h = optimize.lsq_linear(A, b, bounds=bounds-p1, method='bvls').x

            p2 = p1 + h
            f2 = self.call(func, p2)

            # Actual versus predicted chi2 reduction
            actred = 1 - chi2(f2)/chi2(f1)
            prered = 1 - chi2(f1 + J1.dot(h))/chi2(f1)
            ratio = actred/prered

            status = self.check_conv(lam, f2, p2, h, D, actred, prered)

            if status != -1:
                break

            # eq.(5.2.7) of Fletcher (1987)
            # Algorithm 4.1 in Nocedal & Wright (2006)
            if ratio < 0.25:
                lam *= 4
            elif ratio > 0.75:
                lam /= 2

            if ratio > 0.0001:  # Successful step: move on
                J2 = self.fdjac(func, p2, f2)
                p1, f1, J1 = p2, f2, J2

        self.x = p2
        self.cost = 0.5*chi2(f2)  # as in least_squares()
        self.fun = f2
        self.jac = J1
        self.grad = J1.T.dot(f2)
        self.status = status
        self.success = status > 0
        self.cov, self.x_err = cov_err(J1)

################################################################################

    def fdjac(self, func, pars, f):

        self.njev += 1
        jac = np.zeros([f.size, pars.size])
        if self.abs_step is None:
            h = self.diff_step*np.maximum(1.0, np.abs(pars))  # as in least_squares()
        else:
            h = self.abs_step

        # Compute derivative for free parameters
        w = np.flatnonzero(self.free)
        for j in w:
            pars1 = pars.copy()
            pars1[j] += h[j]
            f1 = self.call(func, pars1)
            jac[:, j] = (f1 - f)/h[j]

        return jac

################################################################################

    def call(self, func, p):

        self.nfev += 1
        w = np.flatnonzero(self.tied != '')
        for j in w:   # loop can be empty
            exec('p[' + str(j) + ']=' + self.tied[j])
        resid = func(p, *self.args, **self.kwargs)

        return resid

################################################################################

    def check_conv(self, lam, f, p, h, D, actred, prered):

        status = -1
        if self.nfev > self.max_nfev:
            status = 0   # Terminating on function evaluations count
        if prered < self.ftol and abs(actred) < self.ftol and actred <= 2*prered:
            status = 2   # Terminating for small function variation (More'+80)
        if linalg.norm(h) < self.xtol*(self.xtol + linalg.norm(D*p)):
            if status == 2:
                status = 4   # Both ftol and xtol convergence tests are satisfied
            else:
                status = 3   # Terminating for small step

        if self.verbose == 2:
            print('\niter:', self.njev, ' lam:', lam, ' chi2:', chi2(f),
                  ' ratio:', actred/prered, '\np:', p, '\nh:', h)

        if status != -1 and self.verbose != 0:
            print('\nFinal iter:', self.njev, ' Func calls:', self.nfev,
                  ' chi2:', chi2(f), ' Status:', status, '\nFinal p:', p)

        return status

################################################################################
