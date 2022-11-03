import COMFIXfuncts as CF
import numpy as np
import math

#Constants
Re = 6378.137
J2 = 1.08263 *10**(-3)
MU = 398600.5

def J2DragPert(inc0, ecc0, n0, ndot2):

    #Find Semimajor Axis
    a = (MU / (n0**2.0))**(1.0/3.0)

    #Find Semilatus Rectum
    p0 = a*(1.0-(ecc0**2.0))

    #Find Nbar
    nbar = n0 * (1.0 + (3.0/2.0)*J2*((Re/p0)**2.0)*math.sqrt(1.0-(ecc0**2.0))*(1.0-(3.0/2.0)*(np.sin(inc0)**2)))

    #Find raandot
    raandot = nbar*((-3.0/2.0)*J2*((Re/p0)**2.0)*np.cos(inc0))

    #Find argpdot
    argpdot = nbar*((3.0/2.0)*J2*((Re/p0)**2.0)*(2.0 -(5.0/2.0)*(np.sin(inc0)**2.0)))

    #Find eccdot
    eccdot = (-2.0*(1.0-ecc0)*ndot2*2.0)/(3*n0)

    return [raandot, argpdot, eccdot]


def newton(M, e):

    return E


def coeupdate(deltat, n0, ndot2, ecc0, eccdot, raan0, raandot, argp0, argpdot, mean0):

    return [n, ecc, raan, argp, nu]


def coes2rpqw(n, ecc, nu):

    return R_pqw


def pqw2ijk(vec_pqw, argp, inc, raan):

    return vec_ijk


def visible(R_ijk, R_site, sitlat, lst, jd):

    return Vis


def rhoazel(R_ijk, R_site, sitlat, lst):

    return [rho, az, el]


def ijk2sez(vec_ijk, lst, sitlat):

    return vec_sez








