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

    #Initialize
    E_g = 0
    if M < math.pi:
        E_n = M + e
    elif M > math.pi:
        E_n = M - e
    else:
        E_n = M

    #i = 0

    #Until error is very small
    while abs(E_n-E_g) > 10.0**(-9.0):

        #Memory Allocation tricks
        TempVal = E_n
        E_n = -1
        E_g = TempVal

        #Newton Raphson Method of finding Eccentric Anomaly
        E_n = E_g + ((M-(E_g -e*np.sin(E_g)))/(1-e*np.cos(E_g)))
        #i += 1
        #print(i)

    E = E_n

    return E


def coeupdate(deltat, n0, ndot2, ecc0, eccdot, raan0, raandot, argp0, argpdot, mean0):

    #Find new eccentricity
    ecc = ecc0 +eccdot*deltat

    #Find new RAAN
    raan = raan0 + raandot*deltat

    #Find new Argument of Perapsis
    argp = argp0 +argpdot*deltat

    #Find ndot
    ndot = ndot2*2.0

    #Find new mean motion
    n = n0 +ndot*deltat

    #Find new Mean Anomaly
    M = mean0 + n0*deltat + ndot2*(deltat**2.0)
    M = CF.revcheck(M, 2*math.pi)

    #Find new Eccentric Anomaly
    E = newton(M, ecc0)

    #Find new True Anomaly
    nu = np.arccos((np.cos(E)-ecc0)/(1.0-ecc0*np.cos(E)))
    if E > math.pi and nu < math.pi:
        nu = 2*math.pi - nu
    if E < math.pi and nu > math.pi:
        nu = 2*math.pi - nu
    nu = CF.revcheck(nu, 2*math.pi)

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








