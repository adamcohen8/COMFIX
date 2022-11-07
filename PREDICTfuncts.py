import COMFIXfuncts as CF
import numpy as np
import math

#Constants
Re = 6378.137
J2 = 1.08263 *10**(-3)
MU = 398600.5

#Tested
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

#Tested
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

#Tested, minor error
def coeupdate(deltat, n0, ndot2, ecc0, eccdot, raan0, raandot, argp0, argpdot, mean0):

    #Find new eccentricity
    ecc = ecc0 + eccdot*deltat
    if ecc < 0.0:
        ecc = 0.0

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

#Not Tested Yet
def coes2rpqw(n, ecc, nu):

    # Find Semimajor Axis
    a = (MU / (n ** 2.0)) ** (1.0 / 3.0)

    # Find Semilatus Rectum
    p0 = a * (1.0 - (ecc ** 2.0))

    #Find Vbar
    PQ = np.array([1.0, 1.0, 0.0])
    PQ[0] = -np.sin(nu)
    PQ[1] = ecc + np.cos(nu)
    temp = math.sqrt(float(MU/p0))
    #print(temp)

    Vbar = CF.scalarMultiply(temp, PQ)

    #Find magnitude of Vbar
    V = CF.mag(Vbar)
    #print(V)

    #Find R
    R = MU/(((V**2.0)/2.0)+(MU/(2.0*a)))
    #print(R)

    #Find R_pqw
    PQ2 = np.array([1.0, 1.0, 0.0])
    PQ2[0] = np.cos(nu)
    PQ2[1] = np.sin(nu)
    R_pqw = CF.scalarMultiply(R,PQ2)

    return R_pqw

#Not Tested Yet
def pqw2ijk(vec_pqw, argp, inc, raan):

    vec_ijk = CF.axisrot(vec_pqw, 3, -argp)
    vec_ijk = CF.axisrot(vec_ijk, 1, -inc)
    vec_ijk = CF.axisrot(vec_ijk, 3, -raan)

    return vec_ijk

#Not Tested
def visible(R_ijk, R_site, sitlat, lst, jd):

    #Find Rho, az, and el
    [rho, az, el] = rhoazel(R_ijk, R_site, sitlat, lst)

    #Find Rsun
    Rsun = CF.Sun(jd)[0]

    #Find Beta, Alpha, and x
    Beta = CF.vecangle(Rsun, R_ijk)
    Alpha = CF.vecangle(R_site, Rsun)
    x = abs(CF.mag(R_ijk)*np.sin(2*math.pi - Beta))

    if Beta > math.pi/2:
        if x < Re:
            return [False, Beta, Alpha, x, 1]

    if Alpha < 102.0*math.pi/180.0:
        return [False, Beta, Alpha, x, 2]

    if rho > 1500.0:
        return [False, Beta, Alpha, x, 3]

    if el < 10.0*math.pi/180.0:
        return [False, Beta, Alpha, x, 4]

    return [True, Beta, Alpha, x]

#Not Tested Yet
def rhoazel(R_ijk, R_site, sitlat, lst):

    #Find Rho in IJK
    rho_ijk = R_ijk - R_site
    #Rotate Rho to SEZ
    rho_sez = ijk2sez(rho_ijk, lst, sitlat)

    #Find Magnitude of Rho
    rho = CF.mag(rho_sez)

    #Find Azimuth Angle
    az = math.atan2(rho_sez[1], -rho_sez[0])
    az = CF.revcheck(az, 2*math.pi)

    #Find Elevation Angle
    el = np.arcsin(rho_sez[2]/rho)
    el = CF.revcheck(el, 2 * math.pi)

    return [rho, az, el]

#Not Tested Yet
def ijk2sez(vec_ijk, lst, sitlat):

    vec_sez = CF.axisrot(vec_ijk, 3, lst)
    vec_sez = CF.axisrot(vec_sez, 2, (math.pi/2.0 - sitlat))

    return vec_sez








