import COMFIXfuncts as CF
import numpy as np
import math

#Constants
Re = 6378.137
J2 = 1.08263 *10**(-3)
MU = 398600.5


def J2DragPert(inc0, ecc0, n0, ndot2):
#########################################################
#
#  Use           : [raandot, argpdot, eccdot] = J2DragPert(inc0, ecc0, n0, ndot2)
#
#   This Function calculates the RAAN rate, Argument of Perigee
#   rate, and eccentricity rate for an orbit due to
#   atmospheric drag and the J2 effect.
#
#   Author       : C2C Adam Cohen, DFAS,      3 Nov 2022
#
#   Input        :
#       inc0     - Initial Inclination (rad)
#       ecc0     - Initial Eccentricity (unitless)
#       n0       - Initial Mean Motion (rad/s)
#       ndot2    - Mean Motion Rate/2 (rad/s^2)
#
#   Output       :
#       raandot  - RAAN Rate (rad/s)
#       argpdot  - Argument of Perigee Rate (rad/s)
#       eccdot   - Eccentricity Rate (1/s)
#
#   Locals       : None.
#
#   Constants    :
#     MU           - Earth's Gravitational Parameter (398600.5 km^3/s^2)
#     J2           - Constant for use in J2 effect equations
#                    (1.08263*10^-3 unitless)
#
#   Coupling     : None.
#
#   References   :
#     Astro Engr 321 Course Handbook PREDICT project description
#
#########################################################
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
#########################################################
#
#  Use           : E = newton(M, e)
#
#   This function finds eccentric anomaly from mean
#   anomaly and eccentricity using the Newton-Raphson
#   method of solving transcendental equations
#
#   Author       : C2C Adam Cohen, DFAS,      3 Nov 2022
#
#   Input        :
#       M     - Mean Anomaly (rad)
#       e     - Eccentricity (unitless)
#
#   Output       :
#       E     - Eccentric Anomaly (rad)
#
#   Locals       : None.
#
#   Constants    :
#     1*10^-9           - Tolerance value for convergence
#
#   Coupling     : None.
#
#   References   :
#     Astro Engr 321 Course Handbook PREDICT project description
#
#########################################################
    #Initialize
    E_g = 0
    # M = or - e depending if M >pi or <pi
    if M < math.pi:
        E_n = M + e
    elif M > math.pi:
        E_n = M - e
    else:
        E_n = M

    #Until error is very small
    while abs(E_n-E_g) > 10.0**(-9.0):

        #Memory Allocation tricks
        TempVal = E_n
        E_n = -1
        E_g = TempVal

        #Newton Raphson Method of finding Eccentric Anomaly
        E_n = E_g + ((M-(E_g -e*np.sin(E_g)))/(1-e*np.cos(E_g)))

    #Assign final value to output
    E = E_n

    return E


def coeupdate(deltat, n0, ndot2, ecc0, eccdot, raan0, raandot, argp0, argpdot, mean0):
#########################################################
#
#  Use           : [n, ecc, raan, argp, nu] = coeupdate(deltat,
#  n0, ndot2, ecc0, eccdot, raan0, raandot, argp0, argpdot, mean0)
#
#   This function updates a set of COEs over a given time period
#   using the time rate of change of the COEs
#
#   Author       : C2C Adam Cohen, DFAS,      3 Nov 2022
#
#   Input        :
#       delatt  - Time to propogate over (sec)
#       n0      - Initial Mean Motion (rad/s)
#       ndot2   - Mean Motion Rate/2 (rad/s^2)
#       ecc0    - Initial Eccentricity (unitless)
#       eccdot  - Eccentricity Rate (1/s)
#       raan0   - Initial RAAN (rad)
#       raandot - RAAN Rate (rad/s)
#       argp0   - Initial Argument of Perigee (rad)
#       argpdot - Argument of Pergigee Rate (rad/s)
#       mean0   - Initial Mean Anomaly (rad)
#
#   Output       :
#       n       - New Mean Motion (rad/s)
#       ecc     - New Eccentricity (uniltess)
#       raan    - New RAAN (rad)
#       argp    - New Argument of Perigee (rad)
#       nu      - New True Anomaly
#
#   Locals       : None.
#
#   Constants    :
#     1*10^-9      - Tolerance value for convergence
#
#   Coupling     :
#       revcheck   - Function that takes an angle in rad or deg
#                    and outputs the equivalent angle between 0
#                    and 360 for degrees and 2 pi for radians
#       newton     - Function that outputs eccentric anomaly
#                    given mean anomaly and eccentricity using
#                    the Newton-Raphson method
#
#   References   :
#     Astro Engr 321 Course Handbook PREDICT project description
#
#########################################################
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

    #Make sure M is between 0 and 2 pi
    M = CF.revcheck(M, 2*math.pi)

    #Find new Eccentric Anomaly
    E = newton(M, ecc0)

    #Find new True Anomaly
    nu = np.arccos((np.cos(E)-ecc0)/(1.0-ecc0*np.cos(E)))

    #Plane Check
    if E > math.pi and nu < math.pi:
        nu = 2*math.pi - nu
    if E < math.pi and nu > math.pi:
        nu = 2*math.pi - nu
    nu = CF.revcheck(nu, 2*math.pi)

    return [n, ecc, raan, argp, nu]


def coes2rpqw(n, ecc, nu):
#########################################################
#
#  Use           : R_pqw = coes2rpqw(n, ecc, nu)
#
#   This function updates a set of COEs over a given time period
#   using the time rate of change of the COEs
#
#   Author       : C2C Adam Cohen, DFAS,      3 Nov 2022
#
#   Input        :
#       n       - Mean Motion (rad/s)
#       ecc     - Eccentricity (untiless)
#       nu      - True Anomaly (rad)
#
#   Output       :
#       R_pqw   - Position Vector in the PQW Frame (km)
#
#   Locals       :
#     p0      - Semilatus Rectum (km)
#     a       - Semimajor Axis (km)
#
#   Constants    :
#     MU      - Earth's Gravitational Paremter (398600.5 (km^3/s^2)
#
#   Coupling     : None.
#
#   References   :
#     Astro Engr 321 Course Handbook PREDICT project description
#
#########################################################
    # Find Semimajor Axis
    a = (MU / (n ** 2.0)) ** (1.0 / 3.0)

    # Find Semilatus Rectum
    p0 = a * (1.0 - (ecc ** 2.0))

    #Find Vbar
    PQ = np.array([1.0, 1.0, 0.0])
    PQ[0] = -np.sin(nu)
    PQ[1] = ecc + np.cos(nu)
    temp = math.sqrt(float(MU/p0))
    Vbar = CF.scalarMultiply(temp, PQ)

    #Find magnitude of Vbar
    V = CF.mag(Vbar)

    #Find R
    R = MU/(((V**2.0)/2.0)+(MU/(2.0*a)))

    #Find R_pqw
    PQ2 = np.array([1.0, 1.0, 0.0])
    PQ2[0] = np.cos(nu)
    PQ2[1] = np.sin(nu)
    R_pqw = CF.scalarMultiply(R,PQ2)

    return R_pqw


def pqw2ijk(vec_pqw, argp, inc, raan):
#########################################################
#
#  Use           : vec_ijk = pqw2ijk(vec_pqw, argp, inc, raan)
#
#   This function rotates a vector from the PQW frame to the
#   IJK frame
#
#   Author       : C2C Adam Cohen, DFAS,      3 Nov 2022
#
#   Input        :
#       vec_pqw - Vector in the PQW Frame (any units)
#       argp    - Argument of Perigee (rad)
#       inc     - Inclination (rad)
#       raan    - RAAN (rad)
#
#   Output       :
#       vec_ijk - Vector in the IJK frame (any units)
#
#   Locals       : None.
#
#   Constants    : None.
#
#   Coupling     :
#       axisrot    - Function that rotates a vector about
#                    a given axis by a given angle
#
#   References   :
#     Astro Engr 321 Course Handbook PREDICT project description
#
#########################################################

    #Rotate 3 times to get from PQW to IJK
    vec_ijk = CF.axisrot(vec_pqw, 3, -argp)
    vec_ijk = CF.axisrot(vec_ijk, 1, -inc)
    vec_ijk = CF.axisrot(vec_ijk, 3, -raan)

    return vec_ijk


def visible(R_ijk, R_site, sitlat, lst, jd):
#########################################################
#
#  Use           : Vis = visible(R_ijk, R_site, sitlat, lst, jd)
#
#   This function rotates a vector from the PQW frame to the
#   IJK frame
#
#   Author       : C2C Adam Cohen, DFAS,      3 Nov 2022
#
#   Input        :
#       R_ijk   - Position Vector in the IJK Frame (km)
#       R_site  - Site Position Vector (km)
#       sitlat  - Site Latitude (rad)
#       lst     - Local Sidereal Time (rad)
#       jd      - Julian day
#
#   Output       :
#       Vis     - Container storing:
#                   True if visibile, False if not visible
#                   Alpha angle
#                   Beta angle
#                   x
#                   1,2,3, or 4 if not visible to indicate
#                   which test was failed
#
#   Locals       :
#       Beta    - Angle between Rsun and R_ijk (rad)
#       Alpha   - Angle between R_site and Rsun (rad)
#       x       - R_ijk*sin(Beta) (km)
#       Rsun    - Vector from the sun to the earth (km)
#       rho     - Range (km)
#       az      - Azimuth Angle (rad)
#       el      - Elevation Angle (rad)
#
#   Constants    :
#       102 (deg)  - Required Alpha angle for Satellite
#                    to be in the dark
#
#   Coupling     :
#       Sun        - Function that finds the vector
#                    From the Sun to the Earth
#       vecangle   - Function that finds the angle
#                    between two vectors
#       mag        - Function that finds the magnitude
#                    of a vector
#
#   References   :
#     Astro Engr 321 Course Handbook PREDICT project description
#
#########################################################
    #Find Rho, az, and el
    [rho, az, el] = rhoazel(R_ijk, R_site, sitlat, lst)

    #Find Rsun
    Rsun = CF.Sun(jd)[0]

    #Find Beta, Alpha, and x
    Beta = CF.vecangle(Rsun, R_ijk)
    Alpha = CF.vecangle(R_site, Rsun)
    x = abs(CF.mag(R_ijk)*np.sin(2*math.pi - Beta))

    #Check if S/C is in Sun
    if Beta > math.pi/2:
        if x < Re:
            Vis = [False, Beta, Alpha, x, 1]
            return Vis
    #Check if S/C is in the dark
    if Alpha < 102.0*math.pi/180.0:
        Vis = [False, Beta, Alpha, x, 2]
        return Vis
    #Check if S/C is close enough to see
    if rho > 1500.0:
        Vis = [False, Beta, Alpha, x, 3]
        return Vis
    #See if elevation angle is large enough to see S/C
    if el < 10.0*math.pi/180.0:
        Vis = [False, Beta, Alpha, x, 4]
        return Vis
    #Return True if all tests pass
    Vis = [True, Beta, Alpha, x]
    return Vis


def rhoazel(R_ijk, R_site, sitlat, lst):
#########################################################
#
#  Use           : [rho, az, el] = rhoazel(R_ijk, R_site, sitlat, lst)
#
#   This function finds the Range, Azimuth angle, and Elevation
#   angle of a satellite from the site vector,site latitude,
#   satellite position vector, and lst
#
#   Author       : C2C Adam Cohen, DFAS,      3 Nov 2022
#
#   Input        :
#       R_ijk   - Position Vector in the IJK Frame (km)
#       R_site  - Site Position Vector (km)
#       sitlat  - Site Latitude (rad)
#       lst     - Local Sidereal Time (rad)
#
#   Output       :
#       rho     - Range (km)
#       az      - Azimuth Angle (rad)
#       el      - Elevation Angle (rad)
#
#   Locals       :
#       rho_ijk - Range in IJK (km)
#       rho_sez - Range in SEZ (km)
#
#   Constants    :
#       102 (deg)  - Required Alpha angle for Satellite
#                    to be in the dark
#
#   Coupling     :
#       mag        - Function that finds the magnitude
#                    of a vector
#       revcheck   - Function that takes an angle and returns
#                    the equivalent angle between 0 and 2 pi
#                    for radians and 360 for degrees
#       ijk2sez    - Function that takes a vector in the IJK
#                    frame and returns the vector in the SEZ
#                    frame
#
#   References   :
#     Astro Engr 321 Course Handbook PREDICT project description
#
#########################################################
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


def ijk2sez(vec_ijk, lst, sitlat):
#########################################################
#
#  Use           : vec_sez = ijk2sez(vec_ijk, lst, sitlat)
#
#   This function takes a vector in the IJK frame and returns the
#   vector in the SEZ frame
#
#   Author       : C2C Adam Cohen, DFAS,      3 Nov 2022
#
#   Input        :
#       vec_ijk - Vector in the IJK Frame (any units)
#       sitlat  - Site Latitude (rad)
#       lst     - Local Sidereal Time (rad)
#
#   Output       :
#       vec_sez - Vector in the SEZ Frame (any units)
#
#   Locals       : None.
#
#   Constants    : None.
#
#   Coupling     :
#       axisrot    - Function that rotates a vector about
#                    a given axis by a given angle
#
#   References   :
#     Astro Engr 321 Course Handbook PREDICT project description
#
#########################################################

    #Rotate Twice to get from IJK to SEZ
    vec_sez = CF.axisrot(vec_ijk, 3, lst)
    vec_sez = CF.axisrot(vec_sez, 2, (math.pi/2.0 - sitlat))

    return vec_sez








