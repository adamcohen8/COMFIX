#Functions for COMFIX
import numpy as np
import math

########################################################
#Constants
MU = 398600.5
Deg = 180.0/math.pi
Rad = math.pi/180.0
Zero_IE = 0.015
Small = 1 * 10**(-6)
Ae = 6378.137
Ee = 0.0818191908426
w = np.array([0, 0, 0.00007292115])

########################################################
#Given Functions
def axisrot(A, axis, alpha):

    B = np.array([0.0,0.0,0.0])

    if axis == 1:
        B[0] = A[0]
        B[1] = A[1] * np.cos(alpha) + A[2] * np.sin(alpha)
        B[2] = -A[1] * np.sin(alpha) + A[2] * np.cos(alpha)

    if axis == 2:
        B[0] = A[0] * np.cos(alpha) - A[2] * np.sin(alpha)
        B[1] = A[1]
        B[2] = A[0] * np.sin(alpha) + A[2] * np.cos(alpha)

    if axis == 3:
        B[0] = A[0] * np.cos(alpha) + A[1] * np.sin(alpha)
        B[1] = -A[0] * np.sin(alpha) + A[1] * np.cos(alpha)
        B[2] = A[2]


    return B


def cuberoot(Xval):

    Temp = Xval/abs(Xval) * abs(Xval) ** (1/3)

    return Temp


def dayofyr2mdhms(Yr, Days):

    DayOfYr = int(math.floor(Days))

    LMonth = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]
    if (Yr-1900) % 4 == 0:
        LMonth[1] = 29

    i = 0
    IntTemp = 0

    while DayOfYr > IntTemp + LMonth[i] and i<11:
        IntTemp = IntTemp + LMonth[i]
        i = i + 1

    Mon = i + 1
    D = DayOfYr - IntTemp
    Temp = (Days - DayOfYr) * 24.0
    H = int(math.floor(Temp))
    Temp = (Temp - H) * 60.0
    M = int(math.floor(Temp))
    Sec = (Temp - M) * 60.0


    return [Mon, D, H, M, Sec]


def elOrb(R, V):

    magr = mag(R)
    magv = mag(V)

    A = -(MU/2) * 1/((magv**2)/2 -(MU/magr))

    h = np.cross(R,V)

    Incl = vecangle(h, [0.0,0.0,1.0])

    n = np.cross([0.0,0.0,1.0], h)

    RAAN = vecangle(n, [1.0,0.0,0.0])
    if n[1] < 0.0:
        RAAN = 2*math.pi - RAAN

    Ev = scalarMultiply((1/MU),np.subtract(scalarMultiply((magv**2 - MU/magr),R),scalarMultiply(np.dot(R,V),V)))

    Ecc = mag(Ev)

    Argp = vecangle(Ev, n)
    if Ev[2] < 0.0:
        Argp = 2*math.pi - Argp

    Nu = vecangle(Ev, R)
    if np.dot(R,V) < 0.0:
        Nu = 2*math.pi - Nu


    # Hbar = np.cross(R, V)
    # RdotV = np.dot(R, V)
    #
    # if mag(Hbar) > Small:
    #     Nbar = [-Hbar[1], Hbar[0], 0.0]
    #     Ebar = np.subtract((np.cross(V, Hbar)/MU), R/mag(R))
    #
    #     SME = mag(V)**2 * 0.5 - MU/mag(R)
    #     if abs(SME) > Small:
    #         A = -MU / (2.0*SME)
    #     else:
    #         A = math.inf
    #
    #     Ecc = mag(Ebar)
    #     P = mag(Hbar)**2 / MU
    #
    #
    #     Hk = Hbar[2]/mag(Hbar)
    #
    #     if abs(abs(Hk) - 1.0) < Zero_IE:
    #         if abs(Hbar[2])> 0.0:
    #             Hk = Hbar[2]/abs(Hbar[2])
    #
    #     Incl = np.arccos(Hk)
    #
    #     TypeOrbit = ["E", "I"]
    #     if Incl < Zero_IE*Rad or abs(math.pi - Incl) < Zero_IE*Rad:
    #         TypeOrbit[1] = 'E'
    #     if Ecc < Zero_IE:
    #         TypeOrbit[0] = "C"
    #     "".join(TypeOrbit)
    #
    #     RAAN = float("NaN")
    #     Argp = float("NaN")
    #     Nu = float("NaN")
    #     U = float("NaN")
    #     CapPi = float("NaN")
    #     L = float("NaN")
    #     M = float("NaN")
    #
    #
    #
    #     if TypeOrbit == "EI":
    #         RAAN = vecangle(Nbar, [1,0,0])
    #         if (not math.isinf(RAAN)) and Nbar[1] < 0.0:
    #             RAAN = 2*math.pi - RAAN
    #
    #         Argp = vecangle(Nbar, Ebar)
    #         if not math.isinf(Argp) and Ebar[2] < 0.0:
    #             Argp = 2*math.pi - Argp
    #
    #         Nu = vecangle(Ebar, R)
    #         if not math.isinf(Nu) and np.dot(R,V) < 0.0:
    #             Nu = 2*math.pi - Nu
    #
    #     if TypeOrbit == "EE":
    #         Nu = vecangle(Ebar, R)
    #         if not math.isinf(Nu) and np.dot(R, V) < 0.0:
    #             Nu = 2 * math.pi - Nu
    #
    #         CapPi = vecangle(Ebar, [1, 0, 0])
    #         if not math.isinf(CapPi) and Ebar[1] < 0.0:
    #             U = 2*math.pi - U
    #
    #     if TypeOrbit == "CI":
    #         RAAN = vecangle(Nbar, [1,0,0])
    #         if not math.isinf(RAAN) and Nbar[1] < 0.0:
    #             RAAN = 2*math.pi - RAAN
    #
    #         U = vecangle(Nbar, R)
    #         if not math.isinf(U) and R[2] < 0.0:
    #             U = 2*math.pi - U
    #
    #     if TypeOrbit == "CE":
    #         L = vecangle(R, [1,0,0])
    #         if not math.isinf(L) and R[2] < 0.0:
    #             L = 2*math.pi - L
    #
    #     if TypeOrbit != "EI" or TypeOrbit != "EE" or TypeOrbit != "CE" or TypeOrbit != "CI":
    #         print("Orbit Type Error in ElOrb")
    #
    #
    #     if Ecc-1.0 > Zero_IE:
    #         F = np.acosh((Ecc + np.cos(Nu))/(1.0+Ecc*np.cos(Nu)))
    #         M = Ecc*np.sinh(F)-F
    #         if Nu > math.pi:
    #             M = -M
    #     else:
    #         if abs(Ecc-1.0) < Zero_IE:
    #             D = math.sqrt(P) * np.tan(Nu*0.5)
    #             M = (3.0*P*D + D*D*D)/6.0
    #             if Nu > math.pi:
    #                 M = -M
    #             else:
    #                 if Ecc > Zero_IE:
    #                     Temp = 1.0 + Ecc*np.cos(Nu)
    #                     if abs(Temp) < Small:
    #                         M = 0.0
    #                     else:
    #                         E = 2 * np.atan(math.sqrt((1-Ecc)/(1+Ecc))*np.tan(Nu/2))
    #                         M = E - Ecc * np.sin(E)
    #                 else:
    #                     if Incl < Zero_IE*Rad or abs(Incl - math.pi) < Zero_IE*Rad:
    #                         M = L
    #                     else:
    #                         M = U
    # else:
    #     P = float("NaN")
    #     A = float("NaN")
    #     Ecc = float("NaN")
    #     Incl = float("NaN")
    #     RAAN = float("NaN")
    #     Argp = float("NaN")
    #     Nu = float("NaN")
    #     M = float("NaN")
    #     U = float("NaN")
    #     L = float("NaN")
    #     CapPi = float("NaN")

    return [A, Ecc, Incl, RAAN, Argp, Nu]

def scalarMultiply(scalar, vector):

    newvec = np.array([0.0,0.0,0.0])
    newvec[0] = scalar*vector[0]
    newvec[1] = scalar*vector[1]
    newvec[2] = scalar*vector[2]

    return newvec

def finddays(Year, Month, Day, Hr, Min, Sec):

    LMonth = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]
    if (Year - 1900) % 4 == 0:
        LMonth[1] = 29

    i = 1
    DDays = 0.0

    while i < Month and i < 12:
        DDays = DDays + LMonth[i-1]
        i = i + 1

    DDays = DDays + (Day - 1) + (Hr / 24.0) + (Min / 1440.0) + (Sec / 86400.0)


    return DDays


def fprintvec():

    return 0


def get_tle():

    return 0


def gstim0(Yr):

    JD = 367.0 * Yr - int(math.floor(7.0 * (Yr) *0.25)) + 30 + 1721014.5
    Tu = (int(math.floor(JD)) + 0.5 - 2451545.0) / 36525.0
    GST0 = 1.753368559 + 628.3319705 * Tu + 6.770708127**(-6) * Tu * Tu

    GST0 = revcheck(GST0, 2*math.pi)
    if GST0 < 0.0:
        GST0 = GST0 + 2*math.pi


    return GST0


def invjulianday(JD):

    Temp = JD - 2415019.5
    Tu = Temp / 365.25
    Yr = 1900 + int(math.floor(Tu))
    LeapYrs = int(math.floor((Yr - 1900 - 1) * 0.25))
    Days = Temp - ((Yr - 1900) * 365.0 + LeapYrs)

    if Days < 1.0:
        Yr = Yr - 1
        LeapYrs = int(math.floor((Yr - 1900 - 1) * 0.25))
        Days = Temp - ((Yr - 1900) * 365.0 + LeapYrs)

    x = dayofyr2mdhms(Yr, Days)

    Mon = x[0]
    D = x[1]
    H = x[2]
    M = x[3]
    S = x[4]

    return [Yr, Mon, D, H, M, S]


def julianday(Yr, Mon, D, H, M, S):

    TERM1 = 367.0 * Yr
    TERM2 = int(math.floor((7.0 * (Yr + int(math.floor((Mon + 9.0) / 12.0)))) * 0.25))
    TERM3 = int(math.floor(275.0 * Mon / 9.0))
    UT = ((S / 60.0 + M) / 60.0 + H) / 24.0

    JD = (TERM1 - TERM2 + TERM3) + D + UT + 1721013.5

    return JD


def mag(vec):

    magnitude = (vec[0]**2 + vec[1]**2 + vec[2]**2)**(1/2)

    return magnitude


def revcheck(x, modby):

    if modby != 0:
        y = x - modby*math.floor(x/modby)
    else:
        y = x

    return y


def Sun(JD):

    N = JD - 2451545.0

    MeanLong = 280.461 + 0.9856474 * N
    MeanLong = revcheck(MeanLong, 360.0)

    MeanAnomaly = 357.528 + 0.9856003 * N
    MeanAnomaly = revcheck(MeanAnomaly * Rad, 2*math.pi)
    if MeanAnomaly < 0.0:
        MeanAnomaly = 2*math.pi + MeanAnomaly

    EclpLong = MeanLong + 1.915*np.sin(MeanAnomaly) +0.02*np.sin(2.0*MeanAnomaly)
    Obliquity = 23.439 - 0.0000004 * N

    MeanLong = MeanLong * Rad
    if MeanLong < 0.0:
        MeanLong = 2*math.pi +MeanLong

    EclpLong = EclpLong * Rad
    Obliquity = Obliquity * Rad

    RtAsc = np.arctan(np.cos(Obliquity)*np.tan(EclpLong))

    if EclpLong < 0:
        EclpLong = EclpLong + 2*math.pi

    if abs(EclpLong-RtAsc) > 0.5*math.pi:
        RtAsc = RtAsc + 0.5*math.pi* round((EclpLong-RtAsc)/(0.5*math.pi))

    Decl = np.arcsin(np.sin(Obliquity)*np.sin(EclpLong))

    RSun4 = 1.00014 -0.01671*np.cos(MeanAnomaly) - 0.00014*np.cos(2.0*MeanAnomaly)

    RSun = np.array([RSun4*np.cos(EclpLong), RSun4*np.cos(Obliquity)*np.sin(EclpLong), RSun4*np.sin(Obliquity)*np.cos(EclpLong)])

    return [RSun, RtAsc, Decl]


def vecangle(A, B):

    tolerance = 0.000001

    if mag(A) * mag(B) >= tolerance:
        temp = np.dot(A,B)/mag(A)/mag(B)
        if abs(temp) > 1.0:
            temp = temp/abs(temp)
        Theta = np.arccos(temp)
    else:
        Theta = "NaN"

    return Theta

def writeOutput(file_out, sitlat, sitlon, sitalt, a, e, Incl, Raan, Argp, Nu, rho, az, el, drho, daz, Del, rho_sez, rho_ijk, drho_sez, drho_ijk, jd, satID, Yr, UT, GST, LST, R, V, Rsite):
    file_out.write(f"*********************Comfix  Satellite: {satID}*******************\n")
    file_out.write("---------------------------------------------------------------\n")
    file_out.write("             ####     Input Data      ####\n\n")
    file_out.write("   LAT      LON     ALT      YEAR     DAY          UT     \n")
    file_out.write("  (deg)    (deg)    (m)              Number   (hr:min:sec)\n")
    file_out.write(
        f"  {format(sitlat * 180 / math.pi, '.2f')}   {format(sitlon * 180 / math.pi, '.2f')}   {format(sitalt * 1000, '.2f')}     {2022}    {220}       {'19:45:10.00'}\n\n")
    file_out.write("   RHO         DRHO        EL           DEL         AZ          DAZ  \n")
    file_out.write("   (km)       (km/s)     (deg)        (deg/s)     (deg)       (deg/s) \n")
    file_out.write(
        f"   {round(rho, 4)}    {round(drho, 5)}    {round((el * 180 / math.pi), 5)}     {round((Del * 180 / math.pi), 5)}    {round((az * 180 / math.pi), 5)}    {round((daz * 180 / math.pi), 5)}\n")
    file_out.write("------------------------------------------------------------\n")
    file_out.write("             ####     Working Data    ####\n\n")
    file_out.write("     LAT        LON         ALT         Julian Date\n")
    file_out.write("    (rad)      (rad)       (km)           (days)    \n")
    file_out.write(f"   {format(sitlat, '.4f')}     {format(sitlat, '.4f')}      {format(sitalt, '.4f')}      {jd}\n\n")
    file_out.write("   RHO         DRHO        EL           DEL         AZ          DAZ  \n")
    file_out.write("   (km)       (km/s)     (rad)        (rad/s)     (rad)       (rad/s) \n")
    file_out.write(
        f"{format(rho, '.4f')}    {format(drho, '.5f')}   {format(el, '.5f')}       {format(Del, '.5f')}     {format(az, '.5f')}    {format(daz, '.5f')}\n\n")
    file_out.write(f"GST  =           {round(GST, 10)} rad     &        {round(GST * 180 / math.pi, 10)} degs\n")
    file_out.write(f"GST  =           {round(LST, 10)} rad     &        {round(LST * 180 / math.pi, 10)} degs\n\n")
    file_out.write("             ####       Vectors        ####\n\n")
    file_out.write(
        f"Rho_sez  =   {format(rho_sez[0], '5.4f')} S    {format(rho_sez[1], '5.4f')} E   {format(rho_sez[2], '5.4f')} Z Mag =   {format(mag(rho_sez), '5.4f')} km \n")
    file_out.write(
        f"Drho_sez =   {format(drho_sez[0], '5.4f')} S     {format(drho_sez[1], '5.4f')} E      {format(drho_sez[2], '5.4f')} Z Mag =     {format(mag(drho_sez), '5.4f')} km/s \n")
    file_out.write(
        f"R_site   =   {format(Rsite[0], '5.4f')} I   {format(Rsite[1], '5.4f')} J   {format(Rsite[2], '5.4f')} K Mag =   {format(mag(Rsite), '5.4f')} km \n")
    file_out.write(
        f"Rho_ijk  =   {format(rho_ijk[0], '5.4f')} I   {format(rho_ijk[1], '5.4f')} J   {format(rho_ijk[2], '5.4f')} K Mag =   {format(mag(rho_ijk), '5.4f')} km \n")
    file_out.write(
        f"Drho_ijk =   {format(drho_ijk[0], '5.4f')} I      {format(drho_ijk[1], '5.4f')} J      {format(drho_ijk[2], '5.4f')} K Mag =     {format(mag(drho_ijk), '5.4f')} km/s \n")
    file_out.write(
        f"R_ijk    =   {format(R[0], '5.4f')} I   {format(R[1], '5.4f')} J   {format(R[2], '5.4f')} K Mag =  {format(mag(R), '5.4f')} km \n")
    file_out.write(
        f"V_ijk    =   {format(V[0], '5.4f')} I      {format(V[1], '5.4f')} J      {format(V[2], '5.4f')} K Mag =      {format(mag(V), '5.4f')} km/s \n")
    file_out.write("------------------------------------------------------------------\n\n")
    file_out.write("             ####     Output Data      ####\n\n")
    file_out.write("CLASSIC ORBITAL ELEMENTS: \n\n")
    file_out.write(f"Semimajor Axis                     (km) =     {format(a, '.4f')}\n")
    file_out.write(f"Eccentricity                            =         {format(e, '.4f')}\n")
    file_out.write(f"Inclination                       (deg) =        {format(Incl, '.4f')}\n")
    file_out.write(f"Right Ascension of Ascending Node (deg) =       {format(Raan, '.4f')}\n")
    file_out.write(f"Argument of Perigee               (deg) =        {format(Argp, '.4f')}\n")
    file_out.write(f"True Anomaly                      (deg) =        {format(Nu, '.4f')}\n\n\n")


    return 0


#####################################################################
#My Functions

def gstlst(jd, sitlon):

    jd = str(jd)

    # Yr = jd[0] +jd[1]
    # Yr = int(Yr)

    day = jd[2] + jd[3] + jd[4]
    day = float(day)

    hr = jd[5] + jd[6]
    hr = float(hr)
    min = jd[7] + jd[8]
    min = float(min)
    sec = jd[9] + jd[10] + jd[11] + jd[12] + jd[13]
    sec = float(sec)

    D = day - 1.0 + (hr/24.0) + (min/1440.0) + (sec/86400.0)

    theta_g0 = 100.63004655

    theta_g = theta_g0 + 1.002737791737697*360*D

    gst = revcheck(theta_g, 360.0)
    gst = gst * math.pi /180.0
    lst = gst + sitlon

    return [gst, lst]


def site(sitlat, sitalt, lst):

    x = abs((Ae / (math.sqrt(1-(Ee**2)*np.sin(sitlat))))+sitalt)*np.cos(sitlat)
    z = abs((Ae*(1-Ee**2))/(math.sqrt(1-(Ee**2)*(math.sin(sitlat)**2)))+sitalt)*np.sin(sitlat)

    R_site = np.array([x*np.cos(lst), x*np.sin(lst), z])

    return R_site


def SEZ2IJK(vec_sez, sitlat, LST):

    vec_ijk = np.array([0.0,0.0,0.0])

    Colat = 0.5*math.pi - sitlat

    vec_ijk = axisrot(axisrot(vec_sez, 2, -Colat),3,-LST)

    return vec_ijk


def OBS2RANGERANGERATE(rho, az, el, drho, daz, Del):

    rho_sez = np.array([0.0,0.0,0.0])
    drho_sez = np.array([0.0, 0.0, 0.0])

    rho_sez[0] = -rho*np.cos(az)*np.cos(el)
    rho_sez[1] = rho*np.sin(az)*np.cos(el)
    rho_sez[2] = rho*np.sin(el)

    drho_sez[0] = -drho*np.cos(az)*np.cos(el) + rho*daz*np.sin(az)*np.cos(el) + rho*Del*np.cos(az)*np.sin(el)
    drho_sez[1] = drho*np.sin(az)*np.cos(el) + rho*daz*np.cos(az)*np.cos(el) - rho*Del*np.sin(az)*np.sin(el)
    drho_sez[2] = drho*np.sin(el) + rho*Del*np.cos(el)

    return [rho_sez, drho_sez]


def ijktorv(rho_ijk, drho_ijk, R_site):

    R = np.array([0.0, 0.0, 0.0])
    V = np.array([0.0, 0.0, 0.0])

    R[0] = R_site[0] + rho_ijk[0]
    R[1] = R_site[1] + rho_ijk[1]
    R[2] = R_site[2] + rho_ijk[2]

    b = np.cross(w, R)

    V[0] = drho_ijk[0] + b[0]
    V[1] = drho_ijk[1] + b[1]
    V[2] = drho_ijk[2] + b[2]

    return [R, V]


