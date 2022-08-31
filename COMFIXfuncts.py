#Functions for COMFIX
import numpy as np
import math

########################################################
#Constants
MU = 398600.5



########################################################
#Given Functions
def axisrot():

    return 0


def cuberoot(Xval):

    Temp = Xval ** (1/3)

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
    Sec = (Temp - M) * 60.0;


    return [Mon, D, H, M, Sec]


def elorb():

    return 0


def finddays():

    return 0


def fprintvec():

    return 0


def get_tle():

    return 0


def gstim0():

    return 0


def invjulianday():

    return 0


def julianday():

    return 0


def mag(vec):

    magnitude = (vec[0]**2 + vec[1]**2 + vec[2]**2)**(1/2)

    return magnitude


def revcheck(x, modby):

    if modby != 0:
        y = x - modby*math.floor(x/modby)
    else:
        y = x

    return y


def sun():

    return 0


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


#####################################################################
#My Functions

def gstlst(jd, sitlon):

    jd = str(jd)

    day = jd[2] + jd[3] + jd[4]
    day = int(day)

    hr = jd[5] + jd[6]
    hr = int(hr)
    min = jd[7] + jd[8]
    min = int(min)
    sec = jd[9] + jd[10] + jd[11] + jd[12] + jd[13]
    sec = float(sec)

    D = day - 1 + (hr/24) + (min/1440) + (sec/86400)

    theta_g0 = 100.63004655

    theta_g = theta_g0 + 1.002737791737697*360*D

    gst = revcheck(theta_g, 360)
    lst = gst + sitlon

    return [gst, lst]


def site(sitlat, sitlon, lst):

    return 0


def SEZ2IJK(vec_sez):

    return 0


def OBS2RANGERANGERATE(rho, az, el, drho, daz, Del):

    return 0


def ijktorv(rho_ijk, drho_ijk, R_site):

    return 0


