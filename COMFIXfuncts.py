#Functions for COMFIX
import numpy
import math
########################################################
#Given Functions
def axisrot():

    return 0


def cuberoot():

    return 0


def dayofyr2mdhms():

    return 0


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


def mag():

    return 0


def revcheck(x, modby):

    if modby != 0:
        y = x - modby*math.floor(x/modby)
    else:
        y = x

    return y


def sun():

    return 0


def vecangle():

    return 0


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

    gst = revcheck(theta_g)
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


