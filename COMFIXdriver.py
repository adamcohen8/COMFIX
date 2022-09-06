#COMFIX driver
#Import Modules
import numpy as np
import COMFIXfuncts as CF
import math


#Place Holders
rho = 0
az = 0
el = 0
drho = 0
daz = 0
Del = 0
jd = 0
sitlon = 0
sitlat = 0
sitalt = 0


x = True

while x:

    #Get next set of Observation Data Here

    #Convert to working units
    az = az * math.pi / 180.0
    el = el * math.pi / 180.0
    daz = daz * math.pi / 180.0
    Del = Del * math.pi / 180.0


    #Find GST and LST
    GST_LST = CF.gstlst(jd, sitlon)
    GST = GST_LST[0]
    LST = GST_LST[1]

    #Find R site
    Rsite = CF.site(sitlat, sitalt, LST)

    #Find Rho and Drho in SEZ
    rho_drho = CF.OBS2RANGERANGERATE(rho, az, el, drho, daz, Del)
    rho_sez = rho_drho[0]
    drho_sez = rho_drho[1]

    #Convert Rho and Drho
    rho_ijk = CF.SEZ2IJK(rho_sez, sitlat, LST)
    drho_ijk = CF.SEZ2IJK(drho_sez, sitlat, LST)

    #Find R and V
    RV = CF.ijktorv(rho_ijk, drho_ijk, Rsite)
    R = RV[0]
    V = RV[1]

    #Find COEs
    COEs = CF.elOrb(R, V)
    a = COEs[1]
    e = COEs[2]
    Incl = COEs[3]
    Raan = COEs[4]
    Argp = COEs[5]
    Nu = COEs[6]

    #Convert to display units
    Incl = Incl * 180 / math.pi
    Raan = Raan * 180 / math.pi
    Argp = Argp * 180 / math.pi
    Nu = Nu * 180 / math.pi


    #Write values to files





