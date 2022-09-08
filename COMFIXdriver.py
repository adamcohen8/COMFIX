#COMFIX driver
#Import Modules
import numpy as np
import COMFIXfuncts as CF
import math


#Place Holders
# rho = 0
# az = 0
# el = 0
# drho = 0
# daz = 0
# Del = 0
# jd = 0
# sitlon = 0
# sitlat = 0
# sitalt = 0


#Open input file
file_in = open("comfix.in", 'r')
contents = file_in.read()
lines = contents.split('\n')

#Open output file
file_out = open("Comfix.out", 'w')

#Change this to when end of file is reached

for line in lines:
    #Get next set of Observation Data Here
    obsData = line.split(' ')
    #print(obsData)
    sitlat = float(obsData[0])
    sitlon = float(obsData[1])
    sitalt = float(obsData[2])
    jd = str(obsData[3])
    satID = obsData[4]
    rho = float(obsData[5])
    az = float(obsData[6])
    el = float(obsData[7])
    drho = float(obsData[8])
    daz = float(obsData[9])
    Del = float(obsData[10])

    #Convert to working units
    az = az * math.pi / 180.0
    el = el * math.pi / 180.0
    daz = daz * math.pi / 180.0
    Del = Del * math.pi / 180.0
    sitalt = sitalt/1000.0
    sitlat = sitlat * math.pi / 180.0
    sitlon = sitlon * math.pi / 180.0


    #Find GST and LST
    GST_LST = CF.gstlst(jd, sitlon)
    GST = GST_LST[0]
    LST = GST_LST[1]


    print(GST, LST)

    #Find R site
    Rsite = CF.site(sitlat, sitalt, LST)

    print(Rsite)


    #Find Rho and Drho in SEZ
    rho_drho = CF.OBS2RANGERANGERATE(rho, az, el, drho, daz, Del)
    rho_sez = rho_drho[0]
    drho_sez = rho_drho[1]

    print(rho_sez, drho_sez)

    #Convert Rho and Drho
    rho_ijk = CF.SEZ2IJK(rho_sez, sitlat, LST)
    drho_ijk = CF.SEZ2IJK(drho_sez, sitlat, LST)

    print(rho_ijk, drho_ijk)

    #Find R and V
    RV = CF.ijktorv(rho_ijk, drho_ijk, Rsite)
    R = RV[0]
    V = RV[1]

    print(R,V)
    #Find COEs
    COEs = CF.elOrb(R, V)
    a = COEs[0]
    e = COEs[1]
    Incl = COEs[2]
    Raan = COEs[3]
    Argp = COEs[4]
    Nu = COEs[5]

    #Convert to display units
    Incl = Incl * 180 / math.pi
    Raan = Raan * 180 / math.pi
    Argp = Argp * 180 / math.pi
    Nu = Nu * 180 / math.pi


    #Write values to files

    print(a,e, Incl, Raan, Argp, Nu)



