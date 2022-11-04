import PREDICTfuncts as PF
import COMFIXfuncts as CF
import math
import numpy as np


#Get user input here
sitlat = 39.006
sitlon = -104.883
sitalt = 2.184
Year = 2022
StartDay = 280
StopDay = 300

jd = 0

#Convert User input to working units
sitlat = 39.006 *math.pi/180.0
sitlon = -104.883 *math.pi/180.0
sitalt = 2.184 *math.pi/180.0

StartJD = CF.julianday(Year, 10, 7, 0, 0, 0)
StopJD = CF.julianday(Year, 10, 27, 0, 0, 0)


#Open TLE File


#Open Output File

#For Each TLE
for i in range(1):
    #Get next TLE

    inc0 = 51.64050000 #deg
    ecc0 = 0.00031520
    n0 = 15.49768145 #rev/day
    ndot2 = 0.00100985 #rev/day^2
    raan0 = 148.07640000 #deg
    argp0 = 245.20250000 #deg
    mean0 = 227.77290000 #deg
    EpochDay = 278.90089709 #days

    #Convert to working units

    inc0 = inc0 *math.pi/180.0 #rad
    n0 = n0 *2*math.pi/86400 #rad/s
    ndot2 = ndot2 *2*math.pi/(86400**2) #rad/s^2
    raan0 = raan0 *math.pi/180.0 #rad
    argp0 = argp0 *math.pi/180.0 #rad
    mean0 = mean0 *math.pi/180.0 #rad

    [Mon, D, H, M, Sec] = CF.dayofyr2mdhms(Year, EpochDay)
    EpochJD = CF.julianday(Year, Mon, D, H, M, Sec)


    #J2DragPert
    [raandot, argpdot, eccdot] = PF.J2DragPert(inc0, ecc0, n0, ndot2)


    #For Each time step
    for t in range(0,20*24*30):

        ObsDay = EpochDay + (t*120.0)/86400.0
        [Mon, D, H, M, Sec] = CF.dayofyr2mdhms(Year, ObsDay)
        jd = CF.julianday(Year, Mon, D, H, M, Sec)

        deltat = t*120.0

        [GST, LST] = CF.gstlst(jd, sitlon)
        Rsite = CF.site(sitlat, sitalt, LST)

        #Find New COEs
        newCOEs = PF.coeupdate(deltat, n0, ndot2, ecc0, eccdot, raan0, raandot, argp0, argpdot, mean0)

        n = newCOEs[0]
        ecc = newCOEs[1]
        raan = newCOEs[2]
        argp = newCOEs[3]
        nu = newCOEs[4]

        #Find Rpqw
        Rpqw = PF.coes2rpqw(n, ecc, nu)

        #Rotate Rpqw into IJK
        Rijk = PF.pqw2ijk(Rpqw, argp, inc0, raan)

        #Find Rho in IJK
        Rho_ijk = Rijk - Rsite

        #Rotate Rho_ijk to SEZ
        Rho_sez = PF.ijk2sez(Rho_ijk, LST, sitlat)

        #Find Rho az el
        [Rho, az, el] = PF.rhoazel(Rijk, Rsite, sitlat, LST)
        print(Rho, az, el)

        #Determin Visibility
        Vis = PF.visible(Rijk, Rsite, sitlat, LST, jd)

        if Vis:
            #Write data to output file
            print("Visible")




