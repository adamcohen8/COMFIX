import PREDICTfuncts as PF
import COMFIXfuncts as CF
import math
import numpy as np


#Get user input here
sitlat = 0
sitlon = 0
sitalt = 0
jd = 0

#Open TLE File


#Open Output File


#For Each TLE

    #Get next TLE

    inc0 = 0
    ecc0 = 0
    n0 = 0
    ndot2 = 0
    raan0 = 0
    argp0 = 0
    mean0 = 0


    #Convert to working units


    #J2DragPert
    [raandot, argpdot, eccdot] = PF.J2DragPert(inc0, ecc0, n0, ndot2)


    #For Each time step
        deltat = 120.0


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

        #Determin Visibility
        Vis = PF.visible(Rijk, Rsite, sitlat, LST, jd)

        if Vis:
            #Write data to output file





