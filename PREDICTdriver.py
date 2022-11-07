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

StartJD = CF.julianday(Year, 10, 7, 0, 0, 0)
StopJD = CF.julianday(Year, 10, 27, 0, 0, 0)

#Open TLE File
file_in = open("predictin.txt", 'r')
contents = file_in.read()
lines = contents.split('\n')

#Open output file
file_out = open("predict.out", 'w')

#Write user inputs to output file
file_out.write("PREDICT Observation Results By: C2C Adam Cohen\n\n")
file_out.write(" Observer Data:\n\n")
file_out.write(f"    Site Latitude  =    {sitlat} deg\n")
file_out.write(f"    Site Longitude =  {sitlon} deg\n")
file_out.write(f"    Site Altitude  =     {sitalt} km\n\n")
file_out.write(f"    Year           =  {Year}\n")
file_out.write(f"    Start Day      =  {StartDay}\n")
file_out.write(f"    End Day        =  {StopDay}\n\n")
file_out.write(" Echo Check Working Units\n\n")
#Convert User input to working units
sitlat = sitlat *math.pi/180.0
sitlon = sitlon *math.pi/180.0
file_out.write(f"    Site Latitude  =  {sitlat}\n")
file_out.write(f"    Site Longitude =  {sitlon}\n\n")
file_out.write(f"    Start Julian Day = {StartJD}\n")
file_out.write(f"    Stop Julian Day  = {StartJD}\n\n\n\n")
file_out.write("*********************************************************************\n")

#For Each TLE
for line in lines:
    #Get Data from Input File
    Data = line.split(' ')
    #Get next TLE
    SatName = Data[0]
    inc0 = float(Data[12]) #51.64050000 #deg
    raan0 = float(Data[13])  # 148.07640000 #deg
    ecc0 = float(f"0.{Data[14]}") #0.00031520
    argp0 = float(Data[15]) #245.20250000  # deg
    mean0 = float(Data[16]) #227.77290000  # deg
    n0 = float(Data[17]) #15.49768145 #rev/day
    ndot2 = float(Data[6]) #0.00100985 #rev/day^2
    Epochjd = Data[4]
    EpochDay = float(f"{Epochjd[2]}{Epochjd[3]}{Epochjd[4]}{Epochjd[5]}{Epochjd[6]}{Epochjd[7]}{Epochjd[8]}{Epochjd[9]}{Epochjd[10]}{Epochjd[11]}{Epochjd[12]}{Epochjd[13]}")# 278.90089709 #days
    Yr = int(f"{Epochjd[0]}{Epochjd[1]}")#22

    #Write input data to output file
    file_out.write(f" Echo Check Input Data for: {SatName}\n")
    file_out.write("*********************************************************************\n\n")
    file_out.write(f"    Epoch Year     =         20{Yr}\n")
    file_out.write(f"    Epoch Day      = {EpochDay}\n")
    file_out.write(f"    N Dot/2        = {ndot2} rev/day^2\n")
    file_out.write(f"    Inclination    = {inc0} deg\n")
    file_out.write(f"    RAAN           = {raan0} deg\n")
    file_out.write(f"    Eccentricity   = {ecc0}\n")
    file_out.write(f"    Arg of Perigee = {argp0} deg\n")
    file_out.write(f"    Mean Anomaly   = {mean0} deg\n")
    file_out.write(f"    Mean Motion    = {n0} rev/day\n\n")


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

    # Write Working Units to output file
    file_out.write(" Echo Check Working Units\n\n")
    file_out.write(f"    Epoch Time (month/day, hr:min:sec) =  {Mon}/ {D},  {H}:{M}:{format(Sec, '.2f')}\n\n")
    file_out.write(f"    Epoch Julian Day (days) =  {EpochJD}\n\n")
    file_out.write(f"    N    = {format(n0, '.11f')} rad/s Ndot2     =  +{ndot2} rad/s^2\n")
    file_out.write(f"    Ecc  = {format(ecc0, '.11f')}       Ndot2     =  {eccdot} 1/s\n")
    file_out.write(f"    Inc  = {format(inc0, '.11f')} rad   Mean Anom =  +{mean0} rad\n")
    file_out.write(f"    RAAN = {format(raan0, '.11f')} rad   RAAN dot  =  {raandot} rad/s\n")
    file_out.write(f"    Argp = {format(argp0, '.11f')} rad   Argp dot  =  +{argpdot} rad/s\n\n")

    file_out.write("SATELLITE NAME             MON/DAY  HR:MIN(UT)    RHO(KM)       AZ(DEG)   EL(DEG)\n")
    file_out.write("---------------------------------------------------------------------------------\n")


    #print(raandot, argpdot, eccdot)
    deltat = ((StartDay-EpochDay)*24*60*60)

    #For Each time step
    for t in range(0,20*24*30):

        ObsDay = EpochDay + (StartDay-EpochDay) + (t*120.0)/86400.0
        [Mon, D, H, M, Sec] = CF.dayofyr2mdhms(Year, ObsDay)
        jd = CF.julianday(Year, Mon, D, H, M, Sec)
        JD = f"{Yr}{int(math.floor(ObsDay))}{H:02d}{M:02d}{Sec:02f}"

        #print(JD)
        deltat = deltat + 120.0
        #print(deltat)

        [GST, LST] = CF.gstlst(JD, sitlon)
        Rsite = CF.site(sitlat, sitalt, LST)

        #Find New COEs
        newCOEs = PF.coeupdate(deltat, n0, ndot2, ecc0, eccdot, raan0, raandot, argp0, argpdot, mean0)
        #print(deltat)
        #print(newCOEs)
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

        if 100480 < deltat < 100490:
            print(Rho, az*180.0/math.pi, el*180.0/math.pi)
            print(Rijk, CF.mag(Rijk))
            #print(ecc)
            #print(raan*180.0/math.pi)
            #print(argp*180.0/math.pi)
            print(deltat)
            #print(Rsite, CF.mag(Rsite))
            #print(Rho_ijk, CF.mag(Rho_ijk))
            #print(Rho_sez, CF.mag(Rho_sez))
            #print(GST, LST)
            print(Vis)
            #print(jd)

        if Vis[0]:
            #Write data to output file
            #print("Visible")
            file_out.write(f"{SatName}                 {Mon}/ {D:02d}     {H:02d}:{M:02d}        {format(Rho, '.3f')}     {format(az*180.0/math.pi,'.3f')}    {format(el*180.0/math.pi, '.3f')}\n")




