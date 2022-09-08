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

    file_out.write(f"*********************Comfix  Satellite: {satID}*******************\n")
    file_out.write("---------------------------------------------------------------\n")
    file_out.write("             ####     Input Data      ####\n\n")
    file_out.write("   LAT      LON     ALT      YEAR     DAY          UT     \n")
    file_out.write("  (deg)    (deg)    (m)              Number   (hr:min:sec)\n")
    file_out.write("\n\n")
    file_out.write("   RHO         DRHO        EL           DEL         AZ          DAZ  \n")
    file_out.write("   (km)       (km/s)     (deg)        (deg/s)     (deg)       (deg/s) \n")
    file_out.write(f"   {round(rho,4)}    {round(drho, 5)}    {round((el*180/math.pi),5)}     {round((Del*180/math.pi),5)}    {round((az*180/math.pi),5)}    {round((daz*180/math.pi),5)}\n")
    file_out.write("------------------------------------------------------------\n")
    file_out.write("             ####     Working Data    ####\n\n")
    file_out.write("     LAT        LON         ALT         Julian Date\n")
    file_out.write("    (rad)      (rad)       (km)           (days)    \n")
    file_out.write("\n\n")
    file_out.write("   RHO         DRHO        EL           DEL         AZ          DAZ  \n")
    file_out.write("   (km)       (km/s)     (rad)        (rad/s)     (rad)       (rad/s) \n")
    file_out.write("\n\n")
    file_out.write(f"GST  =           {round(GST, 10)} rad     &        {round(GST * 180 / math.pi, 10)} degs\n")
    file_out.write(f"GST  =           {round(LST, 10)} rad     &        {round(LST * 180 / math.pi, 10)} degs\n\n")
    file_out.write("             ####       Vectors        ####\n\n")
    file_out.write(f"Rho_sez  =   {round(rho_sez[0], 4)} S    {round(rho_sez[1], 4)} E   {round(rho_sez[2], 4)} Z Mag =   {round(CF.mag(rho_sez),4)} km \n")
    file_out.write("Drho_sez =       {5.9367} S     {-4.4833} E      {4.2295} Z Mag =      {8.5577} km/s \n")
    file_out.write("R_site   =   {-1357.7548} I   {-120.7082} J   {6209.9350} K Mag =   {6357.7796} km \n")
    file_out.write("Rho_ijk  =    {1381.8950} I   {-441.5631} J   {1547.8381} K Mag =   {2121.4180} km \n")
    file_out.write("Drho_ijk =      {-7.0721} I      {3.8723} J      {2.8678} K Mag =      {8.5577} km/s \n")
    file_out.write("R_ijk    =      {24.1401} I   {-562.2713} J   {7757.7731} K Mag =   {7778.1601} km \n")
    file_out.write("V_ijk    =      {-7.0311} I      {3.8740} J      {2.8678} K Mag =      {8.5246} km/s \n")
    file_out.write("------------------------------------------------------------------\n\n")
    file_out.write("             ####     Output Data      ####\n\n")
    file_out.write("CLASSIC ORBITAL ELEMENTS: \n\n")
    file_out.write("Semimajor Axis                     (km) =     {13365.4112}")
    file_out.write("Eccentricity                            =         {0.4991}")
    file_out.write("Inclination                       (deg) =        {93.4987}")
    file_out.write("Right Ascension of Ascending Node (deg) =       {329.8944}")
    file_out.write("Argument of Perigee               (deg) =        {33.3380}")
    file_out.write("True Anomaly                      (deg) =        {54.4301}")

