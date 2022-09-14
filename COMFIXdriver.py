#COMFIX driver
#Import Modules
import numpy as np
import COMFIXfuncts as CF
import math


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


    #print(GST, LST)

    #Find R site
    #print(sitlat, sitalt, LST)
    Rsite = CF.site(sitlat, sitalt, LST)
    #Rsite = np.array([-1357.7548, -120.7082, 6209.9350])

    #print(Rsite)


    #Find Rho and Drho in SEZ
    rho_drho = CF.OBS2RANGERANGERATE(rho, az, el, drho, daz, Del)
    rho_sez = rho_drho[0]
    drho_sez = rho_drho[1]

    #print(rho_sez, drho_sez)

    #Convert Rho and Drho
    rho_ijk = CF.SEZ2IJK(rho_sez, sitlat, LST)
    drho_ijk = CF.SEZ2IJK(drho_sez, sitlat, LST)

    #print(rho_ijk, drho_ijk)

    #Find R and V
    RV = CF.ijktorv(rho_ijk, drho_ijk, Rsite)
    R = RV[0]
    V = RV[1]

    #print(R,V)
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

    Yr = 0
    UT = 0

    #Write values to files

    #print(a,e, Incl, Raan, Argp, Nu)

    CF.writeOutput(file_out, sitlat, sitlon, sitalt, a, e, Incl, Raan, Argp, Nu, rho, az, el, drho, daz, Del, rho_sez,rho_ijk, drho_sez, drho_ijk, jd, satID, GST, LST, R, V, Rsite)

    # file_out.write(f"*********************Comfix  Satellite: {satID}*******************\n")
    # file_out.write("---------------------------------------------------------------\n")
    # file_out.write("             ####     Input Data      ####\n\n")
    # file_out.write("   LAT      LON     ALT      YEAR     DAY          UT     \n")
    # file_out.write("  (deg)    (deg)    (m)              Number   (hr:min:sec)\n")
    # file_out.write(f"  {format(sitlat * 180 / math.pi, '.2f')}   {format(sitlon * 180 / math.pi, '.2f')}   {format(sitalt * 1000, '.2f')}     {2022}    {220}       {'19:45:10.00'}\n\n")
    # file_out.write("   RHO         DRHO        EL           DEL         AZ          DAZ  \n")
    # file_out.write("   (km)       (km/s)     (deg)        (deg/s)     (deg)       (deg/s) \n")
    # file_out.write(f"   {round(rho,4)}    {round(drho, 5)}    {round((el*180/math.pi),5)}     {round((Del*180/math.pi),5)}    {round((az*180/math.pi),5)}    {round((daz*180/math.pi),5)}\n")
    # file_out.write("------------------------------------------------------------\n")
    # file_out.write("             ####     Working Data    ####\n\n")
    # file_out.write("     LAT        LON         ALT         Julian Date\n")
    # file_out.write("    (rad)      (rad)       (km)           (days)    \n")
    # file_out.write(f"   {format(sitlat, '.4f')}     {format(sitlat, '.4f')}      {format(sitalt, '.4f')}      {jd}\n\n")
    # file_out.write("   RHO         DRHO        EL           DEL         AZ          DAZ  \n")
    # file_out.write("   (km)       (km/s)     (rad)        (rad/s)     (rad)       (rad/s) \n")
    # file_out.write(f"{format(rho, '.4f')}    {format(drho, '.5f')}   {format(el, '.5f')}       {format(Del, '.5f')}     {format(az, '.5f')}    {format(daz, '.5f')}\n\n")
    # file_out.write(f"GST  =           {round(GST, 10)} rad     &        {round(GST * 180 / math.pi, 10)} degs\n")
    # file_out.write(f"GST  =           {round(LST, 10)} rad     &        {round(LST * 180 / math.pi, 10)} degs\n\n")
    # file_out.write("             ####       Vectors        ####\n\n")
    # file_out.write(f"Rho_sez  =   {format(rho_sez[0], '5.4f')} S    {format(rho_sez[1], '5.4f')} E   {format(rho_sez[2], '5.4f')} Z Mag =   {format(CF.mag(rho_sez),'5.4f')} km \n")
    # file_out.write(f"Drho_sez =   {format(drho_sez[0], '5.4f')} S     {format(drho_sez[1], '5.4f')} E      {format(drho_sez[2], '5.4f')} Z Mag =     {format(CF.mag(drho_sez),'5.4f')} km/s \n")
    # file_out.write(f"R_site   =   {format(Rsite[0], '5.4f')} I   {format(Rsite[1], '5.4f')} J   {format(Rsite[2], '5.4f')} K Mag =   {format(CF.mag(Rsite),'5.4f')} km \n")
    # file_out.write(f"Rho_ijk  =   {format(rho_ijk[0], '5.4f')} I   {format(rho_ijk[1], '5.4f')} J   {format(rho_ijk[2], '5.4f')} K Mag =   {format(CF.mag(rho_ijk),'5.4f')} km \n")
    # file_out.write(f"Drho_ijk =   {format(drho_ijk[0], '5.4f')} I      {format(drho_ijk[1], '5.4f')} J      {format(drho_ijk[2], '5.4f')} K Mag =     {format(CF.mag(drho_ijk),'5.4f')} km/s \n")
    # file_out.write(f"R_ijk    =   {format(R[0], '5.4f')} I   {format(R[1], '5.4f')} J   {format(R[2], '5.4f')} K Mag =  {format(CF.mag(R),'5.4f')} km \n")
    # file_out.write(f"V_ijk    =   {format(V[0], '5.4f')} I      {format(V[1], '5.4f')} J      {format(V[2], '5.4f')} K Mag =      {format(CF.mag(V),'5.4f')} km/s \n")
    # file_out.write("------------------------------------------------------------------\n\n")
    # file_out.write("             ####     Output Data      ####\n\n")
    # file_out.write("CLASSIC ORBITAL ELEMENTS: \n\n")
    # file_out.write(f"Semimajor Axis                     (km) =     {format(a,'.4f')}\n")
    # file_out.write(f"Eccentricity                            =         {format(e,'.4f')}\n")
    # file_out.write(f"Inclination                       (deg) =        {format(Incl,'.4f')}\n")
    # file_out.write(f"Right Ascension of Ascending Node (deg) =       {format(Raan,'.4f')}\n")
    # file_out.write(f"Argument of Perigee               (deg) =        {format(Argp,'.4f')}\n")
    # file_out.write(f"True Anomaly                      (deg) =        {format(Nu,'.4f')}\n\n\n")

