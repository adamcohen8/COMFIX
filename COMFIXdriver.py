#COMFIX driver
####################################################################################
#
#  Use           : Comfix.out = COMFIXdriver(Comfix.in)
#
#    This program will take observation data of satellites from an input
#    file and converts it into R and V vectors and Classic Orbital Elements.
#    The observation data includes the range, azimuth, elevation, range rate,
#    time derivative of azimuth, and time derivative of elevation.
#
#   Algorithm     : Open input file
#                   Open output file
#                   Enter mainloop
#                       Read in data from input file
#                       Convert to working units
#                       Find GST and LST
#                       Find R site
#                       Find Rho and Drho in SEZ
#                       Convert Rho and Drho to IJK
#                       Find R and V
#                       Find COEs
#                       Write data to output file
#
#
#  Author        : C2C Adam Cohen  USAFA/DFAS  954-980-7732  19 Sep 2022

#
#  Inputs        :
#       Comfix.in  - An input file containing 5 sets of the following Information:
#                       Rho - Magnitude of the range (km)
#                       az  - Azimuth angle (deg)
#                       el  - Elevation angle (deg)
#                       Drho- Magnitude of the range rate (km/s)
#                       daz - Time derivative of Azimuth angle (deg/s)
#                       Del - Time derivative of Elevation angle (deg/s)
#                       Julian Date
#
#
#  Outputs       :
#       Comfix.out - An output file containing 5 sets of the following information:
#                       Lat - latitude (deg and rad)
#                       Lon - Longitude (deg and rad)
#                       Alt - Altitude (km and m)
#                       Rho - Range (Magnitude, SEZ, IJK) (km)
#                       Drho- Range Rate (Magnitude, SEZ, IJK) (km/s)
#                       daz - Time derivative of Azimuth angle (deg/s and rad/s)
#                       Del - Time derivative of Elevation angle (deg/s and rad/s)
#                       Time- HR:MIN:SEC.MILISECOND
#                       Year
#                       GST - Greenwich Sidereal Time (deg and rad)
#                       LST - Local Sidereal Time (deg and rad)
#                       R   - Position Vector (Magnitude and IJK) (km)
#                       V   - Velocity Vector (Magnitude and IJK) (km)
#                       A   - Semimajor Axis (km)
#                       E   - Eccentricity
#                       Incl- Inclination (deg)
#                       Argp- Argument of Perigee (deg)
#                       Raan- Right Ascension of the Ascending Node (deg)
#                       Nu  - True Anomaly (deg)
#                       Julian Date (days)
#
#
#  Locals        :
#
#
#  Constants     :
#       MU    - Mass of the Earth times the gravity constant = 398600.5 (km^3/s^2)
#       Deg   - Conversion factor for degrees from radians = 180/pi (deg/rad)
#       Rad   - Conversion factor for radians from degrees = pi/180 (rad/deg)
#       Zero_IE -Tolerance factor for zero = 0.015
#       Samll - Tolerance factor for zero = 1 * 10^-6
#       Ae    - Semimajor axis of the earth = 6378.137 (km)
#       Ee    - Eccentricity of the earth = 0.0818191908426
#       w     - Anglular velocity of the rotation of the earth = 0.00007292115 K (rad/s)
#
#
#  Coupling      :
#       COMFIXfuncts.py   - A python file containing all of the sub functions
#                           needed for the driver
#       axisrot           - Function that rotates a vector around a given axis
#                           by a given angle
#       dayofyr2mdhms     - Function that takes the year and day of the year and outputs the
#                           month, day, hour, minute, hour, and second
#       elOrb             - Function that takes R and V vectors and outputs COEs
#       finddays          - Function that takes year, month, day, hr, min, sec and outputs dday
#       gstim0            - Function that finds theta g0 for a given year
#       invjulianday      - Function that takes the julian date and outputs Yr, Mon, Day,
#                           Hr, Min, Sec
#       julainday         - Function that takes Yr, Mon, Day, Hr, Min, Sec and outputs Julian Date
#       mag               - Function that takes a 3D vector and outputs its magnitude
#       revcheck          - Function that takes an angle and a modby value and outputs an angle
#                           between 0 and 360 degrees or 0 and 2pi radians
#       vecangle          - Function that takes two vectors and outputs the angle between
#                           the two vectors
#       writeOutput       - Function that takes all output data and writes it to the output file
#       gstlst            - Function that takes julian day and site longitude and outputs
#                           GST and LST
#       site              - Function that takes Site latitude, Site longitude, and LST and
#                           finds R site
#       SEZ2IJK           - Function that rotates vectors from SEZ to IJK
#       OBS2RANGERANGERATE- Function that takes input data and finds Rho and Drho in SEZ
#       ijktorv           - Function that takes Rho, Drho, and R site and outputs R and V
#       numpy             - Python module for advanced math
#       math              - Python module for basic math
#
#
#  References    :
#       Astro Engr 321 Course Handbook COMFIX project description
#
###################################################################################
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


for line in lines:
    #Get next set of Observation Data Here
    obsData = line.split(' ')
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
    CF.writeOutput(file_out, sitlat, sitlon, sitalt, a, e, Incl, Raan, Argp, Nu, rho, az, el, drho, daz, Del, rho_sez,rho_ijk, drho_sez, drho_ijk, jd, satID, GST, LST, R, V, Rsite)


