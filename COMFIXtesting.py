import numpy as np
import math
import COMFIXfuncts as CF

# R = np.array([8840.0, 646.0, 5455.0])
# V = np.array([-0.6950, 5.2500, -1.6500])
# print(CF.elOrb(R,V))

#print(CF.gstim0(2022)*180/math.pi)

#print(CF.julianday(2022,8,8,19,45,10.0))
#print(CF.invjulianday(2459800.3230324076))

#print(CF.dayofyr2mdhms(2022, 220))

#print(CF.Sun(2459823.0788194444))

#print(CF.axisrot(CF.axisrot(np.array([-1636.4022,562.2001,1227.4091]),2,-0.214675498),3,-3.2302623391))

#print(CF.site(1.3561, 0.049, 3.2302623391))

#print(CF.SEZ2IJK(np.array([-1636.4022,562.2001,1227.4091]), 1.3561, 3.2302623391))


#print(CF.ijktorv(np.array([1381.8950, -441.5631, 1547.8381]), np.array([-7.0721, 3.8723, 2.8678]), np.array([-1357.7548, -120.7082, 6209.9350])))


# print(CF.site(1.356120828799594, 0.049, 3.2302623381265736))
# print(CF.site(1.356120828799594, 0.049, 3.2302623391265736))

print(CF.finddays(2022, 8,8,19,45,10.0))