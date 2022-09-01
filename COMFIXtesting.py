import numpy as np
import math
import COMFIXfuncts as CF

# R = np.array([8840.0, 646.0, 5455.0])
# V = np.array([-0.6950, 5.2500, -1.6500])
# print(CF.elorb(R,V))

#print(CF.gstim0(2022)*180/math.pi)

# print(CF.julianday(2022,8,31,13,53,30))
# print(CF.invjulianday(2459823.0788194444))

#print(CF.Sun(2459823.0788194444))

print(CF.axisrot(CF.axisrot(np.array([-1636.4022,562.2001,1227.4091]),2,-0.214675498),3,-3.2302623391))

#print(CF.site(1.3561, 0.049, 3.2302623391))

print(CF.SEZ2IJK(np.array([-1636.4022,562.2001,1227.4091]), 1.3561, 3.2302623391))
