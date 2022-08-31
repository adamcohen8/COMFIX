import numpy as np
import math
import COMFIXfuncts as CF

R = np.array([8840.0, 646.0, 5455.0])
V = np.array([-0.6950, 5.2500, -1.6500])
print(CF.elorb(R,V))
