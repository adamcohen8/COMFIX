import COMFIXfuncts as CF
import PREDICTfuncts as PF
import numpy as np
import math

#Test J2DragPert
#print(PF.J2DragPert(0.90129675238, 0.00031520000, 0.00112702320, (8.49981069216*10.0**-13.0)))

#Test Newton
#print(4.355546-0.2*np.sin(4.355546))
#print(PF.newton(4.54294685891266,0.2))

#Test coeupdate
#print(PF.coeupdate(94962.49143630, 0.00112702320, (8.49981069216*10.0**-13.0), 0.00031520000, (-1.00525958512*10.0**(-9.0)), 2.58442072450, (-1.00068980006*10.0**(-6.0)), 4.27959095912, (7.46318463048*10.0**(-7.0)), 3.97538705182))


#Test COEs2rPQW
#print(PF.coes2rpqw())

print(CF.Sun(2459858.40089709))