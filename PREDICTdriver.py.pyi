import PREDICTfuncts as PF
import COMFIXfuncts as CF
import math
import numpy as np


#Get user input here


#Open TLE File


#Open Output File


#For Each TLE

    #Get next TLE

    inc0 = 0
    ecc0 = 0
    n0 = 0
    ndot2 = 0

    #Convert to working units


    #J2DragPert
    [raandot, argpdot, eccdot] = PF.J2DragPert(inc0, ecc0, n0, ndot2)

