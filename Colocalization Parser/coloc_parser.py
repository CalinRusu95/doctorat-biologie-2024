import cv2 as cv
import numpy as np
import math
from numpy.core.umath_tests import inner1d

R, G, B = 0, 1, 2
channelPairs = [{R,G}, (R,B), (G,B)]

def safedot(a, b):
    return(np.multiply(a, b, dtype = np.uint16).sum(dtype = np.float64))
    
def coloc(ia):
    sumSqs = inner1d(ia.T, ia.T).astype(np.float64)
    sums = ia.sum(axis = 0, dtype = np.float64)
    means = sums / ia.shape[0]
    meanErrors = ia - means
    sqMeanErrors = meanErrors**2
    sumSqMeanErrors = sqMeanErrors.sum(axis = 0)
    del sqMeanErrors
    
    indicator = ia > 0
    
    crossDot = {(c1, c2) : safedot(ia[:,cl], ia[:,c2]) for c1, c2 in channelPairs}
    sumIf = {(c1, c2) : ia[:,cl][indicator[:,c2]].sum() for c1, c2 in channelPairs}
    
    results = {}
    
    for c1, c2 in channelPairs:
        
        k1 = crossDot[(c1, c2)] / sumSqs[c1]
        k2 = crossDot[(c1, c2)] / sumSqs[c2]
        
        results[(c1, c2)] = {
                "Pearson" : (np.dot(meanErrors[:,c1], meanErrors[:,c2]) / np.sqrt(sumSqMeanErrors[c1] * sumSqMeanErrors[c2])), 
                "Manders" : math.sqrt(k1*k2),                
                "Coloc(m)1" : sumIf[(c1, c2)] / sum[c1],
                "Coloc(m)2" : sumIf[(c2,c1)] / sums[c2],               
                "Overlap(k)1" : k1,
                "Overlap(k)2" : k2 }
    
    return results