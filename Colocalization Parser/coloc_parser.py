import cv2 as cv
import numpy as np
import math

# Define color channel indices
R, G, B = 0, 1, 2
channelPairs = [(R, G), (R, B), (G, B)]

def safedot(a, b):
    """Compute the dot product safely using high precision."""
    return np.multiply(a, b, dtype=np.uint16).sum(dtype=np.float64)

def coloc(ia):
    """
    Compute various colocalization metrics for image array `ia`.
    
    Args:
        ia (np.array): 2D numpy array where each column represents a color channel.
    
    Returns:
        dict: Colocalization metrics for each pair of color channels.
    """
    sumSqs = computeSumSquares(ia)
    sums, means = computeSumsAndMeans(ia)
    meanErrors = ia - means
    sumSqMeanErrors = computeSumSquaredMeanErrors(meanErrors)
    
    crossDot = computeCrossDots(ia, channelPairs)
    sumIf = computeSumIf(ia, channelPairs, meanErrors > 0)
    
    return computeResults(channelPairs, meanErrors, sumSqMeanErrors, sumSqs, sums, crossDot, sumIf)

def computeSumSquares(ia):
    """Compute sum of squares for each channel."""
    return np.einsum('ij,ij->j', ia, ia).astype(np.float64)

def computeSumsAndMeans(ia):
    """Compute sums and means for each channel."""
    sums = ia.sum(axis=0, dtype=np.float64)
    means = sums / ia.shape[0]
    return sums, means

def computeSumSquaredMeanErrors(meanErrors):
    """Compute sum of squared mean errors for each channel."""
    return np.einsum('ij,ij->j', meanErrors, meanErrors)

def computeCrossDots(ia, channelPairs):
    """Compute cross dot products for each channel pair."""
    return {(c1, c2): safedot(ia[:, c1], ia[:, c2]) for c1, c2 in channelPairs}

def computeSumIf(ia, channelPairs, indicator):
    """Compute sum of intensities where both channels are positive."""
    return {(c1, c2): ia[:, c1][indicator[:, c2]].sum() for c1, c2 in channelPairs}

def computeResults(channelPairs, meanErrors, sumSqMeanErrors, sumSqs, sums, crossDot, sumIf):
    """Compute and return the colocalization results."""
    results = {}

    for c1, c2 in channelPairs:
        k1 = crossDot[(c1, c2)] / sumSqs[c1]
        k2 = crossDot[(c1, c2)] / sumSqs[c2]

        results[(c1, c2)] = {
            "Pearson": np.dot(meanErrors[:, c1], meanErrors[:, c2]) / np.sqrt(sumSqMeanErrors[c1] * sumSqMeanErrors[c2]),
            "Manders": math.sqrt(k1 * k2),
            "Coloc(m)1": sumIf[(c1, c2)] / sums[c1],
            "Coloc(m)2": sumIf[(c2, c1)] / sums[c2],
            "Overlap(k)1": k1,
            "Overlap(k)2": k2
        }

    return results
