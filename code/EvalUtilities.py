#Evaluation utilities
import os
import sys
import math
import numpy as np
from munkres import Munkres

def VIScore():
    """
    """
    return


def matchScore(truecomp2domins,truecomp2scale,estcomp2domins,estcomp2scale,domains,freqmat):
    """bipartite matching based score
    Args:
       truecomp2domins:
       truecomp2scale:
       estcomp2domins:
       estcomp2scale:
       domains:
       freqmat:
    Returns:
       score:
    """
    compcount = len(truecomp2doins.keys()) 
    scoremat = np.zeros((compcount,compcount),dtype=np.float)    
    for truecomp in xrange(compcount):
        truedomins = truecomp2domins[truecomp]
        truemat = np.zeros(np.shape(freqmat),dtype=np.float)
        for domin in truedomins:
            start,end = domains[domin]
            truemat[start:end+1,start:end+1] += truecomp2scale[truecomp]
        for estcomp in xrange(compcount):
            estdomins = estcomp2domins[estcomp]
            estmat = np.zeros(np.shape(freqmat),dtype=np.float)
            for domin in estdomins:
                start,end = domains[domin]
                estmat[start:end+1,start:end+1] += estcomp2scale[estcomp]
            toterr = np.sum(np.fabs(truemat-estmat))
            scoremat[truecomp,estcomp] = toterr
    m = Munkres()
    indexes = m.compute(scoremat)
    score = sum([scoremat[row,column] for row, column in indexes])
    return score
