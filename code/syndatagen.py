#Synthetic Data Generator
import networkx as nx
import numpy as np
import os
import sys
import math
import random
import EmbedUtilities
import deconvolverun
import itertools



class SynDataTest:

    @staticmethod
    def testSynData(freqmat,comp2domain,comp2scale):
        """tests the synthetic generated data
        Args:
          freqmat:
          comp2domain:
          comp2scale:
        Returns:
          bool: true fale
        """
        maxfreq = np.amax(freqmat)
        maxindices = set((in1,in2) for in1 in xrange(np.shape(freqmat)[0]) for in2 in xrange(np.shape(freqmat)[1]) if freqmat[in1,in2] == maxfreq)
        #each motif can only contribute as its scale
        maxtotcount = 0
        for in1,in2 in maxindices:
            totcount = 0
            for comp in comp2domain.keys():
                count = 0
                for dom in comp2domain[comp]:
                    if in1 >= dom[0] and in2 >= dom[0] and in1 <= dom[1] and in2 <= dom[1]:
                       count += 1
                       totcount += comp2scale[comp]
                assert count <= 1
            if totcount > maxtotcount:
               maxtotcount = totcount
        intersect = lambda (s1,e1), (s2,e2): False if (e1 < s2 or e2 < s1) else True
        #generated domains must not intersect 
        for comp in comp2domain.keys():
            curdoms = list(comp2domain[comp])
            for dom1,dom2 in list(itertools.combinations(curdoms,2)):
                assert not intersect(dom1,dom2)
        assert np.amax(freqmat) == sum(comp2scale.values())
        return True
     
    @staticmethod
    def testReadWrite(freqfile,freqmat,nodecount):
        """
        Args:
          freqfile:
          freqmat:
          nodecount:
        Returns:
          bool:
        """
        freqmat2, allnodes = EmbedUtilities.readFreqFile(freqfile)
        return len(allnodes) == nodecount and (freqmat == freqmat2).all()       

    
def genSynData(nodecount,compcount,avgintlen,intlenstd,avgscale,scalestd,p):
    """generates synthetic 3C data
    Args:
       nodecount: number of nodes
       compcount: component count
       avgintlen: avg interval length
       intlenstd:
       avgscale: avg scale
       scalestd:
       p: non interval proabability
    Returns:
       freqmat: frequency matrix
       comp2domain:
       comp2scale:
    """
    comp2domain = {}
    comp2scale = {}
    freqmat = np.zeros((nodecount,nodecount),dtype=np.int)
    for comp in xrange(compcount):
        curscale = max(1,int(round(np.random.normal(avgscale,scalestd))))
        comp2scale[comp] = curscale
        comp2domain[comp] = set()
        startlen = 0
        while(startlen<nodecount):
           intlen = int(round(np.random.normal(avgintlen,intlenstd)))
           if intlen <= 0:
              while(random.random() < p):
                 startlen += 1
           else:
              intlen = min(intlen,nodecount-startlen) 
              comp2domain[comp].add((startlen+1,startlen+intlen))
              freqmat[startlen:startlen+intlen, startlen:startlen+intlen] += curscale
              startlen += intlen
    return freqmat,comp2domain,comp2scale
    
TESTMODE = True 
def main():
    """
    """
    outfolder = "synout"
    nodecount = 100
    compcount = 20
    avgintlen = 10
    intlenstd = 8
    avgscale = 5
    scalestd = 3
    p = 0.8
    freqmat,comp2domain,comp2scale = genSynData(nodecount,compcount,avgintlen,intlenstd,avgscale,scalestd,p)
    domains = list(set(pair for pairset in comp2domain.values() for pair in pairset))
    freqfile = "freqdeneme.gz"
    domainfile = "mydomain.domain"
    EmbedUtilities.writeFreqData(freqmat,freqfile)
    EmbedUtilities.writeDomainFile(domainfile,domains)
    if TESTMODE:
       assert SynDataTest.testSynData(freqmat,comp2domain,comp2scale) 
       assert SynDataTest.testReadWrite(freqfile,freqmat,nodecount)
     
    method = "fd"        
    objtype = "frobenius"
    deconvolverun.runDeconvolution(method,objtype,freqfile,domainfile,compcount,outfolder)        

    
def geetrun():
    """
    """
    outfolder = "synsol"
    method = "fd"        
    objtype = "frobenius"
    indices = [1,2,3,4]
    for ind in indices:
        synfolder = "synthetic{0}-test".format(ind)
        freqfile = "{0}/countmatrix.gz".format(synfolder)
        domainfile = "{0}/domains.txt".format(synfolder)
        countfile = "{0}/numclasses.txt".format(synfolder)
        sentoutfolder = "{0}-synthetic{1}-test".format(outfolder,ind)
        with open(countfile,"r") as infile:
            for line in infile:
                compcount = int(line.rstrip())
        deconvolverun.runDeconvolution(method,objtype,freqfile,domainfile,compcount,sentoutfolder)        


if  __name__ =='__main__':
    geetrun()
    #main()



    

    

    

    
