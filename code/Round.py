#Rounding related methods
import networkx as nx
import numpy as np
import scipy as sp
import math
import random
import operator
import os
import sys
import itertools
import EmbedUtilities
from copy import deepcopy


def roundSimple():
    """
    """
    maxval = 0                
    comp2domset = {comp:set() for comp in xrange(compcount)}
    for (domin,comp) in xdict.keys():
        if xdict[(domin,comp)] > 0.02:
           comp2domset[comp].add(domin)
        if xdict[(domin,comp)] > maxval:
           maxval = xdict[(domin,comp)]    
    print comp2domset
    print "here"
    newcomp2domset = {comp:set() for comp in xrange(compcount)}
    for comp in comp2domset.keys():
        for dom in comp2domset[comp]:
            if xdict[(dom,comp)] >= 0.3:
               newcomp2domset[comp].add(dom)
    print newcomp2domset
    print maxval      
    comp2dominds = {}
    for comp in comp2domset.keys():
        domsdict= {domin: xdict[(domin,comp)] for domin in comp2domset[comp]}
        sorteddict = sorted(domsdict.iteritems(), key=operator.itemgetter(1),reverse=True)
        curdomins = set()
        for domin,val in sorteddict:
            intflag = False
            for curdomin in curdomins:
                if intersect(domains[curdomin],domains[domin]):
                   intflag = True
                   break
            if not intflag and random.random() <= val:
               curdomins.add(domin)
        comp2dominds[comp] = set(curdomins)

    modxdict = {(domin,comp,1): xdict[(domin,comp)] for domin,comp in xdict.keys()}
    node2dom = EmbedUtilities.getnode2dom(freqmat,domains)
    fracobj = EmbedUtilities.estFracObjective(modxdict,freqmat,node2dom,set([1]),compcount)
    newxdict = {(domin,comp,1): 0.0 for domin in xrange(len(domains)) for comp in xrange(compcount)}
    for comp in comp2dominds:
        for domin in comp2dominds[comp]:
            newxdict[(domin,comp,1)] = 1.0
    fracobj2 = EmbedUtilities.estFracObjective(newxdict,freqmat,node2dom,set([1]),compcount)         
    print "pre rounding ",fracobj
    print "after rounding ",fracobj2
    print "fsum: ",sum([freqmat[ind1,ind2]**2 for ind1 in xrange(np.shape(freqmat)[0]) for ind2 in xrange(np.shape(freqmat)[1])])
    return comp2dominds   


def roundCR(xdict,sideparams):
    """round based on contention resolution 1, 1/e
    Args:
       xdict:
       sideparams:
    Returns:
       comp2dominds:
    """
    intersect = lambda (s1,e1), (s2,e2): False if (e1 < s2 or e2 < s1) else True
    [compcount,freqmat,interdom,domains] = sideparams
    comp2R = {comp:set() for comp in xrange(compcount)}
    for (domin,comp) in xdict.keys():
        if random.random() < 1.0-math.exp(-1.0*xdict[(domin,comp)]):   
           comp2R[comp].add(domin)
    comp2mark = {comp:set() for comp in xrange(compcount)}    
    for (domin,comp) in xdict.keys():
        for domin2 in comp2R[comp]:
            if domin2 == domin:
               continue 
            if domains[domin2][0] <= domains[domin][0] and domains[domin2][1] >= domains[domin][0]:
               comp2mark[comp].add(domin)
               break
    for comp in comp2mark.keys():
        comp2R[comp] -= comp2mark[comp]
    for comp in comp2R.keys():
        alldoms = list(comp2R[comp])
        for dom in alldoms:
            assert dom not in comp2mark[comp]   
    return comp2R

   
class Round():

    @staticmethod
    def greedy3Rand(xdict,sideparams):
        """greedy rounding iterative loop
        Args:
           xdict:
           sideparams:
        Returns:
           comp2dominds: domain indices of each component
           comp2scale: component to scale map
        """
        intersect = lambda (s1,e1), (s2,e2): False if (e1 < s2 or e2 < s1) else True
        comp2dominds = {comp: set() for comp in xrange(compcount)}
        comp2scale = {comp: None for comp in xrange(compcount)}
        scales,compcount,freqmat,interdom,domains = sideparams
        node2dom = EmbedUtilities.getnode2dom(freqmat,domains)
        sorted_x = sorted(xdict.iteritems(), key=operator.itemgetter(1),reverse=True)
        cursol = {(domin,comp,scale): 0 for domin in xrange(len(domains)) for comp in xrange(compcount) for scale in scales}
        curobjval = EmbedUtilities.estFracObjective(cursol,freqmat,node2dom,scales,compcount)
        for (domin,comp,scale),val in sorted_x:
            cursol[(domin,comp,scale)] = 1
            if violateCons(cursol,compcount,domains):
               cursol[(domin,comp,scale)] = 0
               continue
            estobjval = EmbedUtilities.estFracObjective(cursol,freqmat,node2dom,scales,compcount)
            if estobjval <= curobjval:
               curobjval = estobjval
            else:
               cursol[(domin,comp,scale)] = 0
        assert violateCons(cursol,compcount,domains)
        for (domin,comp,scale) in cursol.keys():
            comp2dominds[comp].add(domin)
            comp2scale[comp].add(scale)
        return comp2dominds,comp2scale
      

    @staticmethod
    def greedy2Rand(xdict,sideparams):
        """greey rounding 2
        Args:
           xdict:
           sideparams:
        Returns:
           comp2dominds: domain indices of each component
           comp2scale: component to scale map
        """
        scales,compcount,freqmat,interdom,domains = sideparams
        comp2dominds = {comp: set() for comp in xrange(compcount)}
        comp2scale = {comp: None for comp in xrange(compcount)}
        #divxdict = {(domin,comp,scale): xdict[(domin,comp,scale)]/ydict[(comp,scale)] for domin,comp,scale in xdict.keys()}
        divxdict = {(domin,comp,scale): xdict[(domin,comp,scale)] for domin,comp,scale in xdict.keys()}
        intersect = lambda (s1,e1), (s2,e2): False if (e1 < s2 or e2 < s1) else True
         
        node2dom = EmbedUtilities.getnode2dom(freqmat,domains)
        soldict = {(domin,comp,scale): 0.0 for (domin,comp,scale) in xdict.keys()}  
        curobjval = EmbedUtilities.estFracObjective(soldict,freqmat,node2dom,scales,compcount)
        allcomps = range(compcount)
        random.shuffle(allcomps)
        for comp in allcomps:
            print "before {0} {1}".format(comp,curobjval)
            bestscaleobj = curobjval
            bestscaledict = {}
            bestscale = None
            bestdominds = None
            for scale in scales:
                scaledict = deepcopy(soldict)
                
                x = {domin: divxdict[(domin,comp,scale)] for domin in xrange(len(domains))}
                curdoms = set()
                remdomset = set(range(len(domains)))
                imprflag = True
                while imprflag and len(remdomset) > 0:
                   imprflag = False
                   remdoms = list(remdomset) 
                   random.shuffle(remdoms)
                   p = random.random()
                   found, prob = remdoms[0], 0.0
                   divsum = sum(x.values())
                   #if divsum == 0:
                   #   break 
                   normprob = {domin: x[domin]/divsum for domin in remdomset}
                   for dom in remdoms:
                       prob += normprob[dom]
                       if prob >= p:
                          found = dom
                          imprflag = True
                          break
                   if imprflag:
                      delset = set(domin for domin in remdomset if intersect(domains[domin],domains[found]))
                      remdomset -= delset
                      for deldom in delset:
                          del x[deldom]
                      curdoms.add(found)
            
                #x = {domin: divxdict[(domin,comp,scale)] for domin in xrange(len(domains))}     
                #sorted_x = sorted(x.iteritems(), key=operator.itemgetter(1),reverse=True)
                #curdoms = set()
                #while len(sorted_x) > 0:
                #   found,foundval = sorted_x[0]
                #   intflag = False
                #   for curdom in curdoms:
                #       if intersect(domains[curdom],domains[found]):
                #          intflag = True
                #          break
                #   del sorted_x[0]   
                #   if not intflag:
                #      curdoms.add(found)

                      
                for curdomin in curdoms:
                    scaledict[(curdomin,comp,scale)] = 1.0
                scaleobjval = EmbedUtilities.estFracObjective(scaledict,freqmat,node2dom,scales,compcount)     
                if scaleobjval <= bestscaleobj:
                   bestscaleobj = scaleobjval
                   bestscaledict = deepcopy(scaledict)
                   bestscale = scale
                   bestdominds = set(curdoms)
            if bestscale != None:
               soldict = deepcopy(bestscaledict)
               curobjval = bestscaleobj
               comp2scale[comp] = bestscale
               comp2dominds[comp] = set(bestdominds)
            else:
               break
        print "lastsol ",EmbedUtilities.estFracObjective(soldict,freqmat,node2dom,scales,compcount)  
        emdict = {(domin,comp,scale): 0.0 for (domin,comp,scale) in xdict.keys()}  
        print "emptysol ",EmbedUtilities.estFracObjective(emdict,freqmat,node2dom,scales,compcount)

        print np.shape(freqmat)[0]*np.shape(freqmat)[1]
        print math.sqrt(EmbedUtilities.estFracObjective(emdict,freqmat,node2dom,scales,compcount) / float(np.shape(freqmat)[0]*np.shape(freqmat)[1]))
        print "max: ", np.amax(freqmat)
        print "min: ", np.amin(freqmat)
        print "info"
        for comp in comp2scale.keys():
            print comp,comp2scale[comp]
        print allcomps    
        exit(1)
    
        for comp in xrange(compcount):
            for domin,comp2,scale in newxdict.keys():
                if comp2 == comp:
                   print domin,comp,scale,newxdict[(domin,comp,scale)] 
            break
    
    @staticmethod
    def roundRand(xdict,ydict,sideparams):
        """randomized rounding
        Args:
           xdict:
           ydict:
           sideparams:
        Returns:
           comp2dominds: domain indices of each component
           comp2scale: component to scale map
        """
        intersect = lambda (s1,e1), (s2,e2): False if (e1 < s2 or e2 < s1) else True
        compcount,scales,domains = sideparams
        comp2dominds = {comp: set() for comp in xrange(compcount)}
        comp2scale = {comp: None for comp in xrange(compcount)}
        #comp scale assignment by randomized round
        for comp in xrange(compcount):
            p = random.random()
            tot = 0.0
            curscales = list(scales)
            random.shuffle(curscales)
            selscale = curscales[0]
            for scale in curscales:
                tot += ydict[(comp,scale)]
                if tot >= p:
                   selscale = scale 
                   break 
            comp2scale[comp] = selscale
    
        domcom2val = {(domin,comp): sum([xdict[(dom,tcomp,scale)] for dom,tcomp,scale in xdict.keys() if tcomp==comp and dom==domin]+[0.0]) for domin in xrange(len(domains)) for comp in xrange(compcount)}
        com2dom2val = {comp:{} for comp in xrange(compcount)}
        for domin,comp in domcom2val.keys():
            com2dom2val[comp][domin] = domcom2val[(domin,comp)]
        #comp domain assignment by modified randomized rounding
        for comp in xrange(compcount):
            domset = set()
            remdomset = set(com2dom2val[comp].keys())
            imprflag = True
            while imprflag and len(remdomset) > 0:
               imprflag = False
               curdoms = list(remdomset) 
               random.shuffle(curdoms)
               p = random.random()
               found, prob = curdoms[0], 0.0
               divsum = float(sum(com2dom2val[comp][remdom] for remdom in remdomset))
               normprob = {dom: com2dom2val[comp][dom]/divsum for dom in remdomset}
               for dom in curdoms:
                   prob += normprob[dom]
                   if prob >= p:
                      found = dom
                      imprflag = True
                      break
               if imprflag:
                  delset = set(domin for domin in remdomset if intersect(domains[domin],domains[found]))
                  remdomset -= delset
                  domset.add(found)
            comp2dominds[comp] = set(domset)
        return comp2dominds, comp2scale

    @staticmethod 
    def greedyRand(xdict,ydict,sideparams):
        """randomized rounding
        Args:
           xdict:
           ydict:
           sideparams:
        Returns:
           compdict: domains of each component
           comp2s: component to scale map
        """
        intersect = lambda (s1,e1), (s2,e2): False if (e1 < s2 or e2 < s1) else True
        compcount,scales,domains = sideparams
        comp2dominds = {comp: set() for comp in xrange(compcount)}
        comp2scale = {comp: None for comp in xrange(compcount)}
    
        #greedy comp scale assignment
        #maxvals = {comp: 0.0 for comp in xrange(compcount)}
        #for comp,scale in ydict.keys():
        #    if ydict[(comp,scale)] >= maxvals[comp]:
        #       maxvals[comp] = ydict[(comp,scale)]
        #       comp2scale[comp] = scale
                                 
        domcom2val = {(domin,comp): sum([xdict[(dom,tcomp,scale)] for dom,tcomp,scale in xdict.keys() if tcomp==comp and dom==domin]+[0.0]) for domin in xrange(len(domains)) for comp in xrange(compcount)}
        com2dom2val = {comp:{} for comp in xrange(compcount)}
        for domin,comp in domcom2val.keys():
            com2dom2val[comp][domin] = domcom2val[(domin,comp)]
        #comp domain assignment by modified randomized rounding
        for comp in xrange(compcount):
            domset = set()
            remdomset = set(com2dom2val[comp].keys())
            imprflag = True
            while imprflag and len(remdomset) > 0:
               imprflag = False
               curdoms = list(remdomset) 
               random.shuffle(curdoms)
               p = random.random()
               found, prob = curdoms[0], 0.0
               divsum = float(sum(com2dom2val[comp][remdom] for remdom in remdomset))
               normprob = {dom: com2dom2val[comp][dom]/divsum for dom in remdomset}
               for dom in curdoms:
                   prob += normprob[dom]
                   if prob >= p:
                     found = dom
                     imprflag = True
                     break
               if imprflag:
                  delset = set(domin for domin in remdomset if intersect(domains[domin],domains[found]))
                  remdomset -= delset
                  domset.add(found)
            comp2dominds[comp] = set(domset)
        return comp2dominds, comp2scale

