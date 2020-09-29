#Local Search for Frequency Deconvolution 
import networkx as nx
import numpy as np
import scipy as sp
import math
import random
import os
import sys
import itertools
from copy import deepcopy
import EmbedUtilities

modsign = lambda mod: -1 if mod == "add" else 1

def findAddDelItem(sentitemset,cursol,objval,sideparams,mod):
    """find improving addition item 
    Args:
       sentitemset: items with 0 for add, with 1 for del
       cursol:
       objval:
       sideparams:
       mod:
    Returns:
       founditem:   
    """
    store = deepcopy(cursol)
    itemset = list(sentitemset)
    freqmat,domains,scales,compcount,interG,dom2index = sideparams
    EPSILON = 0.1         
    random.shuffle(itemset)
    founditem = None
    for (domin,comp,scale) in itemset:
        if mod == "add":
           assert not cursol.has_key(domin) or not cursol[domin].has_key(comp) or not cursol[domin][comp].has_key(scale)
        elif mod == "del":
           assert cursol.has_key(domin) and cursol[domin].has_key(comp) and cursol[domin][comp].has_key(scale)
        if not cursol.has_key(domin):
           cursol[domin] = {}
        if not cursol[domin].has_key(comp):
           cursol[domin][comp] = {}         
        cursol[domin][comp][scale] = (-1*modsign(mod)+1)/2
        if violateCons(cursol,compcount,domains):
           continue 
        curobjval = solModSingle(cursol,(domin,comp,scale),mod,sideparams)
        if objval - curobjval > EPSILON:      
           founditem = (domin,comp,scale)
           break
        #del cursol[domin][comp][scale] = (modsign(mod)+1)/2
    for key in store.keys():
        if key != founditem:
           assert store[key] == cursol[key]
        else:
           assert store[key] != cursol[key]     
    return founditem,objval

def formatchange(cursol,domains,compcount,scales):
    """
    """
    return {(domin,comp,scale): cursol[domin][comp][scale] for domin in cursol.keys() for comp in cursol[domin].keys() for scale in cursol[domin][comp].keys()}



def findInitSol(domains,interG,dom2index,compcount,scales):
    """
    Args:
       domains:
       interG:
       dom2index:
       compcount:
       scales:
    Returns:
       initsol:
    """
    intersect = lambda (s1,e1), (s2,e2): False if (e1 < s2 or e2 < s1) else True
    print interG.number_of_nodes()
    #manode, maxitem = 0,0
    #for node in interG.nodes():
    #    if 
    #    max([len(interG.neighbors(node)) for node in interG.nodes()])
    exit(1) 
    initsol = {}
    compinterG = nx.complement(interG)
    cliques = EmbedUtilities.findCliqueDecomp(compinterG)
    print len(cliques)
    for clique in cliques:
        for dom1 in clique:
            for dom2 in clique:
                if dom1 != dom2:
                   assert not intersect(dom1,dom2)
    print "clique info"               
    for clique in cliques:
        print len(clique)                
    for compin in xrange(compcount):
        clique = cliques[compin]
        for dom in clique:
            domin = dom2index[dom]
            if not initsol.has_key(domin):
               initsol[domin] = {}
            if not initsol[domin].has_key(compin):
               initsol[domin][compin] = {}   
            scale = random.choice(list(scales))
            initsol[domin][compin][scale] = 1
    print initsol
    exit(1)
    print violateCons(initsol,compcount,domains)
    exit(1)        
    return initsol


def findSwap():
    """finds swap
    Args:
    Returns:
    """
    return


def findScaleMod(itemset,cursol,comp2scale,objval,sideparams):
    """finds scale modification
    Args:
       itemset:
       cursol:
       comp2scale:
       objval:
       sideparams:
    Returns:  
       cursol:
    """
    for comp in comp2scale.keys():
        for comp in comp2scale.keys():
            pass         
    return    


    

def runLocalSearch(sideparams,node2dom):
    """runs local search for nonmonotone submodular maximization under exchange system
    Args:
       sideparams:
       node2dom:
    Returns:
       comp2dominds:
       comp2scale:
    """
    [freqmat,domains,scales,compcount,interG,dom2index] = sideparams
    comp2dominds = {comp: set() for comp in xrange(compcount)}
    comp2scale = {comp: None for comp in xrange(compcount)}
    basesum = sum([freqmat[ind1,ind2]**2 for ind1 in xrange(np.shape(freqmat)[0]) for ind2 in xrange(np.shape(freqmat)[1])])
    impflag = True
    zeroset = [(domin,comp,scale) for domin in xrange(len(domains)) for comp in xrange(compcount) for scale in scales]  #items that can be added
    oneset = [] #items to be deleted
    cursol = {}
    #initsol = findInitSol(domains,interG,dom2index,compcount,scales)
    #cursol = {domin: {comp: {scale: 0 for scale in scales} for comp in xrange(compcount)} for domin in xrange(len(domains))}
    #cursol2 = formatchange(cursol,domains,compcount,scales)
    #objval = EmbedUtilities.estFracObjective(cursol,freqmat,node2dom,scales,compcount)
    objval = basesum
    while impflag:
        impflag = Falseemre
        opers = ["add","del","scalemod"]
        random.shuffle(opers)
        for oper in opers:
            if oper in ["add","del"]:
               founditem,objval = findAddDelItem((lambda mode: zeroset if mode == "add" else oneset)(oper),cursol,objval,sideparams,oper) #cursol modified inside
            else:
               founditem,objval = findScaleMod((lambda mode: zeroset if mode == "add" else oneset)(oper),cursol,objval,sideparams,oper) #cursol modified inside    
            print founditem
            if founditem != None:
               impflag = True
               if oper == "add":
                  zeroset.remove(founditem)
                  oneset.append(founditem)
               elif oper == "del":
                  oneset.remove(founditem)
                  zeroset.append(founditem)
               break
    print "finished"
    exit(1)       
    assert len(set(zeroset).intersection(oneset)) == 0 and len(zeroset) + len(oneset) == len(scales)*compcount*len(domains)
    return cursol


def solModSingle(cursol,item,mod,sideparams):
    """estimates the solution difference by item addition/deletion
    Args:
       cursol: current solution dictionary dict of dict of dict
       item:
       mod: add or delete
       sideparams:
    Returns:
       objdif: objective value difference
    """
    modsign = lambda mod: -1 if mod == "add" else 1
    freqmat,domains,scales,compcount,interG,dom2index = sideparams
    objdif = 0.0
    itdomin,itcomp,itscale = item
    assert cursol[itdomin][itcomp][itscale] == (modsign(mod)+1)/2
    itdom = domains[itdomin]
    for dom2 in interG[itdom]:
        dom2in = dom2index[dom2]
        interlen = min(itdom[1],dom2[1])-max(itdom[0],dom2[0]) + 1
        scalesum = sum([scale for comp in cursol[dom2in].keys() for scale in cursol[dom2in][comp].keys() if cursol[dom2in][comp][scale] == 1])
        objdif -= modsign(mod)*2*(interlen**2)*scalesum*itscale
    qcoef = (itdom[1]-itdom[0]+1)**2
    fsum = np.sum(freqmat[itdom[0]:itdom[1]+1,itdom[0]:itdom[1]+1])
    scalesum = sum([scale for comp in cursol[itdomin].keys() for scale in cursol[itdomin][comp].keys() if cursol[itdomin][comp][scale] == 1])
    if mod == "add":
       objdif -= (modsign(mod)*qcoef*itscale*((2*scalesum)+itscale))
    elif mod == "del":
       objdif -= (modsign(mod)*qcoef*itscale*((2*scalesum)-itscale))  
    objdif += (modsign(mod)*2.0*itscale*fsum)
    return objdif

    
def violateCons(cursol,compcount,domains):
    """returns if there is any violating constraints
    Args:
       cursol:
       compcount:
       domains:
    Returns:
       bool: 
    """
    intersect = lambda (s1,e1), (s2,e2): False if (e1 < s2 or e2 < s1) else True
    comp2scale = {comp: set() for comp in xrange(compcount)}
    comp2doms = {comp: set() for comp in xrange(compcount)}
    for domin in cursol.keys():
        for comp in cursol[domin].keys():
            if len(cursol[domin][comp]) > 0:
               for existdomin in comp2doms[comp]:
                   if intersect(domains[existdomin],domains[domin]):
                      print "inter",comp,domains[existdomin],domains[domin] 
                      return False       
               comp2doms[comp].add(domin)
               for scale in cursol[domin][comp].keys():
                   comp2scale[comp].add(scale)
                   if len(comp2scale[comp]) > 1:
                      print "scale mis", comp2scale[comp]
                      return False
    return True
         
            
TESTMODE = True
def runLocal(freqmat,domains,scales,compcount):
    """runs local search for sub maximization under matroid constraints
    Args:
       freqmat: frequency matrix
       domains: list of domains(list of set of nodes)
       scales: set of scales
       compcount: number of components
    Returns:
    """
    interdom = EmbedUtilities.getInterDomain(domains)
    interG = EmbedUtilities.inter2Graph(interdom,domains)
    node2dom = EmbedUtilities.getnode2dom(freqmat,domains)
    dom2index = {domains[index]:index for index in xrange(len(domains))}
    sideparams = [freqmat,domains,scales,compcount,interG,dom2index]
    comp2dominds,comp2scale = runLocalSearch(sideparams,node2dom)
    exit(1)

    if TESTMODE:
       TestLocal.testAddDel(sideparams)
       #TestLocal.testMethod()
    
    return comp2dominds,comp2scale

class TestLocal():

    @staticmethod
    def testAddOns(sideparams):
        """tests fast addition/deletion
        Args:
           sideparams:
        Returns:
           bool: true or false
        """
        return
    
    @staticmethod
    def testAddDel(sideparams):
        """tests fast addition/deletion
        Args:
           sideparams:
        Returns:
           bool: true or false
        """
        modsign = lambda mod: -1 if mod == "add" else 1
        freqmat,domains,scales,compcount,interG,dom2index = sideparams
        node2dom = EmbedUtilities.getnode2dom(freqmat,domains)
        cursol = {domin: {comp: {scale: 0 for scale in scales} for comp in xrange(compcount)} for domin in xrange(len(domains))}
        objval = EmbedUtilities.estFracObjective(cursol,freqmat,node2dom,scales,compcount)
        assert objval == sum([freqmat[in1,in2]**2 for in1 in xrange(np.shape(freqmat)[0]) for in2 in xrange(np.shape(freqmat)[1])])
        existitems = []
        upperlimit = len(domains)*len(scales)*compcount
        count = 500
        for index in xrange(count):
            if random.random() < 0.5:
               mode = "add"
            else:
               mode = "del"
            if len(existitems) == 0:
               mode = "add"
            if len(existitems) == upperlimit:
               mode = "del"
            if index <= 10:
               mode = "add"
            randdomin, randcomp, randscale = None, None, None
            if mode == "add":         
               while True:
                  randdomin = random.choice(range(len(domains)))
                  randcomp = random.choice(range(compcount)) 
                  randscale = random.choice(list(scales))
                  if cursol[randdomin][randcomp][randscale] == 0:
                     break
               assert (randdomin,randcomp,randscale) not in existitems  
               existitems.append((randdomin,randcomp,randscale))
            elif mode == "del":
               randdomin, randcomp, randscale = random.choice(existitems)
               existitems.remove((randdomin,randcomp,randscale))
            objdif = solModSingle(cursol,(randdomin,randcomp,randscale),mode,sideparams)
            objval += objdif
            cursol[randdomin][randcomp][randscale] = (-1*modsign(mode)+1)/2
            sentsol = {(domin,comp,scale): cursol[domin][comp][scale] for domin in cursol.keys() for comp in cursol[domin].keys() for scale in cursol[domin][comp].keys()}
            newobj = EmbedUtilities.estFracObjective(sentsol,freqmat,node2dom,scales,compcount)
            assert newobj == objval
        return True 
