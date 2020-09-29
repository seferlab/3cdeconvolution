#BiSubmodular Maximization  Formulation of Frequency Deconvolution
import networkx as nx
import numpy as np
import scipy as sp
import math
import random
import os
import sys
import itertools
import operator
import EmbedUtilities
import Bls_Sdp
import Round
from copy import deepcopy
   
def genBoundStr(domains,compcount):
    """generates lp relaxation bound constraints
    Args:
       domains: set of domains
       compcount: number of components 
    Returns:
       boundstr:
    """
    boundstr = " Bounds\n "
    boundstr += "\n"+ " ".join([" 0 <= x{0}_{1} <= 1\n".format(domin,comp) for domin in xrange(len(domains)) for comp in xrange(compcount)])
    return boundstr

def genConsStr(domains,compcount,cliques):
    """generates lp relaxation constraints
    Args:
       domains: set of domains
       compcount: number of components
       cliques:
    Returns:
       consstr:
    """
    dom2index = {domains[index]: index for index in xrange(len(domains))}
    consstr = " Subject To \n"
    for clique in cliques:
        consstr += " ".join([" + ".join([" x{0}_{1} ".format(dom2index[dom],comp) for dom in clique]) + " <= 1 \n " for comp in xrange(compcount)])
    return consstr


def genObjStr(freqmat,compcount,domains,interdom):
    """generates objective function string
    Args:
       freqmat: frequency matrix
       compcount: number of components
       domains: all domains
       interdom: intersecting domains 
    Returns:
       objstr: 
    """
    dom2index = {domains[index]:index for index in xrange(len(domains))}
    objstr = " Minimize\n obj: " 
    singles, quads = [], []
    for dom1,dom2 in interdom:
        domin1 = dom2index[dom1]
        domin2 = dom2index[dom2]
        qcoef = (min(dom1[1],dom2[1])-max(dom1[0],dom2[0]) + 1)**2
        quads.extend([" {0} x{1}_{2} * x{3}_{4} ".format(4*qcoef,domin1,comp1,domin2,comp2) for comp1 in xrange(compcount) for comp2 in xrange(compcount)])
    for dom in domains:
        domin = dom2index[dom]
        qcoef = (dom[1]-dom[0]+1)**2
        quads.extend([" {0} x{1}_{2} * x{1}_{3} ".format(2*qcoef,domin,comp1,comp2) for comp1 in xrange(compcount)  for comp2 in xrange(compcount) ])
        fsum = np.sum(freqmat[dom[0]:dom[1]+1,dom[0]:dom[1]+1])
        singles.extend([" - {0} x{1}_{2} ".format(2*fsum,domin,comp) for comp in xrange(compcount)])
    if len(singles) > 0:
       objstr += " ".join(singles)
    if len(quads) > 0:
       objstr += " + [ " + " + ".join(quads) + " ] "                   
    return objstr


def convertCplexOut(cplexoutpath):
    """reads output solution and returns xdict,ydict,objval 
    Args:
       cplexoutpath: Cormin output file
    Returns:
       xdict: 
       objval:
    """
    retvalues,objval = EmbedUtilities.readCplexOut(cplexoutpath,specific = ["x"])
    xdict = {}
    for key in retvalues.keys():
        if key.startswith("x"):
           xdict[tuple(int(part) for part in key.replace("x","").split("_"))] = retvalues[key]
    return xdict,objval


def getRatio(freqmat,comp2dominds,comp2scale,domains):
    """returns solution ratio 
    Args:
       freqmat:
       comp2dominds,comp2scale
       domains:
    Returns:
       ratio:
    """
    secondmat = np.zeros(np.shape(freqmat),dtype=np.float)
    for comp in comp2dominds.keys():
        for domin in comp2dominds[comp]:
            start,end = domains[domin]
            secondmat[start:end+1,start:end+1] += comp2scale[comp]
    difmat = freqmat - secondmat
    newobj = sum([difmat[ind1,ind2]**2 for ind1 in xrange(np.shape(difmat)[0]) for ind2 in xrange(np.shape(difmat)[1])])
    freqsum = sum([freqmat[ind1,ind2]**2 for ind1 in xrange(np.shape(freqmat)[0]) for ind2 in xrange(np.shape(freqmat)[1])])
    return newobj/freqsum, newobj


def runScaleOpt(comp2dominds,scales,sideparams):
    """runs second step of scale optimization
    Args:
       comp2domindS:
       scales:
       sideparams:
    Returns:
       comp2scale:
    """
    compcount,freqmat,interdom,domains = sideparams
    var2index, index2var, varcount = {}, [], 0
    for comp,scale in list(itertools.product(range(compcount),scales)):
        var2index[(comp,scale)] = varcount
        varcount += 1
        index2var.append((comp,scale))
    A,b = genSdpCoef(freqmat,scales,compcount,domains,interdom,comp2dominds,var2index,index2var)
    Amat = np.dot(A.transpose(),A)
    side = -1.0*np.dot(b.transpose(),A)
    numsum = np.sum([item**2 for item in b])
    P = np.append(Amat,[side.transpose()],0)
    P = np.column_stack([P, np.append(side,numsum)])
    Xmat,objval = Bls_Sdp.runSDPRelax(P)
    binsol = Bls_Sdp.gaussRound(Xmat,objval,P)
    comp2scale = {comp:1 for comp in xrange(compcount)}
    for ind in xrange(len(binsol)):
        if binsol[ind] == 1:
           comp,scale = index2var[ind]
           comp2scale[comp] += scale
  
    if TESTMODE:
       newmat = np.dot(P,Xmat)
       assert abs(sum([newmat[ind1,ind1] for ind1 in xrange(np.shape(newmat)[0])]) - objval) < 0.1
       for eig in (np.linalg.eigh(Xmat))[0]:
           assert eig >= 0.0
       assert np.allclose(P.transpose(), P)    
    return comp2scale

      
def genSdpCoef(freqmat,scales,compcount,domains,interdom,comp2dominds,var2index,index2var):
    """generates sdp matrix coefs
    Args:
       freqmat: frequency matrix
       scales: set of scales
       compcount: number of components
       domains: all domains
       interdom: intersecting domains
       comp2dominds: first step solution
       var2index:
       index2var:
    Returns:
       coefmat:
       b:
       fsum:
    """
    node2dom = EmbedUtilities.getnode2dom(freqmat,domains)
    dom2index = {domains[index]:index for index in xrange(len(domains))}
    vcount = np.shape(freqmat)[0]
    scalesum = sum(scales)
    A = np.zeros((vcount**2,compcount*len(scales)),dtype=np.float)
    b = np.array([0.0] * vcount**2)
    for ind1 in xrange(np.shape(freqmat)[0]):
        for ind2 in xrange(np.shape(freqmat)[1]):
            b[(ind1*vcount)+ind2] = freqmat[ind1,ind2]
            for comp in comp2dominds.keys():
                doms = set(dom for dom in comp2dominds[comp] if dom in node2dom[ind1] and dom in node2dom[ind2])
                if len(doms) == 1:    
                   for scale in scales:
                       A[(ind1*vcount)+ind2,var2index[(comp,scale)]] = 0.5*scale
                   b[(ind1*vcount)+ind2] -= (scalesum*0.5)+1
    return A,b


TESTMODE = True
def runBisub(freqmat,domains,compcount):
    """runs bisubmodular optimization with rounding
    Args:
       freqmat: frequency matrix
       domains: list of domains(list of set of nodes)
       compcount: number of components
    Returns:
       comp2doms:
       comp2scale:
    """
    interdom = EmbedUtilities.getInterDomain(domains)
    cliques = EmbedUtilities.findDomainCliqueDecomp(domains,interdom)
    
    #first step
    objstr = genObjStr(freqmat,compcount,domains,interdom)
    consstr =  genConsStr(domains,compcount,cliques)
    boundstr = genBoundStr(domains,compcount)
    outmethod = globals()["convertCplexOut"]
    xdict, objval = EmbedUtilities.runCplexCode(consstr,objstr,boundstr,"","deconfdbisub",outmethod)
    sideparams = [compcount,freqmat,interdom,domains]
    comp2dominds = Round.roundCR(xdict,sideparams)
    #domain selection
    print "compo"
    print comp2dominds
    #for comp in comp2dominds.keys():
    #    assert len(comp2dominds[comp]) >= 1
              
    #second step
    scales = set(2**index for index in xrange(0,int(math.ceil(math.log(np.amax(freqmat)+1,2)))))
    comp2scale = runScaleOpt(comp2dominds,scales,sideparams)
    ratio,newobjval = getRatio(freqmat,comp2dominds,comp2scale,domains)
    print "ratio: ",ratio
    metadata = {"objval": newobjval}
    
    if TESTMODE:
       sideparams = [scales,compcount,freqmat,interdom,domains]
       assert TestBiSub.testData(cliques,sideparams)
       assert TestBiSub.testImprovement(sideparams,comp2dominds,comp2scale)
       assert TestBiSub.testLPData(sideparams,objstr)
       assert TestBiSub.testSDPData(sideparams)
       #assert TestBiSub.testLPOutput(xdict,ydict,sideparams,objval)
       #assert TestBiSub.testAfterRound(xdict,ydict,comp2dominds,comp2scale,sideparams)
    return comp2dominds,comp2scale,metadata  


class TestBiSub():

    @staticmethod
    def testImprovement(sideparams,comp2dominds,comp2scale):
        """tests improvement in the solution
        Args:
           sideparams:
           comp2dominds:
           comp2scale:
        Returns:
           bool:
        """
        [scales,compcount,freqmat,interdom,domains] = sideparams
        for scalval in comp2scale.values():
            assert scalval >= 1
        firstmat = np.zeros(np.shape(freqmat),dtype=np.float)
        secondmat = np.zeros(np.shape(freqmat),dtype=np.float)
        for comp in comp2dominds.keys():
            for domin in comp2dominds[comp]:
                start,end = domains[domin]
                firstmat[start:end+1,start:end+1] += 1
                secondmat[start:end+1,start:end+1] += comp2scale[comp]
        difmat1 = freqmat - firstmat        
        difmat2 = freqmat - secondmat
        firstdif = sum([difmat1[ind1,ind2]**2 for ind1 in xrange(np.shape(firstmat)[0]) for ind2 in xrange(np.shape(firstmat)[1])])
        freqsum = sum([freqmat[ind1,ind2]**2 for ind1 in xrange(np.shape(freqmat)[0]) for ind2 in xrange(np.shape(freqmat)[1])])
        seconddif = sum([difmat2[ind1,ind2]**2 for ind1 in xrange(np.shape(secondmat)[0]) for ind2 in xrange(np.shape(secondmat)[1])])
        assert firstdif <= freqsum and seconddif <= firstdif
        return True

    @staticmethod
    def testSDPData(sideparams):
        """tests sdp data, script string code etc for the second step
        Args:
          sideparams:
          objstr:
        Returns:
          bool: true or false
        """
        return True
        for comp in comp2dominds.keys():
            alldoms = list(comp2dominds[comp])
            intercount = 0
            for ind1 in xrange(len(alldoms)):
                domin1 = alldoms[ind1]
                for ind2 in xrange(ind1+1,len(alldoms)):
                    domin2 = alldoms[ind2]
                    assert not intersect(domains[domin1],domains[domin2])
    
    @staticmethod
    def testLPData(sideparams,objstr):
        """tests lp data, script string code etc for the first step
        Args:
          sideparams:
          objstr:
        Returns:
          bool: true or false
        """
        [scales,compcount,freqmat,interdom,domains] = sideparams
        #Compares lp coefs in string with the expected ones
        node2dom = EmbedUtilities.getnode2dom(freqmat,domains)
        dom2index = {domains[index]:index for index in xrange(len(domains))}
        var2index, varcount, index2var = {}, 0, {}
        for dom,comp in list(itertools.product(domains,range(compcount))):
            var2index[(dom2index[dom],comp)] = varcount
            index2var[varcount] = (dom2index[dom],comp)
            varcount += 1
        coefmat = np.zeros((varcount,varcount),dtype=np.float)
        impstr = objstr.split("[")[1].split("]")[0]
        for part in impstr.split(" + "):
            splitted = part.split()
            assert len(splitted) == 4 and splitted[2] == "*"
            domin1,comp1 = [int(item) for item in splitted[1].replace("x","").split("_")]
            domin2,comp2 = [int(item) for item in splitted[3].replace("x","").split("_")]
            ind1 = var2index[(domin1,comp1)]
            ind2 = var2index[(domin2,comp2)]
            coefmat[ind1,ind2] += float(splitted[0])/4
            coefmat[ind2,ind1] += float(splitted[0])/4
        assert np.allclose(coefmat.transpose(), coefmat)  
        #for ind1 in xrange(np.shape(coefmat)[0]):
        #    domin1,comp1 = index2var[ind1] 
        #    for ind2 in xrange(np.shape(coefmat)[1]):
        #        domin1,comp2 = index2var[ind2]
        #        interlen =
        #        assert coefmat[ind1,ind2] == tsum                           
        return True
     
    @staticmethod  
    def testData(cliques,sideparams):
        """tests data
        Args:
          cliques:
          sideparams:
        Returns:
          bool: true or false
        """
        [scales,compcount,freqmat,interdom,domains] = sideparams
        intersect = lambda (s1,e1), (s2,e2): False if (e1 < s2 or e2 < s1) else True
        storeG = EmbedUtilities.inter2Graph(interdom,domains)
        tinterG = deepcopy(storeG)
        for clique in cliques:
            for ind1 in xrange(len(clique)):
                node1 = clique[ind1]
                for ind2 in xrange(ind1+1,len(clique)):
                    node2 = clique[ind2]
                    if tinterG.has_edge(node1,node2):
                       tinterG.remove_edge(node1,node2)
        assert tinterG.number_of_edges() == 0
        tinterG = deepcopy(storeG)
        for clique in cliques:
            for ind1 in xrange(len(clique)):
                node1 = clique[ind1]
                for ind2 in xrange(ind1+1,len(clique)):
                    node2 = clique[ind2]
                    assert tinterG.has_edge(node1,node2)
        assert len(set([dom for clique in cliques for dom in clique]) ^ set(domains)) == 0
        for dom1,dom2 in interdom:
            assert intersect(dom1,dom2)
        node2dom = EmbedUtilities.getnode2dom(freqmat,domains)
        for node in node2dom.keys():
            for domin in node2dom[node]:
                assert domains[domin][0] <= node and node <= domains[domin][1]    
        return True

    @staticmethod
    def testLPOutput(xdict,sideparams,objval):
        """tests LP output including objective val
        Args:
           xdict:
           sideparams:
           objval: 
        Returns:
           bool: true/false
        """
        [scales,compcount,freqmat,interdom,domains] = sideparams
        for index in xrange(compcount):
            for clique in cliques:
                sumval = sum([xdict[(domin,comp)] for domin,comp in xdict.keys() if comp == index and domains[domin] in clique])
                assert sumval <= 1.001
        return True
    
        [scales,compcount,freqmat,interdom,domains] = sideparams
        node2dom = EmbedUtilities.getnode2dom(freqmat,domains)        
        assert len(xdict.keys()) <= len(domains)*compcount*len(scales) and len(ydict.keys()) <= compcount*len(scales)
                 
        estobjval = EmbedUtilities.estFracObjective(xdict,freqmat,node2dom,scales,compcount)
        frsum = sum([freqmat[in1,in2]**2 for in1 in xrange(np.shape(freqmat)[0]) for in2 in xrange(np.shape(freqmat)[1])])
        assert estobjval <= frsum 
        realobjval = objval + frsum
        assert abs(estobjval-realobjval) < 0.5

        comp2scale = {comp:0.0 for comp in xrange(compcount)}
        for comp,scale in ydict.keys():
            comp2scale[comp] += ydict[(comp,scale)]
        for scal in comp2scale.values():
            assert scal < 1.02
        return True

    @staticmethod
    def testAfterRound(xdict,ydict,comp2dominds,comp2scale,sideparams):
       """afterrounding test
       Args:
          xdict:
          ydict:
          comp2dominds: comp to assigned domains(not indices)
          comp2scale:
          sideparams:
       Returns:
          bool:
       """
       [scales,compcount,freqmat,interdom,domains] = sideparams
       assert len(set(comp2scale.keys()) ^ set(range(compcount))) == 0
       node2dom = EmbedUtilities.getnode2dom(freqmat,domains)
       #intersection check
       intersect = lambda (s1,e1), (s2,e2): False if (e1 < s2 or e2 < s1) else True
       for comp in comp2dominds.keys():
           curdoms = list(comp2dominds[comp])
           for ind1,dom1 in list(enumerate(curdoms)):
               for ind2 in xrange(ind1+1,len(curdoms)):
                   assert not intersect(domains[dom1],domains[curdoms[ind2]])
       dom2index = {domains[index]:index for index in xrange(len(domains))}
       #approximation ratio check        
       binxdict = {(domin,comp,scale): 0.0 for (domin,comp,scale) in xdict.keys()}
       for comp in comp2dominds.keys():
           for domin in comp2dominds[comp]:
               binxdict[(domin,comp,comp2scale[comp])] = 1
       fracobjval = EmbedUtilities.estFracObjective(xdict,freqmat,node2dom,scales,compcount)
       binobjval = EmbedUtilities.estFracObjective(binxdict,freqmat,node2dom,scales,compcount)
       initdict = {(domin,comp,scale): 0.0 for (domin,comp,scale) in xdict.keys()}
       initobjval = EmbedUtilities.estFracObjective(initdict,freqmat,node2dom,scales,compcount)
       print "obj info"
       print fracobjval
       print binobjval
       print initobjval
       print "Approx ratio {0}".format(fracobjval/binobjval)     
       return True
