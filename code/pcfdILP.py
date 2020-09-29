#ILP Formulation of Prize Collecting Frequency Deconvolution
import networkx as nx
import numpy as np
import scipy as sp
import scipy.linalg
import scipy.sparse
import math
import os
import sys
import itertools
import EmbedUtilities
              
def genBoundStr(domains,scales,compcount):
    """generates lp relaxation bound constraints
    Args:
       domains: set of domains
       scales: set of scales
       compcount: number of components 
    Returns:
       boundstr:
    """
    boundstr = " Bounds\n "
    boundstr += "\n"+ " ".join([" 0 <= x{0}_{1}_{2} <= 1\n".format(domin,comp,scale) for domin in xrange(len(domains)) for comp in xrange(compcount) for scale in scales])
    boundstr += " ".join([" 0 <= y{0}_{1} <= 1\n".format(comp,scale) for comp in xrange(compcount) for scale in scales])           
    return boundstr

def genConsStrNOTUSED(domains,interdom,scales,compcount,seendoubles,cliques):
    """generates lp relaxation constraints
    Args:
       domains: set of domains
       interdom: interacting domain pairs
       scales: set of scales
       compcount: number of components
       seendoubles:
       cliques:
    Returns:
       objstr:
    """
    dom2index = {domains[index]: index for index in xrange(len(domains))}
    consstr = " Subject To \n" 
    #for zvar in seendoubles:
    #    u1,u2,u3,v1,v2,v3 = zvar.replace("z","").split("_")
    #    consstr += " {0} - x{1}_{2}_{3} <= 0\n ".format(zvar,u1,u2,u3)
    #    consstr += " {0} + x{1}_{2}_{3} <= 1\n ".format(zvar,v1,v2,v3)
    consstr += " \n ".join([" + ".join([" y{0}_{1} ".format(comp,scale) for scale in scales]) + " <= 1" for comp in xrange(compcount)])  
    consstr += "\n" + " ".join([" x{0}_{1}_{2} - y{1}_{2} <= 0\n".format(domind,comp,scale) for domind in xrange(len(domains)) for comp in xrange(compcount) for scale in scales])
    for clique in cliques:
        consstr += " ".join([" + ".join([" x{0}_{1}_{2} ".format(dom2index[dom],comp,scale) for dom in clique]) + " - y{0}_{1} <= 0 \n".format(comp,scale) for comp in xrange(compcount) for scale in scales])
    #for dom1,dom2 in interdom:
    #    consstr += "\n" + " ".join([" x{0}_{1}_{2} + x{3}_{1}_{2} - y{1}_{2} <= 0\n".format(dom2index[dom1],comp,scale,dom2index[dom2]) for comp in xrange(compcount) for scale in scales])
    return consstr


def genConsStr(domains,interdom,scales,compcount,cliques):
    """generates lp relaxation constraints
    Args:
       domains: set of domains
       interdom: interacting domain pairs
       scales: set of scales
       compcount: number of components
       cliques:
    Returns:
       objstr:
    """
    dom2index = {domains[index]: index for index in xrange(len(domains))}
    consstr = " Subject To \n" 
    consstr += " \n ".join([" + ".join([" y{0}_{1} ".format(comp,scale) for scale in scales]) + " <= 1" for comp in xrange(compcount)])
    for clique in cliques:
        consstr += " ".join([" + ".join([" x{0}_{1}_{2} ".format(dom2index[dom],comp,scale) for dom in clique]) + " - y{0}_{1} <= 0 \n".format(comp,scale) for comp in xrange(compcount) for scale in scales])
    return consstr


def genObjStrNOTUSED(freqmat,node2dom,scales,compcount):
    """generates objective function string for maximization
    Args:
       freqmat: frequency matrix
       node2dom: node to domain index mapping
       scales: set of scales
       compcount: number of components 
    Returns:
       objstr:
    """
    seendoubles = set()
    objstr = " Maximize\n obj: "
    isPlus = lambda x: "+" if x >= 0 else " "
    domains = set(dom for doms in node2dom.values() for dom in doms)
    lincoefs = {(dom,comp,scale): 0 for dom in domains for comp in xrange(compcount) for scale in scales}  
    for index1 in xrange(np.shape(freqmat)[0]):
        for index2 in xrange(index1,np.shape(freqmat)[1]):
            domset = node2dom[index1].intersection(node2dom[index2])
            seenset = [(dom,comp,scale) for dom in domset for scale in scales for comp in xrange(compcount)]
            quads = [" + {6} z{0}_{1}_{2}_{3}_{4}_{5} ".format(dom1,comp1,scale1,dom2,comp2,scale2,scale1*scale2) for dom1,comp1,scale1 in seenset for dom2,comp2,scale2 in seenset] 
            seendoubles |= set(["z{0}_{1}_{2}_{3}_{4}_{5}".format(dom1,comp1,scale1,dom2,comp2,scale2) for dom1,comp1,scale1 in seenset for dom2,comp2,scale2 in seenset])  
            objstr += " ".join(quads)
            scalesum = sum([scale for dom,comp,scale in seenset])
            for dom1,comp1,scale1 in seenset:
                lincoefs[(dom1,comp1,scale1)] -= scale1*(scalesum - 2*freqmat[index1,index2]) 
    objstr += " ".join([" {0} {1} x{2}_{3}_{4} ".format(isPlus(lincoefs[(dom,comp,scale)]),lincoefs[(dom,comp,scale)],dom,comp,scale) for (dom,comp,scale) in lincoefs.keys() if lincoefs[(dom,comp,scale)] != 0])
    return objstr,seendoubles


def genObjStr(freqmat,scales,compcount,domains,interdom):
    """generates objective function string
    Args:
       freqmat: frequency matrix
       scales: set of scales
       compcount: number of components
       domains: all domains
       interdom: intersecting domains 
    Returns:
       objstr: 
    """
    dom2index = {domains[index]:index for index in xrange(len(domains))}
    objstr = " Minimize\n obj: " 
    singles, quads = [], []
    pairs = [(comp,scale) for comp in xrange(compcount) for scale in scales]
    for dom1,dom2 in interdom:
        domin1 = dom2index[dom1]
        domin2 = dom2index[dom2]
        qcoef = len(set(range(dom1[0],dom1[1]+1)).intersection(set(range(dom2[0],dom2[1]+1))))**2
        quads.extend([" {0} x{1}_{2}_{3} * x{4}_{5}_{6} ".format(4*qcoef*scale1*scale2,domin1,comp1,scale1,domin2,comp2,scale2) for comp1,scale1 in pairs for comp2,scale2 in pairs])      
    for dom in domains:
        domin = dom2index[dom]
        qcoef = (dom[1]-dom[0]+1)**2
        quads.extend([" {0} x{1}_{2}_{3} * x{1}_{4}_{5} ".format(2*qcoef*scale1*scale2,domin,comp1,scale1,comp2,scale2) for comp1,scale1 in pairs for comp2,scale2 in pairs])
        fsum = np.sum(freqmat[dom[0]:dom[1]+1,dom[0]:dom[1]+1])
        singles.extend([" - {0} x{1}_{2}_{3} ".format(2*scale*fsum,domin,comp,scale) for comp,scale in pairs])
    if len(singles) > 0:
       objstr += " ".join(singles)
    if len(quads) > 0:
       objstr += " + [ " + " + ".join(quads) + " ] "                   
    return objstr


def convertCplexOut(cplexoutpath):
    """reads output solution and returns values plus objval
    Args:
       cplexoutpath: Cormin output file
    Returns:
       xdict:
       ydict:
       objval:
    """
    retvalues,objval = EmbedUtilities.readCplexOut(cplexoutpath,specific = ["x","y"])
    xdict, ydict = {}, {}
    for key in retvalues.keys():
        if key.startswith("x"):
           xdict[tuple(int(part) for part in key.replace("x","").split("_"))] = retvalues[key]
        elif key.startswith("y"):
           ydict[tuple(int(part) for part in key.replace("y","").split("_"))] = retvalues[key]
    return xdict,ydict,objval

def subIntervalScheduling():
    """nonmonotone submodular interval scheduling
    Args:
       
    Returns:
    """
    return


def genSdpCoef(freqmat,scales,compcount,domains,interdom,overpen,scalepen):
    """generates sdp matrix coefs
    Args:
       freqmat: frequency matrix
       scales: set of scales
       compcount: number of components
       domains: all domains
       interdom: intersecting domains
       overpen: overlap penalty
       scalepen: scale unmatch penalty
    Returns:
       objstr: 
    """
    dom2index = {domains[index]:index for index in xrange(len(domains))}
    var2index, index2var, varcount = {}, [], 0
    for dom,comp,scale in list(itertools.product(domains,range(compcount),scales)):
        var2index[(dom2index[dom],comp,scale)] = varcount
        varcount += 1
        index2var.append((dom2index[dom],comp,scale))
    coefmat = np.zeros((varcount,varcount),dtype=np.float) #coefmat = scipy.sparse.bsr_matrix((varcount,varcount), dtype=np.int)
    bvec = [0] * varcount
    pairs = [(comp,scale) for comp in xrange(compcount) for scale in scales]
    for dom1,dom2 in interdom:
        domin1 = dom2index[dom1]
        domin2 = dom2index[dom2]
        interlen = min(dom1[1],dom2[1])-max(dom1[0],dom2[0]) + 1
        interpen = overpen * interlen
        for (comp1,scale1),(comp2,scale2) in itertools.product(pairs,pairs):
            coefmat[var2index[(domin1,comp1,scale1)],var2index[(domin2,comp2,scale2)]] += (interlen**2)*scale1*scale2
            coefmat[var2index[(domin2,comp2,scale2)],var2index[(domin1,comp1,scale1)]] += (interlen**2)*scale1*scale2
        for comp in xrange(compcount):
            for scale1,scale2 in itertools.product(scales,scales):
                coefmat[var2index[(domin1,comp,scale1)],var2index[(domin2,comp,scale2)]] += interpen
                coefmat[var2index[(domin2,comp,scale2)],var2index[(domin1,comp,scale1)]] += interpen
    for dom in domains:
        domin = dom2index[dom]
        qcoef = (dom[1]-dom[0]+1)**2
        fsum = np.sum(freqmat[dom[0]:dom[1]+1,dom[0]:dom[1]+1])
        for (comp1,scale1),(comp2,scale2) in itertools.product(pairs,pairs):
            coefmat[var2index[(domin,comp1,scale1)],var2index[(domin,comp2,scale2)]] += 0.5*qcoef*scale1*scale2
            coefmat[var2index[(domin,comp2,scale2)],var2index[(domin,comp1,scale1)]] += 0.5*qcoef*scale1*scale2
        for comp,scale in pairs:
            #coefmat[var2index[(domin,comp,scale)],var2index[(domin,comp,scale)]] += -2.0*scale*fsum 
            bvec[var2index[(domin,comp,scale)]] += -2.0*scale*fsum
    return coefmat, index2var


def convertCoefMat(coefmat):
    """converts matrix for 0,1 to -1,1 case
    Args:
       coefmat:
    Returns:
       modmat: modified matrix
    """
    modmat = np.zeros(np.shape(coefmat),dtype=np.float64)
    return coefmat /4.0
    #for in1 in xrange(np.shape(modmat)[0]):
    #    for in2 in xrange(np.shape(modmat)[1]):
    #        modmat[in1,in2] += coefmat[in1,in2] / 4.0
    #return modmat


def genSDPProg(modmat,index2var):
    """generates sdp program
    Args:
       modmat:
       index2var:
    Returns:
       fracsol:
    """
    fracsol = {}
    return fracsol
   

def HyperRound(fracsol):
    """hyperplane based rounding
    Args:
       fracsol:
    Returns:
       binsol:
    """
    binsol = {}
    return binsol
            
    
TESTMODE = True
def runSDP(freqmat,domains,scales,compcount,params):
    """runs SDP optimization relaxation
    Args:
       freqmat: frequency matrix
       domains: list of domains(list of set of nodes)
       scales: set of scales
       compcount: number of components
       params: overpen, scalepen, roundmethod
    Returns:
       comp2dom2scale:   
    """
    interdom = EmbedUtilities.getInterDomain(domains)
    node2dom = EmbedUtilities.getnode2dom(freqmat,domains)
    
    for node1 in xrange(np.shape(freqmat)[0]):
        doms1 = node2dom[node1]
        for node2 in xrange(node1,np.shape(freqmat)[1]):
            doms2 = node2dom[node2]
            idomset = node2dom[node1].intersection(node2dom[node2])
            for idomin in idomset:
                idom = domains[idomin]
                assert idom[0] <= node1 and idom[1] >= node1 and idom[0] <= node2 and idom[1] >= node2
    
    matcoef, index2var = genSdpCoef(freqmat,scales,compcount,domains,interdom,params["overpen"],params["scalepen"])
    print "coef gen"
    modmat = convertCoefMat(matcoef)
    print "coef converted"
    fracsol = genSDPProg(modmat,index2var)
    comp2dom2scale = HyperRound(fracsol)    
     
    if TESTMODE:
       assert SDPTest.testCoefs(matcoef,modmat,params,scales,domains,compcount,interdom,freqmat)
       assert SDPTest.testRound(matcoef,modmat,params,scales,domains,compcount,interdom)
       assert SDPTest.testSDPProg(matcoef,modmat,params,scales,domains,compcount,interdom)
    exit(1)   
    return comp2dom2scale 


class SDPTest():
        
    @staticmethod
    def testSDPProg(matcoef,modmat,params,scales,domains,compcount,interdom):
        """tests sdp running
        Args:
           matcoef: matrix coef
           modmat: modified matrix
           params:
           scales:
        Returns:
           bool:   
        """
        EmbedUtilities.estFracObjective()
        return True
    
    @staticmethod
    def testRound(matcoef,modmat,params,scales,domains,compcount,interdom):
        """tests rounding
        Args:
           matcoef: matrix coef
           modmat: modified matrix
           params:
           scales:
        Returns:
           bool:   
        """
        EmbedUtilities.estBinaryObjective()
        return True
   
    @staticmethod
    def testCoefs(matcoef,modmat,params,scales,domains,compcount,interdom,freqmat):
        """tests sdp coefs 
        Args:
           matcoef: matrix coef
           modmat: modified matrix
           params:
           scales:
           domains:
           compcount:
           interdom:
           freqmat:
        Returns:
           bool:   
        """
        node2dom = EmbedUtilities.getnode2dom(freqmat,domains)
        #test quadratic maxtrix coefs
        matcoefsum = np.sum(matcoef)
        scalesum = sum(scales)
        coefsum = sum([(scalesum*compcount*len(node2dom[in1].intersection(node2dom[in2])))**2 for in1 in xrange(np.shape(freqmat)[0]) for in2 in xrange(np.shape(freqmat)[1])])
        coefsum2 = sum([sum([scale for dom in node2dom[in1].intersection(node2dom[in2]) for scale in scales for comp in xrange(compcount)])**2 for in1 in xrange(np.shape(freqmat)[0]) for in2 in xrange(np.shape(freqmat)[1])])
        coefsum3 = 0.0
        for in1 in xrange(np.shape(freqmat)[0]):
            for in2 in xrange(np.shape(freqmat)[1]):
                domset = set(domin for domin in node2dom[in1].union(node2dom[in2]) if domains[domin][0] <= in1 and domains[domin][1] >= in1 and domains[domin][0] <= in2 and domains[domin][1] >= in2)
                coefsum3 += sum([scale for dom in domset for scale in scales for comp in xrange(compcount)])**2       
        assert abs(coefsum - matcoefsum) <= 0.0001 and abs(coefsum2 - matcoefsum) <= 0.0001 and abs(coefsum3 - coefsum2) <= 0.0001

        #tests quadratic matrix's positive semidefiniteness
        assert np.allclose(matcoef.transpose(), matcoef)
        if params["overpen"] == 0:
           E,V = scipy.linalg.eigh(matcoef)
           for eig in E:
               assert eig >= -0.1
           E,V = scipy.linalg.eigh(modmat)
           for eig in E:
               assert eig >= -0.1
        return True         
