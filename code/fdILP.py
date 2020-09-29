#ILP Formulation of Frequency Deconvolution
import networkx as nx
import numpy as np
import scipy as sp
import math
import random
import os
import sys
import itertools
import EmbedUtilities
from Round import Round
from copy import deepcopy


def roundSolution(xdict,ydict,sideparams,roundmethod):
    """rounds solution
    Args:
       xdict:
       ydict:
       sideparams:
       roundmethod:
    Returns:
       comp2doms:
    """
    assert roundmethod in ["random","greedy","greedy2"]
    if roundmethod == "random":
       return Round.roundRand(xdict,ydict,sideparams)
    elif roundmethod == "greedy":
       return Round.greedyRand(xdict,ydict,sideparams)
    elif roundmethod == "greedy2":
       return Round.greedy2Rand(xdict,sideparams)  

              
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
    #boundstr += " ".join([" 0 <= y{0}_{1} <= 1\n".format(comp,scale) for comp in xrange(compcount) for scale in scales])           
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
       cliques2:
    Returns:
       objstr:
    """
    dom2index = {domains[index]: index for index in xrange(len(domains))}
    consstr = " Subject To \n"
    for dom1 in domains:
        for dom2 in domains:
            if (dom1,dom2) not in interdom and (dom2,dom1) not in interdom and dom1 != dom2:
               for comp in xrange(compcount):
                   parts = []
                   for scale in scales:
                       parts.append(" x{0}_{1}_{2} ".format(dom2index[dom1],comp,scale))
                       for scale2 in scales - set([scale]):
                           parts.append(" x{0}_{1}_{2} ".format(dom2index[dom2],comp,scale2))
                   consstr += " + ".join(parts) + " <= 1 \n "
                   parts = []
                   for scale in scales:
                       parts.append(" x{0}_{1}_{2} ".format(dom2index[dom2],comp,scale))
                       for scale2 in scales - set([scale]):
                           parts.append(" x{0}_{1}_{2} ".format(dom2index[dom1],comp,scale2))
                   consstr += " + ".join(parts) + " <= 1 \n "
    for clique in cliques:
        consstr += " ".join([" + ".join([" x{0}_{1}_{2} ".format(dom2index[dom],comp,scale) for dom in clique for scale in scales]) + " <= 1 \n " for comp in xrange(compcount)])
    return consstr


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
        qcoef = (min(dom1[1],dom2[1])-max(dom1[0],dom2[0]) + 1)**2
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
    """reads output solution and returns xdict,ydict,objval 
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


def findCliqueDecomp(domains,interdom):
    """finds maximal clique decomposition 
    Args:
       domains:
       interdom:
    Returns:
       cliques:
    """
    G = nx.Graph()
    for dom in domains:
        G.add_node(dom)
    for dom1,dom2 in interdom:
        G.add_edge(dom1,dom2)
    cliques = []
    while G.number_of_nodes() > 0:
        maxclique,maxsize = None, None
        allcliques = list(nx.find_cliques(G))
        for tclique in allcliques:
            if len(tclique) >= maxsize:
               maxclique = tclique
        for ind1 in xrange(len(maxclique)):
            node1 = maxclique[ind1]
            for ind2 in xrange(ind1+1,len(maxclique)):
                node2 = maxclique[ind2]
                G.remove_edge(node1,node2)
        for node in G.nodes():
            if len(G.neighbors(node)) == 0:
               G.remove_node(node)         
        cliques.append(maxclique)
    return cliques


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
    for domin,comp,scale in cursol.keys():
        if cursol[(domin,comp,scale)] == 1:
           comp2doms[comp].add(domin)
           comp2scale[comp].add(scale)     
    for comp in comp2doms.keys():          
        alldoms = list(comp2doms[comp])
        for ind1 in xrange(len(alldoms)):
            domin1 = alldoms[ind1]
            for ind2 in xrange(ind1+1,len(alldoms)):
                domin2 = alldoms[ind2]
                if intersect(domains[domin1],domains[domin2]):
                   return True
        if len(comp2scale[comp]) > 1:
           return True      
    return False


TESTMODE = True
def runILP(freqmat,domains,scales,compcount,params):
    """runs ILP optimization and then rounds
    Args:
       freqmat: frequency matrix
       domains: list of domains(list of set of nodes)
       scales: set of scales
       compcount: number of components
       params: parameters
    Returns:
       comp2doms:
       comp2scale:
    """
    interdom = EmbedUtilities.getInterDomain(domains)
    cliques = EmbedUtilities.findDomainCliqueDecomp(domains,interdom)
    #interdom2 = []
    #for ind1 in xrange(len(domains)):
    #    dom1 = domains[ind1]
    #    for ind2 in xrange(ind1+1,len(domains)):
    #        dom2 = domains[ind2]
    #        if (dom1,dom2) not in interdom and (dom2,dom1) not in interdom:
    #           interdom2.append((dom1,dom2))           
    #compcliques2 = findCliqueDecomp(domains,interdom2)
         
    objstr = genObjStr(freqmat,scales,compcount,domains,interdom)
    consstr = genConsStr(domains,interdom,scales,compcount,cliques)
    boundstr = genBoundStr(domains,scales,compcount)   
    outmethod = globals()["convertCplexOut"]
    xdict, ydict, objval = EmbedUtilities.runCplexCode(consstr,objstr,boundstr,"","deconilpfolder",outmethod)
     
    sideparams = [scales,compcount,freqmat,interdom,domains]       
    comp2dominds,comp2scale = roundSolution(xdict,ydict,sideparams,params["roundmethod"])
    print "done rounding"
    exit(1)
    
    if TESTMODE:
       assert TestILP.testData(cliques,sideparams)
       assert TestILP.testLPData(sideparams,objstr)
       assert TestILP.testLPOutput(xdict,ydict,sideparams,objval)
       assert TestILP.testAfterRound(xdict,ydict,comp2dominds,comp2scale,sideparams)
    return comp2dominds,comp2scale   


class TestILP():

    @staticmethod
    def testLPData(sideparams,objstr):
        """tests lp data, script string code etc
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
        var2index, varcount = {}, 0
        for dom,comp,scale in list(itertools.product(domains,range(compcount),scales)):
            var2index[(dom2index[dom],comp,scale)] = varcount
            varcount += 1
        coefmat = np.zeros((varcount,varcount),dtype=np.float)
        impstr = objstr.split("[")[1].split("]")[0]
        for part in impstr.split(" + "):
            splitted = part.split()
            assert len(splitted) == 4 and splitted[2] == "*"
            domin1,comp1,scale1 = [int(item) for item in splitted[1].replace("x","").split("_")]
            domin2,comp2,scale2 = [int(item) for item in splitted[3].replace("x","").split("_")]
            ind1 = var2index[(domin1,comp1,scale1)]
            ind2 = var2index[(domin2,comp2,scale2)]
            coefmat[ind1,ind2] += float(splitted[0])/4
            coefmat[ind2,ind1] += float(splitted[0])/4
        scalesum = sum(scales)
        coefsum = sum([(scalesum*compcount*len(node2dom[in1].intersection(node2dom[in2])))**2 for in1 in xrange(np.shape(freqmat)[0]) for in2 in xrange(np.shape(freqmat)[1])])
        assert abs(coefsum - np.sum(coefmat)) <= 0.01 
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
        assert len(set([dom for clique in cliques for dom in clique]) ^ set(domains)) == 0
        for dom1,dom2 in interdom:
            assert intersect(dom1,dom2)
        return True

    @staticmethod
    def testLPOutput(xdict,ydict,sideparams,objval):
        """tests LP output including objective val
        Args:
           xdict:
           ydict:
           sideparams:
           objval: 
        Returns:
           bool: true/false
        """
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
