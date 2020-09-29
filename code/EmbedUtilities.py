#Methods used in Optimization
import os
import sys
import time
import numpy as np
import networkx as nx
import gzip

def findDomainCliqueDecomp(domains,interdom):
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
    return findCliqueDecomp(G)

    
def findCliqueDecomp(interG):
    """finds maximal clique decomposition 
    Args:
       interG
    Returns:
       cliques:
    """
    tinterG = nx.Graph(interG)
    return list(map(list,nx.find_cliques(tinterG)))


def estFracObjective(xdict,freqmat,node2dom,scales,compcount):
    """estimates objective function when xdict may be fractional
    Args:
       xdict:
       freqmat:
       node2dom:
       scales:
       compcount:
    Returns:
       objval:
    """
    objval = 0.0
    for in1 in xrange(np.shape(freqmat)[0]):
        for in2 in xrange(np.shape(freqmat)[1]):
            cursum = freqmat[in1,in2]
            domset = node2dom[in1].intersection(node2dom[in2])
            cursum -= sum([scale*xdict[(dom,comp,scale)] for dom in domset for comp in xrange(compcount) for scale in scales if xdict.has_key((dom,comp,scale))])
            objval += cursum**2                      
    return objval


def getnode2dom(freqmat,domains):
    """return node 2 domain mapping
    Args:
       freqmat:
       domains:
    Returns:
       node2dom:
    """
    node2dom = {node:set() for node in xrange(np.shape(freqmat)[0])}
    for ind,dom in list(enumerate(domains)):
        for node in xrange(dom[0],dom[1]+1):
            node2dom[node].add(ind)
    return node2dom


def inter2Graph(interdom,domains):
    """returns inter domain graph
    Args:
       interdom: list of domains
       domains:
    Returns:
       G: 
    """
    G = nx.Graph()
    for dom in domains:
        G.add_node(dom)
    for dom1,dom2 in interdom:
        G.add_edge(dom1,dom2)
    return G
 
def getInterDomain(domains):
    """returns intersecting domain set
    Args:
       domains: list of domains
    Returns:
       interset: 
    """
    intersect = lambda (s1,e1), (s2,e2): False if (e1 < s2 or e2 < s1) else True
    enumdomains = list(enumerate(domains))
    return set((dom1,enumdomains[in2][1]) for in1,dom1 in enumdomains for in2 in xrange(in1+1,len(enumdomains)) if intersect(dom1,enumdomains[in2][1]))


def writeDomainFile(domoutfile,domains):
    """writes domain file
    Args:
       domoutfile:
       domains:
       comp2scale:
    Returns:
    """
    with open(domoutfile,"w") as outfile:
        for start,end in domains:
            outfile.write("{0},{1}\n".format(start,end))


def readDomainFile(domfile):
    """reads domain file
    Args:
       domfile:
    Returns:
       domains:
    """
    domains = []
    with open(domfile,"r") as infile:
        for line in infile:
            splitted = line.rstrip().split(",")
            domains.append((int(splitted[0]),int(splitted[1])))
    return domains


def readFreqFile(freqfile):
    """reads input matrix file in gz format
    Args:
       freqfile:
    Returns:
       freqmat:
       nodenames:
    """    
    nodenames = []
    with gzip.open(freqfile,"r") as infile:
        for line in infile:
            nodenames.append(line.rstrip().split("\t")[0])
    freqmat = np.zeros((len(nodenames),len(nodenames)),dtype=np.int)
    index = 0       
    with gzip.open(freqfile,"r") as infile:
        for line in infile:
            parts = line.rstrip().split("\t")
            for index2 in xrange(len(nodenames)):
                freqmat[index,index2] = int(float(parts[3+index2]))
            index+=1
    return freqmat,nodenames


def writeFreqData(freqmat,freqfile):
    """writes frequency data to file
    Args:
       freqmat:
       freqfile:
    Returns:
    """
    with gzip.open(freqfile,"w") as outfile:
        for in1 in xrange(np.shape(freqmat)[0]):
            linestr = "{0}\tba\tca\t".format(in1) + "\t".join(["{0}".format(freqmat[in1,in2]) for in2 in xrange(np.shape(freqmat)[1])])
            outfile.write(linestr+"\n")
            

def writeDeconMeta(comp2dominds,comp2scale,domains,metadata,objoutfile):
    """writes deconvolution meta file
    Args:
       comp2dominds:
       comp2scale:
       domains:
       metadata:
       objoutfile:
    Returns:
    """
    with open(objoutfile,"w") as outfile:
        outfile.write("{0}\n".format(metadata["objval"]))

def writeDeconOut(comp2dominds,comp2scale,domains,deconoutfile):
    """writes deconvolution output
    Args:
       comp2dominds: list of tuples
       comp2scale:
       domains:
       outensfile: output ensemble file
    Returns:
    """
    with open(deconoutfile,"w") as outfile:
        for comp in comp2dominds.keys():
            linestr = "\t".join([str(comp2scale[comp])] + ["{0},{1}".format(domains[domin][0],domains[domin][1]) for domin in comp2dominds[comp]])
            outfile.write(linestr+"\n") 
     
            
def readCplexOut(outfile,specific=[]):
    """reads CPLEX output file and returns only SPECIFIC variable values as dictionary
    Args:
       outfile: CPLEX output file
       specific: specific variable prefixes such as x
    Returns:
       retvalues: variable-value dictionary
    """
    retvalues = {}
    varflag = False
    objval = None
    with open(outfile,"r") as file:
        for line in file:
            line = line.rstrip()
            if not varflag and line.find("- Optimal:")!=-1 and line.find("Objective")!=-1:
               objval = float(line.split("=")[1])
               continue    
            if not varflag and line.find("Variable Name")!=-1 and line.find("Solution Value")!=-1:
               varflag=True
               continue
            if varflag:
               for varname in specific: 
                   if line.startswith(varname):
                      key,value = line.split()
                      retvalues[key] = float(value)
                      break
    return retvalues,objval 


def runCplexCode(consstr,objstr,boundstr,varstr,runfolder,outmethod):
    """Runs cplex code
    Args:
        consstr: constraint string
        objstr: objective function string
        boundstr: boundary string
        varstr: variable string
        runfolder: run folder
        outmethod = function to be run before returning output
    Returns:
        xdict:
        ydict:
    """
    if not os.path.exists(runfolder):
       os.makedirs(runfolder)
    PNUM = 1 #processor count
    filepref = "cplexrun"
    outlppath = "{0}/{1}.lp".format(runfolder,filepref)
    with open(outlppath,"w") as file:
       file.write(objstr+"\n")
       file.write(consstr)
       file.write(boundstr+"\n")
       file.write(varstr+"\n")
       file.write("End\n")
    cplexoutpath = "{0}/{1}.lpout".format(runfolder,filepref)
    cplexscriptpath = "{0}/{1}.script".format(runfolder,filepref)
    with open(cplexscriptpath,"w") as file:
       file.write("read {0}\n".format(outlppath))
       file.write("set threads {0}\n".format(PNUM))
       file.write("optimize\n")
       file.write("display solution objective\n")
       file.write("display solution variables -\n")  
    t1 = time.time()
    code="cplex < {0} > {1}".format(cplexscriptpath,cplexoutpath)
    os.system(code)
    t2 = time.time()
    print "Cplex solved problem in {0} seconds".format(t2-t1)
    returns = outmethod(cplexoutpath)
    os.system("rm -rf {0}".format(runfolder))
    return returns                                 
