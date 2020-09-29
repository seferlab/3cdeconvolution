#Main Deconvolution Executor
import networkx as nx
import numpy as np
import scipy as sp
import math
import os
import EmbedUtilities
import fdILP
import pcfdILP
import fdLocal
import fdBiSub


def runDeconvolution(method,objtype,freqfilename,domainfile,compcount,outfolder):
    """runs deconvolution method
    Args:
       method:
       objtype:
       freqfilename:
       domainfile:
       compcount:
       outfolder:
    Returns:
    """
    assert method in ["fd","pcfd","fdadd"] and objtype in ["frobenius","trace"]
    freqmat,nodenames = EmbedUtilities.readFreqFile(freqfilename)
    tdomains = EmbedUtilities.readDomainFile(domainfile)
    domains = [(start-1,end-1) for start, end in tdomains]
    if method == "fd" and objtype == "frobenius":
       comp2dominds,comp2scale,metadata = fdBiSub.runBisub(freqmat,domains,compcount)
    elif method == "pcfd" and objtype == "frobenius":
       pass
    elif method == "fdadd" and objtype == "frobenius":
       pass
           
    if not os.path.exists(outfolder):
       os.makedirs(outfolder) 
    deconoutfile = "{0}/decon.txt".format(outfolder,method,objtype,freqfilename,domainfile,compcount)
    objoutfile = "{0}/objscore.txt".format(outfolder,method,objtype,freqfilename,domainfile,compcount)
    print "Writing output"
    EmbedUtilities.writeDeconOut(comp2dominds,comp2scale,tdomains,deconoutfile)
    EmbedUtilities.writeDeconMeta(comp2dominds,comp2scale,tdomains,metadata,objoutfile)
     
     
def main():
    """
    """
    method, objtype, freqfilename, domainfile, enscount = sys.argv[1:6]
    runDeconvolution(method,objtype,freqfilename,domainfile,enscount)
        
        
