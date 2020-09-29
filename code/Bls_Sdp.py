#SDP Formulation of Binary Least Squares
import cvxopt
import numpy as np
import os
import sys
import scipy as sp
from scipy.io import loadmat


def gaussRound(Xmat,objval,P):  
    """randomized gauss rounding
    Args:
       Xmat:
       objval:
       P: coef P matrix
    Returns:
       binsol:
    """
    L = 60
    print "info"
    print Xmat
    print type(Xmat)
    print np.shape(Xmat)
    samples = np.random.multivariate_normal([0.0 for ind in xrange(np.shape(Xmat)[0])],Xmat,L)
    binsamples = np.zeros(np.shape(samples),dtype=np.float)
    for samind in xrange(np.shape(samples)[0]):
        binsamples[samind,:] = np.array([1 if item >= 0 else -1 for item in samples[samind,:]])   
    minsol,minind = 1000000000000.0, None
    for ind in xrange(L):
        res = np.dot(binsamples[ind,:],P)
        res = np.dot(res,binsamples[ind,:].transpose())
        if res < minsol:
           minsol = res
           minind = ind
    vec = binsamples[minind,:]
    assert len(vec) == np.shape(binsamples)[1]
    binsol = [vec[ind]*vec[-1] for ind in xrange(len(vec)-1)]
    return binsol

def genSDPAOut(quadmat,outfilename):
    """outputs sdpa format for optimization
    Args:
       quadmat:
       outfilename:
    Returns:
    """
    dim = np.shape(quadmat)[0]
    parts = ["{0}".format(dim), "{0} ".format(1), "{0} ".format(dim), " ".join("1" for index in xrange(dim))]
    parts.extend(["{0} {1} {2} {3} {4}".format(0,1,ind1+1,ind2+1,-1.0*quadmat[ind1,ind2]) for ind1 in xrange(dim) for ind2 in xrange(ind1,dim) if quadmat[ind1,ind2] !=0])
    parts.extend(["{0} {1} {2} {3} {4}".format(ind,1,ind,ind,1.0) for ind in xrange(1,dim+1)])
    with open(outfilename,"w") as outfile:
        outfile.write("\n".join(parts)+"\n")
    
def genMatlabScript(mscriptfile,sdpfile,outdatafile,SDPTDIR,curdir):
    """generates matlab script
    Args:
       mscriptfile: matlab script filename
       sdpfile: sdpfile to be run
       outdatafile: output save file
       SDPTDIR:
       curdir: current directory to be returned
    Returns:
    """
    with open(mscriptfile,"w") as outfile:
        outfile.write("cd('{0}')\n".format(SDPTDIR))
        outfile.write("startup\n") 
        outfile.write("[blk,At,C,b] = read_sdpa(\'{0}/{1}\');\n".format(curdir,sdpfile))
        outfile.write("[obj,X,y,Z] = sqlp(blk,At,C,b);\n")
        outfile.write("Xmat = cell2mat(X);\n")
        outfile.write("cd('{0}')\n".format(curdir))
        outfile.write("save(\'{0}\',\'Xmat\',\'obj\')\n".format(outdatafile))

def runSDPRelax(quadmat):
    """runs sdp relaxation via sdpt3
    Args:
       quadmat:
    Returns:
       X: cholesky decomposed matrix
       obj: obj value
    """
    sdptfilename = "tempsdpt"
    sdptoutfile = "sdptoutput.mat"
    mscriptfile = "sdpt.m"
    genSDPAOut(quadmat,sdptfilename)
    olddir = os.getcwd()
    genMatlabScript(mscriptfile,sdptfilename,sdptoutfile,SDPTDIR,olddir)
    code = "{0} -nodisplay -nosplash -nodesktop -r \" {1}; quit; \" ".format(MATLABPATH,mscriptfile.replace(".m",""))
    os.system(code)
    data = loadmat(sdptoutfile)
    if type(data["Xmat"]) == sp.sparse.csc.csc_matrix:
       return data["Xmat"].todense(), data["obj"][0,0]
    else:
       return data["Xmat"], data["obj"][0,0]  


def coefs2Boolean(coefmat):
    """converts matrix and bfor 0,1 to -1,1 case
    Args:
       coefmat:
    Returns:
       modmat: modified matrix
       val: additional value
    """
    modmat = np.array(coefmat)
    for ind1 in xrange(np.shape(coefmat)[0]):
        for ind2 in xrange(np.shape(coefmat)[1]):
            modmat[ind1,ind2] += coefmat[ind1,ind2]/4.0
            modmat[ind1,ind1] += coefmat[ind1,ind2]/4.0
            modmat[ind2,ind2] += coefmat[ind1,ind2]/4.0
    val = 0.25 * np.shape(coefmat)[0] * np.shape(coefmat)[1]
    return modmat,val


TESTMODE = True
MATLABPATH = "/Applications/MATLAB_R2013a.app/bin/matlab"
SDPTDIR = "/Users/esefer/Downloads/SDPT3-4.0"
def runSDP(comp2dominds,domains,freqmat,compcount):
    """runs sdp relaxation
    Args:
       coefmat: coefficient matrix already squared
       b:
       freqmat:
    Returns:
       sol:
    """
    usecoefmat = deepcopy(coefmat)
    assert len(b) == np.shape(usecoefmat)[0]
    for ind in xrange(len(b)):
        usecoefmat[ind,ind] += b[ind]
    #c = freqmat.reshape(np.shape(freqmat)[0]*np.shape(freqmat)[1],1)
    freqsqsum = sum([freqmat[ind1,ind2]**2 for ind1 in xrange(np.shape(freqmat)[0]) for ind2 in xrange(np.shape(freqmat)[1])])
    modmat, val = coefs2Boolean(usecoefmat)
    freqsum += val                                 
    
    P = np.append(coefmat,modb.transpose(),0)
    P = np.append(P,np.append(modb,[[freqsqsum]],0),1)
    
    fracarr,obj = runSDPRelax(modmat)
    bindict = gaussRound(fracarr)
    if TESTMODE:
       SDPTest() 
    exit(1)

