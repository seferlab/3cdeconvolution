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
    touchmat = np.zeros((varcount,varcount),dtype=np.int)
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
            touchmat[var2index[(domin1,comp1,scale1)],var2index[(domin2,comp2,scale2)]] += 1
            touchmat[var2index[(domin2,comp2,scale2)],var2index[(domin1,comp1,scale1)]] += 1
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
            touchmat[var2index[(domin,comp1,scale1)],var2index[(domin,comp2,scale2)]] += 1
            touchmat[var2index[(domin,comp2,scale2)],var2index[(domin,comp1,scale1)]] += 1
        for comp,scale in pairs:
            #coefmat[var2index[(domin,comp,scale)],var2index[(domin,comp,scale)]] += -2.0*scale*fsum 
            bvec[var2index[(domin,comp,scale)]] += -2.0*scale*fsum
    print "touch info"
    print np.amax(touchmat)
    print np.shape(touchmat)
    touchcount = {0:0,1:0,2:0}
    for ind1 in xrange(np.shape(touchmat)[0]):
        for ind2 in xrange(np.shape(touchmat)[1]):
            #if touchmat[ind1,ind2] == 1:
               #print ind1,ind2
               #exit(1)
            touchcount[touchmat[ind1,ind2]] += 1
    print touchcount
    exit(1)                   
    return coefmat, index2var


 dom2val = {domains[domin]: sum([xdict[(tdom,comp,scale)] for tdom,comp,scale in xdict.keys() if tdom == domin]+[0.0]) for domin in xrange(len(domains))}
       comp2val = {comp: sum([xdict[(dom,tcomp,scale)] for dom,tcomp,scale in xdict.keys() if tcomp == comp]+[0.0]) for comp in xrange(compcount)}
       
       print sum(xdict.values())
       print "dom2val print"
       for dom in dom2val.keys():
           print dom,dom2val[dom]
       print sum(dom2val.values())    
       print "comp val:"    
       for comp in comp2val.keys():
           print comp,comp2val[comp]
       domcom2val = {(domin,comp): sum([xdict[(dom,tcomp,scale)] for dom,tcomp,scale in xdict.keys() if tcomp==comp and dom==domin]+[0.0]) for domin in xrange(len(domains)) for comp in xrange(compcount)}
       print domcom2val
    
       for domin,comp in domcom2val.keys():
           print domin,comp,domcom2val[(domin,comp)]
       com2dom2val = {com:{} for com in xrange(compcount)}
       for domin,com in domcom2val.keys():
           com2dom2val[com][domin] = domcom2val[(domin,com)]
       print "info"
       for com in com2dom2val.keys():
           print "comp: ",com, sum(com2dom2val[com].values()) 
           for domin in com2dom2val[com].keys():
               print domin,com2dom2val[com][domin]
