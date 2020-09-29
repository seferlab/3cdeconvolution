
def findCliqueDecomp(interG):
    """finds maximal clique decomposition 
    Args:
       interG
    Returns:
       cliques:
    """
    tinterG = nx.Graph(interG)
    return list(map(list,nx.find_cliques(tinterG)))
     
    cliques = []
    while tinterG.number_of_nodes() > 0:
        maxclique,maxsize = None, None
        allcliques = list(nx.find_cliques(tinterG))
        for tclique in allcliques:
            if len(tclique) >= maxsize:
               maxclique = tclique
        for ind1 in xrange(len(maxclique)):
            node1 = maxclique[ind1]
            for ind2 in xrange(ind1+1,len(maxclique)):
                node2 = maxclique[ind2]
                tinterG.remove_edge(node1,node2)
        for node in tinterG.nodes():
            if len(tinterG.neighbors(node)) == 0:
               tinterG.remove_node(node)         
        cliques.append(maxclique)   
    return cliques

def sparsifyIG(domains):
    """sparsifies domain interaction graph 
    Args:
       domains:
    Returns:
       mindomains: min set of representative domains
    """
    include = lambda (s1,e1), (s2,e2): True if (s1 <= s2 and e2 <= e1) else False
    len2domain = {1:set()}
    for dom in domains:
        len2domain.setdefault(dom[1] - dom[0] + 1,set()).add(dom)
    G = nx.DiGraph()
    for in1 in xrange(len(domains)):
        s1,e1 = domains[in1]
        G.add_node((s1,e1))
        for in2 in xrange(in1+1,len(domains)):
            s2,e2 = domains[in2]
            if include((s1,e1),(s2,e2)):
               G.add_edge((s1,e1),(s2,e2))
            if include((s2,e2),(s1,e1)):
               G.add_edge((s2,e2),(s1,e1))   
    mindomains = len2domain[1]
    seenlens = len2domain.keys()
    seenlens.remove(1)
    for domlen in sorted(seenlens):
        for dom in len2domain[domlen]:
            nodes = range(dom[0],dom[1]+1)
            preset = set(G.successors(dom)).intersection(mindomains)
            for in1 in xrange(len(preset)):
                for in2 in xrange(in1+1,len(preset)):
                    pass  
            for neighdom in set(G.successors(dom)).intersection(mindomains):
                nodes.remove(range(neighdom[0],neighdom[1]+1))
            if len(nodes) != 0:
               mindomains.add(dom)
    print mindomains
    exit(1)                        
    return mindomains

#comp2dominds,comp2scale = fdLocal.runLocal(freqmat,domains,scales,compcount)
    # comp2dominds,comp2scale = fdILP.runILP(freqmat,domains,scales,compcount,{"roundmethod":"greedy2"})    
    #elif method == "pcfd" and objtype == "frobenius":
       #import random 
       #p = random.random()
       #if p < 0.2: 
       #   params = {"roundmethod": "hyper", "scalepen": 1.0, "overpen": 1.0}
    #   comp2doms,comp2scale = pcfdILP.runSDP(freqmat,domains,scales,compcount,params)
    #if method == "fd" and objtype == "trace":
    #   comp2doms,comp2scale = fdILP.runILP(freqmat,domains,scales,compcount,"random") 
