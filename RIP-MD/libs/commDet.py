from functools import reduce
from random import shuffle, choice
import numpy as np
import networkx as nx
from scipy.cluster import hierarchy


defaultWeights = {"Alpha Carbon Distance":1,"Salt Bridge":1,"Arg-Arg":1,"Hydrogen Bond":1,"Disulfide Bridge":1,"Cation-Pi":1,"Pi-Pi":1,"Coulomb":1,"Van der Waals":1}


######################################################################
##
## Filters the graph according to the specified set of
## interactions and time threshold. It also turns multigraphs
## into graphs and removes loops.
##
## Consider that the time threshold only makes sense for
## graphs from MD trajectories and when it's greater or equal
## than the time threshold used when running RIP-MD, which is
## by default 75% but can be modified by the user.
##
## It add the "weight" edge attribute which is used for weighted
## versions of the algorithms. By default all edges have weight
## equal to 1. weightDict is a dictionary whose keys are the
## interaction names and the values are the corresponding weight
## for that interaction type. These values must be strictly positive.
## If you need a value of 0 then just filter that interaction with
## interactions.
##
######################################################################
def graphFilter(grafo,interactions=None,time=75,removeIsolated=False,weightDict=None):
    if not weightDict:
        weightDict = defaultWeights
    else:
        for key in defaultWeights.keys():
            if key not in weightDict.keys():
                weightDict[key] = defaultWeights[key]

    G = nx.Graph()
    G.add_nodes_from(grafo.nodes(data=False))

    if nx.get_edge_attributes(grafo, "time"):  # If graph comes from an MD trajectory.
        if interactions:
            for edge in grafo.edges(data=True):
                if edge[0] != edge[1] and float(edge[2]['time'][:-1]) >= time and edge[2][
                    'interaction'] in interactions:
                    if (edge[0], edge[1]) in G:
                        G.edges[edge[0], edge[1]]["weight"] += weightDict[edge[2]['interaction']]
                    else:
                        G.add_edge(edge[0], edge[1], weight=weightDict[edge[2]['interaction']])
        else:
            for edge in grafo.edges(data=True):
                if edge[0] != edge[1] and float(edge[2]['time'][:-1]) >= time:
                    if (edge[0], edge[1]) in G:
                        G.edges[edge[0], edge[1]]["weight"] += weightDict[edge[2]['interaction']]
                    else:
                        G.add_edge(edge[0], edge[1], weight=weightDict[edge[2]['interaction']])
    else:  # If graph comes from a static PDB.
        if interactions:
            for edge in grafo.edges(data=True):
                if edge[0] != edge[1] and edge[2]['interaction'] in interactions:
                    if (edge[0], edge[1]) in G:
                        G.edges[edge[0], edge[1]]["weight"] += weightDict[edge[2]['interaction']]
                    else:
                        G.add_edge(edge[0], edge[1], weight=weightDict[edge[2]['interaction']])
        else:
            for edge in grafo.edges(data=True):
                if edge[0] != edge[1]:
                    if (edge[0], edge[1]) in G:
                        G.edges[edge[0], edge[1]]["weight"] += weightDict[edge[2]['interaction']]
                    else:
                        G.add_edge(edge[0], edge[1], weight=weightDict[edge[2]['interaction']])

    if removeIsolated:
        listToRemove = []
        for node in G.nodes():
            if G.degree[node] == 0:
                listToRemove.append(node)
        G.remove_nodes_from(listToRemove)

    return G

##################################################################
##
## Uses a divisive algorithm to detect communities
## as the one introduced by Girvan and Newman in
## https://doi.org/10.1073/pnas.122653799
##
## Returns a list of tuples with each edge and its ecb,
## in order of removal. Also returns nodes (set) of the
## filtered graph which also excludes isolated nodes.
##
## Different weights per interaction type can be set with the
## the weights parameter. Its structure details are specified
## in the graphFilter documentation.
##
## Consider that weights=None is not equivalent to
## weights=defaultWeights. The first works as if the algorithm
## were run in its original unweighted version, while the second
## works with each edge (u,v) having weight equal to the number
## of edges (u,v,c) in the input multigraph.
##
##################################################################
def girvanNewman(grafo,interactions=None,time=75,weights=None):
    # Filtering graph
    G = graphFilter(grafo, interactions=interactions, time=time, removeIsolated=True,weightDict=weights)
    if not weights:
        nx.set_edge_attributes(G, 1, "weight")

    etrDict={} # Dict of (edge, ebc) in order of removal.
    # Each iteration removes the edge with the highest betweenness centrality.
    try:
        while len(G.edges) > 0:
            ebc = nx.edge_betweenness_centrality(G,normalized=True,weight="weight")

            m = 0
            edge_to_remove = None
            for item in ebc.items():
                if item[1] > m:
                    m = item[1]
                    edge_to_remove = item[0]
            etrDict[edge_to_remove]=m
            G.remove_edge(*edge_to_remove)
    except TypeError:
        nx.write_gml(G,"weird.gml")

    return etrDict, set(G.nodes)

######################################################################
##
## Takes the nodes of a graph and a list of its edges in order
## of removal. Generates the components after each removal,
## starting from 0 removals. It expects the output of girvanNewman.
##
######################################################################
def getComponentsPerIteration(etrDict,nodes):
    G = nx.Graph()
    G.add_nodes_from(nodes)
    G.add_edges_from(etrDict.keys())
    yield list(nx.connected_components(G))
    for edge in etrDict.keys():
        G.remove_edge(*edge)
        yield list(nx.connected_components(G))

############################################################
##
## Takes the output of girvanNewman and finds the best
## division of the graph in k communities according to
## such output.
##
## k is an integer between the number of nodes and the
## number of connected components of the graph, inclusive.
## If it's outside that interval, it will be set to the
## closest interval limit.
##
## Returns the corresponding graph whose nodes include
## the label "community", outOfBounds which tells if k
## was too low ("cc") or too high ("nodes") and the final
## value of k.
##
############################################################
def findKCommGN(etrDict,nodes,k):
    # Finding k communities with etrList and nodes
    commGen = getComponentsPerIteration(etrDict, nodes)
    CCs = next(commGen)
    outOfBounds=""
    if k < len(CCs):
        k = len(CCs)
        comms = CCs
        outOfBounds = "cc"
    elif k > len(nodes):
        k = len(nodes)
        comms = nodes
        outOfBounds = "nodes"
    else:
        while len(CCs) < k:
            CCs = next(commGen)
        comms = CCs
    comms = list(comms)

    # Setting up graph
    G = nx.Graph()
    G.add_nodes_from(nodes)
    G.add_edges_from(etrDict.keys())

    # Setting "community" label
    for i in range(len(comms)):
        for node in comms[i]:
            G.nodes[node]["community"] = i

    return G,outOfBounds,k

##########################################################
##
## Given a graph G with a community structure determined
## by attribute label, returns its modularity measure as
## the one defined by Newman and Girvan in
## https://doi.org/10.1103/PhysRevE.69.026113
##
## Note that the graph needs to have at least one edge
## or the returned value will be most likely nan.
##
##########################################################
def modularity(G,label="community"):
    communityDict = {}
    k=0
    for node in G.nodes(data=label):
        if node[1] not in communityDict.keys():
            communityDict[node[1]] = k
            k += 1

    E = np.zeros((k,k))
    for edge in G.edges(data=False):
        i = communityDict[G.nodes[edge[0]][label]]
        j = communityDict[G.nodes[edge[1]][label]]
        if i==j: # fixed after sending final
            E[i, j] += 1
        else:
            E[i,j] += 1
            E[j,i] += 1
    E = E / len(G.edges)

    return np.trace(E) - np.sum(E @ E)

############################################################
##
## Takes the output of setupDendrogram and finds a
## division of the graph with the highest modularity
## measure.
##
## Returns the corresponding graph whose nodes include the
## label "community", the respective modularity value and
## the number of communities in the resulting graph.
##
############################################################
def findBestModGN(etrDict,nodes):
    Z, leafLabels,colors,stopColoring = setupDendrogram(etrDict, nodes, jump=lambda x: 0, colored=None)

    nodeCard = len(leafLabels)
    G = nx.Graph()
    for i in range(nodeCard):
        G.add_node(leafLabels[i],community=i)
    G.add_edges_from(etrDict.keys())

    mod = modularity(G)
    CCs=dict([(i,{i}) for i in range(nodeCard)])
    comms = dict(CCs)
    for i in range(len(Z)):
        CCs[i+nodeCard] = CCs[Z[i][0]] | CCs[Z[i][1]]
        del CCs[Z[i][0]]
        del CCs[Z[i][1]]

        for nodeid in CCs[i+nodeCard]:
            G.nodes[leafLabels[nodeid]]["community"] = i+nodeCard

        nmod = modularity(G)
        if nmod >= mod:
            comms = dict(CCs)
            mod = nmod

    for cc in comms.items():
        for nodeid in cc[1]:
            G.nodes[leafLabels[nodeid]]["community"] = cc[0]

    return G,mod,len(comms)

##############################################################################
##
## Sets up the preliminaries to plot a dendrogram with the output
## of girvanNewman.
##
## jump is a function that transforms the ecb of each edge
## to the 'height jump' of the link generated by that edge.
## Suggested jump functions:
##     lambda x: x # Each jump is the corresponding ecb. (linear)
##     lambda x: np.sqrt(x) # Helps smoothe the jumps. (sqrt)
##     lambda x: 1 # Gives homogeneous jumps. (constant)
##
## threshold can be a float in (0,1) or an integer.
##
## If threshold is an integer between the number of nodes and
## the number of connected components of the graph, inclusive,
## then it will find threshold communities. If it's outside that
## interval, it will be set to the closest interval limit.
##
## If threshold is a float between 0 and 1, it will consider as
## a community any cluster whose leading link has a height lower
## than threshold*highestLink.
##
## If colored is set to 'simple', the clusters determined by the
## above criteria will be highlighted in red, in contrast to the
## rest of the dendrogram which is blue. If it is set to 'colorful',
## it will highlight them in varying colors. If it is set to 'alternate',
## they will be highlighted in red and orange alternatingly.
##
##############################################################################
def setupDendrogram(etrDict,nodes,jump=lambda x: np.sqrt(x),colored=None,threshold=0.7):
    nodeCard = len(nodes)
    etaDict = dict(reversed(etrDict.items()))
    leadClID = {-1: set(nodes)}
    father = {}
    children = {}
    Z = []
    leafLabels = []
    k = 0
    n = 0
    h = 0.0
    for edge in etaDict.keys():
        ID0 = None
        ID1 = None
        for key in leadClID.keys():
            if edge[0] in leadClID[key]:
                ID0 = key
            if edge[1] in leadClID[key]:
                ID1 = key
            if not (ID0 is None or ID1 is None):
                break
        if ID0 == -1:
            ID0 = n
            leafLabels.append(edge[0])
            children[n] = {n}

            leadClID[-1] = leadClID[-1] - {edge[0]}
            leadClID[n] = {edge[0]}
            n += 1

        if ID1 == -1:
            ID1 = n
            leafLabels.append(edge[1])
            children[n] = {n}

            leadClID[-1] = leadClID[-1] - {edge[1]}
            leadClID[n] = {edge[1]}

            n += 1
        if ID0 != ID1:
            leadClID[nodeCard + k] = leadClID[ID0] | leadClID[ID1]
            del leadClID[ID0]
            del leadClID[ID1]
            h = h + float(jump(etaDict[edge]))
            Z.append([ID0, ID1, h, len(leadClID[nodeCard + k])])

            father[ID0] = nodeCard + k
            father[ID1] = nodeCard + k
            children[nodeCard + k] = children[ID0] | children[ID1] | {nodeCard + k}

            k += 1

    # Apparently hierarchy.dendrogram can only make dendrograms for connected graphs. To solve this, this will add
    # 'dummy links' which will be transparent so they don't appear in the plot.
    stopColoring = nodeCard + k
    while len(Z) < nodeCard - 1:
        i = None
        j = None
        keys = iter(leadClID.keys())
        for l in keys:
            if l != -1:
                i = l
                j = next(keys)
                if j != -1:
                    break
                else:
                    j = next(keys)
                    break
        leadClID[nodeCard + k] = leadClID[i] | leadClID[j]
        del leadClID[i]
        del leadClID[j]
        Z.append([i, j, h, len(leadClID[nodeCard + k])])
        k += 1

    # Coloring
    stdColor = 'C0'

    colors = dict((u + nodeCard, stdColor) for u in range(stopColoring - nodeCard))

    if colored:
        H = hierarchy.dendrogram(Z, no_plot=True, get_leaves=True)['leaves']

        cindex = iter(range(stopColoring))
        palette = (stdColor,)
        if colored == 'colorful':
            palette = ('C1','C2','C3','C4','C5','C6','C7','C8','C9')
        elif colored == 'simple':
            palette = ('red',)
        elif colored == 'alternate':
            palette = ('red', 'orange')

        if threshold < 1:
            threshold = float(threshold) * h
            i = 0
            while i < nodeCard:
                grandpa = None
                u = H[i]
                if Z[father[u] - nodeCard][2] < threshold:
                    u = father[u]
                    try:
                        while Z[father[u] - nodeCard][2] < threshold:
                            u = father[u]
                    except KeyError:
                        pass
                    grandpa = u

                if grandpa:
                    color = palette[next(cindex) % len(palette)]
                    for v in children[grandpa]:
                        colors[v] = color
                    i += Z[grandpa - nodeCard][3]
                else:
                    i += 1
        else:
            threshold = int(threshold)
            if threshold < 2 * nodeCard - stopColoring:
                threshold = 2 * nodeCard - stopColoring
            elif threshold > nodeCard:
                threshold = nodeCard
            thresholdGrandpas = reduce(lambda x, y: (x | {Z[y][0], Z[y][1]}) - {y + nodeCard}, range(nodeCard - 2, nodeCard - threshold - 1, -1), {2*nodeCard-2})

            i = 0
            while i < nodeCard:
                u = H[i]
                while u not in thresholdGrandpas:
                    u = father[u]
                grandpa = u

                color = palette[next(cindex) % len(palette)]
                for v in children[grandpa]:
                    colors[v] = color

                if grandpa < nodeCard:
                    i += 1
                else:
                    i += Z[grandpa - nodeCard][3]

    return Z, leafLabels, colors, stopColoring




######################################################################
##
## Runs a label propagation algorithm based on the
## one introduced by Raghavan, Albert and Kumara in
## https://doi.org/10.1103/PhysRevE.76.036106
##
## Which nodes will have initial labels can be specified
## with the centralNodes (set). The other nodes will
## start without a label and the algorithm will first
## ensure all nodes are labeled before checking its normal
## stop criterion. If additionally the force option is set
## then centralNodes will not change their label.
##
## Labels are initially by default each nodes' name. Using
## namesAsLabels=False they can be set to en enumeration instead.
##
## Different weights per interaction type can be set with the
## the weights parameter. Its structure details are specified
## in the graphFilter documentation.
##
## Consider that weights=None is not equivalent to
## weights=defaultWeights. The first works as if the algorithm
## were run in its original unweighted version, while the second
## works with each edge (u,v) having weight equal to the number
## of original edges (u,v,c) in the input multigraph.
##
## Returns a filtered graph whose nodes have the label "community".
##
######################################################################
def labelPropagation(grafo,interactions=None,time=75,centralNodes=None,force=False,namesAsLabels=True,removeIsolatedNodes=False,weights=None):
    # Filtering graph
    G = graphFilter(grafo, interactions=interactions, time=time, removeIsolated=removeIsolatedNodes, weightDict=weights)
    if not weights:
        nx.set_edge_attributes(G,1,"weight")

    if not namesAsLabels:
        nodeID = dict([(node,i) for i,node in enumerate(G.nodes)])

    # Runs algorithm for each connected component to ensure proper working of centralNodes. If a component has no
    # central node, then the algorithm starts labeling uniquely each node of that component, as is the 'normal' way.
    for component in nx.connected_components(G):
        listedNodes = list(component)
        if not centralNodes or len(centralNodes & component)==0:
            prelNodes = component
        else:
            prelNodes = centralNodes & component
            if force:  # Makes nodes in prelNodes static
                listedNodes = list(component - prelNodes)

        if namesAsLabels:
            for node in component:
                if node in prelNodes:
                    G.nodes[node]["community"] = node
                else:
                    G.nodes[node]["community"] = None
        else:
            for node in component:
                if node in prelNodes:
                    G.nodes[node]["community"] = nodeID[node]
                else:
                    G.nodes[node]["community"] = None

        flag=True
        while flag:
            flag=False
            shuffle(listedNodes)
            for u in listedNodes:
                if G.nodes[u]["community"] is None:
                    contDict = {None: 0}
                    for v in G[u]:
                        if G.nodes[v]["community"] in contDict.keys():
                            contDict[G.nodes[v]["community"]] += G.edges[u,v]["weight"]
                        else:
                            contDict[G.nodes[v]["community"]] = G.edges[u,v]["weight"]
                    if len(contDict.keys()) == 1:
                        continue

                    del contDict[None]
                    maxFreq = max(contDict.values())
                    G.nodes[u]["community"] = choice([label for label in contDict.keys() if contDict[label] == maxFreq])
                    flag = True
                else:
                    contDict = {G.nodes[u]["community"]: 0, None: 0}
                    for v in G[u]:
                        if G.nodes[v]["community"] in contDict.keys():
                            contDict[G.nodes[v]["community"]] += G.edges[u,v]["weight"]
                        else:
                            contDict[G.nodes[v]["community"]] = G.edges[u,v]["weight"]

                    del contDict[None]
                    maxFreq = max(contDict.values())
                    if contDict[G.nodes[u]["community"]] < maxFreq:
                        G.nodes[u]["community"] = choice([label for label in contDict.keys() if contDict[label] == maxFreq])
                        flag = True
    return G
