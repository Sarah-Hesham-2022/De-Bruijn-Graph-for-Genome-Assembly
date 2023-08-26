import networkx as nx
import matplotlib.pyplot as plot 

#Exercise 1:
#➔ Write a function kmers that takes a string sequence
#and an integer k as input and returns a list of all the k-mers in the sequence.

def kmers(seq,k):

    kmers =[]
    if(k > len(seq)):
        return kmers
    for i in range(len(seq) - k + 1):
        kmers.append(seq[i:i+k])

    return kmers

#Test Function
#Result : ['AGC', 'GCG', 'CGG', 'GGT', 'GTA', 'TAT']
#print(kmers("AGCGGTAT",3))

#Exercise 2:
#➔ Write a function DeBruijnGraph that takes a list of
#reads and an integer k as input and returns the DeBruijn graph of the reads.

def DeBruijnGraph(reads,k):

    #Step 1 : get all k-mers of reads
    k_mers = []
    temp = []
    for i in reads:
        temp = kmers(i,k)
        k_mers.extend(temp)

    #Step 2 : get k-1-mer prefix nodes and k-1-mer suffix nodes connected by an edge (k-mer)
    edges = {}
    for i in k_mers:
        edges[i] = [i[0:k-1],i[1:]]
    
    #Step 3 : decide which nodes and edges to keep on graph
    graph = {}
    graph['nodes'] = []
    graph['edges'] = []
    graphEdges = {}
    for i in edges:

        if( (edges[i][0] not in graph['nodes']) and (edges[i][1] not in graph['nodes']) ):
           graph['nodes'].extend(edges[i][0])
           graph['nodes'].extend(edges[i][1])
           graph['edges'].append(i)
           graphEdges[i] = [edges[i][0],edges[i][1]]

        elif(edges[i][0] not in graph['nodes']):
            graph['nodes'].extend(edges[i][0])
            graph['edges'].append(i)
            graphEdges[i] = [edges[i][0],edges[i][1]]

        else:
            graph['nodes'].extend(edges[i][1])
            graph['edges'].append(i)
            graphEdges[i] = [edges[i][0],edges[i][1]]

    return graphEdges

#Test Function
#Result : {'GTTAC': ['GTTA', 'TTAC'], 'CCGTT': ['CCGT', 'CGTT'], 'GTTCG': ['GTTC', 'TTCG'], 'CGTTA': ['CGTT', 'GTTA'], 'ACGTT': ['ACGT', 'CGTT'], 'CGTTC': ['CGTT', 'GTTC'], 'TTACG': ['TTAC', 'TACG'], 'TACGT': ['TACG', 'ACGT'], 'TTCGA': ['TTCG', 'TCGA']}
#print(DeBruijnGraph(["TTACGTT", "CCGTTA", "GTTAC", "GTTCGA", "CGTTC"],5))

def Draw_DeBruijnGraph(graphEdges):

    dbGraph=nx.DiGraph()
    for i in graphEdges:
        dbGraph.add_edge(graphEdges[i][0],graphEdges[i][1])

    pos = nx.shell_layout(dbGraph)
    nx.draw_networkx(dbGraph, pos,with_labels=True, node_size=3000)

    for i in graphEdges:
        graphEdges[i] = tuple(graphEdges[i])
    graphNodes = dict(zip(graphEdges.values(), graphEdges.keys()))
    nx.draw_networkx_edge_labels(dbGraph, pos, edge_labels=graphNodes)
    plot.savefig("Draw_DeBruijnGraph.png")
    plot.show()
 
#Test Function    
Draw_DeBruijnGraph(DeBruijnGraph(["TTACGTT", "CCGTTA", "GTTAC", "GTTCGA", "CGTTC"],5))
Draw_DeBruijnGraph(DeBruijnGraph(["ATGG", "TGCC", "TAAT", "CCAT", "GGG", "GGATG", "ATGTT"],3))
Draw_DeBruijnGraph(DeBruijnGraph(["CGATTCTAAGT"],4))

