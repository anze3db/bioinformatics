from IPython.core.display import Math
from collections import defaultdict, Counter
from math import log
from numpy import median
import math
import networkx as nx
import numpy
import pylab
import random
import re
import sys



if __name__ == '__main__':
    
    def draw_graph(cluster=None):
        nodes = G.nodes() if cluster == None else [n for n in G.nodes() if cluster_labels[n] == cluster_labels[cluster]]
        
        
        
        nx.draw_networkx_nodes(G, pos, nodes,
            node_size=[ (ngenes[a]) for a in nodes],
            node_color=[ nclusters[a] for a in nodes],
            linewidths=1,
            alpha=0.4)  # just because the colors are dark
        
        labels = shown_labels if cluster == None else dict((n,n) for n in G.nodes() if cluster_labels[n] == cluster_labels[cluster])
        
        nx.draw_networkx_labels(G, pos, labels)
        
        edges = G.edges() if cluster==None else [e for e in G.edges() if cluster_labels[e[0]] == cluster_labels[cluster] and cluster_labels[e[1]] == cluster_labels[cluster]]
        nx.draw_networkx_edges(G, pos, edges, alpha=0.2, width=0.3)
        pylab.axis("off")
        pylab.show()
        
    database = open('morbidmap').readlines()
    
    # 1.Create the network
    
    G = nx.Graph()
    genes = defaultdict(set)
    for d in database[:2000]:
        fields = d.strip().split('|')
        start = re.search('[A-Za-z]', fields[0]).start()
        end = re.search('((\-[A-Za-z])|([A-Za-z ]))*', fields[0][start:]).end()
        name = fields[0][start:start + end].strip()
        genes[name] = genes[name].union(set(fields[1].split(', ')))
        
        if name not in G:
            G.add_node(name, new_name=name)

        for n in genes:
            if n == name: continue
            if len(genes[name].intersection(genes[n])) > 0:
                G.add_edge(name, n)
        
    # 2. Analyze the network [TODO]
        
    conn_comp = [len(c) for c in nx.connected_component_subgraphs(G)]
    
#    pylab.bar(range(len(conn_comp)), conn_comp)
#    pylab.figure(1)
#    pylab.yscale('log')
#    pylab.xlabel("Distribution of sizes of connected components")
#    pylab.ylabel("TODO")
#    pylab.title("TODO")
    # pylab.show()
        
#    hist = nx.degree_histogram(G)
#    pylab.figure(2)
#    pylab.subplot(121)
#    pylab.bar(range(len(hist)), hist, lw=1)
#    pylab.xlabel("Degree distribution of the network ")
#    pylab.ylabel("TODO")
#    pylab.title("TODO")
#    
#    pylab.subplot(122)
#    pylab.bar(range(len(hist)), [log(h+1) for h in hist])
#
#    pylab.xlabel("Degree distribution of the network (log)")
#    pylab.ylabel("TODO")
#    pylab.title("TODO")
    # pylab.show()
    
    G = nx.connected_component_subgraphs(G)[0]
         
    ngenes = {n:len(genes[n]) for n in genes}  # Normalize number of genes

    # 3. Clustering
    
    cluster_labels = {g:g for g in G.nodes()}
    
    for i in range(100):
        nodes = G.nodes()
        random.shuffle(nodes)
        for n in nodes:
            c = Counter(cluster_labels[n] for n in G.neighbors(n))
            if len(c.most_common(1)) > 0: 
                cluster_labels[n] = c.most_common(1)[0][0] 
    unique_labels = list(set(cluster_labels[l] for l in cluster_labels))
    
    pos = nx.spring_layout(G, iterations=500)
    
    # set cluster colors:
    nclusters = {l:unique_labels.index(cluster_labels[l]) for l in cluster_labels}
    
    # labels to be shown:
    shown_labels = dict((n,n) for n in G.nodes() if ngenes[n] == max(ngenes[g] for g in ngenes if g in cluster_labels and cluster_labels[g] == cluster_labels[n]))

    draw_graph()
    # 4. Extract clusters
    draw_graph('Deafness')
    draw_graph('Breast cancer')
    draw_graph('Diabetes mellitus')
    
