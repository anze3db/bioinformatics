import re
import networkx as nx
import pylab
from collections import defaultdict
if __name__ == '__main__':
    database = open('morbidmap').readlines()
    
    # 1.Create the network
    
    G = nx.Graph()
    ngenes = defaultdict(int)
    nclusters = defaultdict(int)
    for d in database:
        fields = d.strip().split('|')
        start = re.search('[A-Za-z]', fields[0]).start()
        end = re.search('[A-Za-z \-]*', fields[0][start:]).end()
        name = fields[0][start:start+end].strip()
        genes = fields[1].split(', ')
        
        G.add_node(name)
        ngenes[name] += len(genes)
        nclusters[name] += 1 # ???
        
    # 2. Analyze the network [TODO]
        
    # 3. Clustering [TODO]
    
    # 4. Extract clusters [TODO]
    
    G.add_edge("Diabetes", "Obesity")
    G.add_edge("Diabetes", "Leukemia")
    G.add_edge("Leukemia", "Obesity")
    G.add_edge("Leukemia", "Breast cancer")
    pos = nx.spring_layout(G, iterations=1)
    
    nodes = G.nodes() #fix node positions
    nx.draw_networkx_nodes(G, pos, nodes,
        node_size = [ 50*ngenes[a] for a in nodes],
        node_color = [ nclusters[a] for a in nodes ],
        linewidths = 0.1,
        alpha=0.4) #just because the colors are dark
    # nx.draw_networkx_labels(G, pos)
    nx.draw_networkx_edges(G, pos, alpha=0.2)
    pylab.axis("off")
    pylab.show()
    #pylab.savefig("nicer.pdf")