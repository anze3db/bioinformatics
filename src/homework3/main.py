from collections import defaultdict
from scipy.cluster.hierarchy import dendrogram, linkage
from itertools import chain

animals = [{
        "id": "NC_002008.4",
        "name": "Gray Wolf",
    }, {
        "id" : "NC_006580.1",
        "name" : "Goldfish"
    }, {
        "id" : "NC_012420.1",
        "name" : "Veiled Chameleon"
    }, {
        "id" : "NC_011391.1",
        "name" : "Daboia"
    }, {
        "id" : "NC_012061.1",
        "name" : "Dolphin"
    }, {
        "id" : "NC_001640.1",
        "name" : "Horse"
    }, {
        "id" : "NC_001645.1",
        "name" : "Gorilla"
    }, {
        "id" : "NC_012920.1",
        "name" : "Human"
    }, {
        "id" : "NC_011137.1",
        "name" : "Neanderthal"
    }, {
        "id" : "NC_001643.1",
        "name" : "Chimpanzee"
    }, {
        "id" : "NC_002083.1",
        "name" : "Orangutan"
    }, {
        "id" : "NC_001665.2",
        "name" : "Rat"
    }, {
        "id" : "NC_014692.1",
        "name" : "Boar"
    }, {
        "id" : "NC_004299.1",
        "name" : "Pufferfish"
}]

def blosum50():
    """Returns the dict containing BLOSUM50 substitution matrix."""
    
    b50 = open('b50.table').readlines()
    # Two horrible lines that turn the above string into a dict. I am sorry.
    arr = [j.split(' ') for j in [i.strip().replace('  ', ' ') for i in b50]]
    return reduce(lambda a,b: dict(a.items() + b.items()), [{(arr[0][j],arr[i][0]): int(arr[i][j]) for j in xrange(1, len(arr[i]))} for i in xrange(1,len(arr))])

def load_entrez(genome):
    """Reads Entrez data from a .pickle file. If the file does not exist it is created"""
    import os.path
    from cPickle import dump, load
    if not os.path.isfile(genome + '.pickle'):
        from Bio import Entrez
        from Bio import SeqIO
        handle = Entrez.efetch(db="nucleotide", rettype="gb", id=genome, email="smotko@smotko.si")
        rec = SeqIO.read(handle, "gb")
        handle.close()
        dump(rec, file(genome + '.pickle', 'w'))
    
    return load(file(genome + '.pickle'))


def cox3():
    c = []
    for a in animals:
        rec = load_entrez(a["id"])
        for f in rec.features:
            if f.type == "CDS":
                if f.qualifiers['gene'][0] == "COX3":
                    c.append(f.qualifiers["translation"][0])
    return c



def dpt(s,t):
    def cost():
        M[i,j] = max(
            M[i-1,j]+gp,
            M[i,j-1]+gp,
            M[i-1,j-1] + blosum[si, tj]
        )
        
    def printTable():
        for i in range(len(s)):
            print ' '.join(["%5d" % M[i,j] for j in range(len(t))])

    M = defaultdict(int)
    gp = -5
    
    for i in range(len(s)):
        M[i,-1] = M[i-1,-1] + gp
    for j in range(len(t)):
        M[-1,j] = M[-1,j-1] + gp

    [cost() for i,si in enumerate(s) for j,tj in enumerate(t)]
    # printTable()
    return M[len(s)-1, len(t)-1]
    
if __name__ == '__main__':
    cox = cox3()
    blosum = blosum50()
    
    dist = list(chain.from_iterable([[dpt(cox[i], cox[j]) for j in range(i+1, len(cox))] for i in range(len(cox))]))

    import pylab
    dendrogram(linkage(dist, 'average'), 15, labels = [a['name'] for a in animals]) 
    pylab.show()
    