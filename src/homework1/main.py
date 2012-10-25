import pylab
from collections import Counter
from math import log
from matplotlib import cm

def get_seq():
    """Return  mitochondrial sequences in FASTA format"""

    f = open("NC_001416.1.fasta")
    fc = f.read()
    fc = fc.replace('N', '') # TODO: figure out how to handle N (any nucleic acid)
    return list("".join(fc.split("\n")[1:]))

def print_probability(cn):
    """Print ordered probabilities"""
    
    l = list(cn)
    l.sort(key=lambda c: cn.get(c), reverse=True) # Sort by number of occurrences
    s = ""
    for c in l:
        s += "%s - %0.2f%%, " % (c, cn.get(c) * 100 / float(sum(cn.values())))
    print s + "\n"
    
def report_changes(seq, window_size):
    """Report changes of the aggregate frequency of G or C"""
    
    window_frequence = []
    
    c = 0
    cnt = 0
    for s in seq:
        cnt += 1
        if s == "C" or s == "G":
            c += 1
            
        if cnt == window_size:
            window_frequence.append(c)
            cnt = 0
            c   = 0
            
    return window_frequence
    
def graph(window_frequence):        
    #pylab.bar(range(len(window_frequence)), window_frequence)
    pass
def cgr(seq):
    m = {'A': 0, 'T': 1, 'C': 0, 'G': 1}
    n = {'A': 0, 'C': 1, 'T': 0, 'G': 1}
    a,b = 0,0
    k = len(seq)
    for i in range(k):
        a += m[seq[i]] * 2**(k-i-1)
        b += n[seq[i]] * 2**(k-i-1)
    
    return a,b

if __name__ == '__main__':
    
    seq = get_seq()
    N = len(seq)
    
    # Calculate occurrences for 1-mers and 2-mers:
    c1 = Counter(seq)
    c2 = Counter([l1 + l2 for l1, l2 in zip(seq[:-1], seq[1:])])    

    print_probability(c1)
    print_probability(c2)

    # Calculate deviations:
    ls = []
    for c in c2:
        l = list(c)
        ls.append((c, log(N * c2.get(c) / float(c1[l[0]] * c1[l[1]]),2)))

    ls.sort(key=lambda c: c[1], reverse=True)
    s = ""
    for c in ls:
        s += "%s - %0.4f, " % (c[0], c[1])
    print s + "\n"
    
    #Report changes on aggregate frequencies:
    
    [graph(report_changes(seq, i)) for i in [10,100,1000]]
    # pylab.show() # TODO: Normalize
    
    # Fill up CGR table:
    seq.reverse()
    k = 4
    
    cgra = []
    for i in range(16):
        cgra.append([])
        for _ in range(16):
            cgra[i].append(0)
            
    for i in range(0, len(seq)-k+1):
        a,b = cgr(seq[i:i+k])
        cgra[a][b] += 1
    
    pylab.imshow(cgra, interpolation='nearest', cmap = cm.gray)
    pylab.show()
