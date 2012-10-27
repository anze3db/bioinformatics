from collections import Counter
from math import log
#from matplotlib import cm
import pylab

def get_seq(fa):
    """Return  mitochondrial sequences in FASTA format"""

    f = open(fa)
    fc = f.read()
    fc = fc.replace('N', '')  # ignore aNy
    return list("".join(fc.split("\n")[1:]))

def print_probability(cn, N):
    """Print ordered probabilities"""
    
    s = ""
    for c in sorted(list(cn), key=lambda c: cn[1], reverse=True):
        s += "%s - %0.2f%%, " % (c[0], c[1] * 100 / float(N))
    print s[:-2] + "\n"
    
def report_changes(seq, window_size):
    """Report changes of aggregate frequency of G or C"""
    
    return [len([j for j in seq[i:i + window_size] if j == "C" or j == "G"]) for i, _ in enumerate(seq) if i % 50 == 0]
    
def graph(window_frequence):   
    """Plot a graph"""
      
    pylab.plot(range(len(window_frequence)), window_frequence, lw=1)
    
def cgr(seq):
    """Calculate the position in the cgr grid for the given k-mer seq"""
    
    
    m = {'A': 0, 'T': 1, 'C': 0, 'G': 1}
    n = {'A': 0, 'C': 1, 'T': 0, 'G': 1}
    a, b = 0, 0
    k = len(seq)
    for i in range(k):
        a += m[seq[i]] * 2 ** (k - i - 1)
        b += n[seq[i]] * 2 ** (k - i - 1)
    
    return a, b


def get_cgra(k):
    cgra = []
    for i in range(2 ** k):
        cgra.append([])
        for _ in range(2 ** k):
            cgra[i].append(0)
            
    for i in range(0, len(seq) - k + 1):
        a, b = cgr(seq[i:i + k])
        cgra[a][b] += 1
    return cgra

if __name__ == '__main__':
    
    seq = get_seq("NC_001416.1.fasta")
    N = len(seq)
    
    # Calculate occurrences for 1-mers and 2-mers:
    c1 = Counter(seq)
    c2 = Counter(l1 + l2 for l1, l2 in zip(seq[:-1], seq[1:]))

    print_probability(c1.items(), N)
    print_probability(c2.most_common(10), N)

    # Calculate deviations:
    ls = []

    for c in c2:
        l = list(c)
        ls.append((c, abs(log(N * c2.get(c) / float(c1[l[0]] * c1[l[1]]), 2))))

    ls.sort(key=lambda c: c[1], reverse=True)
    s = ""
    for c in ls:
        s += "%s - %0.4f, " % (c[0], c[1])
    print s + "\n"
    
    # Report changes on aggregate frequencies:
    
    [graph(report_changes(seq, i)) for i in [10, 100, 1000]]
    
    pylab.yscale('log')
    pylab.xlabel("Genome sequence")
    pylab.ylabel("G or C frequency")
    pylab.title("GC content across the genome")
    
    pylab.show()
    # Fill up CGR table:
    seq.reverse()
    k = 4
    
    pylab.figure()
    pylab.subplot(121, aspect='equal')
    pylab.title("CGR (k = 4)")
    
    cgra = get_cgra(k)
    
    #pylab.imshow(cgra, interpolation='nearest', cmap=cm.gray_r)
    k = 3
    pylab.subplot(122, aspect='equal')
    pylab.title("CGR (k = 3)")
    
    cgra = get_cgra(k)
    
    #pylab.imshow(cgra, interpolation='nearest', cmap=cm.gray_r)
    pylab.show()
