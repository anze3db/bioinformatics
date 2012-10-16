from collections import Counter

def get_seq():
    """Return  mitochondrial sequences in FASTA format"""

    f = open("NC_001416.1.fasta")
    fc = f.read()
    fc = fc.replace('N', '') # TODO: figure out how to handle N (any nucleic acid)
    return "".join(fc.split("\n")[1:])

def report_kmers(seq):
    """Report k-mers probabilities"""
    
    global cnt, c1, c2
    cnt = len(seq)
    c1 = Counter(seq)
    c2 = Counter([l1 + l2 for l1, l2 in zip(seq[:-1], seq[1:])])
    
    def print_probability(cn):
        """Print ordered probabilities"""
        
        l = list(cn)
        l.sort(key=lambda c: cn.get(c), reverse=True)
        for c in l:
            print "%s - %0.2f%%" % (c, cn.get(c) * 100 / float(cnt))
        print ""

    print_probability(c1)
    print_probability(c2)
    
def calculate_deviations():
    """Calculate deviations from joint probabilities"""
    
    # TODO: Globals, really?
    global cnt, c1, c2
    ls = []
    for c in c2:
        l = list(c)
        
        # TODO: Figure out if k is being calculated correctly 
        k = 100*(c1.get(l[0])/float(cnt))*c1.get(l[1])/float(cnt)
        
        # TODO: Refactor so that c2 percentage isn't being calculated twice
        ls.append((c, abs(c2.get(c)*100/float(cnt) - k)))

    ls.sort(key=lambda c: c[1], reverse=True)
    for c in ls:
        print "%s - %0.2f%%" % (c[0], c[1])
    print ""
    
def report_changes(seq):
    """Report changes of the aggregate frequency of G or C"""
    
    window_size = 1000
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
    print window_frequence

if __name__ == '__main__':
    seq = list(get_seq())
    report_kmers(seq)
    calculate_deviations()
    report_changes(seq)