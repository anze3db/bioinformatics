from main import get_seq

def simple_seq_diff(seq1, seq2):
    """Simple seq diff with two step lookahead"""
    
    diff = 0
    offset = 0
    skipped = 0
    for i in xrange(min(len(seq1), len(seq2))):
        if i + offset == len(seq2)-1:
            break
        if(seq1[i] != seq2[i+offset]):
            if(seq1[i] == seq2[i+offset+1] and seq1[i+1] == seq2[i+offset+2]):
                offset += 1
                skipped += 1
                continue
            diff+=1
    print "skipped ", skipped
    return diff+skipped

if __name__ == '__main__':
    
    
    seqs = [get_seq("NC_001807.1.fasta"), get_seq("NC_001807.2.fasta"), get_seq("NC_001807.3.fasta"), 
            get_seq("NC_001807.4.fasta"), get_seq("NC_012920.1.fasta")]
    lens = []
    for s in seqs[:-1]:
        lens.append(simple_seq_diff(seqs[4], s))
        print "" 
        
    import pylab
    pylab.bar(range(len(lens)), lens)
    pylab.xlabel("NC_001807.1 NC_001807.2 NC_001807.3 NC_001807.4")
    pylab.ylabel("Num differences + num skips")
    pylab.show()
    
