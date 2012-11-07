from cPickle import dump, load

def get_seq(fa):
    """Return  mitochondrial sequences in FASTA format"""

    f = open(fa)
    fc = f.read()
    fc = fc.replace('N', '')  # ignore aNy
    return "".join(fc.split("\n")[1:])

def rev_inv(seq):
    """Reverses and invertise a sequence"""
    
    s = {'A': 'T', 'T':'A', 'C': 'G', 'G': 'C'}
    return "".join([s[i] for i in seq[::-1]])

def codon_walk(seq, offset, length=3):
    """"""
    for i in xrange(offset, len(seq)-(len(seq)-offset)%length, length):
        yield i, seq[i:i+length]

if __name__ == "__main__":
    seq = get_seq('NC_006058.1.fasta')
    ris = rev_inv(seq)
    
    stop_codon = {"TGA"}
    start_codon = {"ATG"}
    # Paramecium tetraurelia. Start codon: ATG. Stop codon: TGA.
    # Emiliania huxleyi virus 86. Start codons: ATG, TTG, CTG. Stop codons: TAA, TAG, TGA.
    orf_start = False
    orfs = []
    for s in (seq, ris):
        for frame in range(3):
            for i, codon in codon_walk(s, frame):
                if not orf_start and codon in start_codon:
                    orf_start = i
                if orf_start and codon in stop_codon:
                    orfs.append(s[orf_start:i+3])
                    orf_start = False
    
    print "%s ORFs code for at least 60 amino acids" % sum(1 for g in orfs if len(g)>=1560)
    
    
#    from Bio import Entrez
#    from Bio import SeqIO
#    handle = Entrez.efetch(db="nucleotide", rettype="gb", id='NC_006058', email="smotko@smotko.si")
#    rec = SeqIO.read(handle, "gb")
#    handle.close()
#    
#    dump(rec, file('NC_006058.pickle', 'w'))
    rec = load(file('NC_006058.pickle'))
    
    for f in rec.features:
        if f.type == "CDS":
            continue
            print len(seq[f.location.start.position:f.location.end.position]), seq[f.location.start.position:f.location.end.position]
            print "start", f.location.start.position
            print "end", f.location.end.position
            print "strand", ["-", "+" ][(f.strand+1)/2]
            print
    
    
    