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
    
    stop_codon  = {"TGA"}
    start_codon = {"ATG"}
    # Paramecium tetraurelia. Start codon: ATG. Stop codon: TGA.
    # Emiliania huxleyi virus 86. Start codons: ATG, TTG, CTG. Stop codons: TAA, TAG, TGA.
    orf_start = False
    orfs = []
    reverse = False
    for s in (seq, ris):
        for frame in range(3):
            orf_start = False
            for i, codon in codon_walk(s, frame):
                if not orf_start and codon in start_codon:
                    orf_start = i
                if orf_start and codon in stop_codon:
                    if reverse:
                        orfs.append({"strand": seq[len(s)-(i+3):len(s)-orf_start], "start": len(s)-(i+3), "end": len(s)-orf_start, "reverse": reverse })
                    else:
                        orfs.append({"strand": s[orf_start:i+3], "start":orf_start, "end": i+3, "reverse": reverse })
                    orf_start = False
        reverse = True
    print "%s without filtering" % len(orfs)
    print "%s ORF codes with at least 60 amino acids" % sum(1 for g in orfs if len(g["strand"])>180)
    
#    from Bio import Entrez
#    from Bio import SeqIO
#    handle = Entrez.efetch(db="nucleotide", rettype="gb", id='NC_006058', email="smotko@smotko.si")
#    rec = SeqIO.read(handle, "gb")
#    handle.close()
#    
#    dump(rec, file('NC_006058.pickle', 'w'))
    rec = load(file('NC_006058.pickle'))
    match = []
    for o in orfs:
        for f in rec.features:
            if f.type == "CDS":
                #if o["strand"] == seq[f.location.start.position:f.location.end.position]:
                if o["start"] == f.location.start.position:
                    match.append({"calculated": o["strand"], "original": seq[f.location.start.position:f.location.end.position]})
                continue
                print len(seq[f.location.start.position:f.location.end.position]), seq[f.location.start.position:f.location.end.position]
                print "start", f.location.start.position
                print "end", f.location.end.position
                print "strand", ["-", "+" ][(f.strand+1)/2]
                print
    #dump(orfs, file('orfs.pickle', 'w'))
    min_len = 9999999999
    for m in match:
        if len(m["original"]) < min_len:
            min_len = len(m["original"])
        #print "or", m["original"]
        #print "ca", m["calculated"]
    print "matched: ", len(match), "60>", sum(1 for g in orfs if len(g["strand"])>=60) 
    print min_len, float(min_len)/sum(1 for g in orfs if len(g["strand"])>=280), float(min_len)/sum(1 for g in orfs if len(g["strand"])>=60) 