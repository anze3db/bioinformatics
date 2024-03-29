from cPickle import dump, load
import pylab

def get_seq(fa):
    """Return sequences in FASTA format"""

    f = open(fa)
    fc = f.read()
    fc = fc.replace('N', '')  # ignore aNy
    return "".join(fc.split("\n")[1:])

def rev_inv(seq):
    """Reverses and complement a sequence"""
    
    s = {'A': 'T', 'T':'A', 'C': 'G', 'G': 'C'}
    return "".join([s[i] for i in seq[::-1]])

def codon_walk(seq, offset, length=3):
    """Walk through all possible windows"""
    for i in xrange(offset, len(seq) - (len(seq) - offset) % length, length):
        yield i, seq[i:i + length]


def load_entrez(genome):
    """Reads Entrez data from a .pickle file. If the file does not exist it is created"""
    import os.path
    if not os.path.isfile(genome + '.pickle'):
        from Bio import Entrez
        from Bio import SeqIO
        handle = Entrez.efetch(db="nucleotide", rettype="gb", id=genome, email="smotko@smotko.si")
        rec = SeqIO.read(handle, "gb")
        handle.close()
        dump(rec, file(genome + '.pickle', 'w'))
    
    return load(file(genome + '.pickle'))


if __name__ == "__main__":

    genome = 'NC_006058'
    #genome = 'NC_007346'
    
    seq = get_seq(genome+'.1.fasta')
    
    ris = rev_inv(seq)

    # Paramecium tetraurelia. Start codon: ATG. Stop codon: TGA.

    stop_codon = {"TGA"}
    start_codon = {"ATG"}
    
    if genome == 'NC_007346':
        # Emiliania huxleyi virus 86. Start codons: ATG, TTG, CTG. Stop codons: TAA, TAG, TGA.
        stop_codon = {"TAA", "TAG", "TGA"}
        start_codon = {"ATG", "TTG", "CTG"}
        
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
                        # We need to get the strand from the original sequence:
                        orfs.append({"strand": seq[len(s) - (i + 3):len(s) - orf_start], "start": len(s) - (i + 3), "end": len(s) - orf_start, "reverse": reverse })
                    else:
                        orfs.append({"strand": s[orf_start:i + 3], "start":orf_start, "end": i + 3, "reverse": reverse })
                    orf_start = False
        reverse = True
        
    print "%s ORFs code for at least 60 amino acids" % sum(1 for g in orfs if len(g["strand"]) > 60*3)
    
    rec = load_entrez(genome)
    original = set(f.location.start.position for f in rec.features if f.type == "CDS")

    scale = range(0, 500, 3)
    precision = []
    recall = []    
    for c in scale:
        predictions = set(g['start'] for g in orfs if len(g['strand']) > c * 3)
        TP = len(predictions.intersection(original))
        precision.append(float(TP) / len(predictions))
        recall.append(float(TP) / (TP + len(original.difference(predictions))))

    pylab.plot(scale, precision, lw=1)
    pylab.plot(scale, recall, lw=1)
    pylab.xlabel("Minimum ORF length")
    pylab.ylabel("Score")
    pylab.title("Precision and recall as a function of minimum ORF length")
    pylab.show()


    # Calculate precision and recall for 125 codons or more
    c = 125
    predictions = set(g['start'] for g in orfs if len(g['strand']) > c * 3)
    TP = len(predictions.intersection(original))
    print "%.4f precision, %.4f recall for 125 codons or more " \
        % ((float(TP) / len(predictions)), (float(TP) / (TP + len(original.difference(predictions)))))
