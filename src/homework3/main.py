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


if __name__ == '__main__':
    
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
        "id" : "NC_018038.1",
        "name" : "Goldfish"
    }, {
        "id" : "NC_004299.1",
        "name" : "Pufferfish"
    }]
    for a in animals:
        rec = load_entrez(a["id"])
        for f in rec.features:
            if f.type == "CDS":
                print "Name", f.qualifiers["gene"][0]
                print "Amino acid sequence", f.qualifiers["translation"]