'''
Created on Nov 25, 2012

@author: smotko
'''
from collections import Counter

if __name__ == '__main__':
    s = """>Rosalind_5098
AGGAATCTAACTCTTTTCATCGATACGTGTGCTACCCCAGCCAACTACGTGGCACGGCCC
CCTGGAAACTGGGTACCCCGGTATTCGACCGTCCGGCCCGCGGATCCGCTTCATTTTAAT
CTCGTCATACGTTGCGACCGAACATGGAGTCCTCGAATCGTTCGGTTTTGAGGTCATCAC
CGCTTTGAAGGGCGTCGCGCAAGCAGACGTTAAGGCCAAGGTGGGGCTGTTGGAGACAGT
TCCCGTGTTAGGAGGTTTCCAGGATGGTCTCATTAGGCATATAAATTATGTAATGAGCCG
AAATAGCCCGGAAATTTGAAAATCCCCGACCTCATTCACCTTCCCCCGCCCTGCACAATG
ACATTATGGACTGTACGGAGTTGGTTCATCGTGATAGCAAAAAACCAAACGTAGCCTTTC
CGAGCAAGGTAGTTCCGACACTTCGAATCATGTCACAACTGCGCCTGAACATTCCTCTAG
GAGCCACCAATCAAGCAACGGGGCGAAACTTCACCGCGGGTAGCGCCCTTTAACAGCTCC
AAAGCCGCTGATCTTCGCCCGCGAGGCACAGCGCGTGCTATATCCGGTTAATCCGACGAT
CGTGACTTGCATATCGTTATGTAGTTTATTCTTCGATCTGTCTCGCGTACCATGACACTA
GATCGCTTGCAGCGACCTAGCTATGGCTGAATTTAGTATTTTCTGCGCAGCGTCATGTGG
TGTATATGCGCGGAAAGTAGAAATCTCTACCCTGAGGGTCACATGCTATTGGTATGGAAA
CTGACCCCGTCCCTATTGAAGTTTACGTACGGGGTGGGTTAAACGGGAACCTTCACTATG
TACCCATATGGAGCTTCGTGGTGGAGTAAGATGCAAACCGTGCATCGGCATGCGGGATTA
GTCGTATTGCGAAACCCGGACCGGGA
>Rosalind_6847
CAATTTGAGTCCGAATTCACAGCAGGCCAATGCGAGCCGTGCGGTCGCGCCCTGATCAAC
GATGCCTACGACCTGAATCGACCTGCTATATATGTGCCGACAAGACTTTACAAGGCACAC
ACAACTCGAAGTATTTTAGGTGCGTCCGACCTATAAATGGTGGGGGCACTCAGGTCAGAG
CTTGAGCCAACACGATCGGACGGCCCTGGATATCAAGTTTGACAAACATCAACACTCAAC
ACCAAGTAAGGACTGTAGTGGAGCTAATCTCAACATCGAGCGAGTACTCAGTTGTCTTAT
AACGAAGAATTGTCGGAAAGGAGAAGAATGTGGGATTGGGTTAAGTGTCTACCTGCAAGC
TTTATCCGGGTACAACTGGCGCCGACGTTAGTAACTGGAAATGTCATGTTCGCCCTCTAG
GATCTTAAATTATATAAAAACAAGCTGTATAGACAAATTTCGAGGGCCTGGCAACCTATG
GAAAACGAGCGAGTACACCGTTGTCGATTAACATTCGTCCTTGTGGTGGGCAGATCGGTC
AATGACGGTTACTAGCATAGTACACCATTGCCAGTGAGAAACCACAGAAAACGGGGTGAT
GGAGCAATGGGAACGCACCAACTTCGGTAGGATCGTGCTAGATGTTAAATGCTTTGCTGC
AGCCCTTCCGCACCACATAAAAGTGTAGGTAAACCTCTTTTATGTATTAGTTCCAAATGA
CCGTCGCTGCAAAGGCGTGTCCGCCACGGTTTGAAAGGTGAGCGGCATCGACCCAGAGCC
TAAGACTCCTTCTGTCTATTGTATAGCCCGTGCGTCTTACGCTCCTCAGGCATCGGGTAC
AGCCCCCCGAGCGTTACCTTATTCAGGGTGAGTATCGCTACCCA
>Rosalind_4435
CTTACCCTTCAGACCCCTCCAATACATTCTGATACCTGCTAAGCCTTCTATGTTCAACAT
AACCCCGGCCCTTTAGCACGAAATCTTCAGTGTAGTCCTGTTCGAACCGGCTGGTTACGT
TAACTATTACTCGGGGTAATTGATGAATACGCGTGAATGTGCCCGCGCGTCCACATTTCA
TATACAAAAGCCTCATGTAAGCATGTGTGACCTAGATCGGAGATTAAAGTTTATAGACTG
TTTTCGTTGCCGACGTACCACAAGTACTGTCCAGACCTATCGACCAAAGGTAAGTCCAGC
GAAACCGGGCGTCGCGGTAAGAGTATATCAACTGTTATACGAATCCCTGGAGGTTCTTTA
CATGGGTGACAGAAAAAGTGCCACGCTGGTCCGCAGCATTAGGACGATGAGTTTGTAGGG
TGCTGTCACCTGGAATGGTACTACAAATAAGCATGTTTATGCGGATGCAAGTGCAAGCAT
CGCGGGCTACTCCTCATCGAAACGCGCTTTCAAATGTGCATATTGATAACTAATTGATAC
ACCATGCGTAGCTATCGGTCGGTTGCTCTCACATGCTGGATAAGGAAAGACTGATGTTAC
CGCTAATCAGCTACGGTGTGTGTTCGATGTGTGGACACAACTCTTCTGAGCGAAGACCCG
ACACCGTCCGCTGATGTCGGGGCGAAAGTCCCGATGCCACCGGCTGTATTACACAACTTT
AGAGACCTGTACGTAGTAACTGGGATAACTCCTGCGTCATTGGCACGTCTAGAAAGTAAC
GTTAGCACTGTGGAGTAAGGAATGACTTATCTCCGCCACCTAAACGCCTTTAGCTTTGCG
CAATTCCGGCTTCGTTGATCCAATCTACTTATGACAGGATTTGCGTTCTCAGTAAGATTA
GA
>Rosalind_1383
ACTCTAGCAAACTATCCGCGATCATAAGCTATAAATTGTTTTTGTGAAACTCCTGATGTG
CCAAGTCCATTCATTACCACCAGCTGAACCGACTACAGTGTGATTTGTGTTAACGGCCGC
ACAGGCTGCTCTCATGCAAGCACTATGGTCTACACTCCCCGGATCGTCAGGCATGCTATT
AGCACTAATCAAAGTATCGCATTGAGGCAGAGTACTAAATACCACGCATCGCCATATCAT
ATCTAAAGGATCTGCGCAACGTCCGGAGGTCCGTGTCCGATAAGACTCGCGACCCAGGCG
GGACGTATAGTCATACGGAAGGAAACCTGTCCTTATCGGGTTTATCATGGTTATCATGGG
GAGCCCCCTAGTCTAGAGGCTCGCCGCAACACGGTTCGTAACATTGCAGCGCACTCACGC
GTGTCGTAGTAATGGGGAAGTCGTGTTCCCCGTGCATCTCCCTCAGTTCGGCCTGAGTCA
GCCGCATGAAAGCAAGCCACTCTTTATCCCAGTTTAAGGATCCGTTCGCAGGCATAATGG
CGAAAAGAAACACACCGAGAGCAAACAAAAATTACGGTTCACCATCGTCCGCCCCCGCAT
CAGTAGCAAGGTACTCCAACTCTAATTTCGAGGGAGCGTGTCAATATTTTCCATCTGGAC
CCAGCGACTGATTAACGCGGCCCGGACCTTCCCTCGTGCCGCAAACTCACGTACCCGATG
TGAAGCGCGGACCTGATCGGTATCAAGTTGATGAATGCGAGGAAACATGCGCGAGGACAG
TCTAGGCCCGGACAGCGACCTTTGTCGATTATTACAAGCAGTCGAGCTACGTGATGGCGC
CGTCGTGGTCGAAGGGGAACCCTTGACTGATGAAGGTAGCCTGCGCGTAACCATCGCATT
TTCCGCGAGTGAAGTGACGGGGA
>Rosalind_6633
TACCTCGTTTGAGTATTAACGCGATCTGAGGGTATGCCTTAGGTTTAGGCCAAAGATAGT
TTTTGAGATGTCCTTCGAACTTAATATTCCCATAGTAATGGGTCGCGATATTCCGAAGGT
ATGCTTCATCTTCCGGGAACCAGACAAGACGTACCGATAAGTCCCTGCATTCTTAGATGT
TAAACTGACTCCAACTAAACGACGCGCCTACTAAATACAGAGTCCAATCTTCGAAGGGAG
CTAACCACGAGATATGAAGGCTGGTCGCGGACCTCCCGCACCTCCTAGAATCCGCCCGAT
GTGGCGGGCGTTATGCCCTTTGACCTGGCCAGTCAAAAGGTTTAAAGGTCTATACTGTGT
ATGGGTTAGCACCTGGATCTTCATTGTCTACCAACTCGCTGATTCCCGATTTATGCAAAT
GAGTTTTAGAAAACCCGTGGAGTGGCATTAGTCACATTAGAAACCACCAGGTTTGTTCCA
ACCTGCCTCTAGGCTATCAATTCTATAACCCTTGCCAGCTCTTGCGTTTCGGTGATTAGG
GTAATCGCCACGTATCCAACTTCCCGTCGGAAGTAGATCGCCGTTTCATTGTATTACCCC
CCGTGACCTCGATGGAACGACCGATAGCGACCCGGTGTACCGAAACTAGGGGAAGAGAAC
TCGGAGATTACAGCAGGTATAATGGACCAGCGCACTCAAAGATTCGCATCAACAGCCCAT
TAGTCTCTGATGCACTCTTCGGTATCCTTTATGGTGAGAGTTACTTTCATAATATACACC
AAGATGCGCCTATAACGAGTCATTACATTGTTCTCTGATCTATTAGGCCGCAATCGCTGC
GGTTCAGATCAAGAACTTTACTAGCGACCCAACTTCAATAGGGTACCGCGCTAAATGTAC
TTACG
>Rosalind_3124
ATCCCGCTCCATACCGAGGTGAAATGCCGGTTTGTTTGCAGTACCCGGCTGATTAAAATG
ATAGTATGCAACGTCGACCATGGTCCTGCGTGCCTAATATGGTTAGCGCTATTATTGCTT
CCATGCGTTTCATATTAGCGGTTCCCTAAACAGCCCATCAGGGACGCACTCTAGACAGTC
CGGAGATAGTCCCGGGTGTTTTTCCAGCATGTTGCGATCGGTCCGACTTATCAGGGTCTT
CGTATTGAAGTTCTTCCAAGTAGGTTGACAACCGTCGATCCGGGCCGCTTAGGGTTTAGA
GGGCTCGCAGATTGGTCATACTGTCATGATGTGAATAATGGAGGTGCACGCGTTACTCAG
TCGATAAATAACGCTATGCCTGTAACACACCACTTGCCTTTACTGTGCAACCACCGCACC
CGGAGGAGCTGAAAGCAGGACGTAATCTGTCTTAGCATGCTCAGGCGGACGCCACGGGCT
TATTGCCATTCGTCAAGACCTAACGCTATGGAAGGAAGCCCAGATCGAATCAGATTTGAG
TGACTGCACCGAGCTTCGTTTACCCGATTGTTAAAGTCTGAGGCTTATTTTGCTAAGATG
ACTATGACGGAATAAAACTATAACCCTGATTTGGTGTTCCTACCCACTAGATCGCCAGAG
CCGCTGATGCGCTTCTCGCTGGATTCTGGACTGTGGGTGATACCGGCTTCATTAATATCG
AATGCGCGGGAGGGCCAGGTGTCAGGATCTACGTGTTAAAGTTCCTGATGCTTAGCGAAA
AAATTAAGTCAAACAATCCTCCCCCGATGGATTAGCGATCTTCCACAGGCCCCGCGCGCT
TAAAAACTTTGCGGGAATGCAACTCAGTCCGATATATTTGGGAATCAGAAGAGTGCTCCT
GTCTTGGGCAGTGCTAATTGCCCTTTGCGATTACCTGACACCAGGCTTGCTGGTGCGAGT
GGCGCAATTGGTGACG
>Rosalind_6527
AGCCTGCGCCCAGTGATGTCATCTTGAGACGGCGGTAATAGGGCCGTACTCTAAAATACC
CCGCTGCTTGGAAAAAATTTTCCTGCATCCAGTAAGATCCGGTGTCAGAGGCAAGTTATG
TCAATTTAAAGAGCGGAAGCGGCACCGTGGAAAGGGACGGACTAGCGGCGTGTGGCGTAA
ACGTGGCTGAGGGACGGAATTGCCTCTTTATATTCCAGAAGTTATCGATATTTAGCCCGA
GAGTGACAATGGCAGTGCAGCTCGCTGACGCTAACAAGGTCCCTGTAGCGTGCATTTCCG
CCTCACGAACGTGGTAAGCTCCGTTGCGCGGACATTATAAATAGACTGCGCCACCCGCTC
GCTTGTTGTCAACCGAGGGGTATTCTCTTATCAAATCGCAACGAGGGCGGTCCCAATCTG
TGCAGGGACATTGACGCGATGCCTAAAAATTTCACGGTCTCGTATTAGCCCCCTTTAGGT
ACTAAGCCTTTTCGTAGTCATTCCATAGTTTAAGATTCTGCACTACGAATCCCCCAATTG
GGGTAAAGGGAGCATCATCCCCCCTCAAGATTAGGGTACACGCCGGTTAAGTCTCCGTAA
CAGGGCATATTCACCTCGTTTCGATTGTGCTCCAGAACTAAAGCTAGCGAAGTCGTTTAG
AAGAAGTTAGGCTTTCGCTGGCCCAATAGCAATCTTGTCACAAGTCCTAGGCTTTGTTCT
TGTGGTAGCAATCGCCGGCACCTTCCGATTACAAGAGTAGGTGCTTTGAACCAAGACGGG
GCGAAACCCCACCCGGTGTGGATTAGGTCGCGACCAAGTCGATGATATGTCTTAGGTCGC
TATATTCGATGAGCCCTCCCAATGACGCCATTAAGATCTTAAGGGAAACCA
>Rosalind_3592
GCCAAACCAGTACATTAGGGGTATCCCCGTTAAGGTCACCACGCGACAAACTAGGCGTGA
GCTACCCTTGTTGTCCTAAGGCCAAGGGATGTTGCTAGCTGCTATCTACATACCTTGGAT
GATCCACTTTTGGTTCATAAATAATTCTGGAATAGCGGGTCTATAGCCTGACGGCCATTT
TGACAAAGTCGGTCAACGTCCTAACCCACTACGGCTGTTTGTCGACGCGGATAATACAGT
GGGGATACCGCCTTCAGAATGCAGTCGACGCGCCTATGTTCATGTTTTTAGCCCACGACA
GTGTCAACGGCCGCTTTAAAGACGGTGGAAGGACGCCGTCGTTTCTTTTCCAAAACTTGT
GAGCGTTGTTCGATTGGTTAGAGGTCTTTGGCATAGTGTCCCGGTGGGTGTGCGAGTGAT
CGTGCATGGTTTGATTACGATGGAACTCTTTGAGACGCGCCAACGCGACTTCCTCCCAGT
TGATCGGTCCTTGTCCTGTGTTAGGCACCTGTCGTAGACCAGAGTCTAGTCTGTATTTCT
CTATGCCCCTACCACACGGTACGATCAATACACCCCCGGAATTGACGAAACGCGTGATCA
ATCGATACCAGGTACAGCCTCTCTCGTACAATAACTTGATGATCGCCCACCCTTCAATAG
TACAACATTACATATGGTCTGGTGAGGCAGCTACTCAAACAATACGTAGTGCTACTCTCT
TTGATTGCGCAATTGATATCGCTCGGTGAACACCGGCTGCCGCCTCCTTGGTTGCGGAGG
CCCAATGCTAAGGCGGATATTGTCTGTACTGACGTGCTTAGAACTTTCACGGGATTCGCA
AACCTGCATAGCGGCTTCTCCAGTCGAACCGTCACCATGCGAACAGAAAGGGATGCTGTT
CCGTCACCCCTGAGCTAATTTCGACACGCCCTCAGCTAACTATCGCAGCGGGCCGCACTT
TCGTGTCA"""

    s = s.replace('\n', '').strip().split('>Rosalind_')
    a = {}
    max = (0,0)
    for i in s[1:]:
        
        c = Counter(i[4:])
        m = float(c['G']+c['C'])/len(i[4:])
        if m > max[1]:
            max = (i[:4], m)
    print "Rosalind_" + max[0]
    print "%0.6f%%" % float(max[1]*100) 