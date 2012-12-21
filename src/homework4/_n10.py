import random

# metanje dveh kovancev (postenega in goljufivega)
# pri postenem (skrito stanje p) sta verjetnost cifre (C) in glave (G) enaki
# pri goljufivem (skrito stanje g) je cifra (C) dosti bolj verjetna

T = {
    0: {'p': 0.5, 'g': 0.5}, # z 0 oznacimo zacetek (stanje pred izbiro prvega skritega)
    'p': {'p': 0.95, 'g': 0.05},
    'g': {'g': 0.9, 'p': 0.1},
}

E = {
    'p': {'C': 0.5, 'G': 0.5},
    'g': {'C': 0.9, 'G': 0.1},
}

HMM = (T, E)

# nekaj primerov zaporedij
"""
>>> h, s = gen_zaporedje(HMM, 20); print "h: %s\ns: %s" % (h, s)
h: ppppppgggggggggggggg
s: GGGCGCGCGGCCCCGGCCGC
>>> h, s = gen_zaporedje(HMM, 20); print "h: %s\ns: %s" % (h, s)
h: ggggggggggggggggpppp
s: CCGGGCCCGCCGCCGCCGCC
>>> h, s = gen_zaporedje(HMM, 20); print "h: %s\ns: %s" % (h, s)
h: gggggggggppppppppppp
s: CCCGCCCCCGGGCCCGGCCG
>>> h, s = gen_zaporedje(HMM, 20); print "h: %s\ns: %s" % (h, s)
h: ggggppgggggggggggggg
s: CGCCGGCCCCCCCGCGCCGC
>>> h, s = gen_zaporedje(HMM, 20); print "h: %s\ns: %s" % (h, s)
h: pppggppggggggggggggg
s: CCGCCGGCCCGGGCCCCGGC
>>> h, s = gen_zaporedje(HMM, 20); print "h: %s\ns: %s" % (h, s)
h: pppppggggggggggggggg
s: CCGGGCGGCCCGCCGGGGCC
>>> h, s = gen_zaporedje(HMM, 20); print "h: %s\ns: %s" % (h, s)
h: ppppppppppggggggggpp
s: GGCGCGGCGGGCCGGCGCGG
"""

def izberi_enega(m):
    s = sum(m.values())
    r = random.random()*s
    for v, f in m.iteritems():
        r -= f
        if r <= 0.0:
            return v

def gen_zaporedje(hmm, dolzina):
    t, e = hmm

    # izberi prvo skrito stanje
    h = izberi_enega(t[0])
    s = izberi_enega(e[h[0]])
    prev_h = h[0]

    # vsako nadaljnje je pogojeno s predhodnim
    for i in range(dolzina-1):
        nh = izberi_enega(t[prev_h])
        h = h + nh
        s = s + izberi_enega(e[nh]) 
        prev_h = nh
    return h, s

def viterbi(s, hmm):
    t, e = hmm

    # seznam skritih stanj
    zh = set()
    for h, tmpd in e.iteritems():
        zh.add(h)

    zh = [0] + list(zh)

    # Create table V
    V = [{} for i in range(len(s)+1)]
    ptr = [{} for i in range(len(s)+1)]

    # Initialize i = 0; V(0, 0) = 1; V(k, 0) = 0 for k > 0
    for k in zh: 
        V[0][k] = 0 #t[0][k]*e[k][s[0]]
    V[0][0] = 1.0

    # for 1 = 1 : n, compute
    for i in range(1, len(s)+1):
        for l in zh:
            vals = [(V[i-1][k]*t[k].get(l, 0.0), k) for k in zh]
            max_val, max_k = max(vals)
            V[i][l] = e.get(l, {}).get(s[i-1], 0.0)*max_val
            ptr[i][l] = max_k

    # trace back
    pi = []
    pi_L = max([(V[-1][k], k) for k in zh])[1]
    pi.append(pi_L)

    for p in ptr[-1:1:-1]:
        pi.append(p[pi[-1]])

    pi.reverse()
    return V, zh, ptr, "".join(pi)


def logmv(a):
    min_val = 0.0000000001
    return math.log(max(a, min_val))

import math
def viterbi_log(s, hmm):
    t, e = hmm

    # seznam skritih stanj
    zh = set()
    for h, tmpd in e.iteritems():
        zh.add(h)

    zh = [0] + list(zh)

    # Create table V
    V = [{} for i in range(len(s)+1)]
    ptr = [{} for i in range(len(s)+1)]

    # Initialize i = 0; V(0, 0) = 1; V(k, 0) = 0 for k > 0
    for k in zh: 
        V[0][k] = logmv(0.0) #t[0][k]*e[k][s[0]]
    V[0][0] = logmv(1.0)

    # for 1 = 1 : n, compute
    for i in range(1, len(s)+1):
        for l in zh:
            vals = [(V[i-1][k]+logmv(t[k].get(l, 0.0)), k) for k in zh]
            max_val, max_k = max(vals)
            V[i][l] = logmv(e.get(l, {}).get(s[i-1], 0.0)) + max_val
            ptr[i][l] = max_k

    # trace back
    pi = []
    pi_L = max([(V[-1][k], k) for k in zh])[1]
    pi.append(pi_L)

    for p in ptr[-1:1:-1]:
        pi.append(p[pi[-1]])

    pi.reverse()
    return V, zh, ptr, "".join(pi)

from Bio import Entrez
from Bio import SeqIO

def get_gb_from_entrez(nid):
    """Iz Entrez baze prebere in vrne zapis tipa genebank za zaporedja z id-jem nid."""
    handle = Entrez.efetch(db="nucleotide", rettype="gb", id=nid)
    rec = SeqIO.read(handle, "gb")
    handle.close()
    return rec

def oceni_ujemanje(z1, z2):
    """Vrne delez ujemajocih polozajev med podanima zaporedjema z1 in z2."""
    assert(len(z1) == len(z2))
    m = sum([e1 == e2 for e1, e2 in zip(z1, z2)])
    return float(m)/len(z1)

def segmentiraj_genom(genom_id):
    """Za podani genom z id-jem genom_id vrne niz znakov
    G ali i, ki dolocajo ali je dano mesto v genomu del
    gena ali medgenska regija. """
    genom = get_gb_from_entrez(genom_id)
    c = set()
    for f in genom.features:
        if f.type not in ['CDS', 'mRNA', 'gene']:
            continue    
        c.update(range(f.location.start.position, f.location.end.position))            
    h = ""
    for p in range(len(str(genom.seq))):
#        if p in c:
#            h = h + "G"
#        else:
#            h = h + "i"
        h = h + "G" if p in c else h + "i"
    return h, str(genom.seq)


def zgradi_hmm(s, h, r=0.0, state_set=None):
    """Na podlagi vidnega zaporedja s in skritega
    zaporedja h vrne skriti Markov model."""
    t = {}
    # {key: skrito stanje h1,
    #  item: {key: skrito stanje h2, item: f(h1->h2)}
    # }
    e = {}
    # {key: skrito stanje h1,
    #  item: {key: vidni simbol s1, item: f(s1|h1)}
    # }
    
    prev_hs = 0

    if state_set == None:
        state_set = set(h)

    #add starting counts
    for a in state_set:
        t[a] = {}
        for b in state_set:
            t[a][b] = r

    for a in state_set:
        e[a] = {}
        for b in set(s):
            e[a][b] = r

    for hs, ss in zip(h, s):
        # e
        tmpd = e.setdefault(hs, {})
        tmpd[ss] = tmpd.get(ss, 0) + 1
        # t
        tmpd = t.setdefault(prev_hs, {})
        tmpd[hs] = tmpd.get(hs, 0) + 1
        prev_hs = hs

    # start to all is equally likely - THIS IS DIFFERENT
    for k in t.keys():
        if k != 0:
            t.setdefault(0, {})[k] = 1.

    # normaliziraj (iz abs. v rel. frek.)
    nt = {}
    for h1, h2_f in t.iteritems():
        n = float(sum(h2_f.values()))
        tmpd = {}
        for h2, f in h2_f.iteritems():
            rf = f / n
            tmpd[h2] = rf
        nt[h1] = tmpd
    ne = {}
    for h1, s1_f in e.iteritems():
        n = float(sum(s1_f.values()))
        tmpd = {}
        for s1, f in s1_f.iteritems():
            rf = f / n
            tmpd[s1] = rf
        ne[h1] = tmpd
    hmm = (nt, ne)

    return hmm


def zgradi_hmm_BUG(s, h):
    """Na podlagi vidnega zaporedja s in skritega
    zaporedja h vrne skriti Markov model."""
    t = {}
    # {key: skrito stanje h1,
    #  item: {key: skrito stanje h2, item: f(h1->h2)}
    # }
    e = {}
    # {key: skrito stanje h1,
    #  item: {key: vidni simbol s1, item: f(s1|h1)}
    # }
    prev_hs = 0
    for hs, ss in zip(h, s):
        # e
        tmpd = e.setdefault(hs, {})
        tmpd[ss] = tmpd.get(ss, 0) + 1
        # t
        tmpd = t.setdefault(prev_hs, {})
        tmpd[hs] = tmpd.get(hs, 0) + 1
        prev_hs = hs

    # start to all is equally likely
    for k in t.keys():
        t.setdefault(0, {})[k] = 1
        for k2 in t.keys():
            t[k].setdefault(k2, 1.0)

    # normaliziraj (iz abs. v rel. frek.)
    nt = {}
    for h1, h2_f in t.iteritems():
        n = float(sum(h2_f.values()))
        tmpd = {}
        for h2, f in h2_f.iteritems():
            rf = f / n
            tmpd[h2] = rf
        nt[h1] = tmpd
    ne = {}
    for h1, s1_f in e.iteritems():
        n = float(sum(s1_f.values()))
        tmpd = {}
        for s1, f in s1_f.iteritems():
            rf = f / n
            tmpd[s1] = rf
        ne[h1] = tmpd
    hmm = (nt, ne)
    return hmm


# napoved s "posterior coding"
# za kar se prej rabimo forward in backward algoritma
def forward_chain(s, hmm):
    t, e = hmm

    # seznam skritih stanj
    zh = set()
    for h, tmpd in e.iteritems():
        zh.add(h)

    zh = [0] + list(zh)

    # Create table
    f = [{} for i in range(len(s)+1)]

    # Initialize i = 0; f_0(0) = 1; f_k(0) = 0 for k > 0
    for k in zh: 
        f[0][k] = 0
    f[0][0] = 1.0

    # Recursion (i=1..L):
    for i in range(1, len(s)+1):
        for l in zh:
            sum_val = sum([f[i-1][k]*t[k].get(l, 0.0) for k in zh])
            f[i][l] = e.get(l, {}).get(s[i-1], 0.0)*sum_val

    # P(x)
    ps = sum([f[len(s)][k] for k in zh]) 
    return f, ps

def log_sum(v1, v2):
    nv1 = max(v1, v2)
    nv2 = min(v1, v2)
    return nv1 + math.log(1.0+math.exp(nv2-nv1))

def sum_log(vals):
    s = vals[0]
    for v in vals[1:]:
        s = log_sum(s, v)
    return s

def forward_chain_log(s, hmm):
    t, e = hmm

    # seznam skritih stanj
    zh = e.keys()
    zh = [0] + list(zh)

    # Create table
    f = [{} for i in range(len(s)+1)]

    # Initialize i = 0; f_0(0) = 1; f_k(0) = 0 for k > 0
    for k in zh: 
        f[0][k] = logmv(0.0)
    f[0][0] = math.log(1.0)

    # Recursion (i=1..L):
    for i in range(1, len(s)+1):
        for l in zh:
            sum_val = sum_log([f[i-1][k] + logmv(t[k].get(l, 0.0)) for k in zh])
            f[i][l] = logmv(e.get(l, {}).get(s[i-1], 0.0)) + sum_val
    # P(x)
    ps = sum_log([f[len(s)][k] for k in zh])
    return f, ps

def backward_chain(s, hmm):
    t, e = hmm

    # seznam skritih stanj
    zh = set()
    for h, tmpd in e.iteritems():
        zh.add(h)

    zh = [0] + list(zh)

    # Create table
    b = [{} for i in range(len(s)+1)]

    # Initialize i = 0; b_k(L) = 1.0 for all k
    for k in zh: 
        b[len(s)][k] = 1.0
    # Recursion (i=1..L):
    for i in range(len(s)-1,0,-1):
        for k in zh:
            sum_val = sum([t[k].get(l, 0.0)*e.get(l, {}).get(s[i], 0.0)*b[i+1][l] for l in zh])
            b[i][k] = sum_val

    # P(x)
    ps = sum([t[0].get(l, 0.0)*e.get(l, {}).get(s[0], 0.0)*b[1][l] for l in zh]) 
    return b, ps

def backward_chain_log(s, hmm):
    t, e = hmm

    # seznam skritih stanj
    zh = e.keys()
    zh = [0] + list(zh)

    # Create table
    b = [{} for i in range(len(s)+1)]
    # Initialize i = 0; b_k(L) = 1.0 for all k
    for k in zh: 
        b[len(s)][k] = math.log(1.0)
    # Recursion (i=1..L):
    for i in range(len(s)-1,0,-1):
        for k in zh:
            sum_val = sum_log([logmv(t[k].get(l, 0.0)) + logmv(e.get(l, {}).get(s[i], 0.0)) + b[i+1][l] for l in zh])
            b[i][k] = sum_val

    # P(x)
    ps = sum_log([logmv(t[0].get(l, 0.0)) + logmv(e.get(l, {}).get(s[0], 0.0)) + b[1][l] for l in zh])
    return b, ps

def post_decoding(s, hmm):
    f, fps = forward_chain(s, hmm)
    b, bps = backward_chain(s, hmm)
    if not (abs(fps - bps) < 0.0000001):
        print fps, bps, fps-bps
        print "ERROR"
        return ""

    # posterior coding
    hn = ""
    for i in range(1, len(s)+1):
        max_p, max_k = max([(fp*b[i][k], k) for k, fp in f[i].iteritems()])
        hn = hn + max_k
    return hn

def post_decoding_log(s, hmm):
    f, fps = forward_chain_log(s, hmm)
    b, bps = backward_chain_log(s, hmm)
    if not (abs(fps - bps) < 0.01):
        print fps, bps, fps-bps
        print "ERROR"
        return ""

    # posterior coding
    hn = ""
    for i in range(1, len(s)+1):
        max_p, max_k = max([(fp+b[i][k], k) for k, fp in f[i].iteritems()])
        hn = hn + max_k
    return hn

if __name__ == "__main__":
    # na primeru konvancev preveri delovanje
    h = 'ggggggpppp'
    s = 'CCCCCCGGGC'

    #random.seed(42)
    h, s = gen_zaporedje(HMM, 20)
    print "dejanski zaporedji"
    print "h:", h
    print "s:", s
    print

    print "napoved z Viterbi"
    V, zh, ptr, hn = viterbi_log(s, HMM)
    print "h:", hn
    ovit = oceni_ujemanje(hn, h)
    print "ujemanje:", ovit
    print

    print "napoved s pd (posterior decoding)"
    #hn = post_decoding(s, HMM)
    hn = post_decoding_log(s, HMM)
    print "h:", hn
    opd = oceni_ujemanje(hn, h)
    print "ujemanje:", opd 
    print

    trials = 2000
    print "primerjava obeh postopkov na %s zaporedjih" % trials
    pd_win_cn = 0
    vit_win_cn = 0
    for r in range(trials):
        h, s = gen_zaporedje(HMM, 100)

        # Viterbi
        V, zh, ptr, hn = viterbi_log(s, HMM)
        ovit = oceni_ujemanje(hn, h)

        # posterior decoding
    #    hn = post_decoding(s, HMM)
        hn = post_decoding_log(s, HMM)

        opd = oceni_ujemanje(hn, h)
     
        pd_win_cn += (opd > ovit)
        vit_win_cn += (ovit > opd)

    print "viterbi zmaga: %s (%0.1f%%)" % (vit_win_cn, 100.0*vit_win_cn/trials)
    print "pd zmaga: %s (%0.1f%%)" % (pd_win_cn, 100.0*pd_win_cn/trials)
    ties = trials - vit_win_cn - pd_win_cn
    print "izenacenje: %s (%0.1f%%)" % (ties, 100.0*ties/trials)

