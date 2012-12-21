from _n10 import * #n10 was also slighty modified (20111222)
import random

def viterbi_training(s, num_hidden):
    zs = set(s)
    zh = [0] + [str(x) for x in range(1, num_hidden+1)]

    # random model
    e = {}
    for k in zh[1:]:
        e[k] = dict((k, random.random()) for k in zs)
    t = {}
    t[0] = {}
    for k in zh[1:]:
        t[0][k] = random.random()
    for k1 in zh[1:]:
        t.setdefault(k1, {})
        for k2 in zh[1:]:
            t[k1][k2] = random.random()
    norm_dd(t)
    norm_dd(e)
    hmm = (t, e)

    for i in range(10):
        _, _, _, hn = viterbi_log(s, hmm)
        hmm = zgradi_hmm(s, hn)
        izpisi_model(hmm)

    return hmm

def norm_dd(d):
    """ normalize a dict of dict structure """
    for a in d:
        sa = float(sum(d[a].values()))
        for b in d[a]:
            d[a][b] = d[a][b]/sa

def baum_welch(s, num_hidden, r=0.0):
    zs = set(s)
    state_set = [str(x) for x in range(1, num_hidden+1)]
    zh = [0] + state_set

    # random model
    e = {}
    for k in zh[1:]:
        e[k] = dict((k, random.random()) for k in zs)
    t = {}
    t[0] = {}
    for k in zh[1:]:
        t[0][k] = random.random()
    for k1 in zh[1:]:
        t.setdefault(k1, {})
        for k2 in zh[1:]:
            t[k1][k2] = random.random()

    norm_dd(t)
    norm_dd(e)

    hmm = (t, e)

    for P in range(10): #make a better stopping criterium
        F,ps = forward_chain_log(s, hmm)
        B,_ = backward_chain_log(s, hmm)
        #we have a single training sequence

        #make a compatible structure, but with pseudocounts
        A = {}
        for a in t:
            A[a] = {}
            for b in t[a]:
                if a == 0: #need something positive for zero state
                    A[a][b] = max(1.,r)
                else:
                    A[a][b] = r
        E = {}
        for a in e:
            E[a] = {}
            for b in e[a]:
                E[a][b] = r
    
        for i in range(0,len(s)):
            for k in state_set:
                for l in state_set:
                    #niz je zamaknjen
                    A[k][l] += math.exp(F[i][k] + logmv(t[k][l]) + logmv(e[l][s[i]]) + B[i+1][l] - ps)

        for i in range(1,len(s)+1):
            for k in state_set:
                E[k][s[i-1]] += math.exp(F[i][k] + B[i][k] - ps)
        
        norm_dd(A)
        norm_dd(E)
        hmm = (A,E)
    
    return hmm

def baum_welch_no_log(s, num_hidden, r=0.0):
    zs = set(s)
    state_set = [str(x) for x in range(1, num_hidden+1)]
    zh = [0] + state_set

    # random model
    e = {}
    for k in zh[1:]:
        e[k] = dict((k, random.random()) for k in zs)
    t = {}
    t[0] = {}
    for k in zh[1:]:
        t[0][k] = random.random()
    for k1 in zh[1:]:
        t.setdefault(k1, {})
        for k2 in zh[1:]:
            t[k1][k2] = random.random()

    norm_dd(t)
    norm_dd(e)

    hmm = (t, e)

    for P in range(10): #make a better stopping criterium
        F,ps = forward_chain(s, hmm)
        B,_ = backward_chain(s, hmm)
        #we have a single training sequence

        #make a compatible structure, but with pseudocounts
        A = {}
        for a in t:
            A[a] = {}
            for b in t[a]:
                if a == 0: #need something positive for zero state
                    A[a][b] = max(1.,r)
                else:
                    A[a][b] = r
        E = {}
        for a in e:
            E[a] = {}
            for b in e[a]:
                E[a][b] = r
    
        for i in range(0,len(s)):
            for k in state_set:
                for l in state_set:
                    #niz je zamaknjen
                    A[k][l] += F[i][k]*t[k][l]*e[l][s[i]]*B[i+1][l]/ps

        for i in range(1,len(s)+1):
            for k in state_set:
                E[k][s[i-1]] += F[i][k]*B[i][k]/ps
        
        norm_dd(A)
        norm_dd(E)
        hmm = (A,E)
    
    return hmm


def izpisi_model(hmm):
    print hmm
    t, e = hmm
    # izpisi verjetnosti prehodov med skritimi stanji (t)

    # izpisi verjetnosti oddaje simbolov v skritih stanjih (e)
    pass

def oceni_ujemanje2(z1, z2):
    # oceni podobnost med zaporedjema, ceprav sestavljena iz razlicnih simbolov
    pass

h = "gggggggggppppppppppp" # dejansko skrito zaporedje, ki ga uporabimo le na
# koncu za ocenjevanje
s = "CCCGCCCCCGGGCCCGGCCG"


# staro: za napoved uporabimo znani HMM
_, _, _, hn_pravi = viterbi_log(s, HMM)

# novo: HMM se naucimo iz vidnega zaporedja
naucen_hmm = viterbi_training(s, 2)
print naucen_hmm

naucen_hmm = baum_welch(s, 2)
print naucen_hmm

# napovej zaporedja skritih stanj
_, _, _, hn_naucen = viterbi_log(s, naucen_hmm)
hn_naucen = post_decoding_log(s, naucen_hmm)


# rezultati
print "PRAVI MODEL"
izpisi_model(HMM)
print "NAUCEN MODEL"
izpisi_model(naucen_hmm)
print

print "dejansko zaporedje h:", h
print "napovedano zaporedje h (pravi HMM):", hn_pravi
print oceni_ujemanje2(h, hn_pravi)
print "napovedano zaporedje h (naucen HMM):", hn_naucen
print oceni_ujemanje2(h, hn_naucen)

