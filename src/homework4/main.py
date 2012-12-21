from collections import Counter, defaultdict
from random import random
from math import log,exp
import textwrap

def log_sum(v1, v2):
    nv1 = max(v1, v2)
    nv2 = min(v1, v2)
    return nv1 + log(1.0+exp(nv2-nv1))

def sum_log(vals):
    s = vals[0]
    for v in vals[1:]:
        s = log_sum(s, v)
    return s
def backward_chain_log(s, hmm):
    t, e = hmm

    # seznam skritih stanj
    zh = e.keys()
    zh = [0] + list(zh)

    # Create table
    b = [{} for i in range(len(s)+1)]
    # Initialize i = 0; b_k(L) = 1.0 for all k
    for k in zh: 
        b[len(s)][k] = log(1.0)
    # Recursion (i=1..L):
    for i in range(len(s)-1,0,-1):
        for k in zh:
            sum_val = sum_log([logmv(t[k].get(l, 0.0)) + logmv(e.get(l, {}).get(s[i], 0.0)) + b[i+1][l] for l in zh])
            b[i][k] = sum_val

    # P(x)
    ps = sum_log([logmv(t[0].get(l, 0.0)) + logmv(e.get(l, {}).get(s[0], 0.0)) + b[1][l] for l in zh])
    return b, ps

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
    f[0][0] = log(1.0)

    # Recursion (i=1..L):
    for i in range(1, len(s)+1):
        for l in zh:
            sum_val = sum_log([f[i-1][k] + logmv(t[k].get(l, 0.0)) for k in zh])
            f[i][l] = logmv(e.get(l, {}).get(s[i-1], 0.0)) + sum_val
    # P(x)
    ps = sum_log([f[len(s)][k] for k in zh])
    return f, ps

def baum_welch(s, num_hidden, r=0.0):
    zs = set(s)
    state_set = [str(x) for x in range(1, num_hidden+1)]
    zh = [0] + state_set

    # random model
    e = {}
    for k in zh[1:]:
        e[k] = dict((k, random()) for k in zs)
    t = {}
    t[0] = {}
    for k in zh[1:]:
        t[0][k] = random()
    for k1 in zh[1:]:
        t.setdefault(k1, {})
        for k2 in zh[1:]:
            t[k1][k2] = random()

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
                    A[k][l] += exp(F[i][k] + logmv(t[k][l]) + logmv(e[l][s[i]]) + B[i+1][l] - ps)

        for i in range(1,len(s)+1):
            for k in state_set:
                E[k][s[i-1]] += exp(F[i][k] + B[i][k] - ps)
        
        norm_dd(A)
        norm_dd(E)
        hmm = (A,E)
    
    return hmm

def norm_dd(d):
    """ normalize a dict of dict structure """
    for a in d:
        sa = float(sum(d[a].values()))
        for b in d[a]:
            d[a][b] = d[a][b]/sa

def logmv(a):
    min_val = 0.0000000001
    return log(max(a, min_val))

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

def viterbi_training(s, num_hidden):
    zs = set(s)
    zh = [0] + [str(x) for x in range(1, num_hidden+1)]

    # random model
    e = {}
    for k in zh[1:]:
        e[k] = dict((k, random()) for k in zs)
    t = {}
    t[0] = {}
    for k in zh[1:]:
        t[0][k] = random()
    for k1 in zh[1:]:
        t.setdefault(k1, {})
        for k2 in zh[1:]:
            t[k1][k2] = random()
    norm_dd(t)
    norm_dd(e)
    hmm = (t, e)

    for i in range(10):
        _, _, _, hn = viterbi_log(s, hmm)
        hmm = zgradi_hmm(s, hn)

    return hmm


        



if __name__ == '__main__':
    
    train = open('train.csv').readlines()
    trains = "".join(t.split(',')[1] for t in train[1:])
    inner = "".join(t.split(',')[0] for t in train[1:])
    currState = train[1].split(',')[0]
    states = defaultdict(int)
    
    naucen_hmm = viterbi_training(inner, 3)
    _, _, _, hn_naucen = viterbi_log(trains, naucen_hmm)
    
    
    print naucen_hmm
    
    print hn_naucen[2000:2100]
    print inner[2000:2100]
    
    
    for t in range(2, len(train)-1):
        l = train[t].split(',')
        lprev = train[t-1].split(',')
        states[lprev[0], l[0]] += 1
        

    for s in states:
        states[s] = states[s]/float(len(train)-1)
        

    
    cnt = Counter()
    cnt.update(t.split(',')[0] for t in train)
    print cnt.most_common(4)
    
    test = open('test.csv').readlines()
    
    
    tests = "".join(t.split(',')[0] for t in test[1:])
    print "tests: " + tests[:100]
    _, _, _, hn_naucen = viterbi_log(tests, naucen_hmm)
    
    open('result_test.csv', 'w').write("\n".join(textwrap.wrap(hn_naucen, width=1)))
    
    
    result = ""
    
    
    for t in range(1, len(test)):
        l = test[t].split(',')
        if currState != l[0]:
            currState = l[0]
        result += "2\n"#test[t].split(',')
    print test[0].split(',')
    open('result.csv', 'w').write(result)
    
    result = open('result.csv').readlines()
    