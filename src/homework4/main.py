from collections import defaultdict
from math import log

def logr(a):
    min_val = 0.000000000001
    return log(max(a, min_val))*-1

def get_hmm(s, h):
    
    def normalize(states):
        for key, value in states.iteritems():
            n = float(sum(value.values()))
            for key2, value2 in value.iteritems():
                states[key][key2] = value2 / n
        return states
    
    states = set(h)
    t = {}; e = {}
    t[0] = defaultdict(float)
    for a in states:
        t[a] = defaultdict(float)
        e[a] = defaultdict(float)
        
    prev_hs = 0
    for hs, es in zip(h, s):
        e[hs][es] += 1
        t[prev_hs][hs] += 1
        prev_hs = hs

    for a in states:
        t[0][a] = h.count(a)
    return (normalize(t), normalize(e))

def forward(obs, states, start_p, trans_p, emit_p):
    V = [{}]
 
    # Initialize base cases (t == 0)
    for y in states:
        V[0][y] = start_p[y] * emit_p[y][obs[0]]
 
    # Run Viterbi for t > 0
    for t in range(1,len(obs)):
        V.append({})
        for y in states:
            prob = sum([V[t-1][y0] * trans_p[y0][y] * emit_p[y][obs[t]] for y0 in states])
            V[t][y] = prob
 
    prob = sum([V[len(obs) - 1][y] for y in states])
    return prob

def viterbi(obs, states, start_p, trans_p, emit_p):
    V = [{}]
    path = {}
 
    # Initialize base cases (t == 0)
    for y in states:
        V[0][y] = logr(start_p[y]) + logr(emit_p[y][obs[0]])
        path[y] = [y]
    # Run Viterbi for t > 0
    for t in range(1,len(obs)):
        V.append({})
        newpath = {}
 
        for y in states:
            (prob, state) = min([(logr(V[t-1][y0]) + logr(trans_p[y0][y]) + logr(emit_p[y][obs[t]]), y0) for y0 in states])
            V[t][y] = prob
            newpath[y] = path[state] + [y]
 
        # Don't need to remember the old paths
        path = newpath
        
    (prob, state) = min([(V[len(obs) - 1][y], y) for y in states])
    return (prob, path[state])

def viterbi_log(obs, states, start_p, trans_p, emit_p):
    V = [{}]
    path = {}
 
    # Initialize base cases (t == 0)
    for y in states:
        V[0][y] = start_p[y] * emit_p[y][obs[0]]
        path[y] = [y]
 
    # Run Viterbi for t > 0
    for t in range(1,len(obs)):
        V.append({})
        newpath = {}
 
        for y in states:
            (prob, state) = max([(V[t-1][y0] * trans_p[y0][y] * emit_p[y][obs[t]], y0) for y0 in states])
            V[t][y] = prob
            newpath[y] = path[state] + [y]
 
        # Don't need to remember the old paths
        path = newpath
 
    (prob, state) = max([(V[len(obs) - 1][y], y) for y in states])
    return (prob, path[state])

if __name__ == '__main__':
    
    train = open('train.csv').readlines()
    seq = "".join(t.split(',')[1] for t in train[1:])
    hid = "".join(t.split(',')[0] for t in train[1:])
    
    transition, emit = get_hmm(seq, hid)
    begin = transition[0]
    test = open('test.csv').readlines()
    tests = tuple(t.split(',')[0] for t in test[1:])
    results = viterbi(tests, ['1','2','3'], begin, transition, emit)
    open('result_test.csv', 'w').write("\n".join(results[1]))
    print "Done"
    