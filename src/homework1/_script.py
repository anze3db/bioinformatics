"""Sample script that observes dinucleotide frequency in a sequence and compares it to
 the expected frequency for sequence generated from multinomial model (log odds ratio).
 We permute sequences to derive null-distribution of this statistics."""

import collections
import random
from matplotlib import pyplot as plt
from math import log

log2 = lambda x: log(x, 2)

def count_kmers(s, k):
    """Count k-mers in string s and return a dictionary with counts."""
    d = collections.defaultdict(int)
    for i in xrange(len(data)-(k-1)):
        d[s[i:i+k]] += 1
    return d

def permutation_test(s, n=100):
    """Assess odds ratio on permuted string and return a list with scores."""
    res = []
    for _ in range(n):
        random.shuffle(s)
        f = count_kmers("".join(s), 2)
        res.append(log2(N * f[a + b] / float(f1[a] * f1[b])))
    return res

from main import get_seq

data = "".join(get_seq())

N = float(len(data))
f2 = count_kmers(data, 2)
f1 = count_kmers(data, 1)

a = "G"; b = "C"

odds_ratio = log2(N * f2[a + b] / float(f1[a] * f1[b]))
print odds_ratio

list_data = list(data)
res = permutation_test(list_data)
p = sum(1 for r in res if r>odds_ratio) / float(len(res))
print "OR=%5.3f, p=%7.5f" % (odds_ratio, p)

plt.close()
plt.hist(res)
plt.title("Distribution of log odds ratios for dinucleotide %s%s\n" % (a,b) +
          "in a sequence generated from a multinomial model")
plt.show()