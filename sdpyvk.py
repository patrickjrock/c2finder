from align import *
from hmm_graph import *
from hmmlearn.hmm import MultinomialHMM
from fasta import *
import numpy as np
import sys
from evaluate import *

ts = transitions('fasta/out.fasta')
es = emissions('fasta/out.fasta')

# initial distribution
I = [1,0,0,0,0,0,0,0]

# generate emission matrix
E = []
for e in es:
  dist = HMM_graph.frequency_distribution(seq2int(e))
  E.append(dist)

# generate transition matrix
A = np.zeros((8,8))
for t in ts:
  A[t[0]-1,t[1]-1] += 1
A = [HMM_graph.norm(row) for row in A]

hmm = MultinomialHMM(8)
hmm.startprob_ = np.array(I)
hmm.emissionprob_ = np.array(E)
hmm.transmat_ = A


print evaluate(hmm)
