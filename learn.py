from align import *
from hmm_graph import *
from hmmlearn.hmm import MultinomialHMM
from fasta import *
import numpy as np
import sys
from evaluate import *

# set up labels for supervised learning
labels = []
for i in range(1,37):
  labels.append(1)

# first beta strand
labels.extend([2,3,4,5,6,7,8,9,10,11])
for i in range(47,53):
  labels.append(12)

# conserved binding area
labels.extend([13,14])

for i in range(55,68):
  labels.append(15)

# second beta strand SDPYVK
labels.extend([16,17,18,19,20,21,22,23,24])
for i in range(77,98):
  labels.append(25)

# third conserved strand
labels.extend([26,27,28,29,30,31,32,33,34,35,36])
for i in range(109,211):
  labels.append(37)
print labels

ts = transitions('fasta/out.fasta', labels)
es = emissions('fasta/out.fasta', labels)

# initial distribution
I = [0 for x in range(0,max(labels))]
I[0] = 1

# generate emission matrix
E = []
for e in es:
  dist = HMM_graph.frequency_distribution(seq2int(e))
  E.append(dist)

# generate transition matrix
A = np.zeros((max(labels),max(labels)))
for t in ts:
  A[t[0]-1,t[1]-1] += 1
A = [HMM_graph.norm(row) for row in A]

hmm = MultinomialHMM(max(labels))
hmm.startprob_ = np.array(I)
hmm.emissionprob_ = np.array(E)
hmm.transmat_ = A

protein = Fasta(sys.argv[1])
query = protein.hmmseq()

print run(hmm, sys.argv[1])
print search(hmm, query)
