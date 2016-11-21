from align import *
from hmm_graph import *
from hmmlearn.hmm import MultinomialHMM
from fasta import *
import numpy as np
import sys
from evaluate import *
from Bio import SeqIO
from unsupervised import *
from utils import *

def remove_dead_states(A):
  """transition matrix can have states with no outgoing edges
     make these point to the end state for now, should reshape 
     matrix and stuff later"""
  for a in A:
    if sum(a) == 0:
      a[len(a)-1] = 1
  return A

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
for i in range(109, 123):
  labels.append(37)
labels.extend(range(38,59))
for i in range(144,211):
  labels.append(59)
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


# Try unsupervised stuff
type1 = [s.seq.tostring() for s in SeqIO.parse(open('fasta/type1.fasta'), 'fasta')]
type2 = [s.seq.tostring() for s in SeqIO.parse(open('fasta/type2.fasta'), 'fasta')]

training_seqs = type1[:len(type1)/2] + type2[:len(type2)/2]
print training_seqs
training_seqs = map(hmmseq, training_seqs)

tA, tE = fit(hmm, training_seqs)

thmm = MultinomialHMM(max(labels))
thmm.startprob_ = np.array(I)
thmm.emissionprob_ = tE
thmm.transmat_ = remove_dead_states(tA)


#query = hmmseq("KECDRKFRVKIRGIDIPVLPRNTDLTVFVEANIQHGQQVLCQRRTSPKPFTEEVLWNVWLEFSIKIKDLPKGALLNLQIYCLLYYVNLLLIDHRFLLRRGEYVLHMWQISGFNADKLTSATNPDKENSMSISILLDN")
#maxi, maxs, pred = search(hmm, query)

#pretty_print_prediction(pred)
#print maxs

query = Fasta(sys.argv[1]).hmmseq()
index, score, pred = search(thmm, query)
pretty_print_prediction(pred)
print score
