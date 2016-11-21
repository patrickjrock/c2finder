from evaluate import *
from hmmlearn import *
import numpy as np
from hmm_graph import *


def pairwise(iterable):
  it = iter(iterable)
  a = next(it, None)

  for b in it:
    yield (a, b)
    a = b


def estimate_emission(predictions, model):
  n = len(model.emissionprob_)
  m = len(model.emissionprob_[0])
  E = [[] for i in range(0, n)]

  for p in predictions:
    for state in p:
      print state[1], state[0][0]
      E[state[1]].append(state[0][0])
  E = [HMM_graph.frequency_distribution(e) for e in E]
  return E  

def estimate_transition(predictions, model):
  n = len(model.transmat_)
  m = len(model.transmat_[0])
  A = np.zeros((n,m))
  
  for p in predictions:
    for first, second in pairwise(p):
      A[first[1], second[1]] += 1
  
  A = [HMM_graph.norm(a) for a in A]
  return A

def fit(model, seqs):
  """use viterbi to refine estimates of model parameters"""

  predictions = [zip(seq, model.predict(seq)) for seq in seqs]
 
  E = estimate_emission(predictions, model)
  A = estimate_transition(predictions, model)
  return A, E
