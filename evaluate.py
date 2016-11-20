# Run cross validation stuff and evaluate bias and variance at some point
from hmmlearn import *
from hmmlearn.hmm import MultinomialHMM
from fasta import *
import random
from numpy import inf

f = file('data/test_c2.txt')
ps = f.readlines()
f.close()

f = file('data/no.txt')
allns = f.readlines()
ns = [allns[i] for i in random.sample(range(0,len(allns)), 78)]
f.close()


def run(model, code):
  protein = Fasta(code)
  querey = protein.hmmseq()
  pred = model.predict(querey)
  pred = zip(protein.seq, pred)
  score = model.score(querey)
  return pred, score

def test(model, code):
  """ return score for model given the rcsb code """
  protein = Fasta(code)
  querey = protein.hmmseq()
  return model.score(querey)

def evaluate(model, positives=ps, negatives=ns):
  """evaluate the model and report number of tp, fp, tn, fn"""
  presult = [(test(model, p), p) for p in positives]
  nresult = [(test(model, n), n) for n in negatives]
  
  print "probability class code"
  for p in presult:
    print str(p[0]) + " c2 " + p[1][:-1]
  for n in nresult:
    print str(n[0]) + " not " + n[1][:-1]

def search(model, query, window_size=150):
  """search the query with a scrolling window"""
  maxi = -1
  maxs = -inf
  maxq = ""

  if len(query) < window_size:
    return -1, model.score(query)

  for i in range(0,len(query)-window_size+1):
    subquery = query[i:i+window_size]
    score = model.score(subquery)
    if score > maxs:
      maxs = score
      maxi = i
      maxq = subquery
  maxp = model.predict(maxq) 
  maxq = [x[0] for x in maxq]
  return maxi, maxs, zip(int2seq(maxq), maxp)


