# Run cross validation stuff and evaluate bias and variance at some point
from hmmlearn import *
from hmmlearn.hmm import MultinomialHMM
from fasta import Fasta
import random

f = file('data/test_c2.txt')
ps = f.readlines()
f.close()

f = file('data/no.txt')
allns = f.readlines()
ns = [allns[i] for i in random.sample(range(0,len(allns)), 78)]
f.close()


def test(model, code):
  """ return score for model given the rcsb code """
  protein = Fasta(code)
  querey = protein.hmmseq()
  pred = model.predict(querey)
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
