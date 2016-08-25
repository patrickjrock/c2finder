import urllib2
from Bio import pairwise2
from sklearn import svm
import sklearn
from sklearn.feature_selection import RFE
from sklearn.linear_model import PassiveAggressiveClassifier as pac
import numpy as np
import sys

from fasta import Fasta

hydrophobicity = {
  'L': 97,
  'I': 99,
  'F': 100,
  'W': 97,
  'V': 76,
  'M': 74,
  'C': 49,
  'Y': 63,
  'A': 41,
  'T': 13,
  'E': -31,
  'G': 0,
  'S': -5,
  'Q': -10,
  'D': -55,
  'R': -14,
  'K': -23,
  'N': -28,
  'H': 8,
  'P': -46,
}

#5fq9
data = [
  '2nsq',
  #'319b',
  '2dmh',
  '2Fk9',
  '2ep6',
  '2cm5',
  '2d8k',
  '2uzp',
  '1rh8',
  '1dsy',
  '1ugk',
  '2chd',
  '1w15',
  '3rpb',
  '2zkm',
  '1e8y',
  '1rlw',
  '1a25',
  '1rsy',
  '1tjx', 
]

target = [1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1]

class Align:
  """computes alignment to a template"""
  def __init__(self, pdb_code):
    self.template = Fasta("4wee")
    self.template.seqs[0] = self.template.seqs[0][26:]

  def get_alignment(self, f1, f2):
    align1 = pairwise2.align.localms(f1.seqs[0], f2.seqs[0], 2, -1, -5, -1)
    return align1

  def get_features(self, pdb_code):
    query = Fasta(pdb_code)
    align = self.get_alignment(self.template, query)
    
    features = []
    for i, c in enumerate(align[0][0]):
      if c == '-':
        pass
      else: 
        amino_acid = align[0][1][i]
        try:
          features.append(hydrophobicity[amino_acid])
        except:
          features.append(0)
    return np.array(features)

class Clf: 



  def __init__(self, d, t, name="passive_aggressive", n=10):
    self.data = d
    self.target = t
    self.align = Align("4wee") 
    self.tdata = np.array(map(lambda x:self.align.get_features(x), self.data))
    

    clf = pac()

    #rfe = RFE(estimator = clf, n_features_to_select=n)
    clf.fit(self.tdata, self.target)
    self.model = clf

  def predict(self, pdb_code):
    features = self.align.get_features(pdb_code)
    if self.model.predict(features.reshape(1,-1)) == np.array([1]):
      return "yes"
    else:
      return "no"
  
def get_data(fname, targ):
  f = open(fname)
  adata = []
  atarget = []
  for i, line in enumerate(f):
      adata.append(line[:-1])
      atarget.append(targ)
  f.close()
  assert len(adata) == len(atarget)
  return adata, atarget

sys.stderr.write("getting data...\n")
rand_data, rand_target = get_data('no.txt', 0)
fdata = data + rand_data
ftarget = target + rand_target

print fdata
print ftarget

sys.stderr.write("building model...\n")
c = Clf(fdata, ftarget)

rcsb = open(sys.argv[1])
for i, line in enumerate(rcsb):
  try:
    print line[:-1] + ' ' + c.predict(line[:-1])
  except IndexError:
    print "no"
rcsb.close()
