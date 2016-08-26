import pickle
import urllib2
from Bio import pairwise2
from sklearn import svm
import sklearn
from sklearn.feature_selection import RFE

from sklearn.ensemble import RandomForestClassifier as rfc
from sklearn.linear_model import PassiveAggressiveClassifier as pac
from sklearn.linear_model import Perceptron as perceptron

import numpy as np
import sys
import random

from fasta import Fasta

hydrophobicity = {
  'A': -0.519,
  'C': -1.343,
  'D': 1.050,
  'E': 1.357,
  'F': -1.006,
  'G': -0.384,
  'H': 0.336,
  'I': -1.239,
  'K': 1.831,
  'L': -1.019,
  'M': -0.663,
  'N': 0.945,
  'P': 0.189,
  'Q': 0.931,
  'R': 1.538,
  'S': -0.228,
  'T': -0.032,
  'V': -1.337,
  'W': -0.595,
  'Y': 0.260
}

structure = {
  'A': -1.302,
  'C': 0.465,
  'D': 0.302,
  'E': -1.453,
  'F': -0.590,
  'G': 1.652,
  'H': -0.417,
  'I': -0.547,
  'K': -0.561,
  'L': -0.987,
  'M': -1.524,
  'N': 0.828,
  'P': 2.081,
  'Q': -0.179,
  'R': -0.055,
  'S': 1.399,
  'T': 0.326,
  'V': -0.279,
  'W': 0.009,
  'Y': 0.830
}

bulk = {
  'A': -0.733,
  'C': -0.862,
  'D': -3.656,
  'E': 1.477,
  'F': 1.891,
  'G': 1.330,
  'H': -1.673,
  'I': 2.131,
  'K': 0.553,
  'L': -1.505,
  'M': 2.219,
  'N': 1.299,
  'P': -1.628,
  'Q': -3.005,
  'R': 1.502,
  'S': -4.760,
  'T': 2.213,
  'V': -0.544,
  'W': 0.672,
  'Y': 3.097
}

composition = {
  'A': 1.570,
  'C': -1.020,
  'D': -0.259,
  'E': 0.113,
  'F': -0.397,
  'G': 1.045,
  'H': -1.474,
  'I': 0.393,
  'K': -0.277,
  'L': 1.266,
  'M': -1.005,
  'N': -0.169,
  'P': 0.421,
  'Q': -0.503,
  'R': 0.440,
  'S': 0.670,
  'T': 0.908,
  'V': 1.242,
  'W': -2.128,
  'Y': -0.838
}


charge = {
  'A': -0.146,
  'C': -0.255,
  'D': -3.242,
  'E': -0.837,
  'F': 0.412,
  'G': 2.064,
  'H': -0.078,
  'I': 0.816,
  'K': 1.648,
  'L': -0.912,
  'M': 1.212,
  'N': 0.933,
  'P': -1.392,
  'Q': -1.853,
  'R': 2.897,
  'S': -2.647,
  'T': 1.313,
  'V': -1.262,
  'W': -0.184,
  'Y': 1.512
}

#5fq9
data = [
  '2nsq',
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

  def get_alignment(self, s1, s2):
    align1 = pairwise2.align.localms(s1, s2, 2, -1, -5, -1)
    return align1

  def get_features_from_string(self, s1, s2):
      features = []
      for i, c in enumerate(s1):
        if c == '-':
          pass
        else: 
          amino_acid = s2[i]
          try:
            features.append(structure[amino_acid])
            features.append(hydrophobicity[amino_acid])
          except:
            features.append(0)
            features.append(0)
      return features

  def get_features_rolling(self, pdb_code):
      query = Fasta(pdb_code)
      qstr = query.seqs[0]
      if len(qstr) < 150:
        align = self.get_alignment(self.template.seqs[0], qstr)
        return [self.get_features_from_string(align[0][0], align[0][1])]
      else:
        fv = []
        for i in range(0, len(qstr)-150): 
          align = self.get_alignment(self.template.seqs[0], qstr[i : i+150])
          fv.append(self.get_features_from_string(align[0][0], align[0][1]))
        return fv

  def get_features(self, pdb_code):
    try:
      query = Fasta(pdb_code)
      align = self.get_alignment(self.template.seqs[0], query.seqs[0])

      return np.array(self.get_features_from_string(align[0][0], align[0][1]))    
    except:
      return np.array([0 for i in range(0,2*118)])

class Classifier: 
  def update_model(self, data, target):
      tdata = np.array(map(lambda x:self.align.get_features(x), data))
      #rfe = RFE(estimator=self.model, n_features_to_select=100, step=1)
      self.model.fit(tdata, target)

  def build_out_of_core(self):
    n = int(len(self.data) / 100)
    split_data = np.array_split(self.data, n)
    split_target = np.array_split(self.data, n)
    for i in range(0,len(split_data)):
      self.update_model(split_data[i], split_target[i])      

  def build_model(self, d, t):
    self.data = d
    self.target = t
    self.align = Align("4wee") 
    self.model = rfc(n_estimators = 100)
    #if len(self.data) > 50000:
    #  self.build_out_of_core()
    #else:    
    self.update_model(d, t)

  def save_model(self, fname):
    s = pickle.dump(self.model, open('pickle.dat','wb'))
 
  def load_model(self, fname):  
    self.model = pickle.load(open('pickle.dat', 'rb'))

  def __init__(self, d=None, t=None):
    if d == None or t == None:
      self.load_model('pickle.data')
    else:
      self.build_model(d,t)
      self.save_model('pickle.data')

  def predict(self, pdb_code):
    align = Align('4wee')
    features = align.get_features(pdb_code)
    return self.model.predict(features.reshape(1,-1))

  def locate(self, pdb_code):
    align = Align('4wee')
    fv =  align.get_features_rolling(pdb_code)
    fv = map(lambda x : np.array(x).reshape(1,-1), fv)
    return map(self.model.predict, fv)    
 
def get_data(fname, targ, n=2):
  f = open(fname)
  adata = []
  atarget = []
  for i, line in enumerate(f):
    if random.randint(1,n) == 1:
      adata.append(line[:-1])
      atarget.append(targ)
  f.close()
  assert len(adata) == len(atarget)
  return adata, atarget

def test_build():
  rand_data, rand_target = get_data('training_no.txt', 0, 500)
  c2_data, c2_target = get_data('training_c2.txt', 1)
  fdata = data + rand_data + c2_data
  ftarget = target + rand_target + c2_target
  
  c = Classifier(fdata, ftarget)
  return c

def test_pickle():
  return Classifier() 

def evaluate(c, f):
  rcsb = open(f)
  yes = 0
  no = 0
  for i, line in enumerate(rcsb):
    try:
      p = c.predict(line[:-1])
      if p == np.array([1]):
        yes = yes + 1
      else:
        no = no + 1
    except IndexError:
      print "IndexError"
  rcsb.close()
  return (yes, no)
