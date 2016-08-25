import pickle
import urllib2
from Bio import pairwise2
from sklearn import svm
import sklearn
from sklearn.feature_selection import RFE

from sklearn.linear_model import PassiveAggressiveClassifier as pac
from sklearn.linear_model import Perceptron as perceptron

import numpy as np
import sys
import random

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
    try:
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
    except:
      return np.array([0 for i in range(0,len(self.template.seqs[0]))] )

class Classifier: 
  def update_model(self, data, target):
      tdata = np.array(map(lambda x:self.align.get_features(x), data))
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
    self.model = svm.SVC()
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
  rand_data, rand_target = get_data('no.txt', 0, 100)
  c2_data, c2_target = get_data('c2.txt', 1)
  fdata = data + rand_data + c2_data
  ftarget = target + rand_target + c2_target
  
  c = Classifier(fdata, ftarget)
  return c

def test_pickle():
  return Classifier() 

def test(c, f):
  rcsb = open(f)
  for i, line in enumerate(rcsb):
    try:
      print c.predict(line[:-1])
    except IndexError:
      print "IndexError"
  rcsb.close()
