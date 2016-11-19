import urllib2
import numpy as np

aadict = { 'A' : 0,
	   'R' : 1,
	   'N' : 2,
	   'D' : 3,
	   'C' : 4,
	   'Q' : 5,
	   'E' : 6,
	   'G' : 7,
	   'H' : 8,
	   'I' : 9,
	   'L' : 10,
	   'K' : 11,
	   'M' : 12,
	   'F' : 13,
	   'P' : 14,
	   'S' : 15,
	   'T' : 16,
	   'W' : 17,
	   'Y' : 18,
	   'V' : 19,
	   'U' : 20,
	   'X' : 21, }

inv_aadict = {v: k for k, v in aadict.iteritems()}

def aa2int(aa):
  return aadict[aa]

def int2aa(i):
  return inv_aadict[i]

def int2seq(i):
  return [int2aa(c) for c in i]

def seq2int(seq):
  return [aa2int(c) for c in seq]

class Fasta:
  """encapsulates FASTA fetching from rcsb and parsing"""
  def parse_record(self, record):
    lines = record.split('\n')
    return "".join(lines[1:])
    
  def __init__(self, pdb_code):
    url = 'http://www.rcsb.org/pdb/files/fasta.txt?structureIdList=' + pdb_code
    response = urllib2.urlopen(url)
    raw = response.read()
    records = raw.split('>')[1:] 
    self.seqs = map(self.parse_record, records)
    self.name = pdb_code
    self.seq = self.seqs[0]

  def seq2int(self):
    return seq2int(self.seqs[0])

  def hmmseq(self):
    seq = self.seq2int()
    return np.array([[x] for x in seq])


