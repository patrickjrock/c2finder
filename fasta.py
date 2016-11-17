import urllib2

aadict = { 'A' : 1,
	   'R' : 2,
	   'N' : 3,
	   'D' : 4,
	   'C' : 5,
	   'Q' : 6,
	   'E' : 7,
	   'G' : 8,
	   'H' : 9,
	   'I' : 10,
	   'L' : 11,
	   'K' : 12,
	   'M' : 13,
	   'F' : 14,
	   'P' : 15,
	   'S' : 16,
	   'T' : 17,
	   'W' : 18,
	   'Y' : 19,
	   'V' : 20, 
	   'B' : 21,
	   'Z' : 22 }

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

  def aa2int(self, aa):
    return aadict[aa]

  def seq2int(self):
    return [self.aa2int(c) for c in self.seqs[0]]

def test():
  p = Fasta('4wee')
  print p.seq2int()
  print p.seqs[0]

