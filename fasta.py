import urllib2

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
