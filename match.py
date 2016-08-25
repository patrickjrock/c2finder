import regex
import sys 
from fasta import Fasta
from numpy.random import shuffle

regions = []


def preprocess(s):
  hydrophobic = '[LFIVCMAGTWPG]'
  hydrophilic = '[KRDNHQYESYGMA]'
  s = s.replace('h', hydrophobic)
  s = s.replace('p', hydrophilic)
  return s

regions.append('phphphpphp.{3,16}')
regions.append('p[DQS][IRLKPM]p.{3,11}')
regions.append('[DKS]hphphph.{12,25}')
regions.append('phphphphphp.{7,35}')
regions.append('hphphp[D]p[ERDN]p[ERDN]p[MLYFV]ppp[RADEP]hhphph.{10,40}')
regions.append('hphphph')



pattern = '(' + ''.join(regions) + ')' 
pattern = preprocess(pattern)
#pattern = '(' + s.replace("H", hydrophobic) + ')'

print ">pattern " + pattern
print ">allowing " + str(int(100*float(sys.argv[1])/42)) + "% wobble\n"

def match(fasta, error):
  fuzzy = '{d<2,s<' + error + '}'
  reg = regex.search(pattern + fuzzy, fasta)
  print reg

f = open(sys.argv[2], 'r')
lines = [l for l in f]
shuffle(lines)
for line in lines:
  fa = Fasta(line[:-1])
  print "Query>" + fa.name
  sys.stderr.write(fa.name + '\n')
  match(fa.seqs[0], sys.argv[1])

f.close()
