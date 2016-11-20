from align import *
from hmm_graph import *
from hmmlearn.hmm import MultinomialHMM
from fasta import *
import numpy as np
import sys
from evaluate import *

# set up labels for supervised learning
labels = []
for i in range(1,53):
  labels.append(1)
labels.append(2)
labels.append(3)
for i in range(55,68):
  labels.append(4)
labels.extend([5,6,7,8,9,10])
for i in range(74,211):
  labels.append(11)

ts = transitions('fasta/out.fasta', labels)
es = emissions('fasta/out.fasta', labels)

# initial distribution
I = [1,0,0,0,0,0,0,0,0,0,0]

# generate emission matrix
E = []
for e in es:
  dist = HMM_graph.frequency_distribution(seq2int(e))
  E.append(dist)

# generate transition matrix
A = np.zeros((11,11))
for t in ts:
  A[t[0]-1,t[1]-1] += 1
A = [HMM_graph.norm(row) for row in A]

hmm = MultinomialHMM(11)
hmm.startprob_ = np.array(I)
hmm.emissionprob_ = np.array(E)
hmm.transmat_ = A

#proetin = Fasta(sys.argv[1])
#query = protein.hmmseq()

query = """MAHHHHHHGTALQLEPPTVVETLRRGSKFIKWDEETSSRNLVTLRVDPNGFFLYWTGPNMEVDTLDISSIRDTRTGRYAR
LPKDPKIREVLGFGGPDARLEEKLMTVVSGPDPVNTVFLNFMAVQDDTAKVWSEELFKLAMNILAQNASRNTFLRKAYTK
LKLQVNQDGRIPVKNILKMFSADKKRVETALESCGLKFNRSESIRPDEFSLEIFERFLNKLCLRPDIDKILLEIGAKGKP
YLTLEQLMDFINQKQRDPRLNEVLYPPLRPSQARLLIEKYEPNQQFLERDQMSMEGFSRYLGGEENGILPLEALDLSTDM
TQPLSAYFINSSHNTYLTAGQLAGTSSVEMYRQALLWGCRCVELDVWKGRPPEEEPFITHGFTMTTEVPLRDVLEAIAET
AFKTSPYPVILSFENHVDSAKQQAKMAEYCRSIFGDALLIEPLDKYPLAPGVPLPSPQDLMGRILVKNKKRPKKPTTDEG
TASSEVNATEEMSTLVNYIEPVKFKSFEAARKRNKCFEMSSFVETKAMEQLTKSPMEFVEYNKQQLSRIYPKGTRVDSSN
YMPQLFWNVGCQLVALNFQTLDVAMQLNAGVFEYNGRSGYLLKPEFMRRPDKSFDPFTEVIVDGIVANALRVKVISGQFL
SDRKVGIYVEVDMFGLPVDTRRKYRTRTSQGNSFNPVWDEEPFDFPKVVLPTLASLRIAAFEEGGKFVGHRILPVSAIRS
GYHYVCLRNEANQPLCLPALLIYTEASDYIPDDHQDYAEALINPIKHVSLMDQRARQLAALIGESEAQAGQET"""

print search(hmm, hmmseq(query.replace('\n','')))
