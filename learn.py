from hmmlearn import *
from fasta import *
from hmm_graph import HMM_graph

training = map(Fasta, ['1rsy','1tjm','1tjx', '1uov'])

hmm = HMM_graph()

dist = hmm.frequency_distribution(seq2int('DDQQDQDSDD'))

hmm.add_linear(1, dist)

model = hmm.get_model()

X, Z = model.sample(200)
print X
