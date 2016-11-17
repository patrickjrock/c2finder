from hmmlearn import *
from fasta import *
from hmm_graph import HMM_graph
import numpy as np

training = map(Fasta, ['1rsy','1tjm','1tjx', '1uov'])

hmm = HMM_graph()

dist = hmm.frequency_distribution(seq2int('DDQQDQDSDD'))

hmm.add_loop(30)
hmm.add_linear(1, dist)
hmm.add_jump(200)

model = hmm.get_model()

training = [t.hmmseq() for t in training]

Xt = np.concatenate(training)
ls = map(len, training)
#model.fit(Xt, lengths=ls)


X, Z = model.sample(200)
print X,Z
