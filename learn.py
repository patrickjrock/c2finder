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


for i in range(0,21):
  print [i] in training[0].hmmseq().tolist()


training = [t.hmmseq() for t in training]


X = np.concatenate(training)
model.fit(X)


X, Z = model.sample(200)
print X,Z
