from hmmlearn import *
from fasta import Fasta
from hmm_graph import HMM_graph

training = map(Fasta, ['1rsy','1tjm','1tjx', '1uov'])

hmm = HMM_graph()
hmm.add_linear(10)
hmm.add_loop(10)
hmm.add_linear(10)
model = hmm.get_model()

print model.sample(40)
