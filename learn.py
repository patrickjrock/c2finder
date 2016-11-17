from hmmlearn import *
from fasta import Fasta
from hmm_graph import HMM_graph

training = map(Fasta, ['1rsy','1tjm','1tjx', '1uov'])

hmm = HMM_graph()
hmm.add_linear(200)
model = hmm.get_model()
print model.sample(100)
