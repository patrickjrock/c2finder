from hmmlearn import hmm
import numpy as np
model = hmm.MultinomialHMM(n_components=3)
model.startprob_ = np.array([0.6, 0.3, 0.1])
model.transmat_ = np.array([[0.7, 0.2, 0.1],
 [0.3, 0.5, 0.2],
 [0.3, 0.3, 0.4]])
model.emissionprob_ = np.array([[.5,.5],[.7,.3],[.2,.8]])
X, Z = model.sample(100)
print X, Z

