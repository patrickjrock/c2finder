import networkx as nx
import numpy as np
from hmmlearn.hmm import MultinomialHMM

class HMM_graph:
  """ Class that represents hmm as a graph
      outputs matrix representations of graph for use in hmm libraries
  """
  def __init__(self):
    self.G = nx.DiGraph()

  @staticmethod
  def norm(raw):
    if sum(raw) == 0:
      return raw
    return [float(i)/sum(raw) for i in raw]

  def uniform_distribution(self):
    """returns a distribution with an equal probability of emitting each aa"""
    v = [1 for i in range(1,23)]
    return self.norm(v)

  @staticmethod
  def frequency_distribution(seq, smooth=0.1):
    """returns a distribution representing the frequency of aas in seq"""
    freq = [0 for i in range(1,23)]
    for aa in seq:
      freq[aa] += 1
    freq = HMM_graph.norm(freq)
    if 0 in freq:
      freq = [(1-smooth)*float(x) for x in freq]
      freq = [smooth/22+float(x) for x in freq] 
    return freq

  def get_max(self):
    return len(self.G.node)

  def add_linear(self, n, dist=None):
    """Add a linear block of size n"""
    if dist == None:
      dist = self.uniform_distribution()

    node_range = range(self.get_max(), self.get_max()+n)
    for i in node_range:
      self.G.add_node(i, e=dist) # abstract out the distribution later
    for i in node_range[:-1]:
      self.G.add_edge(i,i+1)
    if node_range[0] != 0: 
      self.G.add_edge(node_range[0]-1, node_range[0])

  def add_jump(self, n):
    """Add a jump block of size n"""
    node_range = range(self.get_max(), self.get_max()+n)
    for i in node_range:
        self.G.add_node(i, e=self.uniform_distribution()) # using uniform distribution for now
    for i in node_range[:-1]:
        self.G.add_edge(i, i+1)
    for i in node_range[2:]:
        self.G.add_edge(node_range[0], i)
    if node_range[0] != 0: 
      self.G.add_edge(node_range[0]-1, node_range[0])


  def add_loop(self, n):
    """Add a loop block of size n"""
    node_range = range(self.get_max(), self.get_max()+n)
    for i in node_range:
        self.G.add_node(i, e=self.uniform_distribution())
    for i in node_range[:-1]:
        self.G.add_edge(i, i+1) # Connect to forward node
        self.G.add_edge(i, i) # Connect node to itself
    self.G.add_edge(node_range[-1], node_range[-1])
    if node_range[0] != 0: 
      self.G.add_edge(node_range[0]-1, node_range[0])

  def get_transition(self):
    """ converts the adjacency matrix to the transition matrix by normalizing the rows,
        nodes may contain their own transition probability matrix, the i,j entry is the 
	probability from transitioning from i to j """
    A = nx.adjacency_matrix(self.G).todense().getA()
    A = [self.norm(row) for row in A]
    return np.array(A)

  def get_emission(self):
    """ returns the emission matrix for the hmm. If no distribution is present in a node assume uniform """
    # for now just assume everything is uniform
    E =  [self.G.node[i]['e'] for i in range(0, self.get_max())]
    return np.array(E)

  def get_start(self):
    """ returns the start probability matrix, just points to the first node"""
    S = [0 for i in range(0, self.get_max())]
    S[0] = 1
    return np.array(S)

  def get_model(self):
    """ returns a multinomial hmm"""
    model = MultinomialHMM(n_components=self.get_max(), params='e', init_params='')
    model.startprob_ = self.get_start()
    model.transmat_ = self.get_transition()
    model.emissionprob_ = self.get_emission()
    return model


def test():
  hmm = HMM_graph()
  hmm.add_loop(4)
  print hmm.G.node
  print hmm.G.edges()
  print "transition matrix..."
  print hmm.get_transition()
  print hmm.uniform_distribution()
  print hmm.G.node[1]['e']

if __name__ == '__main__':
  test()

#hmm.add_jump(100)
#hmm.add_linear(6)
#hmm.add_jump(100)
#m = hmm.get_transition()
