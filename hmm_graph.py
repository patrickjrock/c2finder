import networkx as nx
from hmmlearn.hmm import MultinomialHMM

class HMM_graph:
  """ Class that represents hmm as a graph
      outputs matrix representations of graph for use in hmm libraries
  """
  def __init__(self):
    self.G = nx.Graph()

  def norm(self, raw):
    return [float(i)/sum(raw) for i in raw]

  def uniform_distribution(self):
    v = [1 for i in range(1,21)]
    return self.norm(v)

  def get_max(self):
    return len(self.G.node)

  def add_linear(self, n):
    node_range = range(self.get_max(), self.get_max()+n)
    for i in node_range:
      self.G.add_node(i, e=self.uniform_distribution()) # abstract out the distribution later
    for i in node_range[:-1]:
      self.G.add_edge(i,i+1)

  def add_jump(self, n):
    pass #Adding something

  def add_loop(self, n):
    pass

  def get_transition(self):
    """ converts the adjacency matrix to the transition matrix by normalizing the rows """
    A = nx.adjacency_matrix(self.G).todense().getA()
    return [self.norm(row) for row in A]

  def get_emission(self):
    """ returns the emission matrix for the hmm. If no distribution is present in a node assume uniform """
    # for now just assume everything is uniform
    return [self.uniform_distribution() for i in range(0, self.get_max())]

  def get_start(self):
    """ returns the start probability matrix, just points to the first node"""
    pass

  def get_model(self):
    """ returns a multinomial hmm"""
    model = MultinomialHMM(n_components=self.get_max())
    model.transmat_ = self.get_transition()
    model.emissionprob_ = self.get_emission()
    return model


def test():
  hmm = HMM_graph()
  hmm.add_linear(10)
  print hmm.G.node
  print hmm.G.edges()
  print hmm.get_transition()
  print hmm.uniform_distribution()
  print hmm.G.node[1]['e']

if __name__ == '__main__':
  test()

#hmm.add_jump(100)
#hmm.add_linear(6)
#hmm.add_jump(100)
#m = hmm.get_transition()
