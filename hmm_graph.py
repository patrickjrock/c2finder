import networkx as nx
import hmmlearn

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
      self.G.add_node(i)
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


def test():
  hmm = HMM_graph()
  hmm.add_linear(10)
  print hmm.G.node
  print hmm.G.edges()
  print hmm.get_transition()
  print hmm.uniform_distribution()

test()

#hmm.add_jump(100)
#hmm.add_linear(6)
#hmm.add_jump(100)
#m = hmm.get_transition()