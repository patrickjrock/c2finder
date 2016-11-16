import networkx as nx
import hmmlearn

class HMM_graph:
  """ Class that represents hmm as a graph
      outputs matrix representations of graph for use in hmm libraries
  """
  def __init__(self):
    self.G = nx.Graph()


  def get_max(self):
    return len(self.G.node)

  def add_linear(self, n):
    node_range = range(self.get_max(), self.get_max()+n)
    for i in node_range:
      self.G.add_node(i)
    for i in node_range[:-1]:
      self.G.add_edge(i,i+1)

  def add_jump(self, n):
    pass

  def add_loop(self, n):
    pass

  def get_transition():
    pass
  

def test():
  hmm = HMM_graph()
  hmm.add_linear(10)
  print hmm.G.node
  print hmm.G.edges()

test()

#hmm.add_jump(100)
#hmm.add_linear(6)
#hmm.add_jump(100)
#m = hmm.get_transition()
