import networkx as nx
import matplotlib.pyplot as plt
import numpy as np

class Graph_ConnectedComponents_Visualizer:
    
	def __init__(self, nodes, edges):
		self._nodes = nodes
		self._edges = edges
		self._graph = self._construct_graph()

	def _construct_graph(self):
		graph = nx.Graph()
		graph.add_nodes_from(self._nodes)
		graph.add_edges_from(self._edges)
		return graph
	
	def visualize_connected_components(self):

		connected_components = list(nx.connected_components(self._graph))

		colours = plt.cm.rainbow(np.linspace(0, 1, len(connected_components)))

		node_colours = {}
		for indx, component in enumerate(connected_components):
			for node in component:
				node_colours[node] = colours[indx]

		list_node_colours = [node_colours[node] for node in list(self._graph.nodes())]
		pos = nx.spring_layout(self._graph, seed=42)
		nx.draw(self._graph, pos, with_labels = True, edge_color = "gray", node_color = list_node_colours)
		plt.show()

nodes = [1, 2, 3, 4]
edges = [(1, 2), (2, 1), (4, 3), (3, 4)]
gccmp_visualizer = Graph_ConnectedComponents_Visualizer(nodes, edges)
gccmp_visualizer.visualize_connected_components()
