class Graph:

	def __init__(self, nodes, edges):
		self._nodes = nodes
		self._edges = edges

		# construct adjacency list (graph)
		self._graph = self._construct_graph()

	def _construct_graph(self):
		
		# adjacency list to store graph for effient neighbour look ups.
		
		graph = {node:[] for node in self._nodes}

		# insert edges
		for (node, adj_node) in self._edges:
			graph[node].append(adj_node)
	
		return graph
	
	def adjacent_nodes(self, node):

		# return the adjacent nodes of a given node.
		if node in self._graph:
			return self._graph[node]
		
		raise KeyError("Node not found in the graph")
	
	@property
	def node_count(self):
		
		# number of nodes in the graph.
		return len(self._nodes)

	@property
	def edge_count(self):
		
		# number of edges in the graph.
		return len(self._edges)

	def get_connected_components(self):
		
		# list to collect the connected components.
		connected_components = []

		if self._graph:

			# set to keep track of remaining nodes after fining connnected comps one after the other.
			remaining_nodes = set(self._nodes)
			
			# loop until all nodes are exausted.
			while remaining_nodes:
				
				# pick one of the remaining nodes as the start node for the next connected component.
				start_node = remaining_nodes.pop()

				# begin with start node and find the connected component.
				conn_comp = self._get_connected_component(start_node)

				# collect the connected component
				connected_components.append(conn_comp)

				# delete the nodes in the connected component from the set of remaining nodes.
				remaining_nodes -= conn_comp

		return connected_components

	def _get_connected_component(self, start_node):
		
		# set to keep track of already visited nodes.
		visited = set()
		
		# add the start node to stack.
		stack = [start_node]

		# mark the start node visited.
		visited.add(start_node) 
		
		# while stack is not empty, continue
		while stack:
			
			# dfs: access the recently added node.
			top_node = stack[-1]
			
			# add to stack and mark visited the unvisited adjacent node of the recently added node.
			unvisted_node_found = False
			adj_nodes = self.adjacent_nodes(top_node)
			for adj_node in adj_nodes:

				# if adj node not visited yet, add it to stack and mark it visited.
				if adj_node not in visited:

					stack.append(adj_node)
					visited.add(adj_node)
					unvisted_node_found = True
					break
			
			# if all the adjacent nodes were already visited, then pop the stack.
			if not unvisted_node_found:
				stack.pop()

		# return the nodes that are connected (visited list)
		return visited
	

# Test case 1: Simple connected graph
nodes = [1, 2, 3, 4]
edges = [(1, 2), (2, 1), (2, 3), (3, 2), (1, 3), (3, 1)]
graph = Graph(nodes, edges)
print(graph.get_connected_components())  # Expected: [{1, 2, 3}, {4}]

# Test case 2: Empty graph
nodes = []
edges = []
graph = Graph(nodes, edges)
print(graph.get_connected_components())  # Expected: []

# Test case 3: Disconnected graph
nodes = [1, 2, 3, 4]
edges = [(1, 2), (2, 1),(3,4),(4,3)]
graph = Graph(nodes, edges)
print(graph.get_connected_components())  # Expected: [{1, 2}, {3}, {4}]

# Test case 4: Graph with a single node
nodes = [1]
edges = []
graph = Graph(nodes, edges)
print(graph.get_connected_components())  # Expected: [{1}]

