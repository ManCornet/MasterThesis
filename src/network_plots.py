#-----------------------------------------------------------------------------
#
#                           - TFE : Bilevel DNEP - 
#                             University of Liege
#
#-----------------------------------------------------------------------------
# Created By  : Manon Cornet
# Created Date: Saturday May 20 2023
#
# network_plots:
#   Python file that draws beautiful plots for the simulations of the optimization 
#   program.
#
# =============================================================================
#                                   Imports
# =============================================================================
import json
from pyvis.network import Network

# =============================================================================
#                                Core of the text
# =============================================================================

f = open('network_topology.json')
topology_data = json.load(f)
f.close()

f = open('network_data.json')
network_data = json.load(f)
f.close()

# print(network_data)

# g = gt.Graph(g=len(topology_data['nodes']), directed=False)

# pos = dict()
# for n in topology_data['nodes']:
#     node_id = int(n['id'])
#     coords  = n['coord']

#     v = g.add_vertex()
#     gt.vertex_id = node_id - 1
#     coords = gt.new_vertex_property("tuple")
#     coords[v] = (float(coords['x']),  float(coords['y']))


# all_edges = []
# for e in topology_data['edges']:
#     from_node = e['from_node']
#     to_node   = e['to_node']
#     g.add_edge(int(from_node['id']), int(to_node['id']))

# gt.graph_draw(g, pos=pos, output="two-nodes.pdf")

net = Network()
net.add_node(1, label='Alex')
net.add_node(2, label='Cathy')

net.show('nodes.html')

# print(topology_data['nodes'])
# 
#gt.graph_draw(g, vertex_text=g.vertex_index, output="two-nodes.pdf")
#gt.graphviz_draw(g, output="lattice-planar.pdf")