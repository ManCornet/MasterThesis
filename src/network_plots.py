# --------- Code to print the network using python to have better plots --
import networkx as nx
import json
import matplotlib.pyplot as plt

f = open('network_topology.json')
data = json.load(f)

# Creato
G = nx.Graph()
for n in data['nodes']:
    G.add_node(n['id'])


subax1 = plt.subplot(121)
nx.draw(G, with_labels=True, font_weight='bold')
plt.show()