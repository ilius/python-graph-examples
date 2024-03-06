#!/usr/bin/env python3

import igraph as ig

g = ig.Graph(1)  # one vertix
g.vs[0]["name"] = "A"

g.add_vertices(2)
g.vs[1]["name"] = "B"
g.vs[2]["name"] = "C"


g.add_vertices(4)

g.vs[g.vcount() - 1]["name"] = "Z"

print(
	g.vs["name"],
)  # read only, any changes (other than assigning the whole list) is ignored


g.add_edges(
	[
		(0, 1),
		(0, 2),
		(1, 2),
		(3, 4),
		(5, 6),
	],
)
print(f"{g.vcount()} vertices")
print(f"{g.ecount()} edges")

cc = g.biconnected_components()
print(f"{len(cc)} connected components")
for i in range(len(cc)):
	gs = cc.subgraph(i)  # a Graph
	print(gs.get_adjacency())
	print(gs.vs["name"])
