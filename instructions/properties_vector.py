import igraph as ig
import numpy as np
from argparse import ArgumentParser

## read argument (algorithm label)
parser = ArgumentParser()
parser.add_argument("-a", "--algo", dest="algo", help="algorithm used", required=True, default="algo")
args = parser.parse_args()
algo = args.algo

## make the vertex and communities 0-based
g = ig.Graph.Read_Ncol('./network.dat', directed=False)
g.simplify() ## should not be required

L = []
L.append(g.transitivity_avglocal_undirected())
L.append(g.transitivity_undirected())
ec = g.evcent()
L.extend([np.mean(ec),np.std(ec)])
sp = []

for i in range(100):
    v = np.random.choice(g.vcount(),size=2,replace=False)
    sp.append(g.shortest_paths_dijkstra(v[0],v[1])[0][0])

L.extend([np.mean(sp),np.std(sp)])
gc = g.coreness()
L.extend([np.min(gc),np.mean(gc),np.std(gc),np.max(gc)])
X = [algo]+[round(x,6) for x in L]
print(', '.join(map(str,X)))
