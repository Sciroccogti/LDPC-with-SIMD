'''
@file tanner.py
@author Yifan Zhang (scirocco_gti@yeah.net)
@brief 
@date 2020-10-29 21:29:24
@modified: 2020-10-31 12:38:14
'''

# Draw tanner from a H matrix

import argparse

import matplotlib.pyplot as plt
import networkx as nx

from alist import alist2sparse2

parser = argparse.ArgumentParser()
parser.add_argument("-H", "--dec-h-path", help="path to H matrix in .alist")
args = parser.parse_args()

H = alist2sparse2(args.dec_h_path)
# H = alist2sparse2("example/CCSDS_ldpc_n128_k64.alist")
N = H.shape[1]
M = H.shape[0]

vNodes = []  # variable nodes
cNodes = []  # parity nodes
colorMap = []
pos = {}
edges = []
for i in range(N):
    vNodes.append("V%d" % i)
    colorMap.append("red")
    pos["V%d" % i] = [i, 5]
for i in range(M):
    cNodes.append("C%d" % i)
    colorMap.append("blue")
    pos["C%d" % i] = [i / M * N, 0]
    for j in range(N):
        if H[i][j] == 1:
            edges.append(("V%d" % j, "C%d" % i))

G = nx.Graph()
G.add_nodes_from(vNodes + cNodes)
G.add_edges_from(edges)
nx.draw(G, pos=pos, with_labels=True, node_color=colorMap, font_color="white")
plt.show()
