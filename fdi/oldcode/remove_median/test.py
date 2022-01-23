import numpy as np
from networkx.generators.community import gaussian_random_partition_graph
from numpy.core.fromnumeric import size
import numpy as np
import math
import pandas as pd
import networkx as nx
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from networkx.drawing.nx_agraph import graphviz_layout
import community as community_louvain
from numpy import genfromtxt

file1 = open('adj.adjlist', 'r')

store = []
for line in file1:
    segline = line.split(" ")
    striped = []
    for x in segline:
        striped.append(x.strip())
    store.append(striped)

new_store = store

holdarray = []
seconedstore = []

for counter in range(5):
    if seconedstore:
        new_store = seconedstore
    for line in new_store:
        for node in line:
            match = False
            for replace in store:
                if line != replace and node in replace:
                    for fuck in replace:
                        if fuck != node:
                            seconedstore.append(fuck)
            if match == False:
                seconedstore.append(node)

        holdarray.append(seconedstore)

lastlist = []
for x in holdarray:
    smalllist = []
    for y in x:
        if y not in smalllist and "mv" not in y:
            smalllist.append(y)
    lastlist.append(smalllist)

file1 = open('array.txt', 'w')
network = []
for x in lastlist:
    for y in x:
        if x[0] != y:
            network.append([x[0],y])
            file1.write(x[0] + " " + y + "\n")


df = pd.read_csv("../data/FinalUSMeta.csv")

names = df["Accession ID"].to_list()
symptoms = df["Symptoms"].to_list()

G = nx.Graph()
for x in network:
    G.add_nodes_from([x[0],x[1]])
    G.add_edge(x[0], x[1])

posa = graphviz_layout(G)
plt.figure(figsize=(64,36))
nx.draw(G, pos = posa, font_weight='bold')
plt.savefig("location.png")