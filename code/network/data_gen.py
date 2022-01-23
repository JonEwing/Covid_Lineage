import pandas as pd
from Bio import pairwise2
from Bio.pairwise2 import format_alignment
import matplotlib.pyplot as plt
import datetime
import matplotlib.pyplot as plt
import numpy
import networkx as nx

f = open("jonfe.out", "r")
array = f.readlines()

hold = False
counter = 0
startcounter = 0
while hold != True:
    if array[counter] == '**************************************************\n':
        hold = True
    if array[counter] == 'Final links : \n':
        startcounter = counter + 2
    counter += 1
counter -= 1
array = array[startcounter:counter]
f.close()

adjacentcy = []
for x in array:
    string = x.split(" ")
    removedempty = []
    for y in string:
        if y != '' and y != "mv":
            removedempty.append(y)
            
    adjacentcy.append(eval(str(removedempty[4])+ ", "+ str(removedempty[7][0:len(removedempty[7])-1])))

G = nx.Graph()
G.add_edges_from(adjacentcy)

nx.draw(G, cmap = plt.get_cmap('jet'))
plt.savefig("test.png")

f = open("1000out.fasta", "r")
array = f.readlines()
keeps = []
for x in array:
    if x[0] == '>':
        keeps.append(x[1:len(x)-1])

removed = []
for x in adjacentcy:
    if x[0] <= len(keeps) and x[1] <= len(keeps):
        removed.append(x)

counter = []
for x in range(len(keeps)):
    counter.append([" ",0,0,0])

for x in removed:
    counter[x[0] - 1][2] += 1
    counter[x[1] - 1][3] += 1

for x in range(len(counter)):
    counter[x-1][0] = keeps[x-1]
    counter[x-1][1] = counter[x-1][2] + counter[x-1][3]

df = pd.DataFrame(counter, columns=["Sample","Total Connections", "Connections to node", "Connections from node"])
df.to_csv('test.csv', index=False)