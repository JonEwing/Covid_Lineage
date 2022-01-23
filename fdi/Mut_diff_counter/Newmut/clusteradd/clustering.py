import pandas as pd

mf = pd.read_csv("mutations.csv", index_col = 0)
gf = pd.read_csv("protein_diff.csv", index_col = 0)
cf = pd.read_csv("communities.csv", index_col = 0)

cf = cf.sort_index()
gf = gf.sort_index()

x = gf.to_string(header=False, index=False,index_names=False).split('\n')
vals = [' - '.join(ele.split()) for ele in x]

totalappend = []

lista = cf.index.to_list()
listb = gf.index.to_list()
listc = cf["Community"].tolist()
for x in range(len(listb)):
    testbool = False
    for y in range(len(lista)):
        if listb[x] == lista[y]:
            testbool = True
            totalappend.append(listc[y])
    if testbool == False:
        gf = gf.drop(listb[x], axis = 0)

gf.insert(0, "Cluster", totalappend)

gf.to_csv('list.csv')

#df.to_csv('list.csv', index=False)
