import pandas as pd

filepath = "directed.adjlist"

nodearray = ["REFERE", "475633", "475585", "677673", "475590", "475674", "ENGLAN", "436946", "463286", "463295", "436944", "940958", "942009", "516771", "900696", 
"mv1", "mv58", "mv5", "mv86", "mv49", "mv74", "mv31", "mv125", "mv24", "mv92", "mv73", "mv16", "mv25", "mv41", "mv19", "mv23", "mv32", "mv13"]

finalnodelist = []
with open(filepath) as fp:
    line = fp.readline()
    while line:
        linesplit = line.split(' ')
        linesplit[1] = linesplit[1].strip()

        for x in nodearray:
            if x == linesplit[0] and 'mv' not in x:
                finalnodelist.append(linesplit[0])
            elif x == linesplit[1] and 'mv' not in x:
                finalnodelist.append(linesplit[1])
            
            elif x == linesplit[0] and 'mv' in x:
                if "mv" not in linesplit[1]:
                    finalnodelist.append(linesplit[1])
            elif x == linesplit[1] and 'mv' in x:
                if "mv" not in linesplit[0]:
                    finalnodelist.append(linesplit[0])



        line = fp.readline()


finalnodelist = list(set(finalnodelist))
percentile_list = pd.DataFrame({'Sample Name': list(finalnodelist)})
percentile_list.to_csv("final_nodes.csv", index = False)