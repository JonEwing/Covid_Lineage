import pandas as pd
from Bio import pairwise2
from Bio.pairwise2 import format_alignment
import matplotlib.pyplot as plt
import datetime


#https://stackoverflow.com/questions/2600775/how-to-get-week-number-in-python
#Month
#Day
#Week
#Country

#look at atcg precents for full group
#then Month
#then week
#then country

df = pd.read_csv('AllMetadata.csv',index_col=0)

names = df["Accession ID"].values

dates = df["Collection date"].values

locations = df["Location"].values

f = open("1stseg.fasta", "r") 

metadata = []
data = []
counter = 0
for x in f:
    if counter % 2 == 0:
        holdx = x[1:].rstrip()
        for y in range(len(names)):
            if names[y] == holdx:
                try:
                    today = dates[y].split('-')
                    a_date = datetime.date(int(today[0]), int(today[1]), int(today[2]))
                    week_number = a_date.isocalendar()[1]
                    month = today[1]
                except:
                    month = -1
                    week_number = -1
                metadata.append([holdx, month, week_number, locations[y]])
                break
    elif counter % 2 == 1:
        data.append(x)
    counter += 1
    if counter % 10000 == 0:
        print("Progress: "+ str(counter) + "/66666")

df = pd.DataFrame(metadata)
df.to_csv("test.csv", index=False)

print("Built Data")
totalarray = []
ntdsgraph = []
ngraph = []
counter = 0
removed = 0
for x in data:
    countA = 0
    countT = 0
    countC = 0
    countG = 0
    countN = 0
    countdash = 0
    tail = 0
    for char in x:
        if countN / len(x) < 0.005 and len(x) >= 29750:
            if char == 'A':
                countA += 1
                tail += 1
            elif char == 'T':
                countT += 1
                tail = 0
            elif char == 'C':
                countC += 1
                tail = 0
            elif char == 'G':
                countG += 1
                tail = 0
            elif char == 'N':
                countN += 1
                tail = 0
            elif char == '-':
                countdash += 1

    if countN / len(x) < 0.005 and countA + countT + countC + countG >= 29750:
        subtotalATCG = countA + countT + countC + countG
        subtotalbad = countN + countdash
        totalval = subtotalATCG + subtotalbad
        totalarray.append([metadata[counter][0], metadata[counter][1], metadata[counter][2], metadata[counter][3], countA, countT, countC, countG, subtotalATCG, countN, countdash, subtotalbad, tail])
    else:
        removed += 1

    if counter % 1000 == 0:
        print("Progress: "+ str(counter) + "/33333")
    ntdsgraph.append(len(x))
    ngraph.append(countN)
    counter += 1

plt.figure(figsize=(20,10))
plt.title("Nucleotides / Sample")
plt.scatter(range(len(ntdsgraph)), ntdsgraph)
plt.savefig("Nucleotides.png")
plt.clf()

plt.figure(figsize=(20,10))
plt.title("Ns / Sample")
plt.scatter(range(len(ngraph)), ngraph)
plt.savefig("Ns.png")
plt.clf()

plt.figure(figsize=(20,10))
plt.title("Ns / Nucleotides")
plt.scatter(ngraph, ntdsgraph)
plt.savefig("NS_nucleotides.png")
plt.clf()

print("Removed:",removed)
df = pd.DataFrame(totalarray,columns=['Sequence ID', 'A', 'T', 'C', 'G','Subtotal', 'N', 'No Call', 'Subtotal', 'PolyA Tail'])
df.to_csv("counts.csv", index=False)

print("Average A:", df["A"].mean())
print("Average T:", df["T"].mean())
print("Average G:", df["G"].mean())
print("Average C:", df["C"].mean())