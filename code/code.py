import pandas as pd
from Bio import pairwise2
from Bio.pairwise2 import format_alignment
import matplotlib.pyplot as plt
import datetime
import matplotlib.pyplot as plt
import numpy

#Change of Ref v other samples
#pull out just us samples month and week

df = pd.read_csv('AllMetadata.csv',index_col=0)

names = df["Accession ID"].values

dates = df["Collection date"].values

locations = df["Location"].values

f = open("FinalAllSequences.fasta", "r") 

metadata = []
data = []
counter = 0
for x in f:
    if counter == 0:
        metadata.append([x, "0", "0", "NAN"])
    elif counter % 2 == 0:
        holdx = x[1:].strip()
        for y in range(len(names)):
            if names[y] == holdx:
                try:
                    today = dates[y].split('-')
                    a_date = datetime.date(int(today[0]), int(today[1]), int(today[2]))
                    week_number = a_date.isocalendar()[1]
                    month = today[1]
                except:
                    month = 0
                    week_number = 0
                try:
                    country = locations[y].split("/")[1]
                except:
                    country = "NAN"
                metadata.append([holdx, month, week_number, country])
                break
    elif counter % 2 == 1:
        data.append(x)
    counter += 1
    if counter % 20000 == 0:
        print("Progress: "+ str(counter) + "/200000")

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
        if char == 'A':
            countA += 1
            tail += 1
        if char == 'T':
            countT += 1
            tail = 0
        if char == 'C':
            countC += 1
            tail = 0
        if char == 'G':
            countG += 1
            tail = 0
        if char == 'N':
            countN += 1
            tail = 0
        if char == '-':
            countdash += 1
            tail = 0
    if counter == 0:
        A = countA - tail
        T = countT
        C = countC
        G = countG
        print('Reference A, T, C, G')
        print(str(A), str(T), str(C), str(G))

    countA -= tail
    if countN < 150 and countA + countT + countC + countG >= 29750:
        totalarray.append([metadata[counter][0], metadata[counter][1], metadata[counter][2], metadata[counter][3], countA - A, countT - T, countC - C, countG - G, countN, tail])
    else: 
        removed += 1

    counter += 1
    if counter % 20000 == 0:
        print("Progress: "+ str(counter) + "/100000")
    ntdsgraph.append(len(x))
    ngraph.append(countN)

print("Removed:",removed)
df = pd.DataFrame(totalarray,columns=['Sequence ID',"month", "week", "location", 'A', 'T', 'C', 'G', 'N', 'PolyA tail'])
df.to_csv("samples/counts.csv", index=False)

################################################################## MONTH
dfvals = ['A', 'T', 'C', 'G', 'N']
allruns = []
for counter in dfvals:
    total = []
    monthnum = df["month"].values
    vals = df[counter].values
    for num in range(10):
        levelcount = []
        for x in range(len(monthnum)):
            if num == 0 and int(monthnum[x]) == 12:
                levelcount.append(vals[x])
            elif int(monthnum[x]) == num:
                levelcount.append(vals[x])
        total.append(levelcount)

    averages = []
    for x in total:
        try:
            averages.append(sum(x)/len(x))
        except:
            averages.append(0)
    
    names = ["Dec", "Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep"]
    plt.figure(figsize=(20,10))
    plt.title(str(counter) + " Difference form reference")
    plt.plot(names, averages)
    plt.xlabel("Month")
    plt.ylabel("Average " + str(counter) + "Difference")
    plt.savefig("samples/months/"+str(counter) + "Month.png")
    plt.clf()

    allruns.append(averages)
allruns.append(names)

tf = pd.DataFrame(numpy.transpose(allruns),columns=['A', 'T', 'G', 'C', 'N', 'Name'])
tf.to_csv("samples/months/avgmonths.csv", index=False)

################################################################### Week
allruns = []
for counter in dfvals:
    total = []
    vals = df[counter].values
    weeknum = df["week"].values
    for num in range(38):
        levelcount = []
        for x in range(len(weeknum)):
            if num == 0 and int(weeknum[x]) == 52:
                levelcount.append(vals[x])
            elif int(weeknum[x]) == num:
                levelcount.append(vals[x])
        total.append(levelcount)

    averages = []
    for x in total:
        try:
            averages.append(sum(x)/len(x))
        except:
            averages.append(0)

    names = ["52", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", '21', "22", '23', '24', "25", '26', "27", "28", "29", "30", "31", "32", "33", "34", "35", "36", "37"]
    plt.figure(figsize=(20,10))
    plt.title(str(counter) + " Difference form reference")
    plt.plot(names, averages)
    plt.xlabel("Week")
    plt.ylabel("Average " + str(counter) + "Difference")
    plt.savefig("samples/week/" + str(counter) + "week.png")
    plt.clf()

    allruns.append(averages)
allruns.append(names)
tf = pd.DataFrame(numpy.transpose(allruns),columns=['A', 'T', 'G', 'C', 'N', 'Names'])
tf.to_csv("samples/week/avgweeks.csv", index=False)

################################################################ Location

uniqloc = []
locs = df["location"].values

for x in locs:
    Hit = False
    x = x.strip()
    if x != 'NAN':
        for y in uniqloc:
            y = y.strip()
            if x == y:
                Hit = True
        if Hit == False:
            uniqloc.append(x)

allruns = []
for counter in dfvals:
    total = []
    monthnum = df["location"].values
    vals = df[counter].values
    for num in uniqloc:
        levelcount = []
        for x in range(len(monthnum)):
            temp = monthnum[x].strip()
            if temp == num:
                levelcount.append(vals[x])
        total.append(levelcount)

    averages = []
    for x in total:
        try:
            averages.append(sum(x)/len(x))
        except:
            averages.append(0)

    names = uniqloc
    plt.figure(figsize=(20,10))
    plt.plot(names, averages)
    plt.title(str(counter) + " Difference form reference")
    plt.xlabel("Country")
    plt.xticks(rotation=90)
    plt.ylabel("Average " + str(counter) + "Difference")
    plt.savefig("samples/location/" + str(counter) + "loc.png")
    plt.clf()

    allruns.append(averages)
allruns.append(uniqloc)
tf = pd.DataFrame(numpy.transpose(allruns),columns=['A', 'T', 'G', 'C', 'N', 'Names'])
tf.to_csv("samples/location/avglocs.csv", index=False)

####################################################################################### USA Month
allruns = []
for counter in dfvals:
    total = []
    vals = df[counter].values
    weeknum = df["month"].values
    locs = df["location"].values
    for num in range(10):
        levelcount = []
        for x in range(len(weeknum)):
            temp = locs[x].strip()
            if temp == "USA":
                if num == 0 and int(weeknum[x]) == 12:
                    levelcount.append(vals[x])
                elif int(weeknum[x]) == num:
                    levelcount.append(vals[x])
        total.append(levelcount)

    averages = []
    for x in total:
        try:
            averages.append(sum(x)/len(x))
        except:
            averages.append(0)

    names = ["Dec", "Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep"]
    plt.figure(figsize=(20,10))
    plt.title(str(counter) + " Difference form reference")
    plt.plot(names, averages)
    plt.xlabel("Week")
    plt.ylabel("Average " + str(counter) + "Difference")
    plt.savefig("samples/usamonths/" + str(counter) + "month.png")
    plt.clf()

    allruns.append(averages)
allruns.append(names)
tf = pd.DataFrame(numpy.transpose(allruns),columns=['A', 'T', 'G', 'C', 'N', 'Names'])
tf.to_csv("samples/usamonths/avgweeks.csv", index=False)

###################################################################USA Week
allruns = []
for counter in dfvals:
    total = []
    vals = df[counter].values
    weeknum = df["week"].values
    locs = df["location"].values
    for num in range(38):
        levelcount = []
        for x in range(len(weeknum)):
            temp = locs[x].strip()
            if temp == "USA":
                if num == 0 and int(weeknum[x]) == 52:
                    levelcount.append(vals[x])
                elif int(weeknum[x]) == num:
                    levelcount.append(vals[x])
        total.append(levelcount)

    averages = []
    for x in total:
        try:
            averages.append(sum(x)/len(x))
        except:
            averages.append(0)

    names = ["52", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", '21', "22", '23', '24', "25", '26', "27", "28", "29", "30", "31", "32", "33", "34", "35", "36", "37"]
    plt.figure(figsize=(20,10))
    plt.plot(names, averages)
    plt.title(str(counter) + " Difference form reference")
    plt.xlabel("Week")
    plt.ylabel("Average " + str(counter) + "Difference")
    plt.savefig("samples/usaweek/" + str(counter) + "week.png")
    plt.clf()

    allruns.append(averages)
allruns.append(names)
tf = pd.DataFrame(numpy.transpose(allruns),columns=['A', 'T', 'G', 'C', 'N', 'Names'])
tf.to_csv("samples/usaweek/avgweeks.csv", index=False)