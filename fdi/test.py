import pandas as pd
import math

df = pd.read_csv("data/FinalUSMeta.csv")
filepath = "data/AlignedFilteredUS.fdi"

############################################################################## Symptoms
symptoms = df["Symptoms"].to_list()
fout = open("results/symptoms.fdi", "w")

with open(filepath) as fp:
    line = fp.readline()
    ref = True
    listcounter = 0
    while line:
        linesplit = line.split(';')
        holdline = line

        print(linesplit[9])
        if linesplit[0] == "TAXON_NAME" and linesplit[1][0] != "m":
            if ref == True:
                linesplit[30] = 16744703
                ref = False
            else:
                if symptoms[listcounter] == "Asymptomatic or Mild Symptoms":
                    linesplit[30] = 16776960
                elif symptoms[listcounter] == "Symptomatic Moderate to Severe Symptoms":
                    linesplit[30] = 255
                else:
                    linesplit[30] = 65280

                listcounter += 1 

            newstring = ""
            for x in linesplit:
                if x != '':
                    newstring = newstring + ";" + str(x)
                else:
                    newstring = newstring + ";"
            holdline = newstring[1:]

        fout.write(holdline)
        line = fp.readline()
fout.close()
fp.close()

############################################################################## Outcome

outcome = df["Outcome"].to_list()
fout = open("results/outcomes.fdi", "w")

with open(filepath) as fp:
    line = fp.readline()
    ref = True
    listcounter = 0
    while line:
        linesplit = line.split(';')
        holdline = line
        if linesplit[0] == "TAXON_NAME" and linesplit[1][0] != "m":
            if ref == True:
                linesplit[30] = 16744703
                ref = False
            else:
                if outcome[listcounter] == "Alive":
                    linesplit[30] = 255
                elif outcome[listcounter] == "Deceased":
                    linesplit[30] = 0
                else:
                    linesplit[30] = 65280

                listcounter += 1 

            newstring = ""
            for x in linesplit:
                if x != '':
                    newstring = newstring + ";" + str(x)
                else:
                    newstring = newstring + ";"
            holdline = newstring[1:]

        fout.write(holdline)
        line = fp.readline()
fout.close()
fp.close()

############################################################################## State

location = df["Location"].to_list()

states = []
for x in location:
    temp = x.split(' / ')[2]
    if temp == "South Carolina":
        temp = "California"
    states.append(temp)

topstates = []
for x in states:
    match = False
    for y in topstates:
        if x == y[0]:
            y[1] += 1
            match = True
    if match == False:
        topstates.append([x,1])

top1 = ["",0]
top2 = ["",0]
top3 = ["",0]
for x in topstates:
    if x[1] > top1[1]:
        top3 = top2
        top2 = top1
        top1 = x
    elif x[1] > top2[1]:
        top3 = top2
        top2 = x
    elif x[1] > top3[1]:
        top3 = x

print(top1)
print(top2)
print(top3)
fout = open("results/states.fdi", "w")

with open(filepath) as fp:
    line = fp.readline()
    ref = True
    listcounter = 0
    while line:
        linesplit = line.split(';')
        holdline = line
        if linesplit[0] == "TAXON_NAME" and linesplit[1][0] != "m":
            if ref == True:
                linesplit[30] = 16744703
                ref = False
            else:
                if states[listcounter] == top1[0]:
                    linesplit[30] = 65535
                elif states[listcounter] == top2[0]:
                    linesplit[30] = 12632256
                elif states[listcounter] == top3[0]:
                    linesplit[30] = 4194432
                else:
                    linesplit[30] = 65280

                listcounter += 1 

            newstring = ""
            for x in linesplit:
                if x != '':
                    newstring = newstring + ";" + str(x)
                else:
                    newstring = newstring + ";"
            holdline = newstring[1:]

        fout.write(holdline)
        line = fp.readline()
fout.close()
fp.close()

############################################################################## Time

datedash = df["Collection date"].to_list()

dates = []
for x in datedash:
    temps = x.split('-')
    m = int(temps[1])
    d = int(temps[2][0]) 

    m = (12 - m) * 10
    d = (3 - d)

    total = 255 - (m + d)
    dates.append(total)

fout = open("results/dates.fdi", "w")

with open(filepath) as fp:
    line = fp.readline()
    ref = True
    listcounter = 0
    while line:
        linesplit = line.split(';')
        holdline = line
        if linesplit[0] == "TAXON_NAME" and linesplit[1][0] != "m":
            if ref == True:
                linesplit[30] = 65535
                ref = False
            else:
                linesplit[30] = dates[listcounter]

                listcounter += 1 

            newstring = ""
            for x in linesplit:
                if x != '':
                    newstring = newstring + ";" + str(x)
                else:
                    newstring = newstring + ";"
            holdline = newstring[1:]

        fout.write(holdline)
        line = fp.readline()
fout.close()
fp.close()

############################################################################## Time and symp

datedash = df["Collection date"].to_list()
symptoms = df["Symptoms"].to_list()

dates = []
for x in datedash:
    temps = x.split('-')
    m = int(temps[1])
    d = int(temps[2][0]) 

    m = (12 - m) * 10
    d = (3 - d)

    total = 255 - (m + d)
    dates.append(total)

fout = open("results/dates_symptom.fdi", "w")

with open(filepath) as fp:
    line = fp.readline()
    ref = True
    listcounter = 0
    while line:
        linesplit = line.split(';')
        holdline = line
        if linesplit[0] == "TAXON_NAME" and linesplit[1][0] != "m":
            if ref == True:
                linesplit[30] = 65535
                ref = False
            else:
                if symptoms[listcounter] == "Asymptomatic or Mild Symptoms":
                    linesplit[30] = dates[listcounter] * 16 * 16 * 16 * 16
                elif symptoms[listcounter] == "Symptomatic Moderate to Severe Symptoms":
                    linesplit[30] = dates[listcounter]
                else:
                    linesplit[30] = dates[listcounter] * 16 * 16
                    
                listcounter += 1 
                
            newstring = ""
            for x in linesplit:
                if x != '':
                    newstring = newstring + ";" + str(x)
                else:
                    newstring = newstring + ";"
            holdline = newstring[1:]

        fout.write(holdline)
        line = fp.readline()
fout.close()
fp.close()