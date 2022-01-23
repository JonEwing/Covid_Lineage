import pandas as pd
lista = ["AlignedFilteredUS.nex"]
import math

aminoacidsnum = [0, 4405, 11501, 12774, 13049, 13124, 13346, 13407, 13528, 13571, 13692, 14111, 14149] 
aminoacidnames = ["ORF1a/ORF1ab", "ORF1ab", "ORF2", "ORF3a", "ORF4", "ORF5", "ORF6", "ORF7a", "ORF7b", "ORF8", "ORF9", "ORF10", ]
nspnums = [0, 179, 817, 2762, 3262, 3568, 3858, 3941, 4139, 4252, 4391, 4404, 5323, 5924, 6451, 6797, 7095]
nspnames = ["Nsp1", "Nsp2", "Nsp3", "Nsp4", "Nsp5", "Nsp6", "Nsp7", "Nsp8", "Nsp9", "Nsp10", "Nsp11", "Nsp12", "Nsp13", "Nsp14", "Nsp15", "Nsp16", ]

#DNA to amino Acids
def translate(seq):
    table = {
        'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M',
        'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
        'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K',
        'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',                 
        'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L',
        'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
        'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q',
        'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',
        'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V',
        'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
        'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E',
        'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
        'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S',
        'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L',
        'TAC':'Y', 'TAT':'Y', 'TAA':'_', 'TAG':'_',
        'TGC':'C', 'TGT':'C', 'TGA':'_', 'TGG':'W',
    }
    protein = ""
    for x in range(len(seq)):
        if x%3 == 0:
            try:
                codon = seq[x].upper() + seq[x+1].upper() + seq[x+2].upper()
                if codon in table:
                    if seq[x].islower() or seq[x + 1].islower() or seq[x + 2].islower():
                        protein += table[codon].lower()
                    else:
                        protein += table[codon]
                else:
                    protein += "-"
            except:
                protein += "-"
    return protein

#Read data
df = pd.read_csv("SARS-CoV-2 Reference Sequence_All Regions.csv")
genename = df["name"]
genlocs = df["start"]
genloce = df["end"]

for x in lista:
    fp = open(x, 'r')

    holdarray = fp.readlines() 

    holdstr = ""
    holdarr = []
    for x in holdarray:
        holdarr.append([x.split(' ')[0],x.split(' ')[1]])

#Find Mutation Diffrences
totalarray = []
for x in range(len(holdarr)):
    if x != 0:
        partalarray = []
        for y in range(len(holdarr[x][1])):
            if holdarr[0][1][y] != holdarr[x][1][y]:
                tempstr = ""
                for z in range(len(genlocs)):
                    if y + 77 >= genlocs[z] and y + 78 <= genloce[z]:
                        tempstr += genename[z] + ", "       
                if tempstr:
                    partalarray.append(tempstr + "at " + str(y + 78) + ": " + holdarr[0][1][y] + "->" + holdarr[x][1][y])
                    temp = list(holdarr[x][1])
                    temp[y] = temp[y].lower()
                    holdarr[x][1] = "".join(temp)
                else:
                    partalarray.append(str(y + 77) + ": " + holdarr[0][1][y] + "->" + holdarr[x][1][y])
                    temp = list(holdarr[x][1])
                    temp[y] = temp[y].lower()
                    holdarr[x][1] = "".join(temp)
        
        partalarray = [holdarr[x][0]] + partalarray
        totalarray.append(partalarray)

df = pd.DataFrame(list(totalarray))
df.to_csv('mutations.csv', index=False)

#Remove bases until ORF1
replaceh = []
for x in holdarr:
    replaceh.append([x[0],x[1][x[1].find('ATGGAGAG'):]])

#Convert DNA to amino acids\

aaaaaaaaaa = []
for x in replaceh:
    totalp = ""
    for y in range(len(genename)):
        x[1] = x[1].replace("\n", "") 
        x[1] = x[1].replace("\r", "")
        if "ORF" in genename[y]:
            if "ab" in genename[y]:
                pro = translate(x[1][genlocs[y] - 266 : 13469 - 266])
                totalp += pro
                pro = translate(x[1][13468 - 266 : genloce[y] - 266])
                totalp += pro

            else:
                pro = translate(x[1][genlocs[y] - 266 : genloce[y] - 266])
                totalp += pro
    aaaaaaaaaa.append([x[0],totalp])

#Read protein data
df = pd.DataFrame(list(aaaaaaaaaa))
df.to_csv('protein.csv', index=False)

proarray = []
for x in replaceh:
    totalpro = []
    for y in range(len(genename)):
        x[1] = x[1].replace("\n", "") 
        x[1] = x[1].replace("\r", "")
        if "ORF" in genename[y]:
            if "ab" in genename[y]:
                tempsty = ""
                pro = translate(x[1][genlocs[y] - 266 : 13469 - 266])
                tempsty += pro
                pro = translate(x[1][13468 - 266 : genloce[y] - 266])
                tempsty += pro
                totalpro.append(tempsty)
            else:
                pro = translate(x[1][genlocs[y] - 266 : genloce[y] - 266])
                totalpro.append(pro)
    proarray.append([x[0],totalpro])

total = []
for a in range(len(proarray)):
    if a != 0:
        storearray = []
        for b in range(len(proarray[a][1])):
            for c in range(len(proarray[a][1][b])):
                if proarray[a][1][b][c].islower():
                    if "ORF1ab" == aminoacidnames[b] and c < 4405:
                        continue
                    elif "ORF1" in aminoacidnames[b]:
                        for x in range(len(nspnums)):
                            try:
                                if c+1 > nspnums[x] and c+1 < nspnums[x+1]:
                                    storearray.append(nspnames[x] + ": " + aminoacidnames[b] + " at " +str(c + 1) + ": " + proarray[0][1][b][c] + "->" + proarray[a][1][b][c].upper())
                            except:
                                storearray.append("IBroke")
                    else:
                        storearray.append(aminoacidnames[b] + " at " +str(c + 1) + ": " + proarray[0][1][b][c] + "->" + proarray[a][1][b][c].upper())
        storearray = [proarray[a][0]] + storearray
        total.append(storearray)

#Output
df = pd.DataFrame(list(total))
df.to_csv('protein_diff.csv', index=False)