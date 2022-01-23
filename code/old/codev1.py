import pandas as pd
from Bio import pairwise2
from Bio.pairwise2 import format_alignment
import matplotlib.pyplot as plt


#https://stackoverflow.com/questions/2600775/how-to-get-week-number-in-python
#Month
#Day
#Week
#Country

#look at atcg precents for full group
#then Month
#then week
#then country

f = open("1stseg.fasta", "r") 

metadata = []
data = []
counter = 0
for x in f:
    if counter % 2 == 0:
        metadata.append(x)
    elif counter % 2 == 1:
        data.append(x)
    counter += 1

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

    FivepUTRs = False
    if x[0:10] == "ATTAAAGGTT":
        FivepUTRs = True
    FivepUTRe = False
    if x[255:265] == "GAAAGGTAAG":
        FivepUTRe = True

    NSP1s = False
    if x[265:275] == "ATGGAGAGCC":
        NSP1s = True
    NSP1e = False
    if x[795:805] == "TAACGGAGGG":
        NSP1e = True

    NSP2s = False
    if x[805:815] == "GCATACACTC":
        NSP2s = True
    NSP2e = False
    if x[2709:2719] == "CAAAGGCGGT":
        NSP2e = True
    
    NSP3s = False
    if x[2719:2729] == "GCACCAACAA":
        NSP3s = True
    NSP3e = False
    if x[8542:8554] == "CTTAAGGGTGGT":
        NSP3e = True

    NSP4s = False
    if x[8554:8564] == "AAAATTGTTA":
        NSP4s = True
    NSP4e = False
    if x[10044:10054] == "TGTTTTGCAG":
        NSP4e = True

    NSP5s = False
    if x[10054:10064] == "AGTGGTTTTA":
        NSP5s = True
    NSP5e = False
    if x[10962:10972] == "TACTTTCCAA":
        NSP5e = True

    NSP6s = False
    if x[10972:10982] == "AGTGCAGTGA":
        NSP6s = True
    NSP6e = False
    if x[11832:11842] == "CACTGTACAG":
        NSP6e = True

    NSP7s = False
    if x[11842:11852] == "TCTAAAATGT":
        NSP7s = True
    NSP7e = False
    if x[12081:12091] == "AACCTTACAA":
        NSP7e = True

    NSP8s = False
    if x[12091:12101] == "GCTATAGCCT":
        NSP8s = True
    NSP8e = False
    if x[12675:12685] == "CAAATTACAG":
        NSP8e = True

    NSP9s = False
    if x[12685:12695] == "AATAATGAGC":
        NSP9s = True
    NSP9e = False
    if x[13014:13024] == "ACGTCTACAA":
        NSP9e = True

    NSP10s = False
    if x[13024:13036] == "GCTGGTAATGCA":
        NSP10s = True
    NSP10e = False
    if x[13431:13441] == "CATGCTTCAG":
        NSP10e = True

    NSP11s = False
    NSP12s = False
    if x[13441:13451] == "TCAGCTGATG":
        NSP11s = True
        NSP12s = True
    NSP11e = False
    if x[13473:13483] == "TGCGGTGTAA":
        NSP11e = True

    NSP12e = False
    if x[16226:16236] == "AGTCTTACAG":
        NSP12e = True

    NSP13s = False
    if x[16236:16246] == "GCTGTTGGGG":
        NSP13s = True
    NSP13e = False
    if x[18029:18039] == "AACTTTACAA":
        NSP13e = True

    NSP14s = False
    if x[18039:18049] == "GCTGAAAATG":
        NSP14s = True
    NSP14e = False
    if x[19610:19620] == "AAGACTTCAG":
        NSP14e = True

    NSP15s = False
    if x[19620:19630] == "AGTTTAGAAA":
        NSP15s = True
    NSP15e = False
    if x[20649:20659] == "AAAATTACAA":
        NSP15e = True

    NSP16s = False
    if x[20659:20669] == "TCTAGTCAAG":
        NSP16s = True
    NSP16e = False
    if x[21546:21556] == "TAACAACTAA":
        NSP16e = True

    UTR2 = False
    if x[21556:21563] == "ACGAACA":
        UTR2 = True

    ORF2s = False
    if x[21563:21573] == "ATGTTTGTTT":
        ORF2s = True
    ORF2e = False
    if x[25375:25385] == "TTACACATAA":
        ORF2e = True

    UTR3 = False
    if x[25385:25393] == "ACGAACTT":
        UTR3 = True

    ORF3s = False
    if x[25393:25405] == "ATGGATTTGTTT":
        ORF3s = True
    ORF3e = False
    if x[26211:26221] == "GCCTTTGTAA":
        ORF3e = True

    UTR4 = False
    if x[26221:26245] == "GCACAAGCTGATGAGTACGAACTT":
        UTR4 = True

    ORF4s = False
    if x[26245:26255] == "ATGTACTCAT":
        ORF4s = True
    ORF4e = False
    if x[26463:26473] == "TCTGGTCTAA":
        ORF4e = True

    UTR5s = False
    if x[26473:26483] == "ACGAACTAAA":
        UTR5s = True
    UTR5e = False
    if x[26513:26523] == "AATTTTAGCC":
        UTR5e = True

    ORF5s = False
    if x[26523:26533] == "ATGGCAGATT":
        ORF5s = True
    ORF5e = False
    if x[27182:27192] == "TGTACAGTAA":
        ORF5e = True

    UTR6 = False
    if x[27192:27202] == "GTGACAACAG":
        UTR6 = True

    ORF6s = False
    if x[27202:27214] == "ATGTTTCATCTC":
        ORF6s = True
    ORF6e = False
    if x[27378:27388] == "GATTGATTAA":
        ORF6e = True

    UTR7 = False
    if x[27388:27394] == "ACGAAC":
        UTR7 = True

    ORF7as = False
    if x[27394:27404] == "ATGAAAATTA":
        ORF7as = True
    ORF7ae = False
    if x[27760:27760] == "GACAGAATGA":
        ORF7ae = True

    ORF7bs = False
    if x[27756:27768] == "ATGATTGAACTT":
        ORF7bs = True
    ORF7be = False
    if x[27878:27888] == "TCACGCCTAA":
        ORF7be = True

    UTR8 = False
    if x[27888:27894] == "ACGAAC":
        UTR8 = True

    ORF8s = False
    if x[27894:27906] == "ATGAAATTTCTT":
        ORF8s = True
    ORF8e = False
    if x[28250:28260] == "TTTCATCTAA":
        ORF8e = True

    UTR9 = False
    if x[28260:28274] == "ACGAACAAACTAAA":
        UTR9 = True

    ORF9s = False
    if x[28274:28286] == "ATGTCTGATAAT":
        ORF9s = True
    ORF9e = False
    if x[29524:29534] == "TCAGGCCTAA":
        ORF9e = True

    UTR10 = False
    if x[29534:29558] == "ACTCATGCAGACCACACAAGGCAG":
        UTR10 = True

    ORF10s = False
    if x[29558:29568] == "ATGGGCTATA":
        ORF10s = True
    ORF10e = False
    if x[29665:29675] == "TCTCACATAG":
        ORF10e = True

    p3UTRs = False
    if x[29675:29685] == "CAATCTTTAA":
        p3UTRs = True
    p3UTRe = False
    if x[29861:29871] == "GGAGAATGAC":
        p3UTRe = True

    if countN / len(x) < 0.005 and countA + countT + countC + countG >= 29750:
        subtotalATCG = countA + countT + countC + countG
        subtotalbad = countN + countdash
        totalval = subtotalATCG + subtotalbad
        totalarray.append([metadata[counter], countA, countT, countC, countG, subtotalATCG, countN, countdash, subtotalbad, tail, FivepUTRs, FivepUTRe, NSP1s, NSP1e, NSP2s, NSP2e, NSP3s, NSP3e, NSP4s, NSP4e, NSP5s, NSP5e, NSP6s, NSP6e, NSP7s, NSP7e, NSP8s, NSP8e, NSP9s, NSP9e, NSP10s, NSP10e, NSP11s, NSP11e, NSP12s, NSP12e, NSP13s, NSP13e, NSP14s, NSP14e, NSP15s, NSP15e, NSP16s, NSP16e, UTR2, ORF2s, ORF2e, UTR3, ORF3s, ORF3e, UTR4, ORF4s, ORF4e, UTR5s, UTR5e, ORF5s, ORF5e, UTR6, ORF6s, ORF6e, UTR7, ORF7as, ORF7ae, ORF7bs, ORF7be, UTR8, ORF8s, ORF8e, UTR9, ORF9s, ORF9e, UTR10, ORF10s, ORF10e, p3UTRs, p3UTRe])
    else:
        removed += 1

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
pd.DataFrame(totalarray,columns=['Sequence ID', 'A', 'T', 'C', 'G','Subtotal', 'N', 'No Call', 'Subtotal', 'PolyA Tail',"5'UTR Start", "5'UTR End","NSP1 Start","NSP1 End","NSP2 Start","NSP2 End", "NSP3 Start", "NSP3 End", "NSP4 Start", "NSP4 End", "NSP5 Start", "NSP5 End", "NSP6 Start", "NSP6 End", "NSP7 Start", "NSP7 End", "NSP8 Start", "NSP8 End", "NSP9 Start", "NSP9 End", "NSP10 Start", "NSP10 End", "NSP11 Start", "NSp11 End", "NSP12 Start", "NSP12 End", "NSP13 Start", "NSP13 End", "NSP14 Start", "NSP14 End", "NSP15 Start", "NSP15 End", "NSP16 Start", "NSP16 End", "UTR2", "ORF2 Start", "ORF2 End", "UTR3", "ORF3 Start", "ORF3 End", "UTR4", "ORF4 Start", "ORF4 End", "UTR5 Start", "UTR5 End", "ORF5 Start", "ORF5 End", "UTR6", "ORF6 Start", "ORF6 End", "UTR7", "ORF7a Start", "ORF7a End", "ORF7b Start", "ORF7b End", "UTR8", "ORF8 Start", "ORF8 End", "UTR9", "ORF9 Start", "ORF9 End", "UTR10", "ORF10 Start", "ORF10 End", "3' UTR Start", "3' UTR End"]).to_csv("counts.csv", index=False)
