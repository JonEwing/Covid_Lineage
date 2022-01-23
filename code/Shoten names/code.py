import pandas as pd
from Bio import pairwise2
from Bio.pairwise2 import format_alignment
import matplotlib.pyplot as plt
import datetime
import matplotlib.pyplot as plt
import numpy

#Change of Ref v other samples
#pull out just us samples month and week

f = open("out.fasta", "r") 
newf = open("shortname.fasta", "w") 

metadata = []
data = []
counter = 0
for x in f:
    if counter == 0:
        x = x
    elif x[0] == '>':
        x = x[9:len(x)]
        x = ">" + x
    counter += 1
    newf.write(x)

f.close()
newf.close()