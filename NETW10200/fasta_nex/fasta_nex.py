lista = ["AlignedFilteredUS.fasta"]

for x in lista:
    fp = open('data/' + x, 'r')
    Nfp = open('out/' + x[:len(x)-6]+ '.nex', 'w')

    holdarray = fp.readlines() 

    holdstr = ""
    holdarr = []
    for x in range(len(holdarray)):
        if x == 0:
            name = holdarray[x][9:len(holdarray[x])-1]
            if name == 'e':
                name = holdarray[x][1:len(holdarray[x])-4]
        elif ">" in holdarray[x]:
            holdarr.append([name,holdstr])
            name = holdarray[x][9:len(holdarray[x])-1]
            holdstr = ""
        else:
            holdstr += holdarray[x][:len(holdarray[x])-1]
    holdarr.append([name,holdstr])

    Nfp.write("#NEXUS\n\n")
    Nfp.write("begin data;\n")
    ntax = len(holdarr)
    nchar = len(holdarr[10][1])
    Nfp.write("dimensions ntax=" + str(ntax) + " nchar=" + str(nchar) + ";\n")
    Nfp.write("format datatype=dna gap=-;\n")
    Nfp.write("matrix\n")

    for x in holdarr:
        Nfp.write(x[0] + " " + x[1] + "\n")

    Nfp.write(";\n")
    Nfp.write("end;\n")
