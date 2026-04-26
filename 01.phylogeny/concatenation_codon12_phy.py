#!/usr/bin/env python
import sys, re, os
from Bio import SeqIO

if len(sys.argv) != 3:
    print("Usage: <Script> <in_root> <outname>")
    exit()

k = re.compile("\s+")

inroot = sys.argv[1]
outroot = sys.argv[2]

dic = {}

n = 0
spset = set()

for i in sorted(os.listdir(inroot)):
    if i.endswith(".fas"):
        inpath = os.path.join(inroot, i)
        for seq_record in SeqIO.parse(inpath, "fasta"):
            name = str(seq_record.id)
            spset.add(name)

for name in spset:
    dic[name] = ""

print(len(spset))
count = 0
for i in sorted(os.listdir(inroot)):
    if i.endswith(".fas"):
        inpath = os.path.join(inroot, i)
        namelist = list(spset)
        for seq_record in SeqIO.parse(inpath, "fasta"):
            name = str(seq_record.id)
            namelist.remove(name)
            sequences = str(seq_record.seq)
            dic[name] += sequences
        gene = i.split(".")[0]
        n0 = n
        n += len(sequences)
        seqlen = len(sequences)
        for individual in namelist:
            dic[individual] += "-" * seqlen

with open(outroot, "w") as outfile:
    #del dic["Lich.1_0"]
        
    for name in sorted(dic):
        nseq = len(dic[name])
        nsp = len(dic)
        break
    
    if nseq % 3 != 0:
        sys.exit("length of sequence is is not divisible by 3")

    nseq = nseq/3
    nsp = nsp

    for i in range(0, 2):
        outfile.write("%d\t%d\n" % (nsp, nseq))
        for name in sorted(dic):
            s = ""
            for j in range(0, len(dic[name]), 3):
                s += dic[name][j+i]
            outline = name + " "*10 + s + "\n"
            outfile.write(outline)
        outfile.write("\n")
