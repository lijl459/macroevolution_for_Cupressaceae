#!/usr/bin/env python
import sys, re, os
from Bio import SeqIO

if len(sys.argv) != 3:
    print("Usage: <Script> <in_root> <outfile>")
    exit()

k = re.compile("\s+")

inroot = sys.argv[1]
outroot = sys.argv[2]

dic = {}

with open(outroot, "w") as f:
    dic = {}
    for i in sorted(os.listdir(inroot)):
        if i.endswith(".fas"):
            inpath = os.path.join(inroot, i)
            lst = []
            gene = i.split(".")[0]
            for seq_record in SeqIO.parse(inpath, "fasta"):
                name = str(seq_record.id)
                sequences = str(seq_record.seq)
                seqlen = len(sequences)
                lst.append(seqlen)
                if name not in dic:
                    dic[name] = ">%s\n" % name
                seq = ""
                for a in range(0, len(sequences), 3):
                    seq += sequences[a]
                    seq += sequences[a+1]
                dic[name] += seq
    for name in dic.keys():
        f.write(dic[name] + "\n")
    

    
