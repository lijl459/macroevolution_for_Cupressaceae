#!/usr/bin/env python
import sys, re, os
from Bio import SeqIO

if len(sys.argv) != 3:
    print("Usage: <Script> <in_root> <outroot>")
    exit()

k = re.compile("\s+")

inroot = sys.argv[1]
outroot = sys.argv[2]

dic = {}

for i in sorted(os.listdir(inroot)):
    if i.endswith(".fas"):
        inpath = os.path.join(inroot, i)
        lst = []
        gene = i.split(".")[0]
        with open(os.path.join(outroot, gene+".codon12.fas"), "w") as f:
            for seq_record in SeqIO.parse(inpath, "fasta"):
                name = str(seq_record.id)
                sequences = str(seq_record.seq)
                seqlen = len(sequences)
                lst.append(seqlen)
                if name not in dic:
                    dic[name] = ""
                seq = ""
                for a in range(0, len(sequences), 3):
                    seq += sequences[a]
                    seq += sequences[a+1]
                f.write(">%s\n%s\n" % (name, seq))    
        

    
