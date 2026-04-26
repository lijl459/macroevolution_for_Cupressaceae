#!/usr/bin/python
import sys, os, re
from ete3 import Tree

if len(sys.argv) != 5:
    print("Usage: <script> <in_MO> <in_RT> <in_out.txt> <outdir>")
    exit()


in_MO, in_RT, in_out, outdir = sys.argv[1:]

inlist = []
outlist = []

#os.chdir(outdir)

with open(in_out) as f:
    for line in f:
        l = re.split("\s+", line)
        if l[0] == "OUT":
            outlist.append(l[1])
        elif  l[0] == "IN":
            inlist.append(l[1])


for name in os.listdir(in_RT):
    if name.endswith(".tre"):
        gene = name.split(".")[0]
        RTpath = os.path.join(in_RT, name)
        MOpath = os.path.join(in_MO, gene + ".ortho.tre")
        SCGpath = os.path.join(in_MO, gene + ".1to1ortho.tre")
        RTtree = Tree(RTpath)
        outleaves = []
        outpath = os.path.join(outdir, name)

        if os.path.exists(SCGpath):
            #print(gene)
            continue

        if os.path.exists(MOpath):
            MOtree = Tree(MOpath)
            leaves = MOtree.get_leaf_names()
            for leaf in leaves:
                sample = leaf.split("@")[0]
                if sample in outlist:
                    outleaves.append(leaf)
            if len(outleaves) > 0:
                outtree = MOtree
                outtree.prune(outleaves)
                t = Tree()
                t.add_child(outtree)
                t.add_child(RTtree)
            else:
                t = RTtree
                
        else:
            t = RTtree
        t.write(outfile=outpath)
