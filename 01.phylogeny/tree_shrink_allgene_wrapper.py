#!/usr/bin/env python
import os
import re
import sys
import shutil
import glob

if len(sys.argv) != 4:
    sys.exit("Usage: <Script> <mm_dir> <quantile> <outdir>")

indir, quantile, outdir = sys.argv[1:]

shrinkin = os.path.join(outdir, "input")
shrinkout = os.path.join(outdir, "output")
finalout = os.path.join(outdir, "final")



if not os.path.exists(shrinkin):
    os.makedirs(shrinkin)

if os.path.exists(shrinkout):
    shutil.rmtree(shrinkout)
os.makedirs(shrinkout)


if not os.path.exists(finalout):
    os.makedirs(finalout)


for name in os.listdir(indir):
    if name.endswith(".mm"):
        gene = name.split(".")[0]
        outpath = os.path.join(shrinkin, gene)
        if os.path.exists(outpath):
            shutil.rmtree(outpath)
        os.makedirs(outpath)
        inname = os.path.join(indir, name)
        outname = os.path.join(outpath, "input.tre")
        os.system("cp {} {}".format(inname, outname))

cmd= str("/data/soft/TreeShrink/TreeShrink-master/run_treeshrink.py -i " + shrinkin +
         " -t input.tre -c -m all-genes " +
         " -q " + quantile +
         " -o " + shrinkout)

try:
    os.system(cmd)
except:
    print("Run treeshrink error")
else:
    outtrees = glob.glob(os.path.join(shrinkout, "*", "output.tre"))
    for tree in outtrees:
        gene = tree.split("/")[-2]
        outname = os.path.join(finalout, gene + ".tt")
        os.system("cp {} {}".format(tree, outname))

        

