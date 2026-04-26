#/usr/bin/env python

import os
import sys

if len(sys.argv) != 4:
    sys.exit("Script <03.absrel> <species_tree_witd_Node> <OUT_prefix>")

indir, intree, outname = sys.argv[1:]


for filename in os.listdir(indir):
    if filename.endswith(".log"):
        with open(filename) as f:
            for line in f:
                
