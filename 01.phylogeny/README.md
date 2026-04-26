# Phylogenomic analyses

This directory contains scripts for sequence processing and phylogenetic reconstruction.

## Main steps

- CDS processing and filtering
- Codon position extraction
- Sequence concatenation
- Alignment (MAFFT)
- Phylogenetic inference (IQ-TREE)

## Key scripts

- `codon12_from_cds.py`  
  Extract codon positions

- `concatenation_codon12_phy.py`  
  Concatenate sequences

- `phyloneny_mafft.cmd`  
  Alignment using MAFFT

- `iqtree.command`  
  Phylogenetic inference

- `tree_shrink_allgene_wrapper.py`  
  Remove long branches

