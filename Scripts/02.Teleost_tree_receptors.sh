#!/bin/bash

# This script is used to create a phylogenetic tree with teleost nonapeptide receptor sequences
# MAFFT version 7.526
# IQTREE version 2.0.3

path_seq="Data/01.Telost_tree_receptors/in/sequences_receptors.fa"
folder_out="Data/01.Telost_tree_receptors/out/"

# Create directory
mkdir ${folder_out}/tree_receptors ; cd ${folder_out}/tree_receptors

# Multiple Sequence Alignment
mafft --auto ${path_seq} > aln_receptors.fa # MAFFT version 7.526

# Phylogenetic Tree
iqtree -s aln_receptors.fa # IQTREE version 2.0.3

exit
