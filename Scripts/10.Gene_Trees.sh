#!/bin/bash
# This script is used to generate gene trees from the alignments produced in script 08, using both nucleotide and amino acid sequences.
# iqtree version 2.0.3

# Paths
path_aln="Data/08.Alignments/out/"
path_trees="Data/10.Gene_Trees/out/"
# Array genes
genes=( $(cat "Data/arrays/transcripts.txt") )

cd "${path_trees}"
for gene in "${genes[@]}"; do

    # Create a directory for each gene and subdirectories for aa and cds trees
    mkdir "${gene}" ; mkdir "${gene}/${gene}_aa_tree" ; mkdir "${gene}/${gene}_cds_tree"
        
    # CDS
    cd "${gene}/${gene}_cds_tree/"
    cp "${path_aln}/${gene}/aln_${gene}_CDS_p2n.fa" . ; iqtree -s "aln_${gene}_CDS_p2n.fa" # version 2.0.3
    
    # aa
    cd "${path_trees}/${gene}/${gene}_aa_tree"
    cp "${path_aln}/${gene}/aln_${gene}_aa.fa" . ; iqtree -s "aln_${gene}_aa.fa"


done
