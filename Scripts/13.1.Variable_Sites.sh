#!/bin/bash

# This script returns a dataset for each variable amino acid site in each transcript. The dataset contains a column with species name
# and a second column with the allele variant: 0 for the major allele, and 1 for the minor allele.
# Ambiguous nucleotides are assigned an 'X'.



source Scripts/functions_unix.sh # To access custom functions defined in functions_unix.sh

# Path
path_aln="Data/08.Alignments/out/"
path_var_sites="Data/13.1.Variable_sites/out/variable_sites/"
# Array
genes=( $(cat "Data/arrays/transcripts.txt") )


cd "${path_var_sites}"

for gen in "${genes[@]}" ; do

    # Create a directory for each gene
    mkdir "${gen}" ; cd "${gen}"

    # Define alignment file
    aln_file="${path_aln}/${gen}/aln_${gen}_aa.fa" # Amino acid alignment

    # Outputs alleles for variable sites for later lm
    python Scripts/functions_python.py allele_site ${aln_file} ${folder}/${gen} ${gen}
    
    cd ${folder}
    
done

exit