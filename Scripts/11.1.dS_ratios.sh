#!/bin/bash 
# This script calculates the dS ratios of each gene and species separately, using the genes of the reference genome.
# It uses Yang and Nielsen method (Yang & Nielsen, 2000).
# mafft v7.526
# biopython 1.80

source Scripts/functions_unix.sh # To access custom functions defined in functions_unix.sh

# Paths
path_aln="Data/08.Alignments/out/"
path_dS="Data/11.dS_ratios/out/dS_genes/"
# Control file yn00
ctl_file="Data/11.dS_ratios/in/yn00.ctl"
# Arrays
genes=( $(cat "Data/arrays/transcripts.txt") )
species=( $(cat "Data/arrays/array_spp_tree.txt") )

cd "${path_dS}"

for gen in "${genes[@]}" ; do

    # Create a directory for each gene
    mkdir "${gen}" ; cd "${gen}"

    # Ensure single-lined fasta alignment
    single_line_fa "${path_aln}/${gen}/aln_${gen}_CDS_p2n.fa"

    # Loop through each species
    for spp in "${species[@]}" ; do

        # Create a directory for each species
        mkdir "${spp}" ; cd "${spp}"

        # Create pairwise sequence file
        grep -A1 "Orenil\|${spp}" "${path_aln}/${gen}/aln_${gen}_CDS_p2n.fa" | grep ">\|a" > seq.fa
        
        # Quick alignment
        quick_align_fa_files seq_out.fa seq.fa # This function is taken from functions_unix.sh
        
        # Yn00 function
        python Scripts/functions_python.py yn00_an seq_out.fa "${ctl_file}" # Biopython version 1.80
        
        # Print dS + SE
        dS=$( awk '{print $(NF-2), $(NF-1), $NF}' rst | sed "s/ +- /\t/g" )
        echo -e "${spp}\t${dS}"

        # Go back to the gene directory
        cd "${path_dS}/${gen}"

    done > dS_all_spp.txt

    # Remove blank lines from output file
    sed '/^$/d' dS_all_spp.txt > dS_"${gen}".txt

    cd "${path_dS}"

done

# The output only includes the dataset with all dS ratio values to avoid redundant information.

exit
