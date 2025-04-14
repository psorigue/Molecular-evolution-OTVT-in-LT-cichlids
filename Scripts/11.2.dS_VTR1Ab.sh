#!/bin/bash
# This script calculates the dS ratio of the two exons of VTR1Ab gene separately, using the VTR1Ab gene of the reference genome.
# mafft v7.526
# biopython 1.80

source Scripts/functions_unix.sh # To access functions defined in functions_unix.sh

# Paths
path_aln="Data/11.dS_ratios/in/VTR1Ab_exons/"
path_dS="Data/11.dS_ratios/out/dS_VTR1Ab/"
# Control file yn00
ctl_file="Data/11.dS_ratios/in/yn00.ctl"
# Arrays
exons=( exon1 exon2 )
species=( $(cat "Data/arrays/array_spp_tree.txt") )

cd "${path_dS}"

for ex in "${exons[@]}" ; do

    # Create a directory for each exon
    mkdir "${ex}" ; cd "${ex}"

    # Ensure single-lined fasta alignment
    single_line_fa "${path_aln}/${ex}/${ex}.fa"

    # Loop through each species
    for spp in "${species[@]}" ; do

        # Create a directory for each species
        mkdir "${spp}" ; cd "${spp}"

        # Create pairwise sequence file
        grep -A1 "Orenil\|${spp}" "${path_aln}/${ex}/${ex}.fa" | grep ">\|a" > seq.fa
        
        # Quick alignment
        quick_align_fa_files seq_out.fa seq.fa # This function is taken from functions_unix.sh
        
        # Yn00 function
        python Scripts/functions_python.py yn00_an seq_out.fa "${ctl_file}" # Biopython version 1.80
        
        # Print dS + SE
        dS=$( awk '{print $(NF-2), $(NF-1), $NF}' rst | sed "s/ +- /\t/g" )
        echo -e "${spp}\t${dS}"

        # Go back to the gene directory
        cd "${path_dS}/${ex}"

    done > dS_all_spp.txt

    # Remove blank lines from output file
    sed '/^$/d' dS_all_spp.txt > dS_"${ex}".txt

    cd "${path_dS}"

done

# The output only includes the dataset with all dS ratio values to avoid redundant information.

exit
