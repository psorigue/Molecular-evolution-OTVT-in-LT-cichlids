#!/bin/bash
# This script calculates the coverage depth for each position of CDS regions for genomes and transcriptomes using the BAM files obtained with script 01
# Samtools depth version 1.21

source Scripts/functions_bash.sh # To access custom functions defined in functions_bash.sh

# Files
bam_list_file_gn="Local/path/to/bam_list_file_genomes.txt" # "\n-separated list of bam files from which to analyze the coverage"
bam_list_file_tr="Local/path/to/bam_list_file_transcriptomes.txt" # "\n-separated list of bam files from which to analyze the coverage"
file_regions="Data/04.Download_Ref_CDS_regions/out/regions_orenil_cds.txt"

# Paths
path_cov_gn="Data/05.Coverage/out/cov_genomes"
path_cov_tr="Data/05.Coverage/out/cov_transcriptomes"

# IDs
genes=( $( cat "Data/arrays/transcripts.txt" ) )


cd "${path_cov_gn}"

for gen in "${genes[@]}"; do
    
    # Genomes
    out_file_gn="${path_cov_gn}/coverage_${gen}.txt"
    
    # Create ";"-separated string with regions to input in coverage function
    regions=$( awk -F'\t' -v gen="${gen}" '$1 == gen { gsub(",", ";", $5); print $5 }' "${file_regions}" )

    ps_coverage "${bam_list_file_gn}" "${regions}" "${out_file_gn}" # Function taken from "Scripts/functions_bash.sh". samtools depth version 1.6
    
    cat "coverage_${gen}.txt" > "cov_${gen}.txt"
    rm "coverage_${gen}.txt"


    # Transcriptomes
    cd "${path_cov_tr}"
    
    out_file_tr="${path_cov_tr}/coverage_${gen}.txt"
    
    # Create ";"-separated string with regions to input in coverage function
    regions=$( awk -F'\t' -v gen="${gen}" '$1 == gen { gsub(",", ";", $5); print $5 }' "${file_regions}" )

    ps_coverage "${bam_list_file_tr}" "${regions}" "${out_file_tr}"
    
    cat "coverage_${gen}.txt" > "cov_${gen}.txt"
    rm "coverage_${gen}.txt"

done

exit