#!/bin/bash
# This script extracts the consensus sequences for each gene and individual from the BAM files and merges them to create a consensus for each species obtained from genomic sequences
# seqtk version 1.4-r122
# bcftools Version: 1.9
# samtools Version: 1.6


source Scripts/functions_unix.sh # To access custom functions defined in functions_unix.sh

# Paths
output_consensus_transcriptomes="Local/path/to/output/consensus_transcriptomes" # Directory to store merged consensus sequences
bam_files_transcriptomes="Local/path/to/bam_files_directory" # Directory containing BAM files
# Files
file_regions="Data/04.Download_Ref_CDS_regions/out/regions_orenil_cds.txt"
ref_genome_file="Local/path/to/reference/genome/fasta_file.fa" # Accession number GCF_001858045.2
# Arrays
species=( $( cat "Data/arrays/array_spp_transcriptomes.txt" ) )
genes=( $( cat "Data/arrays/transcripts.txt" ) ) 


for gen in "${genes[@]}"; do

    cd "${output_consensus_transcriptomes}"

    mkdir "${gen}" ; cd "${gen}"
    for spp in "${species[@]}"; do
        
        bam_file="${bam_files_transcriptomes}/${spp}*.bam" # BAM files from STAR aligner
        thr=2

        mkdir "${spp} ; cd ${spp}"
        get_consensus_gen "${gen}" "${file_regions}" "${spp}" "${bam_file}" "." "${ref_genome_file}" "${thr}" # Function taken from Scripts/functions_unix.sh.
        cp "${spp}"_"${gen}"_cons.fa ../"${spp}"_"${gen}"_cons_tr.fa # Copy outside the species folder

        # Change header fasta file
        sppuc=${spp^} #Species name in upper case
        sed "/${spp}/c\>${sppuc}" "${spp}"*cons_tr.fa > "${sppuc}"_"${gen}"_cons_tr.fa # Change header fasta

        cd "${output_consensus_transcriptomes}/${gen}"

    done

done

exit