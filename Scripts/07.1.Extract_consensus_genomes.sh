#!/bin/bash
# This script extracts the consensus sequences for each gene and individual from the BAM files and merges them to create a consensus for each species obtained from genomic sequences
# To reduce space, the consensus sequences are only shown in the alignment files (Data 08)
# mafft v7.526
# seqtk version 1.4-r122
# bcftools Version: 1.9
# samtools Version: 1.6

source "Scripts/functions_unix.sh"

# Paths
bam_files_directory="Local/path/to/bam_files_directory" # Directory containing BAM files
output_consensus_individual="Local/path/to/output/consensus_individual" # Directory to store individual consensus sequences
output_consensus_merged="Local/path/to/output/consensus_merged" # Directory to store merged consensus sequences
# Files
file_regions="Data/04.Download_Ref_CDS_regions/out/regions_orenil_cds.txt"
ref_genome_file="Local/path/to/reference/genome/fasta_file.fa" # Accession number GCF_001858045.2
# Arrays
speciesID=( $( cat "Data/arrays/array_spp_genomes.txt" ) )
species_red=( $( cat "Data/arrays/array_spp_tree.txt" ) )
genes=( $( cat "Data/arrays/transcripts.txt" ) )



# Genomes
for gen in "${genes[@]}"; do

    cd "${output_consensus_individual}"

    mkdir "${gen}" ; cd "${gen}"
    for spp in "${speciesID[@]}"; do
        
        bam_file="${bam_files_directory}/${spp}_bwa.bam"
        thr=2

        mkdir "${spp} ; cd ${spp}"
        get_consensus_gen "${gen}" "${file_regions}" "${spp}" "${bam_file}" "." "${ref_genome_file}" "${thr}" # Function taken from Scripts/functions_unix.sh.
        cp "${spp}"_"${gen}"_cons.fa ../ # Copy outside the species folder

        cd "${output_consensus_individual}/${gen}"

    done

    # Merge consensus of individuals to generate spp consensus
    cd "${output_consensus_merged}" ; mkdir "${gen}" ; cd "${gen}"

    for spp in "${species_red[@]}"; 

        spplc=${spp,} # Lowercase the first letter of the species name
        num=$( ls "${output_consensus_individual}/${gen}${spplc}*cons.fa" | wc -l ) # Counts how many consensus sequences there are for the species

        if [[ ! $num == 1 ]] ; then # If the species has more than one individual
            
            cat "${output_consensus_individual}/${gen}${spplc}*cons.fa" > "${spplc}.fa" # Concatenates consensus of different individuals in one file
            sed -i 's/>/\n>/g' "${spplc}.fa" # Separates lines in case they are concatenated without \n
            mafft --auto "${spplc}.fa" > "${spplc}_aln.fa" # Quick alignment to put both sequences in the same reading frame. MAFFT v7.526
            split_fasta "${spplc}_aln.fa" "_div" # Splits aligned sequences to make the merge (two different files are required for the merge)
            seqtk mergefa "${spplc}"*_div.fa > "${spplc}_${gen}_cons4h.fa" # Merge individual sequences into species consensus
            single_line_fa "${spplc}_${gen}_cons4h.fa" # Make single-line fasta
            head=$( grep '>' "${spplc}_${gen}_cons4h.fa" )
            sed -i "s/${head}/>${spp}/g" "${spplc}_${gen}_cons4h.fa" # Change fasta header by the species and gene names
            
            rm "${spplc}.fa" "${spplc}_aln.fa" "${spplc}"*_div.fa # Remove intermediate files
                
        else # If the species has only one individual
                
            cp "${path_cons}/${gen}/${spplc}*cons.fa" "${spplc}_${gen}_cons4h.fa"
            single_line_fa "${spplc}_${gen}_cons4h.fa" # Make single-line fasta
            head=$( grep '>' "${spplc}_${gen}_cons4h.fa" )
            sed -i "s/${head}/>${spp}/g" "${spplc}_${gen}_cons4h.fa" # Change fasta header by the species and gene names
                        
        fi

    done

done


exit