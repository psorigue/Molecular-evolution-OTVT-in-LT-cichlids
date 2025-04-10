#!/bin/bash

# This script aligns the species consensus for each gene with the reference genome. Using the coverage information, it selects the best coverage consensus (either genomes or transcriptomes) and, if necessary, replaces genomic consensus by transcriptomic consensus.
# The script makes a nucleotide alignment and a codon-based alignment (pal2nal.pl) for the protein-coding sequences.
# seqkit v2.8.2
# mafft v7.526
# pal2nal.pl v14


source "Scripts/functions_unix.sh"

# Paths
input_directory="Local/path/to/consensus/files/" # The input file is a fasta file with all concatenated cichlid species consensus sequences. The first sequence is the reference sequence. These files contain the same information as the output of this same script (Alignments), so they are not shown in "Data/08.Alignments/in/".
path_alignments="Data/08.Alignments/out/"
path_consensus_transcriptomes="Local/path/to/transcriptome/consensus/"
# Files
best_coverage_file="Data/05.Coverage/out/best_coverage.txt" 
# Array
species=( $( cat "Data/arrays/array_spp_transcriptomes.txt" ) )



cd $path_alignments

for gen in "${genes[@]}"; do # Loop through the genes

    mkdir "${gen}" ; cd "${gen}"

    file_sequences="${input_directory}/${gen}.fa"  # Concatenated fasta file with all species consensus sequences and reference gene
    single_line_fa "${file_sequences}" # Function taken from Scripts/functions_unix.sh

    for spp in ${spp_tr[@]} ; do
        
        cons_tr_file=$path_consensus_transcriptomes/"$spp"_"$gen"_cons_tr.fa
    
        group_best_cov $gen $spp $file_seq $file_sequences $cons_tr_file # Selects best coverage consensus (either genomes or transcriptomes) and, if necessary, replaces genomic consensus by transcriptomic consensus.
        # The function group_best_cov is taken from Scripts/functions_unix.sh.

    done
        
        # Align nucleotides, translate and align amino acids
        aln_trn $gen $file_seq $thr # seqkit v2.8.2. mafft v7.526. The function aln_trn is taken from Scripts/functions_unix.sh.
        
        #Codon-based alignment         
        pal2nal.pl aln_"$gen"_aa.fa aln_"$gen"_CDS.fa -output fasta > aln_"$gen"_CDS_p2n.fa # pal2nal.pl v14
    
done

# Manual corrections were applied before the final alignments output in "Data/08.Alignments".
# These manual corrections are explained in the publication Supplementary Methods

exit


