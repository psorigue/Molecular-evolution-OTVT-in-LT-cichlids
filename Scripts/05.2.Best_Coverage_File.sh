#!/bin/bash
# This script creates a file indicating the best coverage per species (either genome of transcriptome BAM). 

# Paths
folder_gn="Data/05.Coverage/out/cov_genomes/"
folder_tr="Data/05.Coverage/out/cov_transcriptomes/"
# Arrays
genes=( $( cat "Data/arrays/transcripts.txt" ) )
spp_file="Data/arrays/array_spp_transcriptomes.txt" # Array of species that have brain transcriptomes 

suffix=otvt


# Sum 0-covered positions
cd $folder_gn
python Scripts/functions_python.py sum_0cov $file_array "$suffix"_gn
cd $folder_tr
python Scripts/functions_python.py sum_0cov $file_array "$suffix"_tr

# Best_coverage file
cd $folder_gn
python Scripts/functions_python.py best_cov_file $file_array $spp_file $folder/num0cov_"$suffix"_gn.csv $folder_tr/num0cov_"$suffix"_tr.csv

exit
