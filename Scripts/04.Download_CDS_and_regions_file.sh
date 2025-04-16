#!/bin/bash
# This script downloads the reference transcripts from NCBI and creates a regions file that includes the locations of the CDS regions
# datasets version: 14.6.0

source Scripts/functions_unix.sh # To access custom functions defined in functions_unix.sh

file_ids="Data/04.Download_Ref_CDS_regions/in/acc_numbers_genes.txt" # Information taken online. The format must be "GENE NAME \t GENE ID \t TRANSCRIPT ID"
path_out="Data/04.Download_Ref_CDS_regions/out"
ref_genome_gtf="Local/path/to/reference/genome/gtf_file.gtf" # Accession number GCF_001858045.2

prefix="orenil"

cd ${path_out}

down_CDS_loc ${file_ids} ${path_out} ${ref_genome_gtf} ${prefix} # Function taken from functions_unix.sh. datasets version: 14.6.0

exit