#!/bin/bash

# This script is used to create a local database from PacBio assemblies for BLAST searches
# makeblastdb version 2.16.0+
# blastn version 2.16.0+

path_ref_genes="Data/02.Local_Database_Blast/in/Referece_Genes/"
path_database="Local/path/to/databases/"
gene_array=( OT OTRa OTRb VT VTR1Aa VTR1Ab VTR2Aa VTR2Ab VTR2Ba VTR2Bb )
species_array=( simdia_AUE1 neomul_IRF6 cyplep_ISI2 cunlon_IWD7 cphfro_LEI6 batmin_IXA5 )


cd "Data/02.Local_Database_Blast/out/"

for spp in ${species_array[@]}; do
  
  # Create local database
  makeblastdb -in ${spp}_assembly.fa -input_type fasta -dbtype 'nucl' -out ${path_database}/${spp}_PBdb # makeblastdb version 2.16.0+
  
  for gen in ${gene_array[@]}; do
    
    mkdir ${gen}
    
    # Run BLAST on local database
    blastn -query ${path_ref_genes}/orenil_${gen}_gene.fa -db ${path_database}/${spp}_PBdb -outfmt 3 > ${gen}/${gen}_${spp}.txt  # blastn version 2.16.0+

  done

done
