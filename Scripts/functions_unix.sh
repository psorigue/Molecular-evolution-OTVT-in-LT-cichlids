#!/bin/bash
# This script contains functions for different parts of the pipeline


# Function 1: Create single-line fasta in-place
single_line_fa () {
 
    fa_file=$1
    
    sed -i "s/\r//g" $fa_file
    awk -i inplace '/^>/ {print (NR>1?"\n":"")$0;; next} {printf "%s",$0;} END{print "";}' $fa_file
  
}

#Function 2: Create space-separated list from column
###
col_to_space-sep_list () {

    column=$1
    
    cat $column | sed -e "s/\r//g" | sed -e :a -e '/$/N; s/\n/ /; ta;' | awk '{$1=$1;print}' | sed '/^\s*$/d'
    
}

# Function 3: Download CDS and regions file
down_CDS_loc () {

    dataset=$1
    out_path=$2
    ref_gtf=$3
    prefix=$4
    
    sed -i "s/\r//g" $dataset
    list_tr=( $(cat $dataset | cut -f3 | col_to_space-sep_list ) ) # Takes the third column (transcript ID) and turns it into an array
    
    mkdir $out_path/genes/
    
    for tr in ${list_tr[@]} ; do 
    
      #Download CDS and make it single line FASTA
      gen_nam=$( cat $dataset | grep -w $tr | cut -f1 )
      datasets download gene accession --include cds $tr # datasets version: 14.6.0
      unzip -q ncbi_dataset.zip
      single_line_fa ncbi_dataset/data/cds.fna 
      grep $tr -A1 ncbi_dataset/data/cds.fna > $out_path/genes/"$prefix"_"$gen_nam"_cds.fa
      rm -r ncbi_dataset ncbi_dataset.zip README* md5sum*
      
      #Download locations
      ##Gene location (reg_gene)
      gen_id=$( cat $dataset | grep -w $tr | cut -f2 )
      reg_gen=$( cat $ref_gtf | grep -w $gen_id | grep -v "NW_" | grep -P '\tgene\t' | cut -f1,4,5 | awk '{ print $1,":",$2,"-",$3 }' | sed 's/ //g' )
      
      ##Concatenate all regions in one comma separated line
      reg_ex=$( cat $ref_gtf | grep $tr | grep -P '\tCDS\t' | cut -f1,4,5 | awk '{ print $1,":",$2,"-",$3 }' | sed 's/ //g' | sed -e :a -e '/$/N; s/\n/,/; ta;' )
    
      ##Sense
      if [[ $( cat $ref_gtf | grep $gen_id | grep -P '\tgene\t' | cut -f7 ) == "-" ]] ; then sen=-1 ; else sen=1 ; fi
      
      ##Output in regions file
      echo -e "$gen_nam\t$gen_id\t$reg_gen\t$tr\t$reg_ex\t$sen\t."
  
    done > $out_path/regions.txt

    sed -i '/^$/d' $out_path/regions.txt

}

#Function 4: Per-site coverage
###
ps_coverage () {

    bam_file_list=$1 # List of files to assess the coverage
    regions=$2 # Semicolon-separated exon regions
    output=$3 # Output file name
    
    #Create array of regions
    ar_reg=( $( echo "$regions" | sed 's/;/ /g') )
    
    #Coverage function
    for exon in ${ar_reg[@]} ; do

        samtools depth -aa -r $exon -f $bam_file_list # Version 1.6
                
    done > "$output"

}