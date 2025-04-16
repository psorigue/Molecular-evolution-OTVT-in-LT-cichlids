#!/bin/bash
# This script contains functions for different parts of the pipeline

# Function 1: Create single-line fasta in-place
single_line_fa () {
    fa_file=$1
    sed -i "s/\r//g" "$fa_file"
    awk -i inplace '/^>/ {print (NR>1?"\n":"")$0;; next} {printf "%s",$0;} END{print "";}' "$fa_file"
}

# Function 2: Create space-separated list from column
col_to_space_sep_list () {
    column=$1
    cat "$column" | sed -e "s/\r//g" | sed -e :a -e '/$/N; s/\n/ /; ta;' | awk '{$1=$1;print}' | sed '/^\s*$/d'
}

# Function 3: Download CDS and regions file
down_cds_loc () {
    dataset=$1
    out_path=$2
    ref_gtf=$3
    prefix=$4

    sed -i "s/\r//g" "$dataset"
    list_tr=( $(cat "$dataset" | cut -f3 | col_to_space_sep_list) ) # Takes the third column (transcript ID) and turns it into an array

    mkdir -p "$out_path/genes/"

    for tr in "${list_tr[@]}" ; do 
        # Download CDS and make it single-line FASTA
        gen_name=$(cat "$dataset" | grep -w "$tr" | cut -f1)
        datasets download gene accession --include cds "$tr" # datasets version: 14.6.0
        unzip -q ncbi_dataset.zip
        single_line_fa ncbi_dataset/data/cds.fna 
        grep "$tr" -A1 ncbi_dataset/data/cds.fna > "$out_path/genes/${prefix}_${gen_name}_cds.fa"
        rm -r ncbi_dataset ncbi_dataset.zip README* md5sum*

        # Download locations
        ## Gene location (reg_gene)
        gen_id=$(cat "$dataset" | grep -w "$tr" | cut -f2)
        reg_gene=$(cat "$ref_gtf" | grep -w "$gen_id" | grep -v "NW_" | grep -P '\tgene\t' | cut -f1,4,5 | awk '{ print $1,":",$2,"-",$3 }' | sed 's/ //g')

        ## Concatenate all regions in one comma-separated line
        reg_ex=$(cat "$ref_gtf" | grep "$tr" | grep -P '\tCDS\t' | cut -f1,4,5 | awk '{ print $1,":",$2,"-",$3 }' | sed 's/ //g' | sed -e :a -e '/$/N; s/\n/,/; ta;')

        ## Sense
        if [[ $(cat "$ref_gtf" | grep "$gen_id" | grep -P '\tgene\t' | cut -f7) == "-" ]] ; then 
            sense=-1 
        else 
            sense=1 
        fi

        ## Output in regions file
        echo -e "$gen_name\t$gen_id\t$reg_gene\t$tr\t$reg_ex\t$sense\t."
    done > "$out_path/regions.txt"

    sed -i '/^$/d' "$out_path/regions.txt"
}

# Function 4: Per-site coverage
ps_coverage () {
    bam_file_list=$1 # List of files to assess the coverage
    regions=$2 # Semicolon-separated exon regions
    output=$3 # Output file name

    # Create array of regions
    ar_reg=( $(echo "$regions" | sed 's/;/ /g') )

    # Coverage function
    for exon in "${ar_reg[@]}" ; do
        samtools depth -aa -r "$exon" -f "$bam_file_list" # Version 1.6
    done > "$output"
}

# Function 5: Extract consensus sequence
get_consensus_gen () {
    id=$1
    file_regions=$2 
    spp=$3
    bam_file=$4
    out_path=$5
    ref_genome=$6
    threads=$7

    reg_bcf=$(grep -w "$id" "$file_regions" | cut -f5)
    ref_genome=~/fil/Ref_genome/GCF_001858045.2_O_niloticus_UMD_NMBU_genomic.fa

    cd "$out_path" || exit 1

    # BAM to BCF
    bcftools mpileup --threads "$threads" -Ou -o "${spp}_${id}.bcf" -r "$reg_bcf" -f "$ref_genome" "$bam_file" # bcftools Version: 1.9
    # BCF to VCF
    bcftools call --threads "$threads" -vmO v -o "${spp}_${id}.vcf" "${spp}_${id}.bcf"
    # Index VCF (required for consensus)
    bgzip -@ "$threads" "${spp}_${id}.vcf" ; bcftools index --threads "$threads" "${spp}_${id}.vcf.gz" ; gunzip -k "${spp}_${id}.vcf.gz"
    # Create array of regions
    IFS=',' read -ra regions <<< "$reg_bcf"
    # Extract consensus by region
    for reg in "${regions[@]}" ; do
        samtools faidx -@ "$threads" "$ref_genome" "$reg" | bcftools consensus -HA -I "${spp}_${id}.vcf.gz"
    done > "seq_${id}_${spp}_by_region.fa"
    # Join fragments in one single-lined fasta
    only_seq=$(grep -v '>' "seq_${id}_${spp}_by_region.fa")
    echo -e ">${spp}"'\n'"${only_seq}" > "${spp}_${id}_cons_sl.fa"
    single_line_fa "${spp}_${id}_cons_sl.fa"
    # Reverse-complement if needed
    if [ "$(grep -w "$id" "$file_regions" | cut -f6)" == -1 ] ; then
        seqtk seq -r "${spp}_${id}_cons_sl.fa" > "${spp}_${id}_cons.fa" # version 1.4-r122
    elif [ "$(grep -w "$id" "$file_regions" | cut -f6)" == 1 ] ; then
        mv "${spp}_${id}_cons_sl.fa" "${spp}_${id}_cons.fa"
    fi
    # Remove useless files
    rm "${spp}_${id}_cons_sl.fa" "seq_${id}_${spp}_by_region.fa"
}

# Function 6: Split fasta files into separate files
split_fasta () {
    fa_file="$1"
    suffix="$2"
    awk '/^>/{if(x){close(x)}x=substr($0,2) "'"${suffix}"'.fa";print > x;next}{print > x}' "${fa_file}"
}

# Function 7: Pick best coverage per spp (genome or transcriptome) and, in case transcriptome has less 0-covered positions, it replaces genome consensus by transcriptome consensus.
group_best_cov () {
    gen=$1
    spp=$2
    spp_uc=${spp^}
    file_seq=$3 # File with the fasta sequences from genome
    best_cov_file=$4 # Format: CSV with col = list_genes, rows = list_spp
    file_cons_tr=$5

    num_col=$(head -n1 "$best_cov_file" | sed "s/,/\t/g" | tr \\t \\n | grep -w -i -n "$gen" | awk -F: '{print $1}') # Print column where the gene is
    cov=$(grep -e "$spp" "$best_cov_file" | cut -f"$num_col" -d,) # Either gn or tr

    if [[ $cov == "tr" ]] ; then 
        sed -i -e "/$spp_uc/,+1d" "$file_seq" # Delete head and sequence of the spp to remove

        # Add line feed between the alignment file and the new sequence added
        last_char=$(tail -c 1 "$file_seq")

        # Check if the last character is not a newline
        if [ "$last_char" != $'\n' ]; then
            # Append a newline character to the end of $file_seq
            echo >> "$file_seq"
        fi

        cat "$file_cons_tr" >> "$file_seq" # Append the transcriptome consensus
        sed -i "s/$(grep "$spp" "all_${gen}_CDS.fa")/>$spp_uc/g" "$file_seq" # Change the head of the sequence for the spp name
    fi
}

# Function 8: Translate and align
aln_trn () {
    gen=$1
    file=$2
    threads=$3

    echo "$gen"

    sed -i s/-//g "$file" # Delete gaps

    seqkit translate --threads "$threads" -x "$file" > "all_${gen}_aa.fa" # Output aa sequences. seqkit v2.8.2

    mafft --auto --thread "$threads" "$file" > "aln_${gen}_CDS.fa" # Output alignment CDS. mafft v7.526
    mafft --auto --thread "$threads" "all_${gen}_aa.fa" > "aln_${gen}_aa.fa" # Output alignment aa
}

# Function 9: Quick alignment of separate files
quick_align_fa_files () {
    out_file=$1
    sequences=${@:2}

    # Align sequences
    cat $sequences > tmp_seq.fa
    mafft --auto tmp_seq.fa > "$out_file"  # mafft v7.526
    single_line_fa "$out_file" ; rm tmp_seq.fa
}
