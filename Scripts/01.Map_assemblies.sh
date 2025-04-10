# All raw reads are publicly available from the NCBI SRA database under BioProject PRJNA550295 (genomes) and PRJNA552202 (transcriptomes)
# I did not include the reference genome (Oreochromis niloticus) in the data but it is available from NCBI under the accession number GCF_001858045.2
# The following script maps the raw reads to the reference genomes and generates a BAM file for each sample
# BWA-MEM version 0.7.18
# STAR version 2.7.10a
# IGVTools version 2.5.3
# SAMtools version 1.21



# Path raw reads
path_genome_reads="Local/path/to/genome/reads"
path_transcriptome_reads="Local/path/to/transcriptome/reads"
# Path arrays of species
spp_genomes=( $( cat "Data/01.Map_assemblies/in/array_spp_genomes.txt" ) )
spp_transcriptomes=( $( cat "Data/01.Map_assemblies/in/array_spp_transcriptomes.txt" ) )
# Path reference genome files
ref_genome_file="Local/path/to/reference/genome/fasta_file.fa" # Accession number GCF_001858045.2
ref_genome_gtf="Local/path/to/reference/genome/gtf_file.gtf" # Accession number GCF_001858045.2
ref_genome_folder="Local/path/to/reference/genome/folder"
# Path_output
path_out_genomes="Local/path/to/output/directory/genomes"
path_out_transcriptomes="Local/path/to/output/directory/transcriptomes"
# Threads
thr=34 # Optional


# Map genomes
# First step create bwa index for the reference genome
bwa index "${ref_genome_file}"

# Second step map genomes
cd "${path_out_genomes}"
for spp in "${spp_genomes[@]}"; do
  
    mkdir "${spp}"

    # Map genomes
    bwa mem -t "${thr}" -v 1 "${ref_genome_file}" "${path_genome_reads}/${spp}_1.fastq" "${path_genome_reads}/${spp}_2.fastq" > "${spp}/${spp}_bwa.sam"
    # bwa-mem version 0.7.18

    # Convert SAM to BAM
    samtools view -@ "${thr}" -Sb "${spp}/${spp}_bwa.sam" > "${spp}/${spp}_bwa.bam" # SAMtools version 1.21

    # If BAM file is not empty, remove SAM
    if [ -s "${spp}/${spp}_bwa.bam" ]; then rm "${spp}/${spp}_bwa.sam" ; fi

    # Sort BAM
    samtools sort -@ "${thr}" -o "${spp}/${spp}_bwa.sorted.bam" "${spp}/${spp}_bwa.bam"

    # If sorted BAM is not empty, remove BAM
    if [ -s "${spp}/${spp}_bwa.sorted.bam" ]; then rm "${spp}/${spp}_bwa.bam" ; fi

    # Index BAM
    samtools index "${spp}/${spp}_bwa.sorted.bam"

    # Compress fastq files
    gzip "${path_genome_reads}/${spp}"_[12].fastq

done


# Map transcriptomes
# First step generate genome indices for STAR aligner (check Manual)
cd ${ref_genome_folder} ; mkdir STAR_genome
STAR --runThreadN $thr --runMode genomeGenerate --genomeDir ./STAR_genome/ --genomeFastaFiles $ref_genome_file --sjdbGTFfile $ref_gtf --genomeSAindexNbases 13 --sjdbOverhang 99

path_STAR_genome="${ref_genome_folder}/STAR_genome/"

# Second step map transcriptomes
cd "${path_out_transcriptomes}"

for spp in "${spp_transcriptomes[@]}"; do

    mkdir "${spp}"

    # Map transcriptomes
    STAR --genomeDir $path_STAR_genome --sjdbGTFfile "${ref_genome_gtf}" --readFilesIn "${path_transcriptome_reads}/${spp}.fastq" --readFilesCommand zcat --outFileNamePrefix "${spp}"_ --outSAMtype BAM SortedByCoordinate --outFilterMultimapNmax 1 --outFilterMatchNminOverLread 0.4 --outFilterScoreMinOverLread 0.4 --runThreadN 20
    # STAR version=2.7.10a 

    # Index BAM
    igvtools index ./"${spp}"*.bam # IGV Version 2.5.3 

done

exit