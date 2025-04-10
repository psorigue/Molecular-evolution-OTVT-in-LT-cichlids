# This script calculates nucleotide diversity (pi) for each transcript in the dataset.
# It takes the alignment files (output of script 08) and calculates pi using the pegas package.

# Load libraries
library("ape") # version 5.8
library("pegas") # version 1.3

# Paths
path_alignments <- "Data/08.Alignments/out/"
path_out <- "Data/09.Nucleotide_Diversity/out/"
# Genes
file_array <- "Data/arrays/transcripts.txt"
list_tr <- as.vector(scan(file_array, what = "character", sep = " "))

# Create data frame
df <- data.frame()

# Loop through genes and calculate nucleotide diversity
for (gen in list_tr) {

    # Read alignment
    dna <- ape::read.dna(paste(path_alignments, gen, "/aln_", gen, "_CDS_p2n.fa", sep = ""), format = "fasta")
    # Calculate nucleotide diversity
    nd <- nuc.div(dna, variance = F, pairwise.deletion = T)
    # Create vector
    out <- c(gen, nd)
    # Add to data frame
    df <- rbind(df, out)
}


colnames(df) <- c('gene', 'pi')

# Write output
write.table(df, file = paste0(path_out, "nucleotide_diversity_general.txt"),
            row.names = F,
            col.names = T,
            quote = F,
            sep = "\t")
