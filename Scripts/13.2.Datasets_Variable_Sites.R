# This script creates a dataset with three columns: species name, allele, and phenotype.
# It is used to create datasets for the corelation analysis in BayesTraits (BT).
# Allele and phenotype variables are binarized. In the case of ambiguous amino acids,
# the allele is displayes as "-", the symbol for missing data in BT

library(dplyr)

# Paths
path_var_sites <- "Data/13.1.Variable_sites/out/variable_sites/"
path_out_datasets_BT <- "Data/13.Variable_Sites/out/datasets_for_BT_pair_bonding/" # As an example, pair-bonding
# Files
## Template with phenotype data
pheno_template <- read.csv("Data/13.Variable_Sites/in/templates_phenotypes/template_pb.txt", sep = "\t") # As an example, pair-bonding template
# Genes
file_array <- "Data/arrays/transcripts.txt"
list_tr <- as.vector(scan(file_array, what = "character", sep = " "))


for (gen in genes) {
  
  print(gen)
  
  # Read vector of sites
  sites <- as.vector(scan(paste0(path_var_sites, gen, "/", gen, "_all_positions_with_change.txt")))
  
  for (sit in sites) {
    
      print(c(gen,sit))
      
      if (file.exists(paste0(path_var_sites, gen, "/", gen, "_site_", sit, ".txt"))) {
        
          df <- read.delim(paste0(path_var_sites, gen, "/", gen, "_site_", sit, ".txt"), sep = "\t")
          colnames(df) <- c("Spp", "Allele")
          
          df$Allele <- gsub(df$Allele, pattern = "X", replacement = "-")
          
          # Merge phenotype
          merged_df <- merge(df, pheno_template, by = "Spp", all.x = T)
          
          # Omit NAs
          na_df <- na.omit(merged_df)
          
          # Output file
          write.table(na_df, file = paste0(path_out_datasets_BT, gen, "_", sit, ".txt"), quote = F, sep = "\t", col.names = F, row.names = F)
          
      }
  }

}