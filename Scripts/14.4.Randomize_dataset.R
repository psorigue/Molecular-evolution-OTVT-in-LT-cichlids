# This script randomizes the dataset 200 times. It creates a new dataset with the same species and alleles but with randomized phenotypes.

dataset <- "Path/to/dataset_file.txt" # Path to the dataset
# The dataset (output from script 14.2) contains three columns: 1) Species names, 2) Allele variant, 3) Phenotype state


dt <- read.csv(file = dataset, header = F, sep = "\t") # Read dataset
vec_phe <- dt[,3] # Read vector from phenotype column
  
for (i in 1:200) { # Randomize 200 times
    
rdz <- sample(vec_phe) # Create randomized vector from phenotype column
dt_rdz <- cbind(dt[,c(1,2)], rdz) # Join randomized vector to Spp and Allele columns
write.table(dt_rdz, file = paste0(path_rdz, "rdz", i, "_", sit, ".txt"), sep = "\t", col.names = F, row.names = F, quote = F) # Save dataset
    
  }
