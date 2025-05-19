
# This script is used to run a phylogenetic regression with the read counts of the brain transcriptomes.
# This script is used for pair-bonding phenotype, but it was also run for the sex of the caregiver.

library(ape) # version 5.8
library(phylolm) # version 2.6.5
library(dplyr)

# Dataset counts
exp_dt <- read.csv("Data/15.Phylogenetic_Regression/in/dataset_counts_spp.txt", sep = "\t", header = T)
# Pair-bonding information
pb_info <- read.csv("Data/13.Variable_Sites/in/templates_phenotypes/template_pb.txt", sep = "\t") # As an example, pair-bonding template
colnames(pb_info) <- c("spp", "phenotype")
# Species tree
tree_phyl <- "Data/14.BayesTraits/in/spp_tree.trees"
treetan <- read.nexus(tree_phyl)
# Isotope data
iso_data <- read.csv("Data/15.Phylogenetic_Regression/in/06_stable_isotope_data.csv", header = T)
# Path output
path_out <- "Data/15.Phylogenetic_Regression/out/"



# Make median iso_data per species
iso <- iso_data %>%
  group_by(X) %>%  # Group by species and tissue
  summarise(across(3:4, \(x) median(x, na.rm = TRUE)), .groups = "drop")  # Compute median for gene count columns
colnames(iso) <- c("spp", "d15N", "d13C")


# PHYLO REGRESSION

# 2. SEXES TOGETHER
# Median by species
dt_med <- exp_dt %>%
  group_by(spp,tissue) %>%  # Group by species and tissue
  summarise(across(3:12, \(x) median(x, na.rm = TRUE)), .groups = "drop")  # Compute median for gene count columns

# Filter brain
dt_med_br <- dt_med[dt_med$tissue == "brain",]
dt_med_br <- as.data.frame(dt_med_br)

# Introduce phenotype
# Merge
dt_expr_mer <- merge(dt_med_br, pb_info, by = "spp")
dt_expr <- merge(dt_expr_mer, iso, by = "spp")
species <- dt_expr$spp

# Species names to index rows
rownames(dt_expr) <- dt_expr$spp

# Reduce tree to species in use
treetan_red <- keep.tip(treetan, species)

#Phylo regression
genes <- colnames(dt_expr)[3:12]
# Create data frame with 2 cols
dt_fr <- data.frame(matrix(ncol = 2))
colnames(dt_fr) <- c("gene_ID", "pval")
for (gen in genes) {
  df_ind <- dt_expr[c("spp", "phenotype", "d15N", "d13C",  gen)]
  df_ind["std_counts"] <- (df_ind[[gen]] - mean(df_ind[[gen]]))/sd(df_ind[[gen]])
  df_ind["std_d15N"] <- (df_ind[["d15N"]] - mean(df_ind[["d15N"]]))/sd(df_ind[["d15N"]])
  df_ind["std_d13C"] <- (df_ind[["d13C"]] - mean(df_ind[["d13C"]]))/sd(df_ind[["d13C"]])
  
  try(
    {
      pval <- summary(phyloglm(phenotype~std_counts + std_d15N + std_d13C, btol = 500, data = df_ind, phy = treetan_red, method = "logistic_IG10"))$coefficients["std_counts","p.value"]
    }
    , silent = T
  )
  
  vec_coef <- c(gen, pval)
  dt_fr <- rbind(dt_fr, vec_coef)
}  
dt_fr

name <- "pair-bond"
write.table(dt_fr, paste0(path_out, name, ".txt"),
            quote = F)



# 3. SEXES SEPARATED
# Median by species
dt_med <- exp_dt %>%
  group_by(spp, sex,tissue) %>%  # Group by species and tissue
  #group_by(spp,tissue) %>%  # Group by species and tissue
  summarise(across(2:11, \(x) median(x, na.rm = TRUE)), .groups = "drop")  # Compute median for gene count columns

# Filter brain
dt_med_br <- dt_med[dt_med$tissue == "brain",]
dt_med_br <- dt_med_br[dt_med_br$sex == "F",]
dt_med_br <- as.data.frame(dt_med_br)


# Introduce phenotype
# Merge
dt_expr_mer <- merge(dt_med_br, pb_info, by = "spp")
dt_expr <- merge(dt_expr_mer, iso, by = "spp")
species <- dt_expr$spp

# Species names to index rows
rownames(dt_expr) <- dt_expr$spp

# Reduce tree to species in use
treetan_red <- keep.tip(treetan, species)

#Phylo regression
genes <- colnames(dt_expr)[4:13]


# Create data frame with 2 cols
dt_fr <- data.frame(matrix(ncol = 2))
colnames(dt_fr) <- c("gene_ID", "pval")

for (gen in genes) {
  df_ind <- dt_expr[c("spp", "phenotype", "d15N", "d13C",  gen)]
  df_ind["std_counts"] <- (df_ind[[gen]] - mean(df_ind[[gen]]))/sd(df_ind[[gen]])
  df_ind["std_d15N"] <- (df_ind[["d15N"]] - mean(df_ind[["d15N"]]))/sd(df_ind[["d15N"]])
  df_ind["std_d13C"] <- (df_ind[["d13C"]] - mean(df_ind[["d13C"]]))/sd(df_ind[["d13C"]])
  
  try(
    {
      pval <- summary(phyloglm(phenotype~std_counts + std_d15N + std_d13C, btol = 500, data = df_ind, phy = treetan_red, method = "logistic_IG10"))$coefficients["std_counts","p.value"]
    }
    , silent = T
  )
  
  vec_coef <- c(gen, pval)
  dt_fr <- rbind(dt_fr, vec_coef)
}  
dt_fr

name <- "pair-bonding_male"
write.table(dt_fr, paste0(path_out, name, ".txt"),
            quote = F)
