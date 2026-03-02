# PGLMM analysis for gene expression data. This script runs Phylogenetic Generalized Linear Mixed Models (PGLMM) using MCMCglmm for each gene and phenotype, extracts pMCMC values, and saves results to a file.

library(tidyverse) # version 2.0.0
library(ape) # version 5.8
library(MCMCglmm) # version 2.36
library(coda) # version 0.19-4.1

# Read input data
input_file <- "Data/16.PGLMM/inp_exp_scaled.txt"
df <- read.table(input_file, header = TRUE, sep = "\t")

# Read Phylogeny in Newick format
tree_tan <- read.tree("Data/16.PGLMM/tree_tan.nwk")

# Define output file path
output_file <- paste0("Data/16.PGLMM/out/", phenotype, "_results.txt")
  
# Convert factors in initial data frame
df$spp <- factor(df$spp)
df$sex <- factor(df$sex)
df$Sex_caregiver <- factor(df$Sex_caregiver)
df$Pair_bonding <- factor(df$Pair_bonding)

# Convert to long format
gene_cols <- c("OT","VT","OTRa","OTRb",
               "VTR1Aa","VTR1Ab",
               "VTR2Aa","VTR2Ab",
               "VTR2Ba","VTR2Bb")
df_long <- df %>%
  pivot_longer(
    cols = all_of(gene_cols),
    names_to = "GeneID",
    values_to = "log_expression"
  )

# Function to run PGLMM for a given gene and phenotype
run_pglmm <- function(data, gene_name, phenotype, inv.phylo) {

  df_gene <- subset(data, GeneID == gene_name)
  
  df_gene$spp <- factor(df_gene$spp, levels = rownames(inv.phylo$Ainv))
  
  formula_txt <- paste0(
    "log_expression ~ ",
    phenotype,
    " * sex + d15N + d13C"
  )
  
  formula_obj <- as.formula(formula_txt)
  
  model <- MCMCglmm(
    formula_obj,
    random = ~ spp,
    family = "gaussian",
    ginverse = list(spp = inv.phylo$Ainv),
    data = as.data.frame(df_gene),
    nitt = 300000,
    burnin = 10000,
    thin = 100,
    verbose = FALSE
  )
  
  return(model)
}

# Function to extract results
extract_results <- function(model, gene, phenotype) {
  
  sol <- model$Sol
  
  # -----------------------------
  # Convergence diagnostics
  # -----------------------------
  geweke_vals <- geweke.diag(sol)$z # Geweke diagnostic z-scores
  convergence_ok <- all(abs(geweke_vals) < 3) # Return convergence diagnostic result
  
  # Effective sample size (minimum across fixed effects)
  eff_samp_min <- min(effectiveSize(sol)) # Minimum effective sample size across all parameters
  
  # -----------------------------
  # Extract pMCMC and posterior mean
  # -----------------------------
  extract_stats <- function(term) {
    if (!(term %in% colnames(sol))) return(c(NA, NA))
    
    post <- sol[, term]
    mean_post <- mean(post)                  # Posterior mean
    pMCMC <- 2 * min(mean(post > 0), mean(post < 0))  # Two-tailed pMCMC
    
    return(c(mean_post, pMCMC))
  }
  # Define term names dynamically
  pheno_stats <- extract_stats(pheno_term)
  sex_stats   <- extract_stats(sex_term)
  int_stats   <- extract_stats(inter_term)
  d15N_stats  <- extract_stats("d15N")
  d13C_stats  <- extract_stats("d13C")
  
  return(data.frame(
    phenotype = phenotype,
    gene = gene,
    
    eff_samp_min = eff_samp_min,
    convergence_ok = convergence_ok,
    
    mean_pheno = pheno_stats[1],
    p_pheno = pheno_stats[2],
    
    mean_sex = sex_stats[1],
    p_sex = sex_stats[2],
    
    mean_interaction = int_stats[1],
    p_interaction = int_stats[2],
    
    mean_d15N = d15N_stats[1],
    p_d15N = d15N_stats[2],
    
    mean_d13C = d13C_stats[1],
    p_d13C = d13C_stats[2]
  ))
}

# ----------------------------------
# Run models for both phenotypes
# ----------------------------------
phenotypes <- c("Sex_caregiver", "Pair_bonding")

all_genes <- unique(df_long$GeneID)

for (phenotype in phenotypes) {
  #phenotype <- "Pair_bonding"
  cat("\nRunning phenotype:", phenotype, "\n")
  
  # Remove NAs for this phenotype
  df_pheno <- df_long[!is.na(df_long[[phenotype]]), ]

  # Prune tree to only species with the phenotype
  tree <- keep.tip(tree_tan, as.vector(unique(df_pheno$spp)))
  inv.phylo <- inverseA(tree, nodes = "TIPS", scale = FALSE)
  
  # Ensure factor levels are set
  df_pheno[[phenotype]] <- factor(df_pheno[[phenotype]])
  
  # Define term names dynamically
  pheno_levels <- levels(df_pheno[[phenotype]])
  sex_levels <- levels(df_pheno$sex)
  
  pheno_term <- paste0(phenotype, pheno_levels[2])
  sex_term <- paste0("sex", sex_levels[2])
  inter_term <- paste0(pheno_term, ":", sex_term)
  
  results_list <- list()
  
  # -------------------------------
  # Loop over genes
  # -------------------------------
  
  for (g in all_genes) {
    #g <- "OT"
    cat("   Gene:", g, "\n")
    
    model <- try(run_pglmm(df_pheno, g, phenotype, inv.phylo))
    
    if (inherits(model, "try-error")) {
      cat("   Failed:", g, "\n")
      next
    }
    
    res <- extract_results(model, g, phenotype)
    results_list[[g]] <- res
  }
  
  # Combine gene results
  results_df <- bind_rows(results_list)
  
  # ----------------------------------
  # Save phenotype-specific results
  # ----------------------------------

  write.table(
    results_df,
    file = output_file,
    sep = "\t",
    row.names = FALSE,
    col.names = TRUE,
    quote = FALSE
  )
  
  cat("Saved:", output_file, "\n")
}
