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

  # Subset data for the specific gene
  df_gene <- subset(data, GeneID == gene_name)
  
  # Ensure species factor levels match the phylogeny
  df_gene$spp <- factor(df_gene$spp, levels = rownames(inv.phylo$Ainv))
  
  # Define the formula dynamically based on the phenotype
  formula_txt <- paste0(
    "log_expression ~ ",
    phenotype,
    " * sex + d15N + d13C"
  )
  
  # Convert to formula object
  formula_obj <- as.formula(formula_txt)
  
  # Run the MCMCglmm model
  model <- MCMCglmm(
    formula_obj,
    random = ~ spp,
    family = "gaussian",
    ginverse = list(spp = inv.phylo$Ainv),
    data = as.data.frame(df_gene),
    nitt = 1000000,
    burnin = 100000,
    thin = 100,
    verbose = FALSE
  )
  
  return(model)
}

# Function to extract results
extract_results <- extract_results <- function(model, gene, phenotype,
                            pheno_term, sex_term, inter_term) {
  
  sol <- model$Sol
  
  # -----------------------------
  # Convergence diagnostics (ONLY pheno, sex, interaction)
  # -----------------------------
  get_geweke <- function(term) {
    
    z_val <- geweke.diag(as.mcmc(sol[, term]))$z
    conv_ok <- abs(z_val) < 1.95
    
    return(c(z_val, conv_ok))
  }

  pheno_conv <- get_geweke(pheno_term)
  sex_conv   <- get_geweke(sex_term)
  int_conv   <- get_geweke(inter_term)
  d15N_conv  <- get_geweke("d15N")
  d13C_conv  <- get_geweke("d13C")
  
  # -----------------------------
  # Extract posterior stats + HPD
  # -----------------------------
  extract_stats <- function(term) {
    
    if (!(term %in% colnames(sol))) return(c(NA, NA, NA, NA))
    
    # Extract posterior samples for the term
    post <- sol[, term]
    mean_post <- mean(post)
    
    # Two-tailed Bayesian pMCMC
    pMCMC <- 2 * min(mean(post > 0), mean(post < 0))
    
    # 95% HPD interval
    hpd <- HPDinterval(as.mcmc(post), prob = 0.95)
    CI_low  <- hpd[1]
    CI_high <- hpd[2]
    
    return(c(mean_post, pMCMC, CI_low, CI_high))
  }
  
  # Define term names based on phenotype
  pheno_stats <- extract_stats(pheno_term)
  sex_stats   <- extract_stats(sex_term)
  int_stats   <- extract_stats(inter_term)
  d15N_stats  <- extract_stats("d15N")
  d13C_stats  <- extract_stats("d13C")
  
  # Compile results into a data frame
  return(data.frame(
    gene = gene,

    # Convergence diagnostics
    geweke_pheno = pheno_conv[1],
    conv_pheno = pheno_conv[2],
    
    geweke_sex = sex_conv[1],
    conv_sex = sex_conv[2],
    
    geweke_interaction = int_conv[1],
    conv_interaction = int_conv[2],
    
    geweke_d15N = d15N_conv[1],
    conv_d15N = d15N_conv[2],
    
    geweke_d13C = d13C_conv[1],
    conv_d13C = d13C_conv[2],
    
    # Phenotype effect
    mean_pheno = pheno_stats[1],
    p_pheno = pheno_stats[2],
    HPD_low_pheno = pheno_stats[3],
    HPD_high_pheno = pheno_stats[4],
    
    # Sex effect
    mean_sex = sex_stats[1],
    p_sex = sex_stats[2],
    HPD_low_sex = sex_stats[3],
    HPD_high_sex = sex_stats[4],
    
    # Interaction
    mean_interaction = int_stats[1],
    p_interaction = int_stats[2],
    HPD_low_interaction = int_stats[3],
    HPD_high_interaction = int_stats[4],
    
    # d15N
    mean_d15N = d15N_stats[1],
    p_d15N = d15N_stats[2],
    HPD_low_d15N = d15N_stats[3],
    HPD_high_d15N = d15N_stats[4],
    
    # d13C
    mean_d13C = d13C_stats[1],
    p_d13C = d13C_stats[2],
    HPD_low_d13C = d13C_stats[3],
    HPD_high_d13C = d13C_stats[4]
  ))
}

# ----------------------------------
# Run models for both phenotypes
# ----------------------------------
phenotypes <- c("Sex_caregiver", "Pair_bonding")

all_genes <- unique(df_long$GeneID)

for (phenotype in phenotypes) {
  
  cat("\nRunning phenotype:", phenotype, "\n")
  
  # Define output file path
  output_file <- paste0("Data/16.PGLMM/out/", phenotype, "_results.txt")

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
    
    cat("   Gene:", g, "\n")
    
    # Run model with error handling
    model <- try(run_pglmm(df_pheno, g, phenotype, inv.phylo))
    
    if (inherits(model, "try-error")) {
      cat("   Failed:", g, "\n")
      next
    }
    
    # Extract results
    res <- extract_results(model, g, phenotype,
                       pheno_term, sex_term, inter_term)
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
