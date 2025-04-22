# This script checks the convergence of the MCMC runs using the Geweke diagnostic. It creates a dataframe with the z-scores of the Geweke diagnostic for all sites
library(coda) # version 0.19-4.1

folder_datasets="Local/path/to/BayesTraits/datasets/" # Path to output of script 14.2
folder_out="Data/14.BayesTraits/out/" # Path to output of this script

# Assuming the files are named as they are output in script 14.2
files <- as.vector(scan("Data/14.BayesTraits/in/list_sites.txt/", what = "character", sep = "\n"))


df <- data.frame()

for (sit in sites) {
  
  models <- c("ind", "dep")
  
  for (mod in models) {
    
    file <- paste0(folder_datasets, sit, "_", mod, ".txt.Log.txt") # Read the Log file
    
    # Read file
    lines <- readLines(file)
    
    # Locate the start of the table
    start_line <- grep("Iteration\t", lines)  # Line after the "Iteration" header
    
    # Extract header and data directly from the lines
    header <- lines[start_line]
    
    # Read the table starting from the next line after the header
    table_data <- read.delim(text = paste(lines[(start_line + 1):length(lines)], collapse = "\n"),
                             sep = "\t", header = FALSE, stringsAsFactors = FALSE)
    
    # Set the column names
    colnames(table_data) <- unlist(strsplit(header, "\t"))
    
    # Extract Lh values column
    Lh_val <- as.numeric(table_data$Lh)
    
    # Convert the likelihoods into an MCMC object
    Lh_mcmc <- mcmc(Lh_val)
    
    # Apply Geweke diagnostic
    geweke_diag <- geweke.diag(Lh_mcmc)
    
    # Add row in dataframe
    vec <- c(paste(sit, mod, sep= "_"), geweke_diag$z)
    df <- rbind(vec, df)
    
  }
  
  colnames(df) <- c("site_mod", "z-score_geweke")
  df$`z-score_geweke` <- as.numeric(df$`z-score_geweke`)
  
}

write.table(df, file = paste0(folder_out, "convergence_values.txt"),
            quote = F,
            sep = "\t",
            col.names = T,
            row.names = F)

# The output (Data/14.BayesTraits/out/Lh-convergence.txt) is a summary with the Lh of all runs and all Geweke values. 
# The runs were done with different iterations until they reached convergence