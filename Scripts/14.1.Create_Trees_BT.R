# This scripts prunes the species tree so it only keeps the species that are present in the dataset (output of script 13.2).
# It produces a species tree for each site.


library(ape)

# Paths
path_datasets <- "Data/13.Variable_Sites/out/datasets_for_BT_pair_bonding/"
path_trees <- "Data/14.BayesTraits/in/trees/"
# Tree species
tree_phyl <- "Data/14.BayesTraits/in/spp_tree.trees/"
tree_tan <- read.nexus(tree_phyl)
# Sites
Sites <- as.vector(scan("Data/14.BayesTraits/in/list_sites.txt/", what = "character", sep = "\n"))


for (sit in Sites) {

    # Read the dataset
    dataset <- read.csv(file = paste0(path_datasets, sit, ".txt"), header = F, sep = "\t")

    # Take the species names of the dataset
    spp_list <- dataset[,1]

    # Prune phyl tree and save it
    new_tan_tree <- keep.tip(tree_tan, spp_list) 
    write.nexus(new_tan_tree, file = paste0(path_trees, sit, ".nex"), translate = T) # Important to translate to fit BT requirements
  
}