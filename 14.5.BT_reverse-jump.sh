#!/bin/bash
# This script runs BayesTraits Discrete using Reversible-Jump method. Here I show the corelation between phenotypes, 
# but this was also run for the site that showed corelated evolution with pair-bonding.
# The interpretation of the results of the Reversible-jump needs to be done manually (shown in output files)

# Files
dataset="Data/14.BayesTraits/in/reverse-jump/dataset_pb_sc.txt" # Dataset file to be used.
tree="Data/14.BayesTraits/in/reverse-jump/tree_pb_sc.nex" # Tree file to be used.

# Paths
folder_BT="Local/path/to/BayesTraits" # Path where BayesTraitsV4 is installed
# Path reverse-jump
path_rj="${folder_BT}/folder_out/" # Path to run the reverse-jump


# All files and folders should be inside the folder where BayesTraits is installed
cp $dataset $path_rj/dataset.txt ; cp $tree $path_rj/tree.nex # Copy dataset and tree to the folder where BayesTraits is installed

# Follwing developer's advice, we run 3 times for dependent model and 3 times for independent model.
# Make folders for all runs
cd "$path_rj" ; mkdir ML Dep1 Dep2 Dep3 Ind1 Ind2 Ind3
        
# Copy dataset in all folders
for dir in */ ; do cp dataset.txt "$dir"/dataset.txt ; done

# Run ML model
cd ML/ ; cp dataset.txt dataset_dep.txt ; mv dataset.txt dataset_ind.txt

cd $folder_BT        
tree="$path_rj"/"$sit"/tree.nex
dataset_dep="$path_rj"/ML/dataset_dep.txt ; dataset_ind="$path_rj"/ML/dataset_ind.txt

./BayesTraitsV4 $tree $dataset_dep < ./param_files/param_rj_ML_dep.txt
./BayesTraitsV4 $tree $dataset_ind < ./param_files/param_rj_ML_ind.txt

# Run MCMC models
for i in {1..3} ; do

    dataset_dep="$path_rj"/Dep"$i"/dataset.txt ; dataset_ind="$path_rj"/Ind"$i"/dataset.txt
    
    # Reverse-jump with exponential prior
    ./BayesTraitsV4 $tree $dataset_dep < ./param_files/param_rj_exp_dep.txt # Dependent
    ./BayesTraitsV4 $tree $dataset_ind < ./param_files/param_rj_exp_ind.txt # Independent
    # These commands create output files in the folder where the initial dataset was located. 

done

# The output files are used to fill in the template for interpretation of the results

exit