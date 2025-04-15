#!/bin/bash

# This script runs BayesTraits Discrete method on the datasets obtained in script 13.2. It tests for correlated evolution between
# the amino acid variant and the phenotype data (both discrete).
# This script shows the general script to run BayesTraits.
# For space purposes, the output files are not shown. Only the likelihood values and log Bayes factors are shown.

folder_BT="Local/path/to/BayesTraits" # Path where BayesTraitsV4 is installed
folder_datasets="${folder_BT}/datasets" # Where datasets are stored. Should be inside $folder_BT
folder_trees="${folder_BT}/trees" # Where trees are stored. Basename of the nexus file should be the same as the dataset file. Should be inside $folder_BT
folder_Lh="${folder_BT}/Lh" # Where likelihood files are stored. Should be inside $folder_BT
folder_parameters_files="${folder_BT}/param_files" # Where parameter files are stored. Should be inside $folder_BT


cd $folder_datasets # It needs to be run in the folder where BayesTraits is installed

for file in *.txt ; do # Do it for all files in folder
    
    # Get the name of the file without the extension
    name=$( echo $file | sed "s/.txt//g" ) 

    # Create two copies of the file. One for dependent and one for independent models
    cp $file "$name"_ind.txt ; mv $file "$name"_dep.txt

    cd $folder_BT

    # Run BayesTraits
    date

    # Run both models
    ./BayesTraitsV4 "$folder_trees"/"$name".nex $folder_datasets/"$name"_dep.txt < $folder_parameters_files/dep_it7.txt > $folder_datasets/"$name"_dep.txt
    ./BayesTraitsV4 "$folder_trees"/"$name".nex $folder_datasets/"$name"_ind.txt < $folder_parameters_files/ind_it7.txt > $folder_datasets/"$name"_ind.txt
    
    # Stones file to obtain marginal likelihood values
    indep=$( grep 'marginal' $folder_datasets/"$name"_dep*Stones.txt | cut -f2 ) ; dep=$( grep 'marginal' $folder_datasets/"$name"_ind*Stones.txt | cut -f2 )
    
    # Output the likelihood values into Lh folder
    echo -e "$name\tLh_indep: $indep\tLh_dep: $dep" > "$folder_Lh"/"$name"_Lh.txt
            
    date

    cd $folder_datasets # Retrun to the datasets folder

done

exit