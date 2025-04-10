# This script serves as a repository of functions that are used in the pipeline.
# It is not intended to be run directly, but rather to be imported into other scripts.
# To run a function of this script, call this script by 'python3 functions_python.py <function_name> <args>'

import argparse

# Function 1: Sum 0-covered positions. Input need file format of samtools depth
def sum_0cov(genes_file, suffix):
    import pandas as pd
    
    genes = open(genes_file , 'r')
    gen_list = genes.readline().split(" ") # Make sure the delimiter is correct
    all_df = pd.DataFrame()
    
    for gen in gen_list:
    
        try:
            # Read CSV
            df = pd.read_csv("cov_" + gen + ".txt", sep="\t")
            
            # Drop columns 'chr' and 'pos'
            dff = df.loc[:, ~df.columns.isin(['chr', 'pos'])]
            
            # Count number of zeros
            df0 = dff.eq(0).sum().to_frame()
            
            # Rename column
            df0.columns = [gen]
            
            # Add col to all_df
            all_df[gen] = df0[gen]
            
        except FileNotFoundError:
        
            # Fill the column with NaN if the file is missing
            all_df[gen] = float('nan')
            
    all_df.to_csv("num0cov_" + suffix + ".csv", header = True, index = True)
            

# Function 2: Create file with best coverage per gen and species (either gn or tr)
def best_cov_file(genes_file, species_file, file_0cov1, file_0cov2):
    import pandas as pd
    
    species = open( species_file , 'r') # Transcriptome species
    spp_list = species.readline().split(" ")
    genes = open( genes_file , 'r')
    gen_list = genes.readline().split(" ") # Make sure the delimiter is correct
    
    #Read CSV 'number of 0-cov positions'
    gn_df = pd.read_csv(file_0cov1, index_col = 0, header = 0) # Header and index are read
    tr_df = pd.read_csv(file_0cov2, index_col = 0, header = 0)
    
    #Create new data frame
    df_bcov = pd.DataFrame()
    df_bcov.index = spp_list
    df_bcov[gen_list] = ""
    
    #Fill dataframe with condition
    for gen in gen_list:
        
        for spp in spp_list:
        
            spp_lc = spp[0].lower() + spp[1:]
            
            spp_val_gn = gn_df[gen].filter(like = spp_lc, axis = 0).max() # Takes the maximum of 0-cov positions (more stringent)
            spp_val_tr = tr_df[gen].filter(like = spp, axis = 0).max()
            
            if spp_val_gn > spp_val_tr:
                spp_val = 'tr'
            if spp_val_gn <= spp_val_tr:
                spp_val = 'gn'
                
            #Write result
            df_bcov.loc[spp, gen] = spp_val
    
    df_bcov.to_csv( "best_coverage.csv", header = True, index = True)





# Top level parser
parser = argparse.ArgumentParser(description='Parse function arguments')

# Create subparsers for each function
subparsers = parser.add_subparsers(title='subcommands', description='valid subcommands', help='additional help', dest='subcommand')

## Function1. Sum 0covered positions
parser_sum_0cov = subparsers.add_parser('sum_0cov', help='execute function')
parser_sum_0cov.add_argument('genes_file', type=str, help="Provide gene file path")
parser_sum_0cov.add_argument('suffix', type=str, help="Provide suffix for the output file name")
## Function2. Best_coverage file
parser_best_cov_file = subparsers.add_parser('best_cov_file', help='execute function')
parser_best_cov_file.add_argument('genes_file', type=str, help="Provide gene file path")
parser_best_cov_file.add_argument('species_file', type=str, help="Provide species file path")
parser_best_cov_file.add_argument('file_0cov1', type=str, help="Provide file 1 to compare")
parser_best_cov_file.add_argument('file_0cov2', type=str, help="Provide file 2 to compare")




# Parse the arguments
args = parser.parse_args()

# Determine which function to execute based on the sub-command
if args.subcommand == 'sum_0cov':
    sum_0cov(args.genes_file, args.suffix)
elif args.subcommand == 'best_cov_file':
    best_cov_file(args.genes_file, args.species_file, args.file_0cov1, args.file_0cov2)