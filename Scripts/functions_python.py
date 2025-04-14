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


# Function 3: Yn00. dNdS ratio for pairwise comparisons
def yn00_an(file_seq, path_ctl_file):
    from Bio.Phylo.PAML import yn00 # Biopython version 1.80
    
    #Create Yn object
    yn_obj = yn00.Yn00()
    
    #Read control file
    yn_obj.read_ctl_file(path_ctl_file)
    
    #Read sequences
    yn_obj.alignment = file_seq
    
    #Set output and wd
    yn_obj.out_file = "./results.txt" # Path output
    yn_obj.working_dir = "./"
    
    #Run
    yn_obj.run()


# Function 4. Parse JSON output Hyphy FEL and filter
def hyphy_parse_fel(in_path, out_name, gene_name):
    import phyphy
    import argparse
    import pandas as pd

    p = phyphy.Extractor(in_path)
    p.extract_csv(out_name + ".csv")

    # Filter data and write output
    ## Read CSV
    dfp = pd.read_csv(out_name + ".csv")
    ## Add gene name column and change the order of columns
    dfp['gene'] = gene_name
    ## Filter by p-value
    dfp_filt = dfp[dfp['p-value'] < 0.05]
    ## Filter only positive selection
    dfp_filt_pos = dfp_filt[dfp_filt['alpha'] < dfp_filt['beta']]
    dfp_filt_neg = dfp_filt[dfp_filt['alpha'] > dfp_filt['beta']]
    ## write files
    dfp_filt_pos.to_csv(out_name + "_filt_positive.csv", index = False)
    dfp_filt_neg.to_csv(out_name + "_filt_negative.csv", index = False)


# Function 5: Dataset for allele variant    
def allele_site(alignment_file, out_dir, gen):
    from Bio import AlignIO
    import pandas as pd
    import os

    # Read alignment file
    aln = AlignIO.read(alignment_file, "fasta")

    aln_length = aln.get_alignment_length()

    positions_with_change = []

    # Loop through alignment positions
    for position in range(aln_length):
        # Gather unique amino acids at the position
        unique_amino_acids = set(record.seq[position] for record in aln)

        # Skip sites with no variation
        if len(unique_amino_acids) <= 1:
            continue

        positions_with_change.append(position + 1)  # Save the position

        site = position + 1  # Add 1 to position to match alignment numbering
        sites = []

        for record in aln:
            species_name = record.id
            amino_acid = record.seq[position]

            # Skip species named "Orenil"
            if "Orenil" in species_name:
                continue

            # Retain species with amino acid "X"
            allele = "X" if amino_acid == "X" else amino_acid

            sites.append([species_name, allele])

        # Convert the data to a DataFrame
        df = pd.DataFrame(sites, columns=["Species", "Allele"])

        # Separate "X" from allele counts for determining most common alleles
        df_with_x = df[df["Allele"] == "X"]
        df_no_x = df[df["Allele"] != "X"]
        
        # Filter out species with a third allele (not counting "X")
        allele_counts = df_no_x["Allele"].value_counts()
        if len(allele_counts) < 2:
            continue  # Skip if fewer than two alleles remain

        top_two_alleles = allele_counts.index[:2]  # Get the two most common alleles
        df_no_x = df_no_x[df_no_x["Allele"].isin(top_two_alleles)]

        # Combine back rows with "X"
        df = pd.concat([df_no_x, df_with_x]).drop_duplicates()

        # Map the alleles to 0 (most common), 1 (second most common), "X" (unchanged)
        most_common_allele, second_common_allele = top_two_alleles[:2]

        def map_allele(allele):
            if allele == most_common_allele:
                return 0
            elif allele == second_common_allele:
                return 1
            return "X"

        df["Allele"] = df["Allele"].apply(map_allele)

        # Sort the DataFrame by the Species column
        df = df.sort_values(by="Species")

        # Create output file path
        out_file = os.path.join(out_dir, f"{gen}_site_{site}.txt")

        # Write DataFrame to file
        df.to_csv(out_file, index=False, sep="\t")

    # Write list of positions with change
    with open(os.path.join(out_dir, f"{gen}_all_positions_with_change.txt"), "w") as f:
        f.write(" ".join(map(str, positions_with_change)) + "\n")


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
## Function3. Yn00
parser_yn00 = subparsers.add_parser('yn00_an', help='dNdS pairwise')
parser_yn00.add_argument('file_seq', type=str, help="Provide path of the sequence file")
parser_yn00.add_argument('path_ctl_file', type=str, help="Provide path of the control file")
## Function4. Hyphy parse FEL
parser_hyphy_fel_parse = subparsers.add_parser('hyphy_parse_fel', help='parse hyphy_FEL output')
parser_hyphy_fel_parse.add_argument('in_path', type=str, help="json file path")
parser_hyphy_fel_parse.add_argument('out_path', type=str, help="name output")
parser_hyphy_fel_parse.add_argument('gene_name', type=str, help="gene name")
## Function5. Allele site
parser_allele_site = subparsers.add_parser('allele_site', help='spot_snps_in_alignment')
parser_allele_site.add_argument('alignment_file', type=str, help="alignment file path")
parser_allele_site.add_argument('out_dir', type=str, help="positions file directory")
parser_allele_site.add_argument('gen', type=str, help="gene")



# Parse the arguments
args = parser.parse_args()

# Determine which function to execute based on the sub-command
if args.subcommand == 'sum_0cov':
    sum_0cov(args.genes_file, args.suffix)
elif args.subcommand == 'best_cov_file':
    best_cov_file(args.genes_file, args.species_file, args.file_0cov1, args.file_0cov2)
elif args.subcommand == 'yn00_an':
    yn00_an(args.file_seq, args.path_ctl_file)
elif args.subcommand == 'hyphy_parse_fel':
    hyphy_parse_fel(args.in_path, args.out_path, args.gene_name)
elif args.subcommand == 'allele_site':
    allele_site(args.alignment_file, args.out_dir, args.gen)