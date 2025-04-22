# Molecular evolution of the nonapeptide signalling in Lake Tanganyika cichlid radiation
Pol Sorigue, Walter Salzburger, Rui Oliveira (2025)

---

## **Repository Structure**

EDIT

```
├─ Data/                     # Input and output data files
│   ├── 01.Map_assemblies/    # Mapping assemblies
│   ├── 02.Local_Database_Blast/ # BLAST-related data
│   ├── 09.Nucleotide_Diversity/ # Nucleotide diversity analysis
│   ├── 13.Variable_Sites/    # Variable site datasets
│   ├── 14.BayesTraits/       # Phylogenetic tree data
├── Scripts/                  # Bash and R scripts for analysis
│   ├── 01.Map_assemblies.sh  # Maps raw reads to reference genomes
│   ├── 03.Local_Database_Blast.sh # Creates local BLAST databases
│   ├── 09.2.NucDiv_by_domain.R # Calculates nucleotide diversity by domain
│   ├── 14.1.Create_Trees_BT.R # Prunes species trees for BayesTraits
│   ├── functions_bash.sh     # Helper functions for bash-based operations
├── README.md                 # Project documentation
```

---

## **Scripts**
The versions refer to the versions used during the analysis.
### **`01.Map_assemblies.sh`**: 
Maps raw reads to reference genomes and transcriptomes, generating BAM files for each sample.
- BWA-MEM version 0.7.18
- STAR version 2.7.10a
- IGVTools version 2.5.3
- SAMtools version 1.21
### **`02.Teleost_tree_receptors`**
Creates a phylogenetic tree with teleost nonapeptide receptor sequences.
- MAFFT version 7.526
- IQTREE version 2.0.3
### **`03.Local_Database_Blast.sh`**
Creates local BLAST databases from PacBio assemblies and runs BLAST searches for specific genes.
- makeblastdb version 2.16.0+
- blastn version 2.16.0+
### **`04.Download_CDS_and_regions_file.sh`**
Downloads the reference transcripts from NCBI and creates a regions file that includes the locations of the CDS regions.
- datasets version: 14.6.0
### **`05.1.Coverage.sh`**
Calculates the coverage depth for each position of CDS regions for genomes and transcriptomes using the BAM files.
- SAMtools depth version 1.21
### **`05.2.Best_Coverage_File.sh`**
This script creates a file indicating the source with the best coverage per gene and species (either genome of transcriptome BAM).
- Python3
### **`06.Gene_checkup.sh`**
### **`07.1.Extract_consensus_genomes.sh`**
Extracts the consensus sequences for each gene and individual from the BAM files and merges them to create a consensus for each species obtained from genomic sequences.
- MAFFT v7.526
- Seqtk version 1.4-r122
- BCFtools version 1.9
- SAMtools version 1.21
### **`07.2.Extract_consensus_genomes.sh`**
Extracts the consensus sequences for each gene and individual from the BAM files and merges them to create a consensus for each species obtained from transcriptomes.
- Seqtk version 1.4-r122
- BCFtools Version: 1.9
- SAMtools Version: 1.21
### **`08.Alignments.sh`**
Aligns the species consensus for each gene with the reference genome. Using the coverage information, it selects the best coverage consensus (either genomes or transcriptomes) and, if necessary, replaces genomic consensus by transcriptomic consensus.
- Seqkit v2.8.2
- MAFFT v7.526
- pal2nal.pl v14
### **`09.1.Nucleotide_Diversity.R`**
Calculates nucleotide diversity (pi) for each transcript in the dataset.
- R version 4.4.0
- R packages: Ape (version 5.8), Pegas (version 1.3)
### **`09.2.NucDiv_by_domain.R`**
Calculates nucleotide diversity for specific domains of nonapeptides and receptors.
- R version 4.4.0
- R packages: Ape (version 5.8), Pegas (version 1.3)
### **`10.Gene_Trees.sh`**
generates gene trees from the alignments, using both nucleotide and amino acid sequences.
- IQTREE version 2.0.3
### **`11.1.dS_ratios.sh`**
Calculates the dS ratios of each gene and species separately, using the genes of the reference genome using Yang and Nielsen method.
- MAFFT v7.526
- Biopython version 1.80
### **`11.2.dS_VTR1Ab.sh`**
Calculates the dS ratio of the two exons of VTR1Ab gene separately.
- MAFFT v7.526
- Biopython version 1.80
### **`12.FEL_pos_selection.sh`**
Performs positive selection analysis on the aligned sequences, using Hyphy FEL method.
- Hyphy version 2.5
- Python 3
### **`13.1.Variable_Sites.sh`**
Returns a dataset for each variable amino acid site in each transcript. The dataset contains a column with species name.
- Python 3
### **`13.2.Datasets_Variable_Sites.R`**
Creates a dataset with three columns: species name, allele, and phenotype. This dataset will be used in the script 14.2 for corelation analysis.
- R version 4.4.0
### **`14.1.Create_Trees_BT.R`**
Prunes the species tree to include only species present in the dataset and generates a tree for each variable site analyzed.
- R version 4.4.0
- R package Ape (version 5.8)
### **`14.2.BayesTraits_run.sh`**
Runs BayesTraits Discrete method on the datasets obtained in script 13.2. It tests for correlated evolution between the amino acid variant and the phenotype data.
- BayesTraits V4
### **`14.3.Convergence_BT_run.R`**
Tests convergence of the BayesTraits MCMC runs using the Geweke diagnostic. It creates a dataframe with the z-scores of the Geweke diagnostic for all sites.
- R version 4.4.0.
- R package Coda version 0.19-4.1
### **`14.4.Randomize_dataset.R`**
Randomizes the input dataset for BayesTraits 200 times. It creates a new dataset with the same species and alleles but with randomized phenotypes.
- R version 4.4.0
### **`14.5.BT_reverse-jump.sh`**
Runs BayesTraits Discrete using Reversible-Jump method.
- BayesTraits V4
### **`15.Phylogenetic_regression.R`**
Runs a phylogenetic regression with the read counts of the brain transcriptomes.
- R version 4.4.0.
- R packages: Ape (version 5.8), Phylolm (version 2.6.5)
### **`functions_python.py`**
This file contains a repertoire of functions in Python language used during the analysis. The functions are accessed and called from the scripts.
### **`functions_bash.py`**
This file contains a repertoire of functions in bash language used during the analysis. The functions are accessed and called from the scripts.

---


