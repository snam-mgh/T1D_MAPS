# Calculating T1D Multi-Ancestry Polygenic Score (T1D MAPS)

## Overview

**Type 1 Diabetes Multi-Ancestry Polygenic Score (T1D MAPS)** is a novel polygenic risk score designed to identify individuals at high risk for type 1 diabetes (T1D) using a single, universal threshold. T1D MAPS was developed to improve the accuracy of T1D risk prediction across multiple ancestries. This tool enables consistent risk stratification regardless of genetic ancestry.

This repository provides the necessary resources to calculate T1D MAPS.

## Prerequisites 

The following tools are required to calculate T1D MAPS:

- [PLINK v1.9](https://www.cog-genomics.org/plink/)
- [R ≥ 4.4](https://www.r-project.org/)
  - data.table_1.16.4 is required to run `T1D_MAPS.R`

## Usage

1. Perform HLA imputation from input PLINK files using the Michigan Imputation Server

  - Note that the Michigan Imputation Server has a limit on the number of individuals per cohort (currently 30,000). If necessary, split the PLINK .fam file into subsets of 30,000 each. 

    `split -l 30000 BFILE_PREFIX.fam BFILE_PREFIX.fam.split`

    Then, for each subset, generate a new PLINK fileset. 

    `plink --bfile BFILE_PREFIX --keep BFILE_PREFIX.fam.split.aa --out BFILE_PREFIX.aa`
    
  - Extract the MHC region from your input genomic data: 

    `plink --bfile BFILE_PREFIX.aa --chr 6 --from-mb 28 --to-mb 34 --recode vcf --out BFILE_PREFIX.aa`
	
  - Go to the Michigan Imputation Server (https://imputationserver.sph.umich.edu). Select "HLA Imputation". Impute HLA region using the "Four-digit Multi-ethnic HLA reference panel". Note that this may take many days depending on server resources.
	
2. For the non-HLA score, format SNP IDs in the PLINK .bim file so that each SNP ID follows the CHR:POS format (e.g., 1:10583). Additionally, for improved speed and efficiency, we recommend filtering the PLINK file to include only the SNPs of interest using the` --extract <filename>` option.

	`cut -f6 PRScs_EUR_phi1e-05_weights.txt > extract_list.txt`
	
	(Note: for hg19, use cut -f3. For hg38, use cut -f6)
	
	`plink --bfile BFILE_PREFIX --extract extract_list.txt --make-bed --out BFILE_PREFIX.reduced`

3. Calculate scores one of two ways:   

  **OPTION 1: If continuous genetic ancestry is not available (T1D MAPS)**
    
  Run in command line to calculate scores: `Rscript T1D_MAPS.R [HLA_file] [non_HLA_file] [build] [path_to_PLINK] [working_dir]`
 
  -
    - [HLA_file]: Input file in compressed VCF format (.vcf.gz)
    
    - [non_HLA_file]: PLINK binary file prefix (BFILE format)
    
    - [build]: Reference genome build; must be either "hg19" or "hg38"
    
    - [path_to_PLINK]: Full path to the PLINK executable (version 1.9)
    
    - [working_dir]: Path to the directory containing all files from the cloned Git repository. Make sure the path ends with a forward slash (`/`)
                    
    - The file containing the calculated scores will be named following the pattern: `T1DMAPS_<non_HLA_file>.total_score.txt`
      
  **OPTION 2: If continuous genetic ancestry is available (T1D MAPS2)**
    
  1. Calculate continuous probability of genetic ancestry. Refer to: https://github.com/atgu/ukbb_pan_ancestry/
          
  2. Calculate T1D GRS2. Refer to: https://github.com/sethsh7/PRSedm
          
  3. Run in command line to calculate scores: `Rscript T1D_MAPS.R [HLA_file] [non_HLA_file] [build] [path_to_PLINK] [working_dir] [genetic_prob_file] [t1d_grs2_file]`
                  
      - [HLA_file]: Input file in compressed VCF format (.vcf.gz)
    
      - [non_HLA_file]: PLINK binary file prefix (BFILE format)
    
      - [build]: Reference genome build; must be either "hg19" or "hg38"
    
      - [path_to_PLINK]: Full path to the PLINK executable (version 1.9)
      
      - [working_dir]: Path to the directory containing all files from the cloned Git repository. Make sure the path ends with a forward slash (`/`)
      
      - [genetic_prob_file]: Output file from calculating the continuous probability of genetic ancestry
      
      - [t1d_grs2_file]: Output file from calculating T1D GRS2
                  
      - The file containing the calculated scores will be named following the pattern: `T1DMAPS2_<non_HLA_file>.total_score.txt`
            
## Development:

Developed by Stella Nam, Aaron Deutsch, Josep Mercader, and Miriam Udler at Massachusetts General Hospital, Boston, Massachusetts. For questions, please contact Aaron Deutsch at ajdeutsch@mgh.harvard.edu. 

## Citation:

Deutsch AJ, Bell AS, Michalek DA, Burkholder AB, Nam S, Kreienkamp RJ, Sharp SA, Huerta-Chagoya A, Mandla R, Nanjala R, Luo Y, Oram RA, Florez JC, Onengut-Gumuscu S, Rich SS, Motsinger-Reif AA, Manning AK, Mercader JM, Udler MS. Development and Validation of a Type 1 Diabetes Multi-Ancestry Polygenic Score. medRxiv.

## References:

- Tutorial on HLA imputation: 
	Sakaue et al., Nat Protoc 2023 Sep;18(9):2625-2641.  doi: 10.1038/s41596-023-00853-4
- For more information about T1D GRS2 and T1D GRS2x:
	Sharp et al., Diabetes Care 2019 Feb;42(2):200-207.  doi: 10.2337/dc18-1785.
	Luckett et al., Diabetes Care 2025 Jun 1;48(6):e81-e83.  doi: 10.2337/dc25-0142
- For more information about inferring genetic ancestry in the Pan-UK Biobank:
	Karczewski et al., medRxiv 2024.03.13.24303864. doi: 10.1101/2024.03.13.2430386