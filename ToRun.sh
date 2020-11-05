#!/bin/bash

# This script has to be modified accordingly. Read comments

direct="/PATH/TO"
feno="test"
VCFcommon="test.common"
VCFrare="test.rare"
cd "$direct"

mkdir -p vcf
# Include VCF files in vcf folder (e.g. one VCF for rare and one VCF for common variants)
# By default, VCFrare will aggregate rare variants into a gene unit and VCFcommon will consist of single common variants
# vcf files must be gzipped
# Use as test example test.rare.vcf.gz and test.common.vcf.gz (you can also use directly the .gds files)

mkdir -p "$feno"
# Results and phenotype files will be stored here


mkdir -p "$feno"/phenotype
# Include phenotype file 
# at least 2 columns with sample_name and affection status (2 cases, others are not considered)
# covariates can be included in the following columns (e.g. principal components)
# Follow the example below:


#sample_names aff pc1 pc2 pc3
#M-P1313 2 -8.85007e-05 -0.0149749 0.00505057
#M-P1317 2 -0.0108567 -0.0205866 0.00822436

# Use as a test example pheno_test_191cranio.txt


mkdir -p scriptR
#Include Rscripts that can be downloaded here

mkdir -p "$feno"/di/pvalues


# Run the analysis
Rscript scriptR/my_analysis_di_test.r "$feno" "$VCFcommon" "$VCFrare" "$direct"

