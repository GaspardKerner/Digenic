# Digenic

Instructions on how to use Rscripts are available on the Rscripts themselves. 
A double loop parallelisation is incorporated to optimize pairwise interaction analyses. 
However, it can be faster to adopt manual parallelization by splitting genotype input files and running parallel jobs (one per subfile). 
More detail information and updated scripts are expected to be available in the near future.


# General information

Scripts were developed to perform digenic analyses. The script is an example on how to perform gene x gene interaction analyses. However, gene x variant is also straightforward using these scripts. 
For gene x gene interaction tests, the general idea is to aggregate rare (+ low-frequency) variant information at the gene level and assess unexpected high correlation between carriers of variants from two different genes. The test is conducted under a case-only framework which avoids the use of controls. Correction for population stratification (PS) and linkage disequilibrium (LD) are needed. For LD we propose to focus only on genes for which variants are separated by at least 2 Mb and for PS adjustment on principal components is expected to be enough. For more details refer to https://www.pnas.org/content/early/2020/07/23/1920650117.

# Test example

ToRun.sh provides an example of gene vs variants analysis. Please download vcf example files and follow instructions.
