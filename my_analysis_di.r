## Copyright 2020 HGID Institut Imagine, Gaspard Kerner (gakerner@pasteur.fr) 
#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

# Tab- or space-separated input phenotype file: One line per individual and J columns where J is >=2 (Sample_name, affected status (column_name= "aff"), J-2 covariates (one per column)). 
# Individual identifiers are identical to those in the input genotype file.
# NA label in "aff" column automatically excludes individuals from the analysis
feno=args[1]

# Parallelisation is performed either with doParallel package (loop is split using n cores) or by splitting input genotype files into several smaller files.
# The script expects two genotype files: a background file with all genotype information as for monogenic analyses and one that may have been split by genes or variants or may be the same as the background file.
# For example, jobs may account for 20 x 20,000 tests if the background file includes 20,000 genes and the smaller file was split into files consisting of 20 genes.
# The script can be run using Job arrays in a cluster machine.
# Names are expected without suffix. It assumes a .vcf.gz file
vcf=args[2]
vcf.all=args[3]

# If one of the files was split, subfile number can be passed as input to be outputed in the output results file.
num=args[4]

# Working directory to be included accordingly. 
# In the working directory a data folder is expected including phenotype files, VCF files and Rscripts
dir="~"

# R libraries to be used. Foreach and doParallel are important if several cores are to be used. 
# Otherwise the code is to be modified to perform a normal double loop.
library("SeqArray")
library("foreach")
library("parallel")
library("doParallel")

# ~/scriptR/my_Rfunction_di_NY.r stores the digenic function to be run with the script.
source(paste(dir,"/data/my_Rfunction_di.r",sep=""))

# create Output directory
ifelse(!dir.exists(file.path(dir, paste("di/",pheno,"/pvalues",sep=""))), dir.create(file.path(dir, paste("di/",pheno,"/pvalues",sep=""))), FALSE)

#Initialize variables
# Define the name of your VCF files. VCF.all is the background genotype file
VCF=paste(dir,"/data/",vcf,".vcf.gz",sep="")
VCF.all=paste(dir,"/data/",vcf.all,".vcf.gz",sep="")

# Define the name of your phenotype file
phenofile=paste("pheno_",feno,".txt",sep="")

# Read phenotype file
phenotype<-read.table(paste(dir,"/data/",phenofile,sep=""),header=T)
rownames(phenotype)<-phenotype[,1]

# Creates .gds file from VCF file. Uses SeqArray package. Info should be consistent in the header of the VCF file.
if(!file.exists(paste(dir,"/data/",vcf,".gds",sep=""))){
  seqVCF2GDS(VCF,paste(dir,"/data/",vcf,".gds",sep=""),storage.option="ZIP_RA.max",parallel=8L)
}

# Read and keep information from the VCF file to be used for analysis (e.g. MAF from public databases, annotation, gene name, etc.)
gds<-seqOpen(paste(dir,"/data/",vcf,".gds",sep=""))

# e.g. extract gnomAD MAFs if annotated previously on VCF (from genome [GnomadG_max] or exome data [GnomadE_max]). Maximum of all sub-populations
GnomadE_max<-seqGetData(gds,"annotation/info/AF_E_popmax")$data
GnomadG_max<-seqGetData(gds,"annotation/info/AF_G_popmax")$data
# e.g. Use allele counts of gnomAD to perform call-rate
GnomadE_AC<-seqGetData(gds,"annotation/info/AC_E")$data
GnomadE_AN<-seqGetData(gds,"annotation/info/AN_E")
GnomadG_AC<-seqGetData(gds,"annotation/info/AC_G")$data
GnomadG_AN<-seqGetData(gds,"annotation/info/AN_G")
# e.g. Allele frequency from the population under study (assuming that the VCF is restricted to this population)
AF<-seqGetData(gds,"annotation/info/AF")$data
# e.g. Predicted effect (missense, intronic, frameshift, stop_gained, etc.)
EFFECT<-seqGetData(gds,"annotation/info/EFFECT")
# e.g. Gene name
gene<-seqGetData(gds,"annotation/info/GENE")
# e.g. CADD score
Cadd<-seqGetData(gds,"annotation/info/CADD")
WARN<-seqGetData(gds,"annotation/info/WARN")
# e.g. Transcript type (protein coding, etc.)
TRANSCRIPT<-seqGetData(gds,"annotation/info/TRANSCRIPT_TYPE")

# e.g. extract gnomAD MAFs if annotated previously on VCF (from genome [GnomadG_max] or exome data [GnomadE_max]). MAF by sub-population
GnomadE_afr<-seqGetData(gds,"annotation/info/AF_E_afr")$data
GnomadE_amr<-seqGetData(gds,"annotation/info/AF_E_amr")$data
GnomadE_nfe<-seqGetData(gds,"annotation/info/AF_E_nfe")$data
GnomadE_asj<-seqGetData(gds,"annotation/info/AF_E_asj")$data
GnomadE_eas<-seqGetData(gds,"annotation/info/AF_E_eas")$data
GnomadE_fin<-seqGetData(gds,"annotation/info/AF_E_fin")$data
GnomadE_sas<-seqGetData(gds,"annotation/info/AF_E_sas")$data

GnomadG_afr<-seqGetData(gds,"annotation/info/AF_G_afr")$data
GnomadG_amr<-seqGetData(gds,"annotation/info/AF_G_amr")$data
GnomadG_nfe<-seqGetData(gds,"annotation/info/AF_G_nfe")$data
GnomadG_asj<-seqGetData(gds,"annotation/info/AF_G_asj")$data
GnomadG_eas<-seqGetData(gds,"annotation/info/AF_G_eas")$data
GnomadG_fin<-seqGetData(gds,"annotation/info/AF_G_fin")$data


seqClose(gds)

# Repeat the same process for the background genotype file

if(!file.exists(paste(dir,"/data/",vcf.all,".gds",sep=""))){
  seqVCF2GDS(VCF.all,paste(dir,"/data/",vcf.all,".gds",sep=""),storage.option="ZIP_RA.max",parallel=8L)
}

gds.all<-seqOpen(paste(dir,"/data/",vcf.all,".gds",sep=""))

GnomadE_max.all<-seqGetData(gds.all,"annotation/info/AF_E_popmax")$data
GnomadG_max.all<-seqGetData(gds.all,"annotation/info/AF_G_popmax")$data
GnomadE_AC.all<-seqGetData(gds.all,"annotation/info/AC_E")$data
GnomadE_AN.all<-seqGetData(gds.all,"annotation/info/AN_E")
GnomadG_AC.all<-seqGetData(gds.all,"annotation/info/AC_G")$data
GnomadG_AN.all<-seqGetData(gds.all,"annotation/info/AN_G")
AF.all<-seqGetData(gds.all,"annotation/info/AF")$data
EFFECT.all<-seqGetData(gds.all,"annotation/info/EFFECT")
gene.all<-seqGetData(gds.all,"annotation/info/GENE")
Cadd.all<-seqGetData(gds.all,"annotation/info/CADD")
WARN.all<-seqGetData(gds.all,"annotation/info/WARN")
TRANSCRIPT.all<-seqGetData(gds.all,"annotation/info/TRANSCRIPT_TYPE")
AFNOVA.all<-seqGetData(gds.all,"annotation/info/AAF_nova")


GnomadE_afr.all<-seqGetData(gds.all,"annotation/info/AF_E_afr")$data
GnomadE_amr.all<-seqGetData(gds.all,"annotation/info/AF_E_amr")$data
GnomadE_nfe.all<-seqGetData(gds.all,"annotation/info/AF_E_nfe")$data
GnomadE_asj.all<-seqGetData(gds.all,"annotation/info/AF_E_asj")$data
GnomadE_eas.all<-seqGetData(gds.all,"annotation/info/AF_E_eas")$data
GnomadE_fin.all<-seqGetData(gds.all,"annotation/info/AF_E_fin")$data
GnomadE_sas.all<-seqGetData(gds.all,"annotation/info/AF_E_sas")$data

GnomadG_afr.all<-seqGetData(gds.all,"annotation/info/AF_G_afr")$data
GnomadG_amr.all<-seqGetData(gds.all,"annotation/info/AF_G_amr")$data
GnomadG_nfe.all<-seqGetData(gds.all,"annotation/info/AF_G_nfe")$data
GnomadG_asj.all<-seqGetData(gds.all,"annotation/info/AF_G_asj")$data
GnomadG_eas.all<-seqGetData(gds.all,"annotation/info/AF_G_eas")$data
GnomadG_fin.all<-seqGetData(gds.all,"annotation/info/AF_G_fin")$data


seqClose(gds.all)
#Define variant sets

# You can define a set of variant annotations in order to group variants according to their annotation. LOF are high impact annotation while MISSLOF
# includes also missense or inframe annotations.
LOF=c("frameshift","stop_gained","stop_lost","start_lost","splice_donor","splice_acceptor")
MISSLOF=c("inframe","missense","exon_loss","initiator_codon","frameshift","stop_gained","stop_lost","start_lost","splice_donor","splice_acceptor")

#Run the analysis. Function is defined in my_Rfunction_di.r

# Parameters for the digenic function: MAF thresholds for variant aggregation, variant annotation, column for affected status in phenotype file and column for covariates in phenotype file.
# Analysis can also be restricted to a set of genes for each of the two genotype files (set_genes1, set_genes2). Cores are used if parallelisation inside R is applied. Otherwise to be set to 1.
res<-analysis(outfile="Outfile_name",Gnomad_AF1=0.02,Gnomad_AF2=0.02,set=LOF,pheno="aff",
              covariates=c("pc1","pc2","pc3"),set_genes1="ALL",set_genes2="ALL",ncore=48)
