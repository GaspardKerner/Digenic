#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

feno=args[1]
vcf=args[2]
vcf.all=args[3]
num=args[4]
#dir="/rugpfs/fs0/casa_lab/scratch/gkerner/PARIS_NY_YALE_filtered/"
dir="/home/gkerner/COVID"

library("SeqArray")
library("foreach")
library("parallel")
library("doParallel")
source(paste(dir,"/scriptR/my_Rfunction_di_NY.r",sep=""))

setwd(paste(dir,"/",feno,sep=""))
#Initialize variables
VCF=paste(dir,"/vcf/tmp/",vcf,".vcf.gz",sep="")
VCF.all=paste(dir,"/vcf/",vcf.all,".vcf.gz",sep="")
phenofile=paste("pheno_",feno,".txt",sep="")

#READ PHENOTYPE
phenotype<-read.table(paste("phenotype/",phenofile,sep=""),header=T)
rownames(phenotype)<-phenotype[,1]

if(!file.exists(paste(dir,"/vcf/tmp/",vcf,".gds",sep=""))){
  seqVCF2GDS(VCF,paste(dir,"/vcf/tmp/",vcf,".gds",sep=""),storage.option="ZIP_RA.max",parallel=8L)
}

gds<-seqOpen(paste(dir,"/vcf/tmp/",vcf,".gds",sep=""))

Gnomad<-seqGetData(gds,"annotation/info/Gnomad_AF_UF")$data
GnomadAFR<-seqGetData(gds,"annotation/info/Gnomad_AF_AFR_UF")
GnomadNFE<-seqGetData(gds,"annotation/info/Gnomad_AF_NFE_UF")
GnomadAMR<-seqGetData(gds,"annotation/info/Gnomad_AF_AMR_UF")
GnomadASJ<-seqGetData(gds,"annotation/info/Gnomad_AF_ASJ_UF")
GnomadEAS<-seqGetData(gds,"annotation/info/Gnomad_AF_EAS_UF")
GnomadFIN<-seqGetData(gds,"annotation/info/Gnomad_AF_FIN_UF")
GnomadSAS<-seqGetData(gds,"annotation/info/Gnomad_AF_SAS_UF")
KGP<-seqGetData(gds,"annotation/info/AF_KGP3")$data
AF<-seqGetData(gds,"annotation/info/AF")$data
EFFECT<-seqGetData(gds,"annotation/info/EFFECT")
gene<-seqGetData(gds,"annotation/info/GENE")
Cadd<-seqGetData(gds,"annotation/info/CADD")
#MSC<-seqGetData(gds,"annotation/info/MSC")
WARN<-seqGetData(gds,"annotation/info/WARN")

seqClose(gds)


if(!file.exists(paste(dir,"/vcf/",vcf.all,".gds",sep=""))){
  seqVCF2GDS(VCF.all,paste(dir,"/vcf/",vcf.all,".gds",sep=""),storage.option="ZIP_RA.max",parallel=8L)
}

gds.all<-seqOpen(paste(dir,"/vcf/",vcf.all,".gds",sep=""))

Gnomad.all<-seqGetData(gds.all,"annotation/info/Gnomad_AF_UF")$data
GnomadAFR.all<-seqGetData(gds.all,"annotation/info/Gnomad_AF_AFR_UF")
GnomadNFE.all<-seqGetData(gds.all,"annotation/info/Gnomad_AF_NFE_UF")
GnomadAMR.all<-seqGetData(gds.all,"annotation/info/Gnomad_AF_AMR_UF")
GnomadASJ.all<-seqGetData(gds.all,"annotation/info/Gnomad_AF_ASJ_UF")
GnomadEAS.all<-seqGetData(gds.all,"annotation/info/Gnomad_AF_EAS_UF")
GnomadFIN.all<-seqGetData(gds.all,"annotation/info/Gnomad_AF_FIN_UF")
GnomadSAS.all<-seqGetData(gds.all,"annotation/info/Gnomad_AF_SAS_UF")
KGP.all<-seqGetData(gds.all,"annotation/info/AF_KGP3")$data
AF.all<-seqGetData(gds.all,"annotation/info/AF")$data
EFFECT.all<-seqGetData(gds.all,"annotation/info/EFFECT")
gene.all<-seqGetData(gds.all,"annotation/info/GENE")
Cadd.all<-seqGetData(gds.all,"annotation/info/CADD")
#MSC<-seqGetData(gds,"annotation/info/MSC")
WARN.all<-seqGetData(gds.all,"annotation/info/WARN")

seqClose(gds.all)
#Define variant sets

LOF=c("frameshift","stop_gained","stop_lost","start_lost","splice_donor","splice_acceptor")
MISSLOF=c("inframe","missense","exon_loss","initiator_codon","frameshift","stop_gained","stop_lost","start_lost","splice_donor","splice_acceptor")


# CONN = read.table(paste(dir,"/CONN_TB.txt",sep=""),header=F)
# CONN = as.vector(CONN[,1])
# 

#Run the analysis


res<-analysis(outfile="NY_Paris_YALE_mut",Gnomad_AF1=0.01,Gnomad_AF2=0.01,KGPAF1=0.01,KGPAF2=0.01,CADD=TRUE,set=MISSLOF,pheno="aff",covariates=c("pc1","pc2","pc3"),set_genes1="ALL",set_genes2="ALL",ncore=12,feno=feno)
#res<-analysis(outfile="NY_Paris_YALE_mut",Gnomad_AF1=0.01,Gnomad_AF2=0.01,KGPAF1=0.01,KGPAF2=0.01,CADD=FALSE,set=MISSLOF,pheno="aff",covariates=c("pc1","pc2","pc3"),set_genes1="ALL",set_genes2="ALL",ncore=48,feno=feno)
