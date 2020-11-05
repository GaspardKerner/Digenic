#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

feno=args[1]
vcf=args[2]
vcf.rare=args[3]
vcf.number=args[4]
direct=args[5]

feno="test"
vcf="test.common"
vcf.rare="test.rare"
vcf.number=""
direct="/volumes/vol002/Backup_vol001/dbGAP/cranyio2"

library("SeqArray")
library("foreach")
library("parallel")
library("doParallel")
source(paste(direct,"/scriptR/my_Rfunction_di_test.r",sep=""))


setwd(paste(direct,"/",feno,sep=""))
#Initialize variables
VCF=paste(direct,"/vcf/",vcf,".vcf.gz",sep="")
VCF.rare=paste(direct,"/vcf/",vcf.rare,".vcf.gz",sep="")
phenofile=paste("pheno_",feno,"_191cranio.txt",sep="")

#READ PHENOTYPE
phenotype<-read.table(paste("phenotype/",phenofile,sep=""),header=T)
rownames(phenotype)<-phenotype[,1]

if(!file.exists(paste(direct,"/vcf/",vcf,".gds",sep=""))){
  seqVCF2GDS(VCF,paste(direct,"/vcf/",vcf,".gds",sep=""),storage.option="ZIP_RA.max",parallel=8L)
}

gds<-seqOpen(paste(direct,"/vcf/",vcf,".gds",sep=""))

Gnomad<-seqGetData(gds,"annotation/info/Gnomad_AF_UF")$data
AF<-seqGetData(gds,"annotation/info/AF")$data
EFFECT<-seqGetData(gds,"annotation/info/EFFECT")
gene<-seqGetData(gds,"annotation/info/GENE")
WARN<-seqGetData(gds,"annotation/info/WARN")

seqClose(gds)



if(!file.exists(paste(direct,"/vcf/",vcf.rare,".gds",sep=""))){
  seqVCF2GDS(VCF.rare,paste(direct,"/vcf/",vcf.rare,".gds",sep=""),storage.option="ZIP_RA.max",parallel=8L)
}

gds.rare<-seqOpen(paste(direct,"/vcf/",vcf.rare,".gds",sep=""))

Gnomad.rare<-seqGetData(gds.rare,"annotation/info/Gnomad_AF_UF")$data
AF.rare<-seqGetData(gds.rare,"annotation/info/AF")$data
EFFECT.rare<-seqGetData(gds.rare,"annotation/info/EFFECT")
gene.rare<-seqGetData(gds.rare,"annotation/info/GENE")
WARN.rare<-seqGetData(gds.rare,"annotation/info/WARN")

seqClose(gds.rare)
#Define variant sets

LOF=c("frameshift","stop_gained","stop_lost","start_lost","splice_donor","splice_acceptor")
MISSLOF=c("missense","exon_loss","initiator_codon","frameshift","stop_gained","stop_lost","start_lost","splice_donor","splice_acceptor")

#Run the analysis


res<-analysis(Gnomad_AF1_low=0.1,Gnomad_AF1_up=0.9,Gnomad_AF2=0.0001,AAF1=1,AAF2=1,set1="ALL",set2=MISSLOF,pheno="aff",covariates=c("pc1","pc2","pc3"),set_genes1="ALL",set_genes2="ALL",ncore=24,feno=feno)
