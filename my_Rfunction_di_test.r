CAST<-function(i,j){
  #initialize variables
  pv.dom.dom=NA # Likelihood Ratio Test (LRT) p-value for a (Gene/variant 1) dominant vs (Gene/variant 2) dominant digenic analysis adjusted with given covariates (e.g. principal components, see pheno file)
  pw.dom.dom.noAdj=NA # Likelihood Ratio Test (LRT) p-value for a (Gene/variant 1) dominant vs (Gene/variant 2) dominant digenic analysis without adjustment
  m1.dom1=NA # Gives the number of carriers under a dominant model at genomic location 1
  n1.dom1=NA # Gives the number of non-carriers under a dominant model at genomic location 1 (i.e. individuals that do not harbor any variant in genomic location 1)
  m1.dom2=NA # Gives the number of carriers under a dominant model at genomic location 2
  n1.dom2=NA # Gives the number of carriers under a dominant model at genomic location 2 (i.e. individuals that do not harbor any variant in genomic location 2)
  m1.dom.dom12=NA # Gives the number of "double carriers" under a dominant vs dominant model, i.e. the number of individuals harboring variants at the 2 tested genomic locations.
  y=NULL # phenotype vector
  nvardom1=NA # Number of variants contributing to genomic location 1 (e.g. number of variants tested within genomic location 1) when working under a dominant model
  nvardom2=NA # Number of variants contributing to genomic location 2 (e.g. number of variants tested within genomic location 2) when working under a dominant model
  dom1=NULL # matrix with carrier status for all cases for the given gene under study. Gives the number of cases tested.
  or_0pc_dom=NA # Odds ratio for the interaction effect between the 2 genomic regions under a dominant vs dominant model without adjusment on covariates
  or_3pc_dom=NA # Odds ratio for the interaction effect between the 2 genomic regions under a dominant vs dominant model with adjusment on covariates
  colnames(phenotype)=c("id","aff","pc1","pc2","pc3")
  #get the data for gene 1 and 2
  gds<-seqOpen(paste(direct,"/vcf/",vcf,".gds",sep=""))
  seqSetFilter(gds,variant.id=i,verbose=FALSE)
  chr = seqGetData(gds,"chromosome")
  pos = seqGetData(gds,"position")
  all = seqGetData(gds,"allele")
  Gnomad<-seqGetData(gds,"annotation/info/Gnomad_AF_UF")$data
  AF<-seqGetData(gds,"annotation/info/AF")$data
  EFFECT<-seqGetData(gds,"annotation/info/EFFECT")
  HGVS_C=seqGetData(gds,"annotation/info/HGVS_C")
  HGVS_P=seqGetData(gds,"annotation/info/HGVS_P")
  genotype1<-2-seqGetData(gds,"$dosage")
  Genotype1=cbind.data.frame(id,genotype1)
  Genotype1=merge(phenotype,Genotype1,by="id")
  phenolist=Genotype1[,"aff"]
  covar=cbind(Genotype1[,c(3,4,5)],1)
  genotype1=as.matrix(Genotype1[,6:ncol(Genotype1)])
  chr1 = seqGetData(gds,"chromosome")[1]
  pos1i = seqGetData(gds,"position")[1]
  pos1f = seqGetData(gds,"position")[length(seqGetData(gds,"position"))]
  seqResetFilter(gds,verbose=FALSE)
  seqClose(gds)
  
  gds.rare<-seqOpen(paste(direct,"/vcf/",vcf.rare,".gds",sep=""))
  seqSetFilter(gds.rare,variant.id=j,verbose=FALSE)
  genotype2<-2-seqGetData(gds.rare,"$dosage")
  Genotype2=cbind.data.frame(id,genotype2)
  Genotype2=merge(phenotype,Genotype2,by="id")
  genotype2=as.matrix(Genotype2[,6:ncol(Genotype2)])
  chr2 = seqGetData(gds.rare,"chromosome")[1]
  pos2i = seqGetData(gds.rare,"position")[1]
  pos2f = seqGetData(gds.rare,"position")[length(seqGetData(gds.rare,"position"))]
  seqResetFilter(gds.rare,verbose=FALSE)
  seqClose(gds.rare)
  ## LD
  if(chr1==chr2 & pos1i<pos2i){
    dist=(pos2i-pos1f)/1e+06
  }else if(chr1==chr2 & pos1i>=pos2i){
    dist=(pos1i-pos2f)/1e+06
  }else{
    dist=NA
  }
  
  
  # inds to exclude
  inds.include.1.cases = !is.na(phenolist) & phenolist==2
  
  #estimate alt all freq
  aaf1<-apply(subset(genotype1,inds.include.1.cases)/2,2,mean,na.rm=T)
  aaf1<-ifelse(is.na(aaf1),0,aaf1)
  
  aaf2<-apply(subset(genotype2,inds.include.1.cases)/2,2,mean,na.rm=T)
  aaf2<-ifelse(is.na(aaf2),0,aaf2)
  
  ## non covered individuals for any variant of the gene
  non.cov.ind1=apply(subset(genotype1,T,select=aaf1<=AAF1&aaf1>0&aaf1<1),1,function(i) all(is.na(i)))
  non.cov.ind2=apply(subset(genotype2,T,select=aaf2<=AAF2&aaf2>0&aaf2<1),1,function(i) all(is.na(i)))
  
  inds.include.2.cases = !is.na(phenolist) & !non.cov.ind1 & !non.cov.ind2 & phenolist==2
  
  if(sum(inds.include.2.cases)>0 & gene[i][1]!=gene.rare[j][1]){
    
    
    #Recode in dom mode
    genoDom1<-subset(genotype1,inds.include.2.cases,select=aaf1<=AAF1&aaf1>0&aaf1<1)
    dom1 = as.numeric(apply(genoDom1, 1, function(i) sum(i, na.rm=T)) > 0)
    nvardom1=ncol(genoDom1)
    
    genoDom2<-subset(genotype2,inds.include.2.cases,select=aaf2<=AAF2&aaf2>0&aaf2<1)
    dom2 = as.numeric(apply(genoDom2, 1, function(i) sum(i, na.rm=T)) > 0)
    nvardom2=ncol(genoDom2)
    
 
    if(nvardom1>0& length(table(dom1))>1 & nvardom2>0& length(table(dom2))>0){
      dom1 = as.numeric(apply(genoDom1, 1, function(i) sum(i, na.rm=T)) > 0)
      m1.dom1 = sum(dom1)
      n1.dom1 = length(which(dom1 == 0))
      m1.dom2 = sum(dom2)
      n1.dom2 = length(which(dom2 == 0))
      
      # double carriers
      dom.dom12 = dom1*dom2
      m1.dom.dom12 = sum(dom.dom12)
      
      #choose as the first gene the one with most proportion of carriers
      prop.dom2 = m1.dom2/n1.dom2
      prop.dom1 = m1.dom1/n1.dom1
      if(prop.dom2>=prop.dom1){
        first.gene = dom2
        second.gene = dom1
      }else{
        first.gene = dom1
        second.gene = dom2
      }
      
      
      mod.null.dom =  glm(first.gene ~ ., family = binomial,data=covar[inds.include.2.cases,])
      mod.full.dom = glm(first.gene ~ .+second.gene, family = binomial,data=covar[inds.include.2.cases,])
      
      mod.null.dom.noAdj =  glm(first.gene ~ 1, family = binomial,data=covar[inds.include.2.cases,])
      mod.full.dom.noAdj = glm(first.gene ~ second.gene, family = binomial,data=covar[inds.include.2.cases,])
      
      # only genes or variants with at least 5% of carriers among the covered ones
      if(( mod.full.dom$converged ) & ( ! mod.full.dom$boundary ) & m1.dom1/(m1.dom1+n1.dom1)>0.05 &
           ( mod.full.dom.noAdj$converged ) & ( ! mod.full.dom.noAdj$boundary )){
        lrt.dom = mod.null.dom$deviance - mod.full.dom$deviance
        pv.dom.dom = pchisq(lrt.dom, df = 1, lower.tail = F)
        
        lrt.dom.noAdj = mod.null.dom.noAdj$deviance - mod.full.dom.noAdj$deviance
        pw.dom.dom.noAdj = pchisq(lrt.dom.noAdj, df = 1, lower.tail = F)
        or_0pc_dom<-as.vector(exp(coefficients(mod.full.dom.noAdj)["second.gene"])) 
        or_3pc_dom<-as.vector(exp(coefficients(mod.full.dom)["second.gene"]))
      } 
      
    } 
    
    
    
  }else{
    
  }
                              
  # Output. Columns as in the output file. A small description for each column can be found at the begining of the script.
  
  res<-cbind(chr=chr,pos=pos,all=all,Gnomad_AF=Gnomad,cohort_AF=AF,effect=EFFECT,HGVS_C=HGVS_C,HGVS_P=HGVS_P,
             genename1= gene[i][1], genename2= gene.rare[j][1], LD= dist, sample.size = length(dom1),
             nvardom1 = nvardom1, Ncarriers.dom.dom12.cases = m1.dom.dom12 , 
             Nnoncarriers.dom1.cases=n1.dom1, Ncarriers.cases.dom1 = m1.dom1,nvardom2 = nvardom2, 
             Nnoncarriers.dom2.cases=n1.dom2, Ncarriers.cases.dom2 = m1.dom2, pvalueLRT.dom.dom = pv.dom.dom,
             pvalueLRT.dom.dom.noAdj = pw.dom.dom.noAdj, OR.dom.adj = or_3pc_dom, OR.dom = or_0pc_dom)
             
  
  return(res)
  
}



index1<-function(i){
  ind<-intersect(which(gene==i),variant1)
  return(ind)
}

index2<-function(i){
  ind<-intersect(which(gene.rare==i),variant2)
  return(ind)
}



analysis<-function(outfile=NULL,Gnomad_AF1_low=0.1,Gnomad_AF1_up=0.9,Gnomad_AF2=0.0001,AAF1=1,AAF2=1,set1=c("frameshift","stop_gained","stop_lost","start_lost","splice_donor","splice_acceptor"),set2="ALL",
                   pheno=NULL,covariates=NULL,set_genes1=c("IFNG"),set_genes2=c("IL12RB1"),ncore=4,feno=NULL){
  if(is.null(pheno)) stop("phenotype is missing")
  if(is.null(outfile)) outfile<-"Rcast.out"
  if(!exists("Gnomad")|!exists("EFFECT")|!exists("gene")|!exists("WARN")) {
    gds<-seqOpen(gdsfile)
    if(!exists("Gnomad")) Gnomad<-seqGetData(gds,"annotation/info/Gnomad_AF_UF")$data
    if(!exists("EFFECT")) EFFECT<-seqGetData(gds,"annotation/info/EFFECT")
    if(!exists("gene")) gene<-seqGetData(gds,"annotation/info/GENE")
    if(!exists("WARN")) WARN<-seqGetData(gds,"annotation/info/WARN")
    seqClose(gds)
  }
  
  # Filter variants
  variant1 = which(Gnomad<=Gnomad_AF1_up & Gnomad>=Gnomad_AF1_low & !grepl("WARN",WARN))
  variant2 = which(!grepl("WARN",WARN.rare))
  g2test2<-names(table(gene.rare[variant2]))
  
  #Get phenotype
  gds <- seqOpen(paste(direct,"/vcf/",vcf,".gds",sep=""))
  id<-seqGetData(gds,"sample.id")
  phenolist<-phenotype[as.factor(id),pheno]
  seqClose(gds)
  
  #Get covariates
  covar<-cbind(phenotype[as.factor(id),covariates],1)
  
  
  cl<-makeCluster(ncore)
  clusterExport(cl=cl,varlist=c("gene","variant1","variant2","phenotype","phenolist","pheno","id","covar","AAF1","AAF2","dir","vcf","vcf.rare","gene.rare","direct"), envir=environment())
  clusterEvalQ(cl, {
    library("SeqArray")
    source(paste(direct,"/scriptR/my_Rfunction_di_test.r",sep=""))
  })
  registerDoParallel(cl)
  
  ## Create list of variants or genes to test with CAST function
  liste1<-as.list(variant1)
  liste2<-parLapply(cl,g2test2,index2)

  ## Run CAST function for every parwise combination
  resultpar<-foreach(i=liste1,.combine=rbind,.packages="SeqArray") %:% foreach(j=liste2,.combine=rbind,.packages="SeqArray") %dopar% {
  CAST(i,j)}
  effects=ifelse("missense" %in% set2,"MISSLOF","LOF")
  
  ## Output data to pvalues folder
  write.table(resultpar,file=paste(feno,"/di/pvalues/di_Gnomad_AF1_low",Gnomad_AF1_low,"_Gnomad_AF2",Gnomad_AF2,"_AF1",AAF1,"_AF2",AAF2,"_",effects,"_vcf.out",sep=""),row.names=F,col.names=T,quote=F,sep="\t")
  stopCluster(cl)
  return(resultpar)
}
