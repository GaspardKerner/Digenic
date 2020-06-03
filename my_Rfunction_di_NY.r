CAST<-function(i,j){
  #initialize variables
  pv.dom.dom=NA
  pw.dom.dom.noAdj=NA
  m1.dom1=NA
  n1.dom1=NA
  m1.dom2=NA
  n1.dom2=NA
  pv.rec.dom=NA
  pw.rec.dom.noAdj=NA
  pv.dom.rec=NA
  pw.dom.rec.noAdj=NA
  pv.rec.rec=NA
  pw.rec.rec.noAdj=NA
  m1.rec1=NA
  n1.rec1=NA
  m1.rec2=NA
  n1.rec2=NA
  m1.dom.rec12=NA
  m1.rec.dom12=NA
  m1.rec.rec12=NA
  m1.dom.dom12=NA
  y=NULL
  nvardom1=NA
  nvardom2=NA
  nvarrec1=NA
  nvarrec2=NA
  dom1=NULL
  dom2=NULL
  rec1=NULL
  rec2=NULL
  or_0pc_dom.dom=NA
  or_3pc_dom.dom=NA
  or_0pc_rec.rec=NA
  or_3pc_rec.rec=NA
  or_0pc_dom.rec=NA
  or_3pc_dom.rec=NA
  or_0pc_rec.dom=NA
  or_3pc_rec.dom=NA
  colnames(phenotype)=c("id","aff","pc1","pc2","pc3")
  #get the data for gene 1 and 2
  gds<-seqOpen(paste(dir,"/vcf/tmp/",vcf,".gds",sep=""))
  seqSetFilter(gds,variant.id=i,verbose=FALSE)
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
  
  gds.all<-seqOpen(paste(dir,"/vcf/",vcf.all,".gds",sep=""))
  seqSetFilter(gds.all,variant.id=j,verbose=FALSE)
  genotype2<-2-seqGetData(gds.all,"$dosage")
  Genotype2=cbind.data.frame(id,genotype2)
  Genotype2=merge(phenotype,Genotype2,by="id")
  genotype2=as.matrix(Genotype2[,6:ncol(Genotype2)])
  chr2 = seqGetData(gds.all,"chromosome")[1]
  pos2i = seqGetData(gds.all,"position")[1]
  pos2f = seqGetData(gds.all,"position")[length(seqGetData(gds.all,"position"))]
  seqResetFilter(gds.all,verbose=FALSE)
  seqClose(gds.all)
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
  non.cov.ind1=apply(subset(genotype1,T,select=aaf1>0&aaf1<1),1,function(i) all(is.na(i)))
  non.cov.ind2=apply(subset(genotype2,T,select=aaf2>0&aaf2<1),1,function(i) all(is.na(i)))
  
  inds.include.2.cases = !is.na(phenolist) & !non.cov.ind1 & !non.cov.ind2 & phenolist==2
  
  if(sum(inds.include.2.cases)>0 & gene[i][1]!=gene.all[j][1]){
    
    
    #Recode in dom mode
    genoDom1<-subset(genotype1,inds.include.2.cases,select=aaf1>0&aaf1<1)
    dom1 = as.numeric(apply(genoDom1, 1, function(i) sum(i, na.rm=T)) > 0)
    nvardom1=ncol(genoDom1)
    
    genoDom2<-subset(genotype2,inds.include.2.cases,select=aaf2>0&aaf2<1)
    dom2 = as.numeric(apply(genoDom2, 1, function(i) sum(i, na.rm=T)) > 0)
    nvardom2=ncol(genoDom2)
    
    #Recode in rec mode
    xyz1=genoDom1
    xyz1[which(xyz1==1)]=0
    xyz1=xyz1/2
    hmz1 = as.numeric(apply(xyz1, 2, function(i) sum(i, na.rm=T)) > 0)
    genoRec1<-subset(genoDom1,select=hmz1>0)
    xx1=as.matrix(genoRec1)
    xx1[which(xx1==1)]=0
    xx1=xx1/2
    rec1 = as.numeric(apply(xx1, 1, function(i) sum(i, na.rm=T)) > 0)
    nvarrec1=ncol(genoRec1)
    
    xyz2=genoDom2
    xyz2[which(xyz2==1)]=0
    xyz2=xyz2/2
    hmz2 = as.numeric(apply(xyz2, 2, function(i) sum(i, na.rm=T)) > 0)
    genoRec2<-subset(genoDom2,select=hmz2>0)
    xx2=as.matrix(genoRec2)
    xx2[which(xx2==1)]=0
    xx2=xx2/2
    rec2 = as.numeric(apply(xx2, 1, function(i) sum(i, na.rm=T)) > 0)
    nvarrec2=ncol(genoRec2)
    
    ## dom vs rec
    
    if(nvardom1>0 & (1 %in% names(table(dom1))) & nvarrec2>0 & (1 %in% names(table(rec2)))){
      dom1 = as.numeric(apply(genoDom1, 1, function(i) sum(i, na.rm=T)) > 0)
      m1.dom1 = sum(dom1)
      n1.dom1 = length(which(dom1 == 0))
      m1.rec2 = sum(rec2)
      n1.rec2 = length(which(rec2 == 0))
      
      # inter
      dom.rec12 = dom1*rec2
      m1.dom.rec12 = sum(dom.rec12)
      
      #choose as the first gene the one with most proportion of carriers
      prop.rec2 = m1.rec2/n1.rec2
      prop.dom1 = m1.dom1/n1.dom1
      if(prop.rec2>=prop.dom1){
        first.gene = rec2
        second.gene = dom1
      }else{
        first.gene = dom1
        second.gene = rec2
      }
      
      
      mod.null.dom =  glm(first.gene ~ ., family = binomial,data=covar[inds.include.2.cases,])
      mod.full.dom = glm(first.gene ~ .+second.gene, family = binomial,data=covar[inds.include.2.cases,])
      
      mod.null.dom.noAdj =  glm(first.gene ~ 1, family = binomial,data=covar[inds.include.2.cases,])
      mod.full.dom.noAdj = glm(first.gene ~ second.gene, family = binomial,data=covar[inds.include.2.cases,])
      
      
      if(( mod.full.dom$converged ) & ( ! mod.full.dom$boundary ) & m1.dom1>=1 & m1.rec2>=1 &
           ( mod.full.dom.noAdj$converged ) & ( ! mod.full.dom.noAdj$boundary )){
        lrt.dom = mod.null.dom$deviance - mod.full.dom$deviance
        pv.dom.rec = pchisq(lrt.dom, df = 1, lower.tail = F)
        
        lrt.dom.noAdj = mod.null.dom.noAdj$deviance - mod.full.dom.noAdj$deviance
        pw.dom.rec.noAdj = pchisq(lrt.dom.noAdj, df = 1, lower.tail = F)
        or_0pc_dom<-as.vector(exp(coefficients(mod.full.dom.noAdj)["second.gene"])) 
        or_3pc_dom<-as.vector(exp(coefficients(mod.full.dom)["second.gene"]))
      } 
      
    } 
    
    ## dom vs dom
    
    if(nvardom1>0 & (1 %in% names(table(dom1))) & nvardom2>0 & (1 %in% names(table(dom2)))){
      dom1 = as.numeric(apply(genoDom1, 1, function(i) sum(i, na.rm=T)) > 0)
      m1.dom1 = sum(dom1)
      n1.dom1 = length(which(dom1 == 0))
      m1.dom2 = sum(dom2)
      n1.dom2 = length(which(dom2 == 0))
      
      # inter
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
      
      
      if(( mod.full.dom$converged ) & ( ! mod.full.dom$boundary ) & m1.dom1>=1 & m1.dom2>=1 &
           ( mod.full.dom.noAdj$converged ) & ( ! mod.full.dom.noAdj$boundary )){
        lrt.dom = mod.null.dom$deviance - mod.full.dom$deviance
        pv.dom.dom = pchisq(lrt.dom, df = 1, lower.tail = F)
        
        lrt.dom.noAdj = mod.null.dom.noAdj$deviance - mod.full.dom.noAdj$deviance
        pw.dom.dom.noAdj = pchisq(lrt.dom.noAdj, df = 1, lower.tail = F)
        or_0pc_dom.dom<-as.vector(exp(coefficients(mod.full.dom.noAdj)["second.gene"])) 
        or_3pc_dom.dom<-as.vector(exp(coefficients(mod.full.dom)["second.gene"]))
      } 
      
    } 
    
    ### rec vs rec
    
    
    if(nvarrec1>0 & (1 %in% names(table(rec1))) & nvarrec2>0 & (1 %in% names(table(rec2)))){
      
      m1.rec1 = sum(rec1)
      n1.rec1 = length(which(rec1 == 0))
      m1.rec2 = sum(rec2)
      n1.rec2 = length(which(rec2 == 0))
      
      # inter
      rec.rec12 = rec1*rec2
      m1.rec.rec12 = sum(rec.rec12)
      
      #choose as the first gene the one with most proportion of carriers
      prop.rec2 = m1.rec2/n1.rec2
      prop.rec1 = m1.rec1/n1.rec1
      if(prop.rec2>=prop.rec1){
        first.gene = rec2
        second.gene = rec1
      }else{
        first.gene = rec1
        second.gene = rec2
      }
      
      
      mod.null.dom =  glm(first.gene ~ ., family = binomial,data=covar[inds.include.2.cases,])
      mod.full.dom = glm(first.gene ~ .+second.gene, family = binomial,data=covar[inds.include.2.cases,])
      
      mod.null.dom.noAdj =  glm(first.gene ~ 1, family = binomial,data=covar[inds.include.2.cases,])
      mod.full.dom.noAdj = glm(first.gene ~ second.gene, family = binomial,data=covar[inds.include.2.cases,])
      
      
      if(( mod.full.dom$converged ) & ( ! mod.full.dom$boundary ) & m1.rec1>=1 & m1.rec2>=1 &
           ( mod.full.dom.noAdj$converged ) & ( ! mod.full.dom.noAdj$boundary )){
        lrt.dom = mod.null.dom$deviance - mod.full.dom$deviance
        pv.rec.rec = pchisq(lrt.dom, df = 1, lower.tail = F)
        
        lrt.dom.noAdj = mod.null.dom.noAdj$deviance - mod.full.dom.noAdj$deviance
        pw.rec.rec.noAdj = pchisq(lrt.dom.noAdj, df = 1, lower.tail = F)
        or_0pc_rec.rec<-as.vector(exp(coefficients(mod.full.dom.noAdj)["second.gene"])) 
        or_3pc_rec.rec<-as.vector(exp(coefficients(mod.full.dom)["second.gene"]))
      } 
      
    } 
    
    
    
  }else{
    
  }
  
  res<-cbind(genename1= gene[i][1], genename2= gene.all[j][1], LD= dist,
             sample.size1 = length(dom1),sample.size2 = length(rec2),Ncarriers.dom.rec12.cases = m1.dom.rec12 , Nnoncarriers.dom1.cases=n1.dom1, Ncarriers.cases.dom1 = m1.dom1,Nnoncarriers.rec2.cases=n1.rec2, Ncarriers.cases.rec2 = m1.rec2,
             pvalueLRT.dom.rec = pv.dom.rec, pvalueLRT.dom.rec.noAdj = pw.dom.rec.noAdj, OR.dom.rec.adj = or_3pc_dom.rec, OR.dom.rec = or_0pc_dom.rec,
             sample.size1 = length(dom1),sample.size2 = length(dom2),Ncarriers.dom.dom12.cases = m1.dom.dom12 , Nnoncarriers.dom1.cases=n1.dom1, Ncarriers.cases.dom1 = m1.dom1,Nnoncarriers.dom2.cases=n1.dom2, Ncarriers.cases.dom2 = m1.dom2,
             pvalueLRT.dom.dom = pv.dom.dom, pvalueLRT.dom.dom.noAdj = pw.dom.dom.noAdj, OR.dom.dom.adj = or_3pc_dom.dom, OR.dom.dom = or_0pc_dom.dom,
             sample.size1 = length(rec1),sample.size2 = length(rec2),Ncarriers.rec.rec12.cases = m1.rec.rec12 , Nnoncarriers.rec1.cases=n1.rec1, Ncarriers.cases.rec1 = m1.rec1,Nnoncarriers.rec2.cases=n1.rec2, Ncarriers.cases.rec2 = m1.rec2,
             pvalueLRT.rec.rec = pv.rec.rec, pvalueLRT.rec.rec.noAdj = pw.rec.rec.noAdj, OR.rec.rec.adj = or_3pc_rec.rec, OR.rec.rec = or_0pc_rec.rec)
             
  
  return(res)
  
}


# index es una funcion que asocia para cada gen las posiciones de los variantes includios en él
# Las variantes tenidas en cuenta fueron ya filtradas en el paso previo
index1<-function(i){
  ind<-intersect(which(gene==i),variant1)
  return(ind)
}

index2<-function(i){
  ind<-intersect(which(gene.all==i),variant2)
  return(ind)
}


# Esta es la funcion a la que le paso todos los argumentos y que va devolver el archivo final
# con todos los pvalores para modelos dominante y recesivo, ajustado o no.
analysis<-function(outfile=NULL,Gnomad_AF1=0.001,Gnomad_AF2=0.001,KGPAF1=1,KGPAF2=1,CADD=FALSE,set=c("frameshift","stop_gained","stop_lost","start_lost","splice_donor","splice_acceptor"),
                   pheno=NULL,covariates=NULL,set_genes1=c("IFNG"),set_genes2=c("IL12RB1"),ncore=4,feno=NULL){
  if(is.null(pheno)) stop("phenotype is missing")
  if(is.null(outfile)) outfile<-"Rcast.out"
  if(!exists("Gnomad")|!exists("EFFECT")|!exists("gene")|!exists("WARN")) {
    gds<-seqOpen(gdsfile)
    if(!exists("Gnomad")) Gnomad<-seqGetData(gds,"annotation/info/Gnomad_AF_UF")$data
    if(!exists("EFFECT")) EFFECT<-seqGetData(gds,"annotation/info/EFFECT")
    if(!exists("gene")) gene<-seqGetData(gds,"annotation/info/GENE")
    #if(!exists("Cadd")) Cadd<-seqGetData(gds,"annotation/info/CADD")
    if(!exists("WARN")) WARN<-seqGetData(gds,"annotation/info/WARN")
    seqClose(gds)
  }
  
  
  myGenes=data.frame(gene)
  
  msc=read.table(paste(dir,"/vcf/MSC.txt",sep=""),header=T)  
  MSC=apply(myGenes,1,function(i) {ifelse(i[1] %in% msc$Gene,msc[msc$Gene %in% i[1],]$MSC,0)})
  
  
  
  msc.all=read.table(paste(dir,"/vcf/MSC_by_gene.txt",sep=""),header=F)  
  MSC.all=msc.all[,1]
  MSC.all = as.numeric(as.vector(MSC.all))
  
  #Variant selection
  # KGP_AF=0.0001;ExAC_AF=0.00001;AAF=1; set=c("frameshift","stop_gained","stop_lost","start_lost","splice_donor","splice_acceptor")
  if(CADD==TRUE & set_genes1[1]=="ALL" & set_genes2[1]=="ALL"){
    variant1<-which(Gnomad<=Gnomad_AF1 & GnomadAFR<=Gnomad_AF1 & GnomadNFE<=Gnomad_AF1 & GnomadAMR<=Gnomad_AF1
                    & GnomadASJ<=Gnomad_AF1 & GnomadEAS<=Gnomad_AF1 & GnomadFIN<=Gnomad_AF1 & GnomadSAS<=Gnomad_AF1 &
                      KGP<=KGPAF1 &
                      Cadd>=MSC &
                      !grepl("WARN",WARN) &
                      grepl(paste(set,collapse="|"),EFFECT))
    g2test1<-names(table(gene[variant1]))
    
    variant2<-which(Gnomad.all<=Gnomad_AF2 & GnomadAFR.all<=Gnomad_AF2 & GnomadNFE.all<=Gnomad_AF2 & GnomadAMR.all<=Gnomad_AF2
                    & GnomadASJ.all<=Gnomad_AF2 & GnomadEAS.all<=Gnomad_AF2 & GnomadFIN.all<=Gnomad_AF2 & GnomadSAS.all<=Gnomad_AF2 &
                      KGP.all<=KGPAF2 &
                      Cadd.all>=MSC.all &
                      !grepl("WARN",WARN.all) &
                      grepl(paste(set,collapse="|"),EFFECT.all))
    g2test2<-names(table(gene.all[variant2]))
  }else if(CADD==FALSE & set_genes1[1]=="ALL" & set_genes2[1]=="ALL"){
    variant1<-which(Gnomad<=Gnomad_AF1 & GnomadAFR<=Gnomad_AF1 & GnomadNFE<=Gnomad_AF1 & GnomadAMR<=Gnomad_AF1
                    & GnomadASJ<=Gnomad_AF1 & GnomadEAS<=Gnomad_AF1 & GnomadFIN<=Gnomad_AF1 & GnomadSAS<=Gnomad_AF1 &
                      KGP<=KGPAF1 &
                      !grepl("WARN",WARN) &
                      grepl(paste(set,collapse="|"),EFFECT))
    g2test1<-names(table(gene[variant1]))
    
    variant2<-which(Gnomad.all<=Gnomad_AF2 & GnomadAFR.all<=Gnomad_AF2 & GnomadNFE.all<=Gnomad_AF2 & GnomadAMR.all<=Gnomad_AF2
                    & GnomadASJ.all<=Gnomad_AF2 & GnomadEAS.all<=Gnomad_AF2 & GnomadFIN.all<=Gnomad_AF2 & GnomadSAS.all<=Gnomad_AF2 &
                      KGP.all<=KGPAF2 &
                      !grepl("WARN",WARN.all) &
                      grepl(paste(set,collapse="|"),EFFECT.all))
    g2test2<-names(table(gene.all[variant2]))
  }else if(CADD==TRUE & set_genes1[1]!="ALL" & set_genes2[1]=="ALL"){
    variant1<-which(Gnomad<=Gnomad_AF1 & GnomadAFR<=Gnomad_AF1 & GnomadNFE<=Gnomad_AF1 & GnomadAMR<=Gnomad_AF1
                    & GnomadASJ<=Gnomad_AF1 & GnomadEAS<=Gnomad_AF1 & GnomadFIN<=Gnomad_AF1 & GnomadSAS<=Gnomad_AF1 &
                      KGP<=KGPAF1 &
                      Cadd>=MSC &
                      !grepl("WARN",WARN) &
                      grepl(paste(set,collapse="|"),EFFECT) &
                      grepl(paste(paste("^",set_genes1,"$",sep=""),collapse="|"),gene))
    g2test1<-names(table(gene[variant1]))
    
    variant2<-which(Gnomad.all<=Gnomad_AF2 & GnomadAFR.all<=Gnomad_AF2 & GnomadNFE.all<=Gnomad_AF2 & GnomadAMR.all<=Gnomad_AF2
                    & GnomadASJ.all<=Gnomad_AF2 & GnomadEAS.all<=Gnomad_AF2 & GnomadFIN.all<=Gnomad_AF2 & GnomadSAS.all<=Gnomad_AF2 &
                      KGP.all<=KGPAF2 &
                      Cadd.all>=MSC.all &
                      !grepl("WARN",WARN.all) &
                      grepl(paste(set,collapse="|"),EFFECT.all))
    g2test2<-names(table(gene.all[variant2]))
    
  }else if(CADD==FALSE & set_genes1[1]!="ALL" & set_genes2[1]=="ALL"){
    variant1<-which(Gnomad<=Gnomad_AF1 & GnomadAFR<=Gnomad_AF1 & GnomadNFE<=Gnomad_AF1 & GnomadAMR<=Gnomad_AF1
                    & GnomadASJ<=Gnomad_AF1 & GnomadEAS<=Gnomad_AF1 & GnomadFIN<=Gnomad_AF1 & GnomadSAS<=Gnomad_AF1 &
                      KGP<=KGPAF1 &
                      !grepl("WARN",WARN) &
                      grepl(paste(set,collapse="|"),EFFECT) &
                      grepl(paste(paste("^",set_genes1,"$",sep=""),collapse="|"),gene))
    g2test1<-names(table(gene[variant1]))
    
    variant2<-which(Gnomad.all<=Gnomad_AF2 & GnomadAFR.all<=Gnomad_AF2 & GnomadNFE.all<=Gnomad_AF2 & GnomadAMR.all<=Gnomad_AF2
                    & GnomadASJ.all<=Gnomad_AF2 & GnomadEAS.all<=Gnomad_AF2 & GnomadFIN.all<=Gnomad_AF2 & GnomadSAS.all<=Gnomad_AF2 &
                      KGP.all<=KGPAF2 &
                      !grepl("WARN",WARN.all) &
                      grepl(paste(set,collapse="|"),EFFECT.all))
    g2test2<-names(table(gene.all[variant2]))
    
  }else if(CADD==TRUE & set_genes1[1]=="ALL" & set_genes2[1]!="ALL"){
    variant1<-which(Gnomad<=Gnomad_AF1 & GnomadAFR<=Gnomad_AF1 & GnomadNFE<=Gnomad_AF1 & GnomadAMR<=Gnomad_AF1
                    & GnomadASJ<=Gnomad_AF1 & GnomadEAS<=Gnomad_AF1 & GnomadFIN<=Gnomad_AF1 & GnomadSAS<=Gnomad_AF1 &
                      KGP<=KGPAF1 &
                      Cadd>=MSC &
                      !grepl("WARN",WARN) &
                      grepl(paste(set,collapse="|"),EFFECT))
    g2test1<-names(table(gene[variant1]))
    
    variant2<-which(Gnomad.all<=Gnomad_AF2 & GnomadAFR.all<=Gnomad_AF2 & GnomadNFE.all<=Gnomad_AF2 & GnomadAMR.all<=Gnomad_AF2
                    & GnomadASJ.all<=Gnomad_AF2 & GnomadEAS.all<=Gnomad_AF2 & GnomadFIN.all<=Gnomad_AF2 & GnomadSAS.all<=Gnomad_AF2 &
                      KGP.all<=KGPAF2 &
                      Cadd.all>=MSC.all &
                      !grepl("WARN",WARN.all) &
                      grepl(paste(set,collapse="|"),EFFECT.all) &
                      grepl(paste(paste("^",set_genes2,"$",sep=""),collapse="|"),gene.all))
    g2test2<-names(table(gene.all[variant2]))
  }else if(CADD==FALSE & set_genes1[1]=="ALL" & set_genes2[1]!="ALL"){
    variant1<-which(Gnomad<=Gnomad_AF1 & GnomadAFR<=Gnomad_AF1 & GnomadNFE<=Gnomad_AF1 & GnomadAMR<=Gnomad_AF1
                    & GnomadASJ<=Gnomad_AF1 & GnomadEAS<=Gnomad_AF1 & GnomadFIN<=Gnomad_AF1 & GnomadSAS<=Gnomad_AF1 &
                      KGP<=KGPAF1 &
                      !grepl("WARN",WARN) &
                      grepl(paste(set,collapse="|"),EFFECT))
    g2test1<-names(table(gene[variant1]))
    
    variant2<-which(Gnomad.all<=Gnomad_AF2 & GnomadAFR.all<=Gnomad_AF2 & GnomadNFE.all<=Gnomad_AF2 & GnomadAMR.all<=Gnomad_AF2
                    & GnomadASJ.all<=Gnomad_AF2 & GnomadEAS.all<=Gnomad_AF2 & GnomadFIN.all<=Gnomad_AF2 & GnomadSAS.all<=Gnomad_AF2 &
                      KGP.all<=KGPAF2 &
                      !grepl("WARN",WARN.all) &
                      grepl(paste(set,collapse="|"),EFFECT.all) &
                      grepl(paste(paste("^",set_genes2,"$",sep=""),collapse="|"),gene.all))
    g2test2<-names(table(gene.all[variant2]))
  }else if(CADD==TRUE & set_genes1[1]!="ALL" & set_genes2[1]!="ALL"){
    variant1<-which(Gnomad<=Gnomad_AF1 & GnomadAFR<=Gnomad_AF1 & GnomadNFE<=Gnomad_AF1 & GnomadAMR<=Gnomad_AF1
                    & GnomadASJ<=Gnomad_AF1 & GnomadEAS<=Gnomad_AF1 & GnomadFIN<=Gnomad_AF1 & GnomadSAS<=Gnomad_AF1 &
                      KGP<=KGPAF1 &
                      Cadd>=MSC &
                      !grepl("WARN",WARN) &
                      grepl(paste(set,collapse="|"),EFFECT) &
                      grepl(paste(paste("^",set_genes1,"$",sep=""),collapse="|"),gene))
    g2test1<-names(table(gene[variant1]))
    
    variant2<-which(Gnomad.all<=Gnomad_AF2 & GnomadAFR.all<=Gnomad_AF2 & GnomadNFE.all<=Gnomad_AF2 & GnomadAMR.all<=Gnomad_AF2
                    & GnomadASJ.all<=Gnomad_AF2 & GnomadEAS.all<=Gnomad_AF2 & GnomadFIN.all<=Gnomad_AF2 & GnomadSAS.all<=Gnomad_AF2 &
                      KGP.all<=KGPAF2 &
                      Cadd.all>=MSC.all &
                      !grepl("WARN",WARN.all) &
                      grepl(paste(set,collapse="|"),EFFECT.all) &
                      grepl(paste(paste("^",set_genes2,"$",sep=""),collapse="|"),gene.all))
    g2test2<-names(table(gene.all[variant2]))
  }else{
    variant1<-which(Gnomad<=Gnomad_AF1 & GnomadAFR<=Gnomad_AF1 & GnomadNFE<=Gnomad_AF1 & GnomadAMR<=Gnomad_AF1
                    & GnomadASJ<=Gnomad_AF1 & GnomadEAS<=Gnomad_AF1 & GnomadFIN<=Gnomad_AF1 & GnomadSAS<=Gnomad_AF1 &
                      KGP<=KGPAF1 &
                      !grepl("WARN",WARN) &
                      grepl(paste(set,collapse="|"),EFFECT) &
                      grepl(paste(paste("^",set_genes1,"$",sep=""),collapse="|"),gene))
    g2test1<-names(table(gene[variant1]))
    
    variant2<-which(Gnomad.all<=Gnomad_AF2 & GnomadAFR.all<=Gnomad_AF2 & GnomadNFE.all<=Gnomad_AF2 & GnomadAMR.all<=Gnomad_AF2
                    & GnomadASJ.all<=Gnomad_AF2 & GnomadEAS.all<=Gnomad_AF2 & GnomadFIN.all<=Gnomad_AF2 & GnomadSAS.all<=Gnomad_AF2 &
                      KGP.all<=KGPAF2 &
                      !grepl("WARN",WARN.all) &
                      grepl(paste(set,collapse="|"),EFFECT.all) &
                      grepl(paste(paste("^",set_genes2,"$",sep=""),collapse="|"),gene.all))
    g2test2<-names(table(gene.all[variant2]))
  }
  #Get phenotype
  gds <- seqOpen(paste(dir,"/vcf/tmp/",vcf,".gds",sep=""))
  id<-seqGetData(gds,"sample.id")
  phenolist<-phenotype[as.factor(id),pheno]
  seqClose(gds)
  
  #Get covariates
  covar<-cbind(phenotype[as.factor(id),covariates],1)
  
  
  cl<-makeCluster(ncore)
  clusterExport(cl=cl,varlist=c("gene","variant1","variant2","phenotype","phenolist","pheno","id","covar","dir","vcf","vcf.all","gene.all","num"), envir=environment())
  clusterEvalQ(cl, {
    library("SeqArray")
    source(paste(dir,"/scriptR/my_Rfunction_di_NY.r",sep=""))
  })
  registerDoParallel(cl)
  
  ## liste tiene en cada elemento las posiciones de las variantes que pertenecen a cada gen en el orden 
  ## que tienen los genes en g2test
  liste1<-parLapply(cl,g2test1,index1)
  liste2<-parLapply(cl,g2test2,index2)
  
  liste_new=liste2
  
  resultpar<-foreach(i=liste1,.combine=rbind,.packages="SeqArray") %:% foreach(j=liste_new,.combine=rbind,.packages="SeqArray") %dopar% {
    CAST(i,j)}
  cadd=ifelse(CADD==FALSE,"NoCADD","CADD")
  effects=ifelse("missense" %in% set,"MISSLOF","LOF")
  
  write.table(resultpar,file=paste("di/",feno,"/pvalues/di_co_AllvsAll","_Gnomad_AF1",Gnomad_AF1,"_Gnomad_AF2",Gnomad_AF2,"_KGPAF1",KGPAF1,"_KGPAF2",KGPAF2,"_",cadd,"_",effects,"_",num,".out",sep=""),row.names=F,col.names=T,quote=F,sep="\t")
  stopCluster(cl)
  return(resultpar)
}
