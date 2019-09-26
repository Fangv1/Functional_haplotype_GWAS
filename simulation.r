rm(list=ls())
setwd('/qg-10/data/AGR-QG/liuf/genotype/ClimBar/GWAS/simulation/arabidopsis/sample_data/')

############genotype data
genodata<-read.table('genodata.txt',header=T,stringsAsFactors=F)

#load GWAS result to filter some SNPs
gwas_result<-read.table('snp_gwas_out.txt',header=T,stringsAsFactors=F)
rownames(gwas_result)<-gwas_result$Marker
gwas_result<-gwas_result[order(gwas_result$chr,gwas_result$pos),]
maf1<-gwas_result$MAF
ind<-which(maf1>0.5)
maf2<-1-maf1[ind]
maf1[ind]<-maf2
gwas_result$MAF<-maf1
#### create a folder if this folder does not exist
f<-'./simulation/'
if(!dir.exists(f)){system(sprintf('mkdir %s',f))}
f1<-'./simulation/phenotype/'
if(!dir.exists(f)){system(sprintf('mkdir %s',f1))}

###### Define the function 
pick_QTL<-function(m){
  for(p in 1:10000){ # in case 
    snp1<-sample(allsnp,1)
    pos1<-gwas[snp1,'pos']
    chr<-gwas[snp1,'chr']
    pos2<-pos1+50000
    ind_p<-which(gwas$pos<pos2 & gwas$pos>pos1 & gwas$chr==chr)
    Nsnp<-length(ind_p)
    snp_sub<-allsnp[ind_p]
    Fehler <<- FALSE
    if(Nsnp>2){
      for(k in 2:Nsnp){
        snp2<-sample(snp_sub,1)
        snp_sub<-setdiff(snp_sub,snp2)
        snp_g12<-genodata[,c(snp1,snp2)]  # genotype of 3 snp
        con12<-cor(snp_g12)
        ld12<-con12[upper.tri(con12)]^2
        if(ld12>ld_thre1 & ld12<ld_thre2){
          snp3<-sample(snp_sub,1)
          snp_r<-sort(c(snp1,snp2,snp3))
          snp_g<-genodata[,snp_r]
          con<-cor(snp_g)
          ld<-con[upper.tri(con)]
          ld<-ld^2
          if(all(ld>ld_thre1) & all(ld<ld_thre2) ){
            ld_mean<-mean(ld)
            result_ld<-data.frame(t(ld),ld_mean)
            snp<-paste(snp_r,collapse='_')
            maf<-gwas[snp_r,]$MAF
            res<-data.frame(t(snp),t(snp_r),result_ld,t(maf),mean(maf))
            write.table(res,file=outfile,quote=F,row.names=F,col.names=F,sep='\t',append=T)
            #cat(m,'\n')
            cat(m,'\t','MAF:',maf_thre2,'-',maf_thre1,'\t','LD:',ld_thre1,'-',ld_thre2,'\t','rep:',p,'\n')
            Fehler <<- TRUE
            break
          }
        }else{next}
      }
    }else{next}
    if (Fehler){break}
  }
}
## pick_QTL(snp_maf,maf_thre1,maf_thre2,ld_thre1,ld_thre2,10000)


############ parameter
maf_thre<-matrix(c(0.5,0.4,0.3,0.2,0.1,0),3,2,byrow=T)
ld_thre<-matrix(c(0,0.2,0.3,0.6,0.7,1),3,2,byrow=T)


library(parallel)
PARALLELE<-FALSE

._ncore<-as.numeric(system('grep -c processor /proc/cpuinfo',intern=T))
._mem<-as.numeric(system("free -g|sed -n '3p'|sed 's/.* //'",intern=T))/2
cat("Number of threads:\t",N_CPU_CORE<- min(ceiling(._mem), floor(._ncore*.8)),"\n")

##### loop for MAF
for(i in 1:nrow(maf_thre)){
  maf_thre1<-maf_thre[i,1]
  maf_thre2<-maf_thre[i,2]
  ind_maf<-which(gwas_result$MAF<maf_thre1 & gwas_result$MAF>maf_thre2)
  gwas<-gwas_result[ind_maf,]
  print(dim(gwas))
  allsnp<-gwas$Marker
  ######## loop for LD 
  for(j in 1:nrow(ld_thre)){
    ld_thre1<-ld_thre[j,1]
    ld_thre2<-ld_thre[j,2]
    outfile<-sprintf('simulation/pick_3QTL_50kb_2000_maf_%s_%s_ld_%s_%s.txt',10*maf_thre2,10*maf_thre1,10*ld_thre1,10*ld_thre2)
    start<-data.frame('Name','SNP1','SNP2','SNP3','ld12','ld23','ld13','mean','maf1','maf2','maf3','maf_mean')
    write.table(start,file=outfile,quote=F,row.names=F,col.names=F,sep='\t')
    marker_pos<-1:2 ## the number of 2 can be changed
    if(PARALLELE){
      res<-mclapply(marker_pos,mc.cores=N_CPU_CORE,FUN=pick_QTL,mc.preschedule = FALSE)
    }else{
      res<-sapply(marker_pos,FUN=pick_QTL)
    }
  }
}


##################### simulation pheotypic data and association mapping basing on SNP and haplotype ##############

#load GWAS result to filter some SNPs
m<-match(gwas_result$Marker,colnames(genodata))
allgeno<-as.matrix(genodata[,m])
rownames(allgeno)<-genodata$Genotype

### Roger_matrix
RD_all<-as.matrix(read.table('all_chr_1003_RD',header=T,row.names=1,stringsAsFactors=F))
r_ind<-match(genodata$Genotype,rownames(RD_all))
RD<-RD_all[r_ind,r_ind]
kinship<-(1-RD)
library(BGLR)

### haplotype data
hapdata<-read.table('hapdata.txt',header=F,stringsAsFactors=F)
hap_sample<-read.table('hapdata_info.txt',header=T,stringsAsFactors=F)
colnames(hapdata)<-c('chr','SNP','pos','ref','var',rep(hap_sample[-1,1],each=2))
h_ind<-match(rownames(kinship),hap_sample[,1])
hap_ind<-c()
for(i in h_ind){
  j<-2*i-3
  k<-2*i-2
  hap_ind<-append(hap_ind,j)
  hap_ind<-append(hap_ind,k)
}
hap_data<-t(hapdata[,6:ncol(hapdata)])
colnames(hap_data)<-hapdata$SNP
hap_data<-hap_data[hap_ind,]
cat('kinship and haplotype :', all(rownames(kinship)==rownames(hap_data)[2*(1:nrow(kinship))-1]),'\n')

### parameter for simulate phenotypic data
nmar<-ncol(allgeno)
ngen<-nrow(allgeno)
sigma_g<-400
sigma_e<-100
rate_qtl<-0.06
rate_inter2<-0.03
rate_inter3<-0.01
rate_other<-(1- 3*rate_qtl -3*rate_inter2 -rate_inter3)/(nmar-3)
sd_g<-sqrt(sigma_g)
sd_e<-sqrt(sigma_e)
snp_name<-colnames(allgeno)

## function of simulating phenotype and GWAS
simu_phe_gwas<-function(q){
  ######## simulation pheotypic data
  pfile<-gsub('pick_3QTL',paste('phenotype/phe',q,sep=''),infile) ## file for savitng phenotypic data 
  snp_qtl<-QTL[q,]
  snp_g<-allgeno[,as.character(snp_qtl)]  # genotype of 3 QTL
  snp1_snp2<-snp_g[,1] * snp_g[,2]  # genotype of interaction 
  snp1_snp3<-snp_g[,1] * snp_g[,3]
  snp2_snp3<-snp_g[,2] * snp_g[,3]
  snp1_snp2_snp3<-snp_g[,1] * snp_g[,2] * snp_g[,3]
  snp_other<-setdiff(snp_name,snp_qtl)
  snp_all_g<-cbind(snp_g, snp1_snp2, snp1_snp3, snp2_snp3, snp1_snp2_snp3, allgeno[,snp_other]) # genotype of all snp and interaction
  # effect of all snp and interaction
  snp_eff<-c(rnorm(3,0,sqrt(rate_qtl*sigma_g)),
             rnorm(3,0,sqrt(rate_inter2*sigma_g)),
             rnorm(1,0,sqrt(rate_inter3*sigma_g)),
             rnorm((nmar-3),0,sqrt(rate_other*sigma_g))) 
  
  vlue_p<-snp_all_g %*% snp_eff + rnorm(ngen,0,sd_e)
  phenotype<-data.frame(Genotype=rownames(allgeno),Phenotype=vlue_p)
  write.table(phenotype,file=pfile,quote=F,row.names=F,col.names=T,sep='\t')
  
  ## caculate the coefficient to change the mixed linear model to simple linear model      
  traitY <- phenotype$Phenotype
  ETA <- list(list(K=kinship,model="RKHS"))
  fm <- BGLR(y=traitY,
             ETA=ETA,
             nIter=10000,
             burnIn=1000,
             saveAt="BGLR_GBLUP_",
             verbose=FALSE)
  var.gen <- fm$ETA[[1]]$varU
  var.err <- fm$varE
  I<-diag(1,nrow(RD),ncol(RD))
  decomposition<-eigen(kinship)
  D<-diag(decomposition$values)
  
  TU<-t(decomposition$vectors)
  inv_V<-solve(sqrt(D*var.gen + I*var.err))
  coef<-inv_V %*% TU
  #print(var.gen)
  #print(var.err)
  ######### data for association mapping
  onlyDATA <-data.frame(Y=coef %*% traitY,intercept=coef %*% rep(1,ngen))
  ######## SNP data for association mapping
  subDATA<-data.frame(onlyDATA,coef %*% snp_g)                      
  #### association mapping for three SNPs
  marker<-colnames(subDATA)[grep('SNP',colnames(subDATA))]
  result<-q
  for(m in marker){
    ### mixed linear model 
    lm_marker<-lm(as.formula(paste('Y ~ -1 + intercept +',m)), data= subDATA)
    sadj_r<-summary(lm_marker)$adj.r.squared
    P_value<-summary(lm_marker)$coefficients[m,4]
    result<-data.frame(result,
                       marker=m,
                       adj.r=sadj_r,
                       Pvalue=P_value)
    
  }
  
  ######## haplotype data for association mapping
  window<-hap_data[,as.character(snp_qtl)]
  window <- apply(window,1,paste,collapse='')
  SNP_hap<-window[2*(1:ngen)-1]
  window <- as.factor(window)
  level <- levels(window) 
  nhap <- length(level) ### Number of haplotypes
  ### constructe matrix for haplotype
  design_matrix_haplotype <- matrix(0,nrow=ngen,ncol=nhap)
  for (i in 1:ngen){
    chain1 <- window[2*i-1]
    chain2 <- window[2*i]
    pos1 <- which(level == chain1)
    pos2 <- which(level == chain2)
    design_matrix_haplotype[i,pos1] <- design_matrix_haplotype[i,pos1]+1
    design_matrix_haplotype[i,pos2] <- design_matrix_haplotype[i,pos2]+1
  }
  colnames(design_matrix_haplotype) <- paste0('hap',1:nhap)
  
  
  hapDATA <- data.frame(onlyDATA,coef%*%design_matrix_haplotype)
  #### association mapping for haplotype
  hlm1<-lm(Y ~ -1 +intercept ,data=onlyDATA)
  hlm2<-lm(paste0("Y ~ -1 +intercept +",paste(colnames(design_matrix_haplotype),collapse=" +")),data=hapDATA)
  hap_anova<-anova(hlm1,hlm2)
  Pvalue_hap<-hap_anova$Pr[2]
  hadj_r<-summary(hlm2)$adj.r.squared
  result<-data.frame(result,
                     marker=paste(snp_qtl,collapse='_'),
                     adj.r=hadj_r,
                     Pvalue=Pvalue_hap)
  write.table(result,file=outfile_gwas,quote=F,row.names=F,col.names=F,sep='\t',append=T)
  cat('MAF:',maf_thre2,'-',maf_thre1,'\t','LD:',ld_thre1,'-',ld_thre2,'\t','rep:',q,'\n')
}

PARALLELE<-TRUE
##### loop for MAF
for(i in 1:nrow(maf_thre)){
  maf_thre1<-maf_thre[i,1]
  maf_thre2<-maf_thre[i,2]
  
  ######## loop for LD 
  for(j in 1:nrow(ld_thre)){
    ld_thre1<-ld_thre[j,1]
    ld_thre2<-ld_thre[j,2]
    
    ### input QTL file
    infile<-sprintf(sprintf('simulation/pick_3QTL_50kb_2000_maf_%s_%s_ld_%s_%s.txt',10*maf_thre2,10*maf_thre1,10*ld_thre1,10*ld_thre2))
    outfile_gwas<-gsub('pick_3QTL_50kb_2000','simu_phe_gwas_50kb_chr2',infile) ## file for savitng Plaves' result
    qtl<-read.table(infile,stringsAsFactors=F,head=T)
    nm<-which(qtl$ld12==1 | qtl$ld23==1 |qtl$ld13==1)
    if(length(nm)>0){qtl<-qtl[-nm,]}
    print(dim(qtl))
    if(nrow(qtl)<1){next}else{
      name<-unique(qtl$Name)
      n<-match(name,qtl$Name)
      qtl_pos<-qtl[n,]
      cat('total repeat:',nrow(qtl_pos),'\t','MAF:',maf_thre2,'-',maf_thre1,'\t','LD:',ld_thre1,'-',ld_thre2,'\n')
      
      ### 1000 repeat
      QTL<-qtl_pos[,c('SNP1','SNP2','SNP3')]
      q_number<-nrow(QTL)  
      if(q_number<1000){maker_ind<-1:q_number}else{maker_ind<-sample(1:q_number,1000)} # number of repeat
      ### Initialize the result
      start<-data.frame('Number','SNP1','R_square1','P_value1','SNP2','R_square2','P_value2','SNP3','R_square3','P_value3','haplotype','R_square','P_value')
      write.table(start,file=outfile_gwas,quote=F,row.names=F,col.names=F,sep='\t')
      
      ######## use multi Cpu
      if(PARALLELE){
        res<-mclapply(maker_ind,mc.cores=N_CPU_CORE,FUN=simu_phe_gwas,mc.preschedule = FALSE)
      }else{
        res<-lapply(maker_ind,FUN=simu_phe_gwas)
      }
    }
  }
}