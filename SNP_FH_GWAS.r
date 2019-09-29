#----------------------------------------------------------------------------
# haplotype association mapping in every 50kb(using mixed model, kinship matrix) 
#----------------------------------------------------------------------------
# nohup R CMD BATCH 
rm(list=ls())
setwd('/qg-10/data/AGR-QG/liuf/genotype/ClimBar/GWAS/simulation/arabidopsis/sample_data/')
threshold1<-0.05
threshold2<-0.1

#outfile for GWAS results
outfile_snp <- 'snp_gwas_out.txt' ## SNP GWAS results
outfile_haplotype <- 'haplotype_gwas_out.txt' ## haplotype GWAS results

###phenotype data ## first column is genotype and second column is phenotype
pheno_file<-read.table('FT10_phenotype.txt' ,header=T,stringsAsFactors=F)
colnames(pheno_file)<-c("Genotype","Phenotype")
print(head(pheno_file))

############genotype data
genodata<-read.table('genodata.txt',header=T,stringsAsFactors=F)
DATA <- merge(pheno_file, genodata, by.x = "Genotype", by.y = "Genotype",sort=T)
ngen<-nrow(DATA)
genodata_info<-read.table('genodata.txt',header=T,row.names = 1,stringsAsFactors=F)
### Roger_matrix for kinship 
RD_all<-as.matrix(read.table('all_chr_1003_RD',header=T,row.names=1,stringsAsFactors=F))
r_ind<-match(DATA$Genotype,rownames(RD_all))
RD<-RD_all[r_ind,r_ind]
kinship<-(1-RD)

## caculate the coefficient to change the mixed linear model to simple linear model
library(BGLR)
traitY <- DATA$Phenotype
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
print(var.gen)
print(var.err)

############################################# SNP based GWAS ##########################################
### prepare the data for GWAS
onlyDATA <-data.frame(Y=coef %*% traitY,
                      intercept=coef %*% rep(1,ngen))

### Initialize the result
start<-data.frame('Number','Marker','AF','Estimate','P_value','pos','chr')
write.table(start,file=outfile_snp,quote=F,row.names=F,col.names=F,sep='\t')

marker_names<-rownames(genodata_info)
marker_ind<-1:nrow(genodata_info)

LM_gwas<-function(m){
  print(m)
  label<-marker_names[m]
  
  marker<-coef %*% DATA[,m]
  subDATA<-data.frame(onlyDATA,marker)
  lm_marker<-lm(Y ~ -1 +intercept  + marker,data=subDATA)
  P_alle<-mean(DATA[,m])/2
  P_value<-anova(lm_marker)['marker',]$Pr
  Estimate<-summary(lm_marker)$coefficients['marker',1]
  variance<-(summary(lm_marker)$sigma)**2.
  
  gwas<-data.frame(number=m,
                   SNP=label,
                   P_alle=P_alle,
                   Estimate=Estimate,
                   P_value=P_value,
                   #Prgene_value=Prgene_value,
                   #Estimate_rgene=Estimate_rgene
                   pos=genodata_info[m,1],
                   chr=genodata_info[m,4])
  #print(gwas)
  write.table(gwas,file=outfile_snp,quote=F,row.names=F,col.names=F,sep='\t',append=T)
}

## using  multiple CUP
library(parallel)
PARALLELE<-TRUE

._ncore<-as.numeric(system('grep -c processor /proc/cpuinfo',intern=T))
._mem<-as.numeric(system("free -g|sed -n '3p'|sed 's/.* //'",intern=T))/5
cat("Number of threads:\t",N_CPU_CORE<- min(ceiling(._mem), floor(._ncore*.8)),"\n")

## run SNP based GWAS 
if(PARALLELE){
  res<-mclapply(marker_ind,mc.cores=N_CPU_CORE,FUN=LM_gwas)
}else{
  res<-lapply(marker_ind,FUN=LM_gwas)
}

########################################### haplotype based GWAS #######################################
### haplotype data
hapdata<-read.table('hapdata.txt',header=F,stringsAsFactors=F)
hap_sample<-read.table('hapdata_info.txt',header=T,stringsAsFactors=F)
colnames(hapdata)<-c('chr','SNP','pos','ref','var',rep(hap_sample[-1,1],each=2))
h_ind<-match(DATA$Genotype,hap_sample[,1])
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
cat('phenotype and haplotype :', all(DATA$Genotype==rownames(hap_data)[2*(1:nrow(DATA))-1]),'\n')
cat('kinship and haplotype :', all(rownames(kinship)==rownames(hap_data)[2*(1:nrow(DATA))-1]),'\n')

#load GWAS result to filter some SNPs
gwas_result<-read.table(outfile_snp,header=T,stringsAsFactors=F)
dim(gwas_result)
gwas_result<-gwas_result[which(gwas_result$P_value<0.01 & gwas_result$AF>0.05 & gwas_result$AF<0.95),]
gwas_result<-gwas_result[order(gwas_result$chr,gwas_result$pos),]
rownames(gwas_result)<-gwas_result$Marker
dim(gwas_result)

### prepare the data for GWAS
onlyDATA <-data.frame(Y=coef %*% traitY,
                      intercept=coef %*% rep(1,ngen))

### Initialize the results for haplotype
start<-data.frame('Number','SNP1','SNP2','SNP3','Pvalue_HAP','Pvalue_SNP1','Pvalue_SNP2','Pvalue_SNP3','ld12','ld23','ld13','mean','maf1','maf2','maf3','maf_mean','pvalue_snp1','pvalue_snp2','pvalue_snp3','pvalue_snp12','pvalue_snp13','pvalue_snp23','pvalue_snp123')
write.table(start,file=outfile_haplotype,quote=F,row.names=F,col.names=F,sep='\t')

# function for haplotype based GWAS
HAP_LM<-function(m){
  start_pos<-gwas_result$pos[m]-1
  end_pos<-gwas_result$pos[m]+50000
  chr<-gwas_result$chr[m]
  ind_p<-which(gwas_result$pos<end_pos & gwas_result$pos>start_pos & gwas_result$chr==chr)
  snp_number<-length(ind_p)
  
  if(snp_number>2){
    snp_p<-gwas_result$Marker[ind_p]
    comb_m<-apply(combn(snp_p,3),2,function(x){paste(x,collapse = '_')})
    if(m==marker_pos[1]){
      snp_f<-comb_m
    }else{ # delete the same haplotype
      n=m-1
      start_posn<-gwas_result$pos[n]-1
      end_posn<-gwas_result$pos[n]+50000
      chrn<-gwas_result$chr[n]
      ind_pn<-which(gwas_result$pos<end_posn & gwas_result$pos>start_posn & gwas_result$chr==chrn)
      if(length(ind_pn)>2){
        snp_pn<-gwas_result$Marker[ind_pn]
        comb_n<-apply(combn(snp_pn,3),2,function(x){paste(x,collapse = '_')})
        snp_f<-setdiff(comb_m,comb_n)
      }else{
        snp_f<-comb_m
      }
    }
    if(length(snp_f)>0){
      for(s in snp_f){
        snp_r<-unlist(strsplit(s, "_"))        #select 3 SNP
        maf<-gwas_result[snp_r,]$AF
        snp_g<-DATA[,snp_r]  # genotype of 3 snp
        con<-cor(snp_g)
        ld<-con[upper.tri(con)]
        ld<-ld^2
        ld_mean<-mean(ld)
        result_ld<-data.frame(t(ld),ld_mean)
        colnames(snp_g)<-c('SNP1','SNP2','SNP3')
        # genotype of interaction 
        SNP12<-snp_g[,1] * snp_g[,2]  
        SNP13<-snp_g[,1] * snp_g[,3]
        SNP23<-snp_g[,2] * snp_g[,3]
        SNP123<-snp_g[,1] * snp_g[,2] * snp_g[,3]
        snp_all_g<-as.matrix(cbind(snp_g, SNP12, SNP13, SNP23, SNP123))
        # SNP data for GWAS
        subDATA<-data.frame(onlyDATA,coef %*% snp_all_g)
        Fehler <<- FALSE
        tryCatch(
          slm <-lm(Y ~ -1 +intercept  + SNP1 + SNP2 + SNP3 + SNP12 + SNP13 + SNP23 + SNP123, data=subDATA),
          error = function(err) {
            Fehler <<- TRUE
          }
        )
        if(Fehler){
          next
        }else{
          lm_anova<-anova(slm)
          lm_summary<-summary(slm)$coefficients
          name_summary<-rownames(lm_summary)
          P_snp1<-lm_summary[intersect(c('SNP1','SNP2','SNP3'),name_summary),4]
          P_snp2<-lm_summary[intersect(c('SNP12','SNP13','SNP23'),name_summary),4]
          ## judge the effect of snp and epistasis 
          if(sum(P_snp1<threshold1)>1 & sum(P_snp2<threshold2)>1){ 
            ## at least 2 p value of SNP are smaller than 0.05 and at least 2 p value of epistasis are smaller than 0.1
            window<-hap_data[,snp_r]
            window <- apply(window,1,paste,collapse='')
            SNP_hap<-window[2*(1:ngen)-1]
            window <- as.factor(window)
            level <- levels(window) 
            nhap <- length(level) ### Number of haplotypes
            if(nhap<2){
              next
            }else{
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
              
              # haplotype data for GWAS
              hapDATA <- data.frame(onlyDATA,coef%*%design_matrix_haplotype)
              Fehler <<- FALSE
              hlm1<-lm(Y ~ -1 +intercept ,data=hapDATA)
              tryCatch(
                hlm2<-lm(paste0("Y ~ -1 +intercept +",paste(colnames(design_matrix_haplotype),collapse=" +")),data=hapDATA),
                error = function(err) {
                  Fehler <<- TRUE
                }
              )
              if(Fehler){
                next
              }else{
                hap_anova<-anova(hlm1,hlm2)
                Pvalue_hap<-hap_anova$Pr[2]
                Pvalue_snp<-gwas_result[snp_r,]$P_value
                result<-data.frame(m,t(snp_r),t(Pvalue_hap), t(Pvalue_snp),result_ld,t(maf),mean(maf), t(P_snp1), t(P_snp2))
                write.table(result,file=outfile_haplotype,quote=F,row.names=F,col.names=F,sep='\t',append=T)
              }
            }
          }else{
            next
          }
        }
      }
    }
  }
  cat('Marker:',m,'\t',snp_number,'\n')
}

## run haplotype based GWAS 
marker_pos<-1:nrow(gwas_result)
if(PARALLELE){
  res<-mclapply(marker_pos,mc.cores=N_CPU_CORE,FUN=HAP_LM,mc.preschedule = FALSE)
}else{
  res<-lapply(marker_pos,FUN=HAP_LM)
}