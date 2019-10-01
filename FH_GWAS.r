#----------------------------------------------------------------------------
# functional haplotype association mapping(using mixed model, kinship matrix) 
#----------------------------------------------------------------------------
#Rscript FH_GWAS.R /qg-10/data/AGR-QG/liuf/Yongyu 128_6_rowed_GBS_filter_18584_snps.hmp.txt traits.txt 10
##"C:\Program Files\R\R-3.4.3\bin\Rscript.exe"

FH_GWAS<-function(path=NULL,geno=NULL,genoinfo=NULL,pheno=NULL,kin=NULL,out=NULL,dis=NULL,windowsize=50000,len=3,thr_add=0.05,thr_eps=0.1,N_CPU_CORE=1){

  
  library(parallel)
  tryCatch(library(BGLR),error = function(err) {
    cat("There is no package 'BGLR', please install it with 'install.packages('BGLR')' ",'\n')
    q()
  })
  
  if(!is.null(path)){setwd(path)}
  
  if(is.null(pheno) | is.null(geno) | is.null(genoinfo)){
    print('genotypic data or phenotypic data or marker infomation is missing')
    q()
  }else{
    
    ## for the kinship 
    if(is.null(kin)){
      if(is.null(dis)){
        print('There is no kinship or distance matrix');q()
      }else{
        RD<-as.matrix(read.table(dis,header=T,row.names=1,stringsAsFactors=F))
        K <- 2*(1-RD)
      }
    }else{
      K<-as.matrix(read.table(kin,header=T,row.names=1,stringsAsFactors=F))
    }
    K[1:5,1:5]
    ## phenotypic data
    cat ('phenotyppic data is',pheno,'\n','       Load the phenotypic data ......','\n')
    pheno_file<-read.table(pheno,header=T,stringsAsFactors=F)
    ## outfile name
    if(is.null(out)){
      outfile<-colnames(pheno_file)[2]
    }else{
      outfile<-out
    }
    cat('The results will be saved in the file:', sprintf('%s_FH_GWAS.txt',outfile),'\n')
    colnames(pheno_file)[1:2]<-c('Genotype','Phenotype')
    #print(head(pheno_file))
    
    ## genotypica data
    cat('genotypic data is',geno,'\n','       Load the genotypic data ......','\n')
    genodata<-read.table(geno,header=T,stringsAsFactors=F)
    colnames(genodata)[1]<-"Genotype"
    
    ## information of makers
    cat('        Load the marker information ......','\n')
    genodata_info<-read.table(genoinfo,header=T,row.names = 1,stringsAsFactors=F)
    genodata_info<-genodata_info[order(genodata_info[,3],genodata_info[,4]),]
    cat('The number of marker is:',nrow(genodata_info),'\n')
    
    DATA <- merge(pheno_file, genodata, by.x = "Genotype", by.y = "Genotype",sort=T)
    ngen<-nrow(DATA)
    
    if(ngen<1){
      cat('There is no common genotype between phenotypic data and genotypic data, Please check it again.','\n')
      q()
    }
    cat(ngen, 'genotypes and', nrow(genodata_info), 'markers are used for functional haplotype GWAS analysis','\n')
    
    ## match the genotype in kinship and genotypic data
    r_ind<-match(DATA$Genotype,rownames(K))
    if(length(r_ind)<ngen){
      cat('The name of genotypes in kinship matrix is different with those in genotypic data, Please check it again.','\n');q()
      }
    kinship<-K[r_ind,r_ind]
    print(all(rownames(kinship)==DATA$Genotype))
    rownames(kinship)<-colnames(kinship)<-DATA$Genotype
    DATA$Genotype<-factor(DATA$Genotype)
    
    
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
    I<-diag(1,ngen,ngen)
    decomposition<-eigen(kinship)
    D<-diag(decomposition$values)
    
    TU<-t(decomposition$vectors)
    inv_V<-solve(sqrt(D*var.gen + I*var.err))
    coef<-inv_V %*% TU
    print(var.gen)
    print(var.err)
    
    ############################################# functional haplotype based GWAS ##########################################
    print('prepare the data for functional haplotype')
    onlyDATA <-data.frame(Y=coef %*% traitY,
                          intercept=coef %*% rep(1,ngen))
    
    marker_names<-rownames(genodata_info)
    marker_ind<-1:nrow(genodata_info)
    
    ### prepare the data for GWAS
    onlyDATA <-data.frame(Y=coef %*% traitY,
                          intercept=coef %*% rep(1,ngen))
    alleles<-read.table('hap_alleles.txt')
    ### Initialize the results for haplotype
    start<-data.frame('SNP1','SNP2','SNP3','chr','pos','Pvalue_HAP','pvalue_snp1','pvalue_snp2','pvalue_snp3','pvalue_snp12','pvalue_snp13','pvalue_snp23','pvalue_snp123')
    write.table(start,file=sprintf('%s_FH_GWAS.txt',outfile),quote=F,row.names=F,col.names=F,sep='\t')
    start1<-data.frame('SNP','combination_num')
    write.table(start1,file=sprintf('%s_hap_number.txt',outfile),quote=F,row.names=F,col.names=F,sep='\t')
    
    
    
    # function for haplotype based GWAS
    hlm1<-lm(Y ~ -1 +intercept ,data=onlyDATA)
    HAP_MLM<-function(m){
      start_pos<-genodata_info[m,4]-1
      end_pos<-genodata_info[m,4]+windowsize
      chr<-genodata_info[m,3]
      ind_p<-which(genodata_info[,4]<end_pos & genodata_info[,4]>start_pos & genodata_info[,3]==chr)
      snp_number<-length(ind_p)
      
      if(snp_number>(len-1)){
        snp_p<-rownames(genodata_info)[ind_p]
        if(m==1){
          snp_f<-apply(combn(snp_p,3),2,paste,collapse='^')
        }else{ # delete the same haplotype
          n=m-1
          start_posn<-genodata_info[n,4]-1
          end_posn<-genodata_info[n,4]+windowsize
          chrn<-genodata_info[n,3]
          ind_pn<-which(genodata_info[,4]<end_posn & genodata_info[,4]>start_posn & genodata_info[,3]==chrn)
          nn<-setdiff(ind_p,ind_pn)
          if(length(nn)>0){
            if(length(ind_pn)>2){
              comb_m<-apply(combn(snp_p,3),2,paste,collapse='^')
              snp_pn<-rownames(genodata_info)[ind_pn]
              comb_n<-apply(combn(snp_pn,3),2,paste,collapse='^')
              snp_f<-setdiff(comb_m,comb_n)
            }else{
              snp_f<-apply(combn(snp_p,3),2,paste,collapse='^')
            }
          }else{snp_f<-setdiff('a','a')}
        }
        cat('Marker:',m,chr,'\t','snp_number:',snp_number,'\t','hap_number:',length(snp_f),'\n')
        write.table(t(c(rownames(genodata_info)[m],length(snp_f))),file=sprintf('%s_hap_number.txt',outfile),quote=F,row.names=F,col.names=F,sep='\t',append = T)
        
        
        if(length(snp_f)>0){
          for(s in snp_f){
            snp_r<-unlist(strsplit(s, "\\^"))        #select 3 SNP
            snp_g<-DATA[,snp_r]  # genotype of 3 snp
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
              cat('warning message: something is wrong in calculation epistasis for marker',rownames(genodata_info)[m],'\n')
            }else{
              lm_anova<-anova(slm)
              lm_summary<-data.frame(summary(slm)$coefficients)
              name_summary<-rownames(lm_summary)
              P_snp1<-lm_summary[intersect(c('SNP1','SNP2','SNP3'),name_summary),4]
              P_snp2<-lm_summary[intersect(c('SNP12','SNP13','SNP23'),name_summary),4]
              p_eps<-lm_summary[c('SNP1','SNP2','SNP3','SNP12','SNP13','SNP23','SNP123'),4]
              #cat(snp_r,p_eps,'\n')
              ## judge the effect of snp and epistasis 
              if(sum(P_snp1<thr_add)>1 & sum(P_snp2<thr_eps)>1){ 
                ## at least 2 p value of SNP are smaller than 0.05 and at least 2 p value of epistasis are smaller than 0.1
                hap <- apply(snp_g,1,paste,collapse='')
                design_matrix_haplotype<-as.matrix(alleles[hap,])
                # haplotype data for GWAS
                hapDATA <- data.frame(onlyDATA,coef%*%design_matrix_haplotype)
                Fehler <<- FALSE
                
                tryCatch(
                  hlm2<-lm(paste0("Y ~ -1 +intercept +",paste(colnames(design_matrix_haplotype),collapse=" +")),data=hapDATA),
                  error = function(err) {
                    Fehler <<- TRUE
                  }
                )
                if(Fehler){
                  cat('warning message: something is wrong in calculation functional haplotype',rownames(genodata_info)[m],'\n')
                }else{
                  hap_anova<-anova(hlm1,hlm2)
                  Pvalue_hap<-hap_anova$Pr[2]
                  mean_pos<-round(mean(genodata_info[snp_r,4]),0)
                  result<-data.frame(t(snp_r),chr,mean_pos,t(Pvalue_hap), t(p_eps))
                  write.table(result,file=sprintf('%s_FH_GWAS.txt',outfile),quote=F,row.names=F,col.names=F,sep='\t',append=T)
                }
              }
            }
          }
          return(NULL)
        }
        
      }else{
        cat('Marker:',m,chr,'\t','snp_number:',snp_number,'\n')
        write.table(t(c(rownames(genodata_info)[m],0)),file=sprintf('%s_hap_number.txt',outfile),quote=F,row.names=F,col.names=F,sep='\t',append = T)
      }
    }
    
    if(N_CPU_CORE==1){
      res<-sapply(marker_ind,HAP_MLM)
    }else if(N_CPU_CORE>1){
      res<-mclapply(marker_ind,mc.cores=N_CPU_CORE,FUN=HAP_MLM,mc.preschedule = FALSE)
    }
    
  }
}

