# Functional haplotype based GWAS
Selecting closely-linked SNPs based on local epistatic effects for haplotype construction improves power of association mapping

## To run functional haplotype based GWAS

### 1. put all the data in the same folder

For example, firstly, put the `FH_GWAS.r`, `sample_genodata.txt`, `sample_genodata_info.txt`, `sample_pheno.txt`, `sample_kinship.txt` and `hap_alleles.txt` (can be download in the folder 'sample') in the folder `E:study/data/gwas/`

### 2. run FH_GWAS.r  
if all the data saved in the folder of `E:study/data/gwas/`  
  ```R
  setwd('E:study/data/gwas/')  
  source('FH_GWAS.r')
  FH_GWAS(path='E:study/data/gwas/',geno='sample_genodata.txt',genoinfo='sample_genodata_info.txt',pheno='sample_pheno.txt',kin='sample_kinship.txt',out=NULL,dis=NULL,windowsize=50000,len=3,thr_add=0.05,thr_eps=0.1,N_CPU_CORE=1)
  ```
  ```path``` is the working path for FH_GWAS.r   
  ```geno``` is the genotypic data  
  ```genoinfo``` is the iformation of markers  
  ```pheno``` is the phenotypic data  
  ```kin``` is the kinship matrix  
  ```dis``` is the distance matrix  
  ```windowsize``` is the window size for haplotype, the defined value is 50000bp  
  ```len``` is the number of SNP in each haplotype, the defined value is 3  
  ```thr_add``` is the threshod of P value for additive effect, the defined value is 0.05  
  ```thr_eps``` is the threshod of P value for epistasis effect, the defined value is 0.1  
  
  *if R working path is the same with the folder where the data are saved, the path='E:study/data/gwas/' is not necessary.*  
  the defined windowsize is 50kb and 
  `FH_GWAS(geno='sample_genodata.txt',genoinfo='sample_genodata_info.txt',pheno='sample_pheno.txt',kin='sample_kinship.txt' )`  
  #### FH_GWAS(path=NULL,geno=NULL,genoinfo=NULL,pheno=NULL,kin=NULL,out=NULL,dis=NULL,windowsize=50000,len=3,thr_add=0.05,thr_eps=0.1,N_CPU_CORE=1) 


### SNP_FH_GWAS.r
SNP_FH_GWAS.r is the code for SNP_based GWAS and functional_haplotype_based GWAS.

The SNP_based GWAS code is based three data: phenotypic data (FT10_phenotype.txt), genotypic data (genodata.txt), file of information of markers (genodata_info.txt) and kinship matrix (all_chr_1003_RD). 

To run functional_haplotype_based GWAS, another two files are needed for constructing haplotype. The GEN file (hapdata.txt) and SAMPLE file (hapdata_info.txt) can be obtained by software 'SHAPEIT'. The results of SNP_based GWAS (snp_gwas_out.txt) is used for filter the SNPs.
 
Subsets of the original data set (The 1001 Genome Consortium 2016) for running the code can be found in the folder 'data' and the format of all the files are described detailedly in README.md under folder 'data'.


### simulation.r
simulation.r is the code of simulation study.

This code is based on genotypic data (genodata.txt), Rogers' ditance matrix  (all_chr_1003_RD), haplotype phases file (hapdata.txt), haplotype information file (hapdata_info.txt) and the results of SNP_based GWAS (snp_gwas_out.txt).


