# Functional haplotype based GWAS
Selecting closely-linked SNPs based on local epistatic effects for haplotype construction improves power of association mapping

## To run functional haplotype based GWAS


### 1. prepare data  
* genotypic data   
* information of markers   
* phenotypic data   
* kinship matrix or distance matrix   
#### [The detail of data format are shown in the folder 'sample'](https://github.com/Fangv1/Functional_haplotype_GWAS/tree/master/sample) 

* put all af the dataset in the same folder   
For example, put the`FH_GWAS.r`, `sample_genodata.txt`, `sample_genodata_info.txt`, `sample_pheno.txt` and `sample_kinship.txt` in the folder `E:study/data/gwas/`



### 2. run FH_GWAS

*The FH_GWAS is based on the package [BGLR](https://github.com/gdlc/BGLR-R)*

If all the data saved in the folder of `E:study/data/gwas/`  
  ```R
  setwd('E:study/data/gwas/')  
  library(BGLR)
  source('FH_GWAS.r')  ## load the function of FH_GWAS 
  FH_GWAS(path='E:study/data/gwas/',geno='sample_genodata.txt',genoinfo='sample_genodata_info.txt',pheno='sample_pheno.txt',kin='sample_kinship.txt',dis=NULL,out=NULL,windowsize=50000,thr_add=0.05,thr_eps=0.1,N_CPU_CORE=1)
  ```
  ```path``` is the working path for FH_GWAS  
  ```geno``` is the genotypic data  
  ```genoinfo``` is the iformation of markers  
  ```pheno``` is the phenotypic data  
  ```kin``` is the kinship matrix  
  ```dis``` is the distance matrix  
  `out` is the file name for saving the result
  ```windowsize``` is the window size for haplotype, the defined value is 50000bp  
  ```thr_add``` is the threshod of P value for additive effect, the defined value is 0.05  
  ```thr_eps``` is the threshod of P value for epistasis effect, the defined value is 0.1  
  `N_CPU_CORE` is the number of cup used for FH_GWAS, the defined value is 1.   
  *in the windows system, `N_CPU_CORE` can be just 1. in the linux, the `N_CPU_CORE` can be more than 1.*
  
  If R working path is the same with the folder where the data are saved, the `path='E:study/data/gwas/'` is not necessary.   
  If the parameters (`windowsize`,`thr_add`,`thr_eps`,`N_CPU_CORE`) use the defined one, they can be omited, so the code can be wrote as follows:
  ```R
  FH_GWAS(geno='sample_genodata.txt',genoinfo='sample_genodata_info.txt',pheno='sample_pheno.txt',kin='sample_kinship.txt')  
  ```
 
 
 ## 3. results   
 The outputs of FH_GWAS will be saved in two files.  
 When the `out='test'`, then one is `test_FH_GWAS.txt` and another is `test_hap_number.txt`. However, when the defined `out=NULL` is used, then the outfile will use the name in the phenotypic data, taking the sample data for example, in this case, the outfile are  `flowering_time_FH_GWAS.txt` and `flowering_time_hap_number.txt`.
 
 #### *  XX_FH_GWAS.txt  
 ```R
 SNP1	SNP2	SNP3	chr	pos	Pvalue_HAP	pvalue_snp1	pvalue_snp2	pvalue_snp3	pvalue_snp12	pvalue_snp13	pvalue_snp23	pvalue_snp123
SNP2692135chr1	SNP2697161chr1	SNP2712514chr1	chr1	2700603	0.0722437259742261	0.0335043010881686	0.0394849514380679	0.0356463351104863	0.0201170509171901	0.023280751078962	NA	NA
SNP2692135chr1	SNP2697161chr1	SNP2712657chr1	chr1	2700651	0.0722437259742261	0.0335043010881686	0.0394849514380679	0.0356463351104863	0.0201170509171901	0.023280751078962	NA	NA
 ````  
 The name of the three SNPs, chromosome, mean position of the three SNPs, P value of haplotype, three P value of SNPs and three P value of pairwise epistasis and P value of three-term epistasis.
 
  ####  *  XX_hap_number.txt 
  ```R
  SNP	combination_num
  ```
The name of SNP and the number of combination of triple SNPs within the windowsize. The sum of this number is suggested for the Bonferroni correction.

---------------------------
### SNP_FH_GWAS.r
SNP_FH_GWAS.r is the code for SNP_based GWAS and functional_haplotype_based GWAS for Arabidopsis dataset.

The SNP_based GWAS code is based three data: phenotypic data (`FT10_phenotype.txt`), genotypic data (`genodata.txt`), file of information of markers (`genodata_info.txt`) and kinship matrix (`all_chr_1003_RD`). 

To run functional_haplotype_based GWAS, another two files are needed for constructing haplotype. The GEN file (`hapdata.txt`) and SAMPLE file (`hapdata_info.txt`) can be obtained by software [SHAPEIT](https://mathgen.stats.ox.ac.uk/genetics_software/shapeit/shapeit.html). The results of SNP_based GWAS (`snp_gwas_out.txt`) is used for filter the SNPs.
 
#### [Subsets of the original data set (The 1001 Genome Consortium 2016) for running the code can be found in the folder 'data' and the format of all the files are described detailedly in README.md under folder 'data'](https://github.com/Fangv1/Functional_haplotype_GWAS/tree/master/data/).


### simulation.r
simulation.r is the code of simulation study.

This code is based on genotypic data (`genodata.txt`), Rogers' ditance matrix  (`all_chr_1003_RD`), haplotype phases file (`hapdata.txt`), haplotype information file (`hapdata_info.txt`) and the results of SNP_based GWAS (`snp_gwas_out.txt`).


