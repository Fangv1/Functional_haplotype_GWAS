# Functional haplotype based GWAS
Selecting closely-linked SNPs based on local epistatic effects for haplotype construction improves power of association mapping

### SNP_FH_GWAS.r
SNP_FH_GWAS.r is the code for SNP_based GWAS and functional_haplotype_based GWAS. 
SNP_based GWAS code is based on the phenotypic data (FT10_phenotype.txt), genotypic data (genodata.txt) and Rogers' ditance matrix  (all_chr_1003_RD using for caculating the kinship). 

Functional_haplotype_based GWAS is based on the phenotypic data (FT10_phenotype.txt), genotypic data (genodata.txt), Rogers' ditance matrix  (all_chr_1003_RD using for caculating the kinship), haplotype phases file (hapdata.txt), haplotype information file (hapdata_info.txt) and the results of SNP_based GWAS (snp_gwas_out.txt)

Subsets of the original data set (The 1001 Genome Consortium 2016) for running the code can be found in the folder 'data'.


### simulation.r

