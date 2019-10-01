# Functional haplotype based GWAS
Selecting closely-linked SNPs based on local epistatic effects for haplotype construction improves power of association mapping

# To run functional haplotype based GWAS

## 1. data prepare

#### a. genotypic data (sample_genodata.txt)

Genotypic data is file with n+1 rows and m+1 columns, in which n is the number of individuals and m is the number of markers.
the first column is the name of individuals and the first row is the names of markers.
the genotypic value for each marker should be 0, 1 and 2 (0 for the reference allele, 2 for alternative allele and 1 for heterozygote). please look at the sample data 'sample_genodata.txt' in the floder of 'sample'.

#### b. marker information (sample_genodata_info.txt)

Marker information includes 5 columns, which are the name of markers, reference allele, alternative allele, chromosome and posistion.

#### c. phenotypic data (sample_pheno.txt)
phenotypic data have two colums. the first column is the name of individuals and the second column is the phenotypic value for the corresponded individual.

#### d. kinship matrix or distance matrix


### SNP_FH_GWAS.r
SNP_FH_GWAS.r is the code for SNP_based GWAS and functional_haplotype_based GWAS.

The SNP_based GWAS code is based three data: phenotypic data (FT10_phenotype.txt), genotypic data (genodata.txt), file of information of markers (genodata_info.txt) and kinship matrix (all_chr_1003_RD). 

To run functional_haplotype_based GWAS, another two files are needed for constructing haplotype. The GEN file (hapdata.txt) and SAMPLE file (hapdata_info.txt) can be obtained by software 'SHAPEIT'. The results of SNP_based GWAS (snp_gwas_out.txt) is used for filter the SNPs.
 
Subsets of the original data set (The 1001 Genome Consortium 2016) for running the code can be found in the folder 'data' and the format of all the files are described detailedly in README.md under folder 'data'.


### simulation.r
simulation.r is the code of simulation study.

This code is based on genotypic data (genodata.txt), Rogers' ditance matrix  (all_chr_1003_RD), haplotype phases file (hapdata.txt), haplotype information file (hapdata_info.txt) and the results of SNP_based GWAS (snp_gwas_out.txt).


