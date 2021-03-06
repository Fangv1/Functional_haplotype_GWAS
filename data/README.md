# data for the code of SNP_FH_GWAS.r

### 1. phenotypic data (FT10_phenotype.txt)

phenotypic data is a matrix (1003 X 2) with colnames. 
the first column is the name of genotype and second column is the phenotypic value.   
colnames must be 'Genotype' and 'Phenotype'.


### 2. genotypic data (genodata.txt)

genotypic is a matrix (1003 X m+1, m is the nember of marker) with the colnames.   
colnames must be 'Genotype' + the names of marker. The first column is the names of genotype and the other column is the genotypic value of each marker (0 for the reference allele, 2 for alternative allele and 1 for heterozygote)


### 3. Information of markers (genodata_info.txt)

This file includes 5 column, which are the name of marker, the position of marker, reference allele, alternative allele and chromosome.


### 4. Rogers' ditance matrix (all_chr_1003_RD )

A matrix (1003  X 1003) with colnames and rownames.   
Rogers' ditance (RD) matrix is used for caculating the kinship, kinship= 1 - RD


### 5. GEN file (hapdata.txt) and SAMPLE file (hapdata_info.txt) for the haplotype

The two files are obtained by the software [SHAPEIT](https://mathgen.stats.ox.ac.uk/genetics_software/shapeit/shapeit.html)   
The format of the files described at http://mathgen.stats.ox.ac.uk/genetics_software/shapeit/shapeit.html#gensample


### 6. the results of SNP_based GWAS (snp_gwas_out.txt)
The file includes 7 column and m+1 rows (m is the number of markers)

* 1st column: the number of markers.    
* 2nd column: the name of markers.   
* 3rd column: the allele frequency (AF) of markers.   
* 4th column: the the estimated of markers.  
* 5th column:  the p value of markers.  
* 6th column: the posistion of markers.   
* 7th column: the chromosome of markers.   
the fisrt row is the name of the 7 columns and must be 'Number','Marker','AF','Estimate','P_value','pos','chr'.


