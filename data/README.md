# data for the code of SNP_FH_GWAS.r

### phenotypic data (FT10_phenotype.txt)

phenotypic data is a matrix (1003 X 2) with colnames. 
the first column is the name of genotype and second column is the phenotypic value. 
colnames is 'Genotype' and 'Phenotype'.


### genotypic data (genodata.txt)

genotypic is a matrix (1003 X m, m is the nember of marker) with colnames and rownames.

rownames is the name of genotypes and colnames is the name of markers.


### Rogers' ditance matrix (all_chr_1003_RD )

a matrix (1003  X 1003) with colnames and rownames.

Rogers' ditance (RD) matrix is used for caculating the kinship, kinship= 1 - RD



### GEN file (hapdata.txt) and SAMPLE file (hapdata_info.txt) for the haplotype
The two files are obtained by the software 'SHAPEIT' (https://mathgen.stats.ox.ac.uk/genetics_software/shapeit/shapeit.html)

The format of the files described at http://mathgen.stats.ox.ac.uk/genetics_software/shapeit/shapeit.html#gensample



### the results of SNP_based GWAS (snp_gwas_out.txt)
The file includes 7 column and m+1 rows (m is the number of markers)

1st column: the number of markers.

2nd column: the name of markers.

3rd column: the minor allele frequency (MAF) of markers.

4th column: the estimated effect of markers.

5th column: the p value of markers.

6th column: the posistion of markers.

5th column: the chromosome of markers.
