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



### haplotype phases file (hapdata.txt) and haplotype information file (hapdata_info.txt)
The two files are obtained by the software 'SHAPEIT' (https://mathgen.stats.ox.ac.uk/genetics_software/shapeit/shapeit.html)

The format of the files described at http://mathgen.stats.ox.ac.uk/genetics_software/shapeit/shapeit.html#gensample



### the results of SNP_based GWAS (snp_gwas_out.txt)
