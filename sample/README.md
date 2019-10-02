# sample data for functional haplotype

## Data format

### 1.  genotypic data (sample_genodata.txt)

Genotypic data is file with n+1  rows and m+1 columns, in which n is the number of individuals and m is the number of markers.  
The first column is the name of individuals and the first row is the names of markers.    
The genotypic value for each marker should be 0, 1 and 2 (0 for the reference allele, 2 for alternative allele and 1 for heterozygote).   
 *In the sample data, n is 30 and m is 1000.*
 
#### *NOTE:Genotypic data should not include any of missing value.*

### 2.  marker information (sample_genodata_info.txt)

Marker information includes 5 columns, which are the name of markers, reference allele, alternative allele, chromosome and posistion.

### 3.  phenotypic data (sample_pheno.txt)

Phenotypic data have two colums. the first column is the name of individuals and the second column is the phenotypic value for the corresponded individual.

### 4.  kinship matrix or distance matrix (sample_kinship.txt or sample_dis.txt)

Only one of them is needed for GWAS, which means either kinship matrix or distance matrix.   
Kinship or distance is a matrix (n x n, n is the number of individuals) with rownames and colnames. rownames and colnames is the same that is the name of individuals.   
Kinship can be caculated using the method of [VanRaden](https://github.com/Zhiwu-Zhang-Lab/GAPIT/blob/master/GAPIT.kinship.VanRaden.R).   
*Note: the maximum value in the distance matrix should be 1.*

