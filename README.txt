Author: Joanna J. Ilska
Date: 	11/09/2020


The following script carries out a Mendelian inconsistency check on the genotypes. 
Requires:
	a) genotypes in the plink .raw format, i.e. SNPs in columns, individuals in rows, with the (FID IID PAT MAT SEX PHENOTYPE) columns up front. 
	b) pedigree file in csv format, with column headings: Sample.ID,Pedigree,Mother.ID,Father.ID
	
If the genotype file contains only one family, then the script returns a dataframe of genotypes which passed the check. 
If there is more than one family, the Mendelian inconsistency check is run within a family. Genotypes of SNPs which failed within a given family (but passed in others) are assigned as missing for this family. SNPs which failed in all families are removed. 

The script also saves the lists of SNPs which failed within a family. 
