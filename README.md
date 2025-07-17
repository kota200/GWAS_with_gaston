# An Rscript to perform a GWAS with gaston (https://cran.r-project.org/web/packages/gaston)
An Rscript and shell scripts to perform a genome-wide association study (GWAS). The threshold was set with the method developed by Li and Ji (2005).



Prerequisites:

・plink2 

・bcftools

・gaston

・poolr

・vcfR




Input file:

・A genotype file with the VCF file format

・A phtnotype file with the tsv file format


Output file:

・Manhattan plot

・GWAS result files (about SNPs whose p-values are less than 0.01)

