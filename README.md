# An Rscript to perform a GWAS with gaston (https://cran.r-project.org/web/packages/gaston)
An Rscript and shell scripts to perform a genome-wide association study (GWAS). 

Prerequisites:

・gaston (https://cran.r-project.org/web/packages/gaston/gaston.pdf)

・vcfR (https://github.com/knausb/vcfR)



Input file:

・A genotype file with the VCF file format

・A phtnotype file with the tsv file format


Output file:

・Manhattan plots

・GWAS result files (about SNPs whose p-values are less than 0.01)



To get the threshold value:
You can set the threshold value using the method in Li and Ji (2005).
Please use the thresh_calc.sh command.

Prerequisites

・plink2 (https://www.cog-genomics.org/plink/2.0/)

・bcftools (https://samtools.github.io/bcftools/bcftools.html)

.poolr
