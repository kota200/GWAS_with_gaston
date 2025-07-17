# An Rscript to perform a GWAS with gaston (https://cran.r-project.org/web/packages/gaston)
An Rscript and shell scripts to perform a genome-wide association study (GWAS). 

<img width="1200" height="1000" alt="BR_PM_Tift_imputed_20250416_anthocyanin_level_in_shoot_basal_part_gaston_lmm_wo_pca" src="https://github.com/user-attachments/assets/16bce962-1336-457c-9782-7767ce0a4f9b" />


Prerequisites:

・gaston (https://cran.r-project.org/web/packages/gaston/gaston.pdf)

・vcfR (https://github.com/knausb/vcfR)



Input file (please check the example files in the test_file folder):

・A genotypic file with the VCF file format

・A phtnotypic file with the tsv file format

Example

```
Rscript GWAS_gaston.R VCF_file Phenotypic_file p-value_cutoff(value) outputfile_prefix
```

Output file:

・Manhattan plots

・GWAS result files (about SNPs whose p-values are less than 0.01)



To get the threshold value:
You can set the threshold value using the method in Li and Ji (2005)(https://www.nature.com/articles/6800717).
Please use the thresh_calc.sh command (the following link was refferred: http://morotalab.org/apsc5984-2020/day35/day35.html#li-and-ji-2005).

Prerequisites

・plink2 (https://www.cog-genomics.org/plink/2.0/)

・bcftools (https://samtools.github.io/bcftools/bcftools.html)

・poolr（https://cran.r-project.org/web/packages/poolr/index.html）
