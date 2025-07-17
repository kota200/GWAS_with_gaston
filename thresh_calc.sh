#!/bin/sh

#argument 1 is the vcf file (.vcf.gz or .vcf)
#argument 2 is the prefix of output file name
#arugumet 3 is the number of chromosome

file=${1:?}
name=${2:?}
nchr=${3:?}

#minor allele frequency 0.05
echo "plink2 started"
plink2 --make-pgen --new-id-max-allele-len 9999 --vcf ${file} --maf 0.05 --out ${name}_maf0.05 --set-all-var-ids @:#,\$r,\$a --allow-extra-chr
plink2 --pfile ${name}_maf0.05 --recode vcf --out ${name}_maf005

echo "plink2 ended"
bgzip ${name}_maf005.vcf; tabix ${name}_maf005.vcf.gz
echo ""
echo ""

#chromosome0 or other chromosomes was removed
echo "bcftools started"
bcftools view -t ^0 ${name}_maf005.vcf.gz -Oz -o ${name}_maf005_chr.vcf.gz
tabix ${name}_maf005_chr.vcf.gz

for chr in `seq 1 nchr`; do bcftools view ${name}_maf005_chr.vcf.gz -r $chr -Oz -o ${name}_maf005_${chr}.vcf.gz; done
for chr in `seq 1 nchr`; do zcat ${name}_maf005_${chr}.vcf.gz | sed "s/#CHROM/CHROM/g" | cut -f3,10- | grep -v "##" | awk '{
    if (NR == 1) {
        print;
    } else {
        for (i=2; i<=NF; i++) {
      split($i, arr, "|");
            if (arr[1] == 0 && arr[2] == 0) {
                $i = 0;
            } else if (arr[1] != 0 && arr[2] != 0) {
                $i = 2;
            } else {
                $i = 1;
            }
        }
        print;
    }
}' > ${name}_maf005_${chr}_marker_matrix.txt; done

echo ""
echo ""
echo ""
echo "Rscript started"
echo ""
echo "Calculating..."
Rscript liji_cor.R ${name}
