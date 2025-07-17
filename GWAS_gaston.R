#!/usr/bin/env Rscript
#A GWAS will be performed using a linear mixed model.
#arg1：VCF file, arg2: phenotypic file, arg3：p-value cutoff, arg4：output file prefix

args <- commandArgs(trailingOnly = TRUE)

#library
library(gaston)
library(vcfR)


vcf_file <- args[1]      
pheno_file <- args[2]    
pvalue <- args[3]	
pre <- args[4]	

print("----!----")
print("Make sure that the missing values in your phenotypic file are set as "NA".")
print("----!----")


bm <- read.vcf(vcf_file)
pheno <- read.table(pheno_file, header = TRUE, row.names = 1)


for(tra in colnames(pheno)){
	rm(y)
	rm(gwa)
	rm(fdr)
	rm(a)
	rm(bm.wp)
	gc()
	y <- pheno[,tra]
	names(y) <- rownames(pheno)
	y <- na.omit(y)
	bm.wp <- bm[bm@ped$id %in% names(y),]
	y <- y[bm.wp@ped$id]
	standardize(bm.wp) <- "mu_sigma"

	grm <- GRM(bm.wp)
#with PCA result
	gwa <- association.test(bm.wp, Y = y, method = "lmm", response = "quantitative", test = "lrt", eigenK = eigen(grm), p = 10)
	col <- rep("orange3", nrow(gwa))
	col[gwa$chr %% 2 == 1] <- "blue4"
	png(paste0(pre,"_",tra, "_gaston_lmm_w_pca.png"), width = 3000, height = 2500, res = 400)
	manhattan(gwa, pch = 20, col = col)
	abline(h = -log10(as.numeric(pvalue)), col = "red", lwd = 2, lty = 2)
	dev.off()
	a <- gwa[,c("chr", "pos", "A1", "A2", "p")]
	a$"SNP" <- paste(a$A1, a$A2, sep=":")
	a <- subset(a, p <= 0.01)
	colnames(a) <- c("CHR", "BP",  "A1", "A2", "P", "SNP")
	a <- a[,c("CHR", "BP", "SNP", "P")]
	write.csv(a,paste0(pre, "_", tra, "gaston_lmm.gwas"), row.names = FALSE)

#without PCA result
	gwa <- association.test(bm.wp, Y = y, method = "lmm", response = "quantitative", test = "lrt", eigenK = eigen(grm), p = 0)
	col <- rep("orange3", nrow(gwa))
	col[gwa$chr %% 2 == 1] <- "blue4"
	png(paste0(pre,"_",tra, "_gaston_lmm_wo_pca.png"), width = 3000, height = 2500, res = 400)
	manhattan(gwa, pch = 20, col = col)
	abline(h = -log10(as.numeric(pvalue)), col = "red", lwd = 2, lty = 2)
	dev.off()
	a <- gwa[,c("chr", "pos", "A1", "A2", "p")]
	a$"SNP" <- paste(a$A1, a$A2, sep=":")
	a <- subset(a, p <= 0.01)
	colnames(a) <- c("CHR", "BP",  "A1", "A2", "P", "SNP")
	a <- a[,c("CHR", "BP", "SNP", "P")]
	write.csv(a,paste0(pre, "_", tra, "_gaston_lmm_wo_pca.gwas"), row.names = FALSE)
}

