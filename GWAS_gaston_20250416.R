#!/usr/bin/env Rscript
#pcaを上位10グループ、lmmをなし
#第一引数：VCFファイル,第二引数: 表現型ファイル,第三引数：p-value cutoff,第四引数：output file prefix
#コマンドライン引数からvcfファイルと表現型ファイルのパスを取得
args <- commandArgs(trailingOnly = TRUE)

#library
library(gaston)
library(vcfR)

# 引数から値を取得
vcf_file <- args[1]      # 第一引数: VCFファイル
pheno_file <- args[2]    # 第二引数: 表現型ファイル
pvalue <- args[3]	#第三引数：p-value cutoff
pre <- args[4]	#第四引数：output file prefix

print("----!----")
print("表現型データの欠損値はNAにしてください")
print("----!----")

#データの読み込み
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
#GWASの実行
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

#PCA入れないで
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

#最後にコマンドラインで以下を実行するfor i in *.gwas; do sed "s/^/Chr0/g" $i | sed "s/\"//g" | sed "s/Chr0CHR/CHR/g" | sed "s/,/\t/g" > tmp; mv tmp $i; done

