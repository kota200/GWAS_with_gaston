#!/usr/bin/env Rscript

#Get the path of the vcf file
args <- commandArgs(trailingOnly = TRUE)
#get the value from the argument
name <- args[1]

library(poolr)

file_list <- paste0(name, "_maf005_", 1:7, "_marker_matrix.txt")

meff_values <- numeric(length(file_list))

for (i in seq_along(file_list)) {
  marker <- read.table(file_list[i], row.names = 1, header = TRUE)
  rx <- cor(marker)
  meff_values[i] <- meff(rx, method = "liji")
}


Meff <- sum(meff_values)


print(1 - (1 - 0.05)^(1/Meff))


sink("Li_ji_thresh.txt")
print(1 - (1 - 0.05)^(1/Meff))
sink()
