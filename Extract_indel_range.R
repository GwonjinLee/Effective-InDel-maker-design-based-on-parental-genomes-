#!/usr/bin/env Rscript

library(dplyr)

indel <-read.table('Indel.txt',
                   sep='\t',
                   header=FALSE)

marker_ran <- read.csv("marker_region.csv")

#colnames(indel) <- c("chr", "pos", "ID", "ref", "alt", "qual", "filter", "info", "length") # from GSAlign

colnames(indel) <- c("chr", "pos", "ID", "ref", "alt", "qual", "filter", "info", "format", "sample", "length") # from vcftools

mer <- merge(indel, marker_ran, by = "chr")

#Filter for +- 2Mb

mer2 <- mer %>%
  filter(pos > first & pos < first+2000000 | pos > second-2000000 & pos < second)


mer2$start <- with(mer2, pos-500)
mer2$end <-   with(mer2, pos+500)
mer2$name <- paste(mer2$pos,mer2$location, mer2$length)

#bed <- mer2[,c(1, 13,14,15)] # from GSAlign
bed <- mer2[,c(1, 15,16,17)] # from vcftools

write.table(bed, file = "marker_pos.bed", col.names = FALSE, row.names = FALSE, sep = "\t", quote = FALSE)





