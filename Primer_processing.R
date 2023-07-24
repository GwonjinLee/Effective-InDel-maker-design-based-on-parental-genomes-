#!/usr/bin/env Rscript


# Primer processing for finding degenerate sequences

library(dplyr)
library(magrittr)
library(tidyr)
library(misty)
library(stringr)

file.remove("A344_primers_alt.bed") 
file.remove("blast_A344_primers_paired.txt") 
file.remove("A344_primers_altseq.fasta")
file.remove("A344_primers_altseq_alpha.fasta")
file.remove("blast_A344_primers_paired_alpha.fasta")
file.remove("blast_A344_primers_paired.fasta")
file.remove("degen_primers.txt")

# Load blasted primers with the alt genome
alt_primers <- read.delim("blast_A344_alt_primers.txt", header = F, sep = "\t"  )

# Load blatsted and merged primers both with ref and alt genome
merged_primers <- read.delim("blast_A344_primers_merged.txt", header = F)


merged_primers$V1 %<>%
  gsub("_F", " F", .) %>%
  gsub("_R", " R", .)

merged_primers <- str_split_fixed(merged_primers$V1, " ", 2)

merged_primers <- as.data.frame(merged_primers)

dupl <- df.duplicated(merged_primers,V1)

dupl$V3 <- paste(dupl$V1,  dupl$V2, sep= "_")

colnames(dupl) <- c("name", "direction", "V1")

alt_mer <- merge(dupl, alt_primers, by = "V1" )

# delete if chr are different 
alt_mer[c("chr", "name")] <- str_split_fixed(alt_mer$V1, ":", 2)
alt_mer$V2 <- sub("^", "chr", alt_mer$V2 ) # If V2 is numeric add chr before the number

alt_mer <- subset(alt_mer, chr == V2)


# +- 20 bp
alt_mer$start <- ifelse(alt_mer$direction == "F", alt_mer$V9-20, alt_mer$V9+20)
alt_mer$end <- ifelse(alt_mer$direction == "F", alt_mer$V10+20, alt_mer$V10-20)

# change the order of start and end for reverse primers
alt_mer2 <- transform(alt_mer, end = pmax(start, end), start = pmin(start, end))
alt_mer2$num <- 1 
  

#add strand
alt_mer2  <- alt_mer2  %>%mutate(strand=ifelse(grepl("F$",direction),"+","-"))


alt_bed <- alt_mer2[,c(15,16,17,1,18,19)]

write.table(alt_bed, file = "A344_primers_alt.bed", sep = "\t",  col.names = F, row.names = F, quote = F)
write.table(alt_mer$V1, file = "blast_A344_primers_paired.txt", sep = "\t",  col.names = F, row.names = F, quote = F)


