#!/usr/bin/env Rscript


# This script is for data table processing for input of muscle. 
# This will create multiple fasta files for each pair of a primer sequence and a corresponding seq of the alt genome.

library(dplyr)
library(tidyr)
library(stringr)
library(reshape2)
library(purrr)

primer_cand <- read.table("primers_cand.fa", sep = "\t", header = F)

primer <- cbind(primer_cand[c(TRUE, FALSE), ],
                primer_cand[c(FALSE, TRUE), ])
colnames(primer) <- c("V1", "V2", "V3", "V4")

primer$primer <- paste(primer$V1, primer$V3)
primer$seq <- paste(primer$V2, primer$V4)

primer2 <- primer[,c(5,6)]

primer2$V2 <- primer2$primer

primer3 <- cbind(primer2$V2, paste(primer2$primer, primer2$seq))

primer4 <- as.data.frame(primer3) %>% 
  mutate(V2=strsplit(V2, " ")) %>% 
  unnest(V2)

primer5 <- primer4$V2

primer6 <- as.data.frame(primer5)

colnames(primer6) <- "primer_input"

nr <- nrow(primer6)

primer7 <- split(primer6, rep(1:(nr/4), each=4))
              

imap(primer7, function(x, y) imap(x, 
                                  ~write.table(.x, paste0(y, '_', .y , ".fa"), col.names = F,
                                               row.names = FALSE, quote = FALSE, sep = "\t", dec = ",")))


