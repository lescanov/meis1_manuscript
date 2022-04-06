#The purpose of this project is to identify binding of ETS factors to MEIS1 E2 enhancer
#in cells with varying levels of MEIS1 expression.
#According to CCLE, SEM has comparable MEIS1 expression as U937 and KG1, with
#RS4;11 and THP1 having considerably less expression.
#CHiP-seq data I am interested in is:
# > SEM (CTCF, FLI1, ELF1, ERG, SPI1, MYB, H3K27ac (DMSO treatment))
# > THP1 (H3K27ac, H3K4me3, H3K4me1, H3K79me2)
# > RS4;11 data is Hi-C (which also contains SEM)

#According to paper, peaks were called against input controll using HOMER.
#bed files with peaks are uploaded in supplemental info of geo dataset, will extract with GEOquery
#bed will be converted to bedgraph and then into bigwig on bash

library(GEOquery)
library(tidyverse)
library(rtracklayer)
library(BSgenome)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(karyoploteR)
library(org.Hs.eg.db)

#setting wd
wdir <- getwd()

#Extracting bed files from supplemental info
chip <- getGEOSuppFiles(
  GEO = "GSE117864", 
  baseDir = wdir
)



