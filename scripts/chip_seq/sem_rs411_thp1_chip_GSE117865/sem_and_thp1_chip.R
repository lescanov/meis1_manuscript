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
#called peaks for transcription factors as bed, broad histone peaks will be plotted as bigwig

library(GEOquery)
library(tidyverse)
library(rtracklayer)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(karyoploteR)
library(org.Hs.eg.db)
library(BSgenome.Hsapiens.UCSC.hg19)

#setting wd
wdir <- getwd()

#Extracting bed files from supplemental info
#all unzipped using bash:
#for file in *; do gunzip -k $file; done
chip <- getGEOSuppFiles(
  GEO = "GSE117864", 
  baseDir = wdir, 
  makeDirectory = FALSE
)

#### Cleaning bed files and exporting as bedGraph ####
#the problem with the bed files is that they have 40 rows of extra info that is not needed
#they also have extra columns pertaining to statistics in the homer peak call
#furthermore, these are no shorted by chr, chromStart
#must clean these bed files before they can be plotted
#will save as bedGraph so that score is saved as numeric

#creating a function to read in bed files
read_in_bed <- function(bed_file){
  bed <- read.table(
    file = paste(wdir, "/", bed_file, sep = ""), 
    header = FALSE, 
    stringsAsFactors = FALSE,
    quote = "", 
    sep = "\t"
  )
}

#in bash unzipped files

#creating a list of bed files obtained from GEOquery
bed_list <- as.list(c("GSE117864_CTCF_peaks.bed", "GSE117864_ELF1_peaks.bed", "GSE117864_ERG_peaks.bed", "GSE117864_FLI1_peaks.bed", 
                         "GSE117864_H3K27ac_0uM_EPZ_peaks.bed", "GSE117864_H3K79me3_0uM_EPZ_peaks.bed", 
                         "GSE117864_MYB_peaks.bed", "GSE117864_SPI1_peaks.bed", "GSE117864_THP1_H3K27ac_peaks.bed", 
                         "GSE117864_THP1_H3K4me1_peaks.bed", "GSE117864_THP1_H3K4me3_peaks.bed", "GSE117864_THP1_H3K79me2_peaks.bed"))

#creating a list of names, used to name dataframes during export as bedGraph
bed_names <- as.list(c("ctcf", "elf1", "erg", "fli1", "h3k27ac", "h3k79me3", "myb", "spi1", "h3k27ac_thp1", "h3k4me1_thp1", "h3kme3_thp1", "h3k79me_thp1"))

#getting a list of imported bed files as table
test <- lapply(bed_list, read_in_bed)

#editing tables to fit proper bed format of chr, chrstart, chrend, score
for(i in seq_along(test)){
  tmp <- test[[i]]
  tmp <- tmp %>%
    dplyr::select(V2, V3, V4, V8)
  colnames(tmp) <- c("chr", "chromStart", "chromEnd", "score")
  tmp <- tmp %>%
    dplyr::arrange(chr, chromStart)
  test[[i]] <- tmp
}

#renaming elements of list containing pre-bedgraph files
names(test) <- bed_names
list2env(test, .GlobalEnv)

#exporting bed
#i hate this but it's almost 1am and i can't get to automate this, just have to brute force it
#probably better to do it in bash
export.bedGraph(
  h3k79me_thp1, 
  con = paste(wdir,"/h3k79me_thp1.bedgraph", sep = ""),
  format = "bedGraph"
)

#manually importing each of the files...
h3k79me3_thp1_bedgraph <- import.bedGraph(
  con = paste(wdir, "/h3k79me_thp1.bedgraph" ,sep =""), 
  genome = "hg19", 
  format = "bedGraph"
)

#### Plotting called peaks for SEM
#Have called peaks for the following:
# CTCF, ELF1, ERG, FLI1, MYB, SPI1
# Will plot DMSO treated H3k27ac and H3k79me3 on top of this

#first defining E2 enhancer region
MEIS1_E2_region <- toGRanges("chr2:66600000-66799891", genome = "hg19")
#E2 around: chr2:66669878-66675877
#whole meis1 chr2:66662532-66799891"
#defiing plot parameters
pp <- getDefaultPlotParams(plot.type=1)
pp$leftmargin <- 0.15
pp$topmargin <- 15
pp$bottommargin <- 15
pp$ideogramheight <- 5
pp$data1inmargin <- 10
pp$max <- 10

#This function serves as a preliminary troubleshooting tool
kp <- plotKaryotype(zoom = MEIS1_E2_region, plot.params = pp, genome = "hg19")
genes_data <- makeGenesDataFromTxDb(TxDb.Hsapiens.UCSC.hg19.knownGene,
                                    karyoplot=kp,
                                    plot.transcripts = TRUE, 
                                    plot.transcripts.structure = TRUE)
genes_data <- addGeneNames(genes_data)
genes_data <- mergeTranscripts(genes_data)
kpAddBaseNumbers(kp, tick.dist = 1000, minor.tick.dist = 200,
                 add.units = TRUE, digits = 6)
kpPlotGenes(kp, data=genes_data, r0=0, r1=0.1, gene.name.cex = 0.5, 
            avoid.overlapping = TRUE, cex = 0.5)

plotkpregions <- function(bedgraph_file, chr_label, a, b){
  kp <- kpPlotRegions(kp, data = bedgraph_file, r0 = a, r1 = b)
  kpAddLabels(kp, labels = chr_label, r0 = a, r1 = b, label.margin = 0.05, cex = 0.5)
}
#no binding from ELF1, MYB
plotkpregions(ctcf_bedgraph, "CTCF", 0.1, 0.15)
plotkpregions(fli1_bedgraph, "FLI1", 0.2, 0.25)
plotkpregions(erg_bedgraph, "ERG", 0.3, 0.35)
plotkpregions(spi1_bedgraph, "SPI", 0.4, 0.45)
plotkpregions(elf1_bedgraph, "ELF1", 0.5, 0.55)
plotkpregions(myb_bedgraph, "MYB", 0.6, 0.65)
plotkpregions(h3k27ac_bedgraph, "H3K27ac", 0.7, 0.8)
plotkpregions(h3k79me3_bedgraph, "H3K79me3", 0.85, 0.95)

