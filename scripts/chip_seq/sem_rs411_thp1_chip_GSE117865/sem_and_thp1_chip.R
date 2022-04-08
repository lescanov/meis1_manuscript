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

#### Importing bed files as table, cleaning up input, exporting as proper bed ####
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

bed_list <- as.list(c("GSE117864_CTCF_peaks.bed", "GSE117864_ELF1_peaks.bed", "GSE117864_ERG_peaks.bed", "GSE117864_FLI1_peaks.bed", 
                         "GSE117864_H3K27ac_0uM_EPZ_peaks.bed", "GSE117864_H3K79me3_0uM_EPZ_peaks.bed", 
                         "GSE117864_MYB_peaks.bed", "GSE117864_SPI1_peaks.bed", "GSE117864_THP1_H3K27ac_peaks.bed", 
                         "GSE117864_THP1_H3K4me1_peaks.bed", "GSE117864_THP1_H3K4me3_peaks.bed", "GSE117864_THP1_H3K79me2_peaks.bed"))
bed_names <- as.list(c("ctcf", "elf1", "erg", "fli1", "h3k27ac", "h3k79me3", "myb", "spi1", "h3k27ac_thp1", "h3k4me1_thp1", "h3kme3_thp1", "h3k79me_thp1"))

test <- lapply(bed_list, read_in_bed)

for(i in seq_along(test)){
  tmp <- test[[i]]
  tmp <- tmp %>%
    dplyr::select(V2, V3, V4, V8)
  colnames(tmp) <- c("chr", "chromStart", "chromEnd", "score")
  tmp <- tmp %>%
    dplyr::arrange(chr, chromStart)
  test[[i]] <- tmp
}
names(test) <- bed_names
list2env(test, .GlobalEnv)

#exporting as bigwig
load_bigwig_from_input <- function(name){
  
  write.table(name, 
              file = paste(name, ".bed", sep = ""), 
              sep = "\t", 
              row.names = FALSE, 
              col.names = FALSE, 
              quote = FALSE)
  
  test <- import(paste(wdir, "/", name, ".bed", sep = ""), 
                 format = "bedGraph", 
                 genome = "hg19", 
                 sep = c("\t"))
  
  export.bw(object=test, 
            con= paste(wdir, "/", name, ".bw", sep =""), 
            format = "bw")
  
  BigWigFile(paste(name, ".bw", sep = ""))
  
}

lapply(names, load_bigwig_from_input)


### test
write.table(fli1, 
            file = "fli1.bed", 
            sep = "\t", 
            row.names = FALSE, 
            col.names = FALSE, 
            quote = FALSE)

test <- import("fli1.bed", 
               format = "bedGraph", 
               genome = "hg19", 
               sep = c("\t"))

export.bw(object=test, 
          con= paste(wdir, "/fli1.bw", sep =""), 
          format = "bw")
fli1_bw <- BigWigFile("fli1.bw")
fli1_bw <- import.bw(
  con = paste(wdir, "/fli1.bw", sep = ""), 
  format = "bw"
)
#convert bed files to bigwig files on bash
#store bigwig files in folder called bigwig
bigwig_dir <- paste(wdir, "/bigwig/", sep = "")

#importing bigwig files from folder bigwig
bigwig_list <- as.list(c("CTCF_peaks.bw", "ELF1_peaks.bw", "ERG_peaks.bw", "FLI1_peaks.bw", 
                         "H3K27ac_0uM_EPZ_peaks.bw", "H3K79me3_0uM_EPZ_peaks.bw", 
                         "MYB_peaks.bw", "SPI1_peaks.bw", "THP1_H3K27ac_peaks.bw", 
                         "THP1_H3K4me1_peaks.bw", "THP1_H3K4me3_peaks.bw", "THP1_H3K79me2_peaks.bw"))
empty_list <- list()
for (i in seq_along(bigwig_list)){
  empty_list[[i]] <- BigWigFile(paste(bigwig_dir, bigwig_list[[i]], sep = ""))
}
names(empty_list) <- bigwig_list
list2env(empty_list, .GlobalEnv)

#plotting bigwig files with karyoploteR 
#will plot SEM 

#first defining E2 enhancer region
MEIS1_E2_region <- toGRanges("chr2:66662532-66799891", genome = "hg19")
#E2 around: chr2:66669878-66675877
#whole meis1 chr2:66662532-66799891"
#defiing plot paramete
pp <- getDefaultPlotParams(plot.type=1)
pp$leftmargin <- 0.15
pp$topmargin <- 15
pp$bottommargin <- 15
pp$ideogramheight <- 5
pp$data1inmargin <- 10
pp$max <- 10

#This function serves as a preliminary troubleshooting tool
plotbigwig <- function(chip_file, patient_name, ra, rb){
  kp <- kpPlotBigWig(kp, data=chip_file, r0=ra, r1=rb, ymax = "visible.region" , ymin = 0)
  computed.ymax <- kp$latest.plot$computed.values$ymax
  computed.ymin <- kp$latest.plot$computed.values$ymin
  kpAxis(kp, ymin=computed.ymin, ymax=computed.ymax, r0=ra, r1=rb, cex = 0.5)
  kpAddLabels(kp, labels = paste(patient_name), r0=ra, r1=rb, label.margin = 0.07, cex = 0.5)
}

#plotting bigwig file
test <- BigWigFile("fli1_new.bw")

kp <- plotKaryotype(zoom = MEIS1_E2_region, plot.params = pp, genome = "hg19")
genes_data <- makeGenesDataFromTxDb(TxDb.Hsapiens.UCSC.hg19.knownGene,
                                    karyoplot=kp,
                                    plot.transcripts = TRUE, 
                                    plot.transcripts.structure = TRUE)
genes_data <- addGeneNames(genes_data)
genes_data <- mergeTranscripts(genes_data)
kpAddBaseNumbers(kp, tick.dist = 1000, minor.tick.dist = 200,
                 add.units = TRUE, digits = 6)
kpPlotGenes(kp, data=genes_data, r0=0, r1=0.05, gene.name.cex = 0.5, 
            avoid.overlapping = TRUE, cex = 0.5)
kp <- kpPlotRegions(kp, data = fli1_bw, r0 = 0.3, r1 = 0.4)
kp <- kpPlotCoverage(kp, data = fli1_bw, r0 = 0.2, r1 = 0.3)
kp <- kpLines(kp, data = f, ymin = 0, ymax = "visible,region")
kp <- kpPlotBigWig(kp, data = test, r0 = 0.2, r1 = 0.3, ymin = 0, ymax = "visible.region")

plotbigwig(ctcf, 0.1, 0.15, "SPI1")

fli1_bw <- import.bw(
  "fli1_new.bw", 
  format = "bw"
)

