#Purpose of this is to convert txt files into bed format and then convert to bedgraph and then bigwig
library(GEOquery)
library(tidyverse)
library(rtracklayer)
library(BSgenome)
library(TxDb.Hsapiens.UCSC.hg18.knownGene)
library(karyoploteR)
library(org.Hs.eg.db)


#Importing data from 
stemcell <- getGEO("GSE38865",
                   destdir = "/Users/suika-san/Dropbox/miR-210_project/R_projects/GSE38865", 
                   AnnotGPL = TRUE,
                   getGPL = TRUE)

#This creates multiple txt files that contain info for the experiment series
#specifically contains patient info and expression data
#based on these txt files, you can retrieve patient data with expression values
#I have manually edited on the txt files to contain expression values segregated by asscension code
#this is saved as chiparray.txt

#Default of getGEO is to save as a gse matrix
#to extract elements of matrix such as patient data, experiment set up, values, must use Biobase

#This is for ChIP array
histone <- stemcell$`GSE38865-GPL15724_series_matrix.txt.gz`
patient_info <- pData(phenoData(histone))

#the characteristics_ch1.5 column contains patient karyotype
#plan to divide patients into normal karyotype and non-normal karyotype
patient_info <- patient_info %>%
  dplyr::select(title, characteristics_ch1.5, characteristics_ch2.8, characteristics_ch2.9)

#importing chiparray.txt
chip_array <- read.delim("chiparray.txt", header = TRUE)

#importing genome annotation for chip array which was received from here:
#https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GPL15724
#annotation ID are not in order
chip_annotation <- read.delim("chip_annotation.txt")
chip_annotation <- chip_annotation %>%
  dplyr::arrange(CHROMOSOMAL_LOCATION)

#checking if MEIS1 exists in this data
#It does, in the annotation
MEIS1 <- chip_annotation %>%
  dplyr::filter(GENE_SYMBOL %in% "MEIS1")

#Checking if MEIS1 exists in Chip
MEIS1_chip <- chip_array %>%
  dplyr::filter(ID_REF %in% MEIS1$ID)

#### THE FOLLOWING CODE WAS TEST TO CONVERT TXT INTO BIGWIG AND IT PASSED ####

#Will have to follow this format if I want to convert txt to bed:
#https://genome.ucsc.edu/FAQ/FAQformat.html#format1

#test_bed <- data.frame(
#  chrom = as.vector(MEIS1$CHROMOSOME), 
#  chromStart = as.numeric(MEIS1$RANGE_START), 
#  chromEnd = as.numeric(MEIS1$RANGE_END),
#  name = MEIS1$GENE_SYMBOL,
#  score = as.numeric(MEIS1_chip$GSM950948)
#)


#exporting bed file
#write.table(test_bed, "test.bed", sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)

#trying r to do this
#bed_loaded <- import(con = "/Users/suika-san/Dropbox/miR-210_project/R_projects/GSE38865/test.bed", 
#                     format = "bed", 
#                     genome = "hg18", 
#                     sep = c("\t"))

#must be exported as bigwig, wig needs span to be consistent across all genome attributes
#export.bw(object=bed_loaded, 
#          con="/Users/suika-san/Dropbox/miR-210_project/R_projects/GSE38865/test.bw", 
#          format = "bw")

#plotting this

#testing bw
#importing bw just created above
#test_dir <- "/Users/suika-san/Dropbox/miR-210_project/R_projects/GSE38865/test.bw"
#test_bw <- BigWigFile(test_dir)

#### END OF TEST ####


#This actually works
#Next need to write a function to convert txt files to bigwig
#For MEIS1 data

create_txt_ready_for_bed <- function(GSM){
  
  bed <- data.frame(
    chrom = as.vector(MEIS1$CHROMOSOME), 
    chromStart = as.numeric(MEIS1$RANGE_START), 
    chromEnd = as.numeric(MEIS1$RANGE_END),
    name = MEIS1$GENE_SYMBOL,
    score = as.numeric(MEIS1_chip[[GSM]])
  )
  
  return(bed)
  
}

#in this function name is the name of file, and name_character is the name as character
load_bigwig_from_txt <- function(name, name_character){
  
  #export as bed
  write.table(name, 
              paste(name_character, ".bed", sep = ""), 
              sep = "\t", 
              row.names = FALSE, 
              col.names = FALSE, 
              quote = FALSE)
  
  #import as bed
  bed_load <- import(con = paste("/Users/suika-san/Dropbox/miR-210_project/R_projects/GSE38865/",name_character,".bed", 
                                 sep = ""),  
                       format = "bed", 
                       genome = "hg18", 
                       sep = c("\t"))
  
  #export as bw
  export.bw(object=bed_load, 
            con= paste("/Users/suika-san/Dropbox/miR-210_project/R_projects/GSE38865/",name_character,".bw", 
                       sep = ""), 
            format = "bw")
  
  #loading bigwig file
  dir_foo <- paste("/Users/suika-san/Dropbox/miR-210_project/R_projects/GSE38865/",name_character,".bw", 
                   sep = "")
  BigWigFile(dir_foo)
  
}

#from this creating a function to load as many samples as you want
#this function names a list, names, of GSM names for the patient samples which can be retrieved from MEIS1_chip
#name_quoted is a list of names you wish to rename GSM samples, must be list of characters
bigwig_from_txt <- function(GSM, names){
  
  list_of_patient_chip <- lapply(GSM, create_txt_ready_for_bed)
  names(list_of_patient_chip) <- names
  
  empty_list <- list()
  
  for(i in seq_along(names)){
    empty_list[[i]] <- load_bigwig_from_txt(list_of_patient_chip[[i]], names[[i]])
  }
  
  names(empty_list) <- names

  
  return(empty_list)
}

#apply function
#stating variables for AML patients, next will have to be 
GSM_names <- as.list(names(MEIS1_chip[, 2:27]))
sample_names <- gsub(" ", "_", patient_info$title)
patient_names <- sample_names[1:26]
patient_names <- as.list(patient_names)

AML_chip <- bigwig_from_txt(GSM_names, patient_names)
list2env(AML_chip, .GlobalEnv)

#Now stating for CD34+ control
GSM_cd34_names <- as.list(names(MEIS1_chip[, 28:37]))
sample_cd34_names <- sample_names[27:36]
sample_cd34_names <- as.list(sample_cd34_names)

cd34_chip <- bigwig_from_txt(GSM_cd34_names, sample_cd34_names)
list2env(cd34_chip, .GlobalEnv)

#### Plotting karyotype ####

#Plotting Normal Karyotype patients first
#setting E2 region 
MEIS1_E2_region <- toGRanges("chr2:66544400 -66549000", genome = "hg18")
#Whole MEIS1 chr2:66480636-66780636

#Old region: chr2:66528000-66532000

#setting plot parameters
pp <- getDefaultPlotParams(plot.type=1)
pp$leftmargin <- 0.15
pp$topmargin <- 15
pp$bottommargin <- 15
pp$ideogramheight <- 5
pp$data1inmargin <- 10
pp$max <- 10

#This function serves as a preliminary troubleshooting tool
plotbigwig <- function(chip_file, patient_name, ra, rb){
  kp <- kpPlotBigWig(kp, data=chip_file, r0=ra, r1=rb, ymax = 1 , ymin = -1)
  computed.ymax <- kp$latest.plot$computed.values$ymax
  computed.ymin <- kp$latest.plot$computed.values$ymin
  kpAxis(kp, ymin=computed.ymin, ymax=computed.ymax, r0=ra, r1=rb, cex = 0.5)
  kpAddLabels(kp, labels = paste(patient_name), r0=ra, r1=rb, label.margin = 0.07, cex = 0.5)
}

#plotting bigwig file
kp <- plotKaryotype(zoom = MEIS1_E2_region, plot.params = pp, genome = "hg18")
genes_data <- makeGenesDataFromTxDb(TxDb.Hsapiens.UCSC.hg18.knownGene,
                                    karyoplot=kp,
                                    plot.transcripts = TRUE, 
                                    plot.transcripts.structure = TRUE)
genes_data <- addGeneNames(genes_data)
genes_data <- mergeTranscripts(genes_data)
kpAddBaseNumbers(kp, tick.dist = 1000, minor.tick.dist = 200,
                 add.units = TRUE, digits = 6)
kpPlotGenes(kp, data=genes_data, r0=0, r1=0.05, gene.name.cex = 0.5, 
            avoid.overlapping = TRUE, cex = 0.5)

#patient 1
plotbigwig(Patient_1_H3K9ac, "patient 1", 0.1, 0.2)
#patient 7
plotbigwig(Patient_7_H3K9ac, "patient 7", 0.25, 0.35)
#patient 10
plotbigwig(Patient_10_H3K9ac, "patient 10", 0.4, 0.5)
#patient 13
plotbigwig(Patient_13_H3K9ac, "patient 13", 0.55, 0.65)
#patient 15
plotbigwig(Patient_15_H3K9ac, "patient 15", 0.7, 0.8)
#patient 18
plotbigwig(Patient_17_H3K9ac, "patient 17", 0.85, 0.95)

#NK continued
kp <- plotKaryotype(zoom = MEIS1_E2_region, plot.params = pp, genome = "hg18")
kpAddBaseNumbers(kp, tick.dist = 1000, minor.tick.dist = 200,
                 add.units = TRUE, digits = 6)
kpPlotGenes(kp, data=genes_data, r0=0, r1=0.05, gene.name.cex = 0.5, 
            avoid.overlapping = TRUE, cex = 0.5)

#patient 18
plotbigwig(Patient_18_H3K9ac, "patient 18", 0.1, 0.2)
#patient 22
plotbigwig(Patient_22_H3K9ac, "patient 22", 0.25, 0.35)
#patient 23
plotbigwig(Patient_13_H3K9ac, "patient 23", 0.4, 0.5)
#patient 25
plotbigwig(Patient_25_H3K9ac, "patient 25", 0.55, 0.65)
#patient 26
plotbigwig(Patient_26_H3K9ac, "patient 26", 0.7, 0.8)



#Non-normal karyotype
#plotting bigwig file
kp <- plotKaryotype(zoom = MEIS1_E2_region, plot.params = pp, genome = "hg18")
kpAddBaseNumbers(kp, tick.dist = 1000, minor.tick.dist = 200,
                 add.units = TRUE, digits = 6)
kpPlotGenes(kp, data=genes_data, r0=0, r1=0.05, gene.name.cex = 0.5, 
            avoid.overlapping = TRUE, cex = 0.5)

#patient 2
plotbigwig(Patient_2_H3K9ac, "patient 2", 0.1, 0.2)
#patient 3
plotbigwig(Patient_3_H3K9ac, "patient 3", 0.25, 0.35)
#patient 5
plotbigwig(Patient_5_H3K9ac, "patient 5", 0.4, 0.5)
#patient 6
plotbigwig(Patient_6_H3K9ac, "patient 6", 0.55, 0.65)
#patient 16
plotbigwig(Patient_16_H3K9ac, "patient 16", 0.7, 0.8)
#patient 19
plotbigwig(Patient_19_H3K9ac, "patient 19", 0.85, 0.95)

#### CD34+ donors ####
#Starting with H3K9ac
kp <- plotKaryotype(zoom = MEIS1_E2_region, plot.params = pp, genome = "hg18")
kpAddBaseNumbers(kp, tick.dist = 1000, minor.tick.dist = 200,
                 add.units = TRUE, digits = 6)
kpPlotGenes(kp, data=genes_data, r0=0, r1=0.05, gene.name.cex = 0.5, 
            avoid.overlapping = TRUE, cex = 0.5)
plotbigwig(CD34_normal_control_1_H3K9ac, "control 1", 0.1, 0.2)
plotbigwig(CD34_normal_control_3_H3K9ac, "control 3", 0.25, 0.35)
plotbigwig(CD34_normal_control_4_H3K9ac, "control 4", 0.4, 0.5)
plotbigwig(CD34_normal_control_5_H3K9ac, "control 5", 0.55, 0.65)
plotbigwig(CD34_normal_control_7_H3K9ac, "control 7", 0.7, 0.8)

#H3K27me3
kp <- plotKaryotype(zoom = MEIS1_E2_region, plot.params = pp, genome = "hg18")
kpAddBaseNumbers(kp, tick.dist = 1000, minor.tick.dist = 200,
                 add.units = TRUE, digits = 6)
kpPlotGenes(kp, data=genes_data, r0=0, r1=0.05, gene.name.cex = 0.5, 
            avoid.overlapping = TRUE, cex = 0.5)
plotbigwig(CD34_normal_control_1_H3K27me3, "control 1", 0.1, 0.2)
plotbigwig(CD34_normal_control_3_H3K27me3, "control 3", 0.25, 0.35)
plotbigwig(CD34_normal_control_4_H3K27me3, "control 4", 0.4, 0.5)
plotbigwig(CD34_normal_control_5_H3K27me3, "control 5", 0.55, 0.65)
plotbigwig(CD34_normal_control_7_H3K27me3, "control 7", 0.7, 0.8)

#heatmap example
plotheatmap <- function(bed_file, patient_name, ra, rb, value){
  kp <- kpHeatmap(kp, data=bed_file, r0=ra, r1=rb, ymax = 1 , ymin = -1, y = bed_file[[value]])
  kpAxis(kp, ymin=-1, ymax=1, r0=ra, r1=rb, cex = 0.5)
  kpAddLabels(kp, labels = paste(patient_name), r0=ra, r1=rb, label.margin = 0.07, cex = 0.5)
}

#importing bed files 
bed_import <- function(file_name){
  import(con =paste("/Users/suika-san/Dropbox/miR-210_project/R_projects/GSE38865/", file_name, sep = ""),
         format = "bed", 
         genome = "hg18", 
         sep = c("\t"))
}
bed_list <- as.list(c("CD34_normal_control_1_H3K27me3.bed", 
                      "CD34_normal_control_3_H3K27me3.bed",
                      "CD34_normal_control_4_H3K27me3.bed",
                      "CD34_normal_control_5_H3K27me3.bed",
                      "CD34_normal_control_7_H3K27me3.bed"))

bed_test <- lapply(bed_list, bed_import)
names(bed_test) <- bed_list
list2env(bed_test, .GlobalEnv)

#This didn't work out as planned
#Can't manipulate bed files
kp <- plotKaryotype(zoom = MEIS1_E2_region, plot.params = pp, genome = "hg18")
kpAddBaseNumbers(kp, tick.dist = 1000, minor.tick.dist = 200,
                 add.units = TRUE, digits = 6)
kpPlotGenes(kp, data=genes_data, r0=0, r1=0.05, gene.name.cex = 0.5, 
            avoid.overlapping = TRUE, cex = 0.5)
plotheatmap(CD34_normal_control_1_H3K27me3.bed, "control 1", 0.1, 0.2)
plotheatmap(CD34_normal_control_3_H3K27me3.bed, "control 3", 0.25, 0.35)
plotheatmap(CD34_normal_control_4_H3K27me3.bed, "control 4", 0.4, 0.5)
plotheatmap(CD34_normal_control_5_H3K27me3.bed, "control 5", 0.55, 0.65)
plotheatmap(CD34_normal_control_7_H3K27me3.bed, "control 7", 0.7, 0.8)

