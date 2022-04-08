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
#I have manually edited the txt files to contain expression values segregated by asscension code
#this is saved as chiparray.txt (as seen in repo)

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

#standardizing score for MEIS1
#temporary saving ascension ids as rownames
tmp <- MEIS1_chip[1]

#removing names so can standardize
MEIS1_chip <- MEIS1_chip[-1]

#getting number of rows
nrow <- seq(nrow(MEIS1_chip))

#inverting axis
MEIS1_chip <- as.matrix(MEIS1_chip)
#standardizing each row of chip_array
for (i in seq_along(nrow)){
  MEIS1_chip[i,] <- scale(MEIS1_chip[i,])
}
MEIS1_chip <- as.data.frame(MEIS1_chip)
MEIS1_chip <- cbind(tmp, MEIS1_chip)

#This actually works
#Next need to write a function to convert txt files to bigwig
#For MEIS1 data

create_txt_ready_for_bed <- function(GSM){
  
  bed <- data.frame(
    chrom = as.vector(MEIS1$CHROMOSOME), 
    chromStart = as.numeric(MEIS1$RANGE_START), 
    chromEnd = as.numeric(MEIS1$RANGE_END),
    score = as.numeric(MEIS1_chip[[GSM]])
  )
  
  return(bed)
  
}

#setting work directory
wdir <- getwd()

#loading 

#in this function name is the name of file, and name_character is the name as character
load_bedgraph <- function(name, name_character){
  
  #export as bed
  write.table(name, 
              paste(name_character, ".bedGraph", sep = ""), 
              sep = "\t", 
              row.names = FALSE, 
              col.names = FALSE, 
              quote = FALSE)
  
  #importing as a dataranges object
  import.bedGraph(
         con = paste(wdir, "/", name_character, ".bedGraph", sep = ""), 
         format = "bedGraph", 
         genome = "hg18", 
         sep = c("\t"))
}

#from this creating a function to load as many samples as you want
#this function names a list, names, of GSM names for the patient samples which can be retrieved from MEIS1_chip
#name_quoted is a list of names you wish to rename GSM samples, must be list of characters
bedgraph_from_txt <- function(GSM, names){
  
  list_of_patient_chip <- lapply(GSM, create_txt_ready_for_bed)
  names(list_of_patient_chip) <- names
  
  empty_list <- list()
  
  for(i in seq_along(names)){
    empty_list[[i]] <- load_bedgraph(list_of_patient_chip[[i]], names[[i]])
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

histone_bed <- bedgraph_from_txt(GSM_names, patient_names)
list2env(histone_bed, .GlobalEnv)

#Now stating for CD34+ control
GSM_cd34_names <- as.list(names(MEIS1_chip[, 28:37]))
sample_cd34_names <- sample_names[27:36]
sample_cd34_names <- as.list(sample_cd34_names)

cd34_chip <- bigwig_from_txt(GSM_cd34_names, sample_cd34_names)
list2env(cd34_chip, .GlobalEnv)

#determining max and min score for each patient
for(i in seq_along(histone_bed)){
  print(c(max(score(histone_bed[[i]])), min(score(histone_bed[[i]]))))
}

#### Plotting karyotype ####

#Plotting Normal Karyotype patients first
#setting E2 region 
MEIS1_E2_region <- toGRanges("chr2:66544400-66549000", genome = "hg18")
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
kpPlotGenes(kp, data=genes_data, r0=0, r1=0.1, gene.name.cex = 0.5, 
            avoid.overlapping = TRUE, cex = 0.5)

#patient 1
kpHeatmap(kp, data = Patient_1_H3K9ac, r0 = 0.1, r1 = 0.13, y = Patient_1_H3K9ac$score, colors = c("blue", "white", "red"))
#patient 7
kpHeatmap(kp, data = Patient_7_H3K9ac, r0 = 0.15, r1 = 0.18, y = Patient_7_H3K9ac$score, colors = c("blue", "white", "red"))
#patient 10
kpHeatmap(kp, data = Patient_10_H3K9ac, r0 = 0.2, r1 = 0.23, y = Patient_10_H3K9ac$score,colors = c("blue", "white", "red"))
#patient 13
kpHeatmap(kp, data = Patient_13_H3K9ac, r0 = 0.25, r1 = 0.28, y = Patient_13_H3K9ac$score,colors = c("blue", "white", "red"))
#patient 15
kpHeatmap(kp, data = Patient_15_H3K9ac, r0 = 0.3, r1 = 0.33, y = Patient_15_H3K9ac$score,colors = c("blue", "white", "red"))
#patient 18
kpHeatmap(kp, data = Patient_18_H3K9ac, r0 = 0.35, r1 = 0.38, y = Patient_18_H3K9ac$score,colors = c("blue", "white", "red"))
#patient 22
kpHeatmap(kp, data = Patient_22_H3K9ac, r0 = 0.4, r1 = 0.43, y = Patient_22_H3K9ac$score,colors = c("blue", "white", "red"))
#patient 23
kpHeatmap(kp, data = Patient_23_H3K9ac, r0 = 0.45, r1 = 0.48, y = Patient_23_H3K9ac$score,colors = c("blue", "white", "red"))
#patient 25
kpHeatmap(kp, data = Patient_25_H3K9ac, r0 = 0.5, r1 = 0.53, y = Patient_25_H3K9ac$score,colors = c("blue", "white", "red"))
#patient 26
kpHeatmap(kp, data = Patient_10_H3K9ac, r0 = 0.55, r1 = 0.58, y = Patient_26_H3K9ac$score,colors = c("blue", "white", "red"))
#adding label
#kpAddLabels(kp, 
##            labels = "Normal Karyotype", 
#            side = "left",
#            srt = 90)

#kpAddLabels(kp, 
#            labels = "Translocations", 
#            srt = 90, 
#            pos = 4)

#Non-normal karyotype
#patient 2
kpHeatmap(kp, data = Patient_2_H3K9ac, r0 = 0.7, r1 = 0.73, y = Patient_2_H3K9ac$score,colors = c("blue", "white", "red"))
#patient 3
kpHeatmap(kp, data = Patient_3_H3K9ac, r0 = 0.75, r1 = 0.78, y = Patient_3_H3K9ac$score,colors = c("blue", "white", "red"))
#patient 5
kpHeatmap(kp, data = Patient_5_H3K9ac, r0 = 0.8, r1 = 0.83, y = Patient_5_H3K9ac$score,colors = c("blue", "white", "red"))
#patient 6
kpHeatmap(kp, data = Patient_6_H3K9ac, r0 = 0.85, r1 = 0.88, y = Patient_6_H3K9ac$score,colors = c("blue", "white", "red"))
#patient 16
kpHeatmap(kp, data = Patient_16_H3K9ac, r0 = 0.9, r1 = 0.93, y = Patient_16_H3K9ac$score,colors = c("blue", "white", "red"))
#patient 19
kpHeatmap(kp, data = Patient_19_H3K9ac, r0 = 0.95, r1 = 0.98, y = Patient_19_H3K9ac$score,colors = c("blue", "white", "red"))

#Now want to compare H3K9ac to MEIS1 expression
#### Loading RNA microarray data ####
#loading rna microarray to determine correlation of ets factors and meis1 h3k9ac
microarray_dataset <- stemcell$`GSE38865-GPL6884_series_matrix.txt.gz`
rna_microarray <- microarray_dataset@assayData[["exprs"]]

#loading annotation table taken from:
#https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GPL6884
#I trimmed off the extrenuous rows on top and saved as rna_microarray_annotation.csv in repo


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

