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

#cleaning title 
patient_info$title <- gsub(" H3K9ac", "", patient_info$title)

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

cd34_chip <- bedgraph_from_txt(GSM_cd34_names, sample_cd34_names)
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
kpHeatmap(kp, data = Patient_1, r0 = 0.1, r1 = 0.13, y = Patient_1$score, colors = c("blue", "white", "red"))
#patient 7
kpHeatmap(kp, data = Patient_7, r0 = 0.15, r1 = 0.18, y = Patient_7$score, colors = c("blue", "white", "red"))
#patient 10
kpHeatmap(kp, data = Patient_10, r0 = 0.2, r1 = 0.23, y = Patient_10$score,colors = c("blue", "white", "red"))
#patient 13
kpHeatmap(kp, data = Patient_13, r0 = 0.25, r1 = 0.28, y = Patient_13$score,colors = c("blue", "white", "red"))
#patient 15
kpHeatmap(kp, data = Patient_15, r0 = 0.3, r1 = 0.33, y = Patient_15$score,colors = c("blue", "white", "red"))
#patient 18
kpHeatmap(kp, data = Patient_18, r0 = 0.35, r1 = 0.38, y = Patient_18$score,colors = c("blue", "white", "red"))
#patient 22
kpHeatmap(kp, data = Patient_22, r0 = 0.4, r1 = 0.43, y = Patient_22$score,colors = c("blue", "white", "red"))
#patient 23
kpHeatmap(kp, data = Patient_23, r0 = 0.45, r1 = 0.48, y = Patient_23$score,colors = c("blue", "white", "red"))
#patient 25
kpHeatmap(kp, data = Patient_25, r0 = 0.5, r1 = 0.53, y = Patient_25$score,colors = c("blue", "white", "red"))
#patient 26
kpHeatmap(kp, data = Patient_26, r0 = 0.55, r1 = 0.58, y = Patient_26$score,colors = c("blue", "white", "red"))
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
kpHeatmap(kp, data = Patient_2, r0 = 0.7, r1 = 0.73, y = Patient_2$score,colors = c("blue", "white", "red"))
#patient 3
kpHeatmap(kp, data = Patient_3, r0 = 0.75, r1 = 0.78, y = Patient_3$score,colors = c("blue", "white", "red"))
#patient 5
kpHeatmap(kp, data = Patient_5, r0 = 0.8, r1 = 0.83, y = Patient_5$score,colors = c("blue", "white", "red"))
#patient 6
kpHeatmap(kp, data = Patient_6, r0 = 0.85, r1 = 0.88, y = Patient_6$score,colors = c("blue", "white", "red"))
#patient 16
kpHeatmap(kp, data = Patient_16, r0 = 0.9, r1 = 0.93, y = Patient_16$score,colors = c("blue", "white", "red"))
#patient 19
kpHeatmap(kp, data = Patient_19, r0 = 0.95, r1 = 0.98, y = Patient_19$score,colors = c("blue", "white", "red"))

#Now want to compare H3K9ac to MEIS1 expression
#### Loading RNA microarray data ####
#loading rna microarray to determine correlation of ets factors and meis1 h3k9ac
microarray_dataset <- stemcell$`GSE38865-GPL6884_series_matrix.txt.gz`
rna_microarray <- microarray_dataset@assayData[["exprs"]]

#sample information for rna microarray
rna_microarray_sample_data <- pData(phenoData(microarray_dataset))
rna_microarray_sample_data <- rna_microarray_sample_data %>%
  dplyr::arrange(title, geo_accession)

#cleaning title in rna_microarray_sample_data 
rna_microarray_sample_data$title <- gsub("AML ", "", rna_microarray_sample_data$title)

#creating a df that converts terms for patients
#this is spaghetti code, will have to clean this up
patient_ids <- rna_microarray_sample_data %>%
  dplyr::select(title, geo_accession)

#get rid of cd34+ samples here
patient_ids <- patient_ids[1:17, ]

#sorting patients to match row order and adding chip ids to rna ids
temp_var <- patient_info %>%
  dplyr::filter(title %in% patient_ids$title) %>%
  dplyr::arrange(match(title, patient_ids$title)) %>%
  rownames_to_column(var = "sample_id") %>%
  dplyr::select(sample_id)

#combining both to make dictionary
patient_ids <- cbind(patient_ids, temp_var)

#changing colnames
colnames(patient_ids) <- c("patient", "rna", "chip")

#standardizing input of microarray
for(i in seq(nrow(rna_microarray))){
  rna_microarray[i,] <- scale(rna_microarray[i,])
}

#saving as dataframe and converting rownames to colname ID
rna_microarray <- as.data.frame(rna_microarray)
rna_microarray <- rownames_to_column(rna_microarray, var = "ID")

#loading annotation table taken from:
#https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GPL6884
#I trimmed off the extrenuous rows on top and saved as rna_microarray_annotation.csv in repo
rna_microarray_annotation <- read.csv("rna_microarray_annotation.csv", header = TRUE)

#creating a list of ETS factors 
ets <- c("CREB1", "ELF1",  "ELF2",  "ELF3",  "ELF4",  "ELK1",  "ELK3",  "ELK4",  "ERG",   "ETS1",  "ETS2", 
         "ETV3",  "ETV5",  "ETV6",  "FLI1",  "GABPA", "MYB",   "SP3",   "SPI1",  "YY1",   "YY2" )

#extracting probes that correspond only to ets factors
ets_microarray_annotation <- rna_microarray_annotation %>%
  dplyr::filter(Symbol %in% ets)

#now extracting ets factors from rna_microarray
ets_microarray <- rna_microarray %>%
  dplyr::filter(ID %in% ets_microarray_annotation$ID) %>%
  dplyr::arrange(match(ID, ets_microarray_annotation$ID)) %>%
  mutate(symbol = ets_microarray_annotation$Symbol)

#Now to correlate ets factor expression with h3k9me3 expression
#first must find the genomic region of chromosome: 66544400 - 66546800 of chr2
chr2_e2 <- chip_annotation %>%
  dplyr::filter(CHROMOSOME %in% "chr2", RANGE_START > 66544400 , RANGE_END < 66546800)

#fetching chip microarray data corresponding to this region
ch2_e2_chip <- MEIS1_chip %>%
  dplyr::filter(ID_REF %in% chr2_e2$ID)

#selecting patient samples
#first convert list GSM_names to string
chip_sample_ids <- unlist(GSM_names)

patient_e2_chip <- ch2_e2_chip %>%
  dplyr::select(patient_ids$chip)

#getting mean signal intensity across this region of genome
mean_signal_intensity <- sapply(patient_e2_chip, FUN = mean)

#correlating with FLI1
#have to find what sample ids correspond to which patient

rna_and_chip <- intersect(colnames(ets_microarray), patient_ids)

fli1 <- ets_microarray %>%
  dplyr::filter(symbol %in% "FLI1") %>%
  dplyr::select(patient_ids$rna)

fli1 <- as.numeric(fli1)

combined_df <- data.frame(
  chip = mean_signal_intensity,
  rna = fli1
)

#plotting correlation
#looks like there is an outlier here where chip  = 1.5
#chip sample is GSM950962
library(ggpubr)
ggscatter(
  combined_df, 
  x = "chip", 
  y = "rna", 
  cor.method = "spearman", 
  alternative = "two.sided", 
  conf.int = TRUE,
  add = "reg.line", 
  cor.coef = TRUE, 
  ggtheme = theme_pubr()
)

#removing outlier
combined_df_no_outlier <- combined_df[rownames(combined_df)!="GSM950962",]

ggscatter(
  combined_df_no_outlier, 
  x = "chip", 
  y = "rna", 
  cor.method = "spearman", 
  alternative = "two.sided", 
  conf.int = TRUE,
  add = "reg.line", 
  cor.coef = TRUE, 
  ggtheme = theme_pubr()
)
