#purpose of this script is to visualize h3k9ac across the e2 enhances region of MEIS1
#then relate the expression of FLI1 to said acetylation.
#will compare MEIS1-high subtype of Normal Karyotype AMLs against other subtypes

#the overall workflow plan for this script is as follows:
# >1. extract patient data from geodataset to determine patient subgroups
# >2. extract microarray data from geodataset, which comes in form of .txt files
# >3. convert .txt files into bedgraph files
# >4. visualize bedgraph files as heatmap on karyoploteR
# >5. correlate FLI1 expression with H3K9ac of E2 enhancer region in AML patients

library(GEOquery)
library(tidyverse)
library(rtracklayer)
library(BSgenome)
library(TxDb.Hsapiens.UCSC.hg18.knownGene)
library(karyoploteR)
library(org.Hs.eg.db)

#### DOWNLOADING GEODATASET FOR H3K9ac IN AML PATIENTS AND PATEITN CHARACTERISTICS ####
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

#extracting patient data from geo dataset using pData
histone <- stemcell$`GSE38865-GPL15724_series_matrix.txt.gz`
patient_info <- pData(phenoData(histone))

#the characteristics_ch1.5 column contains patient karyotype
#plan to divide patients into normal karyotype and non-normal karyotype
patient_info <- patient_info %>%
  dplyr::select(title, characteristics_ch1.5, characteristics_ch2.8, characteristics_ch2.9)

#cleaning title column of patient_info
patient_info$title <- gsub(" H3K9ac", "", patient_info$title)

#### STANDARDIZING H3K9ac SIGNALS ACROSS PATIENTS ####
#importing chiparray.txt (h3k9ac microarray data)
chip_array <- read.delim("chiparray.txt", header = TRUE)

#importing genome annotation for chip array which was received from here:
#https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GPL15724
#annotation ID are not in order, must be sorted.
chip_annotation <- read.delim("chip_annotation.txt")
chip_annotation <- chip_annotation %>%
  dplyr::arrange(CHROMOSOMAL_LOCATION)

#extracting MEIS1 locus from chip_annotation
MEIS1 <- chip_annotation %>%
  dplyr::filter(GENE_SYMBOL %in% "MEIS1")

#extracting MEIS1 locus from microarray data
MEIS1_chip <- chip_array %>%
  dplyr::filter(ID_REF %in% MEIS1$ID)

#standardizing score for MEIS1
#temporararily saving column (asscension values) to be combined later
#matrix cannot contain rownames, iirc, so this is my solution
tmp <- MEIS1_chip[1]

#removing names so can standardize
MEIS1_chip <- MEIS1_chip[-1]

#getting number of rows
nrow <- seq(nrow(MEIS1_chip))

#inverting axis so that signal intensities are on rows and patients as columns
MEIS1_chip <- as.matrix(MEIS1_chip)

#applying scale across each row of chip_array, standardizing h3k9ac signal intensity across samples
for (i in seq_along(nrow)){
  MEIS1_chip[i,] <- scale(MEIS1_chip[i,])
}

#saving as dataframe and adding back asscension numbers
MEIS1_chip <- as.data.frame(MEIS1_chip)
MEIS1_chip <- cbind(tmp, MEIS1_chip)

#### CONVERT TXT FILES INTO BEDGRAPH FILES FOR VISUALIZATION ####
#writing functions to convert txt files to bedgraph

create_txt_ready_for_bed <- function(GSM){
  
  #this creates a dataframe that fits the bed format as specified by ucsc
  #https://genome.ucsc.edu/FAQ/FAQformat.html#format1.7
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

#load_bedgraph function now creates a bedgraph from the output of create_txt_ready_for_bed, and
#imports it as a granges object bedgraph

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
GSM_names <- as.list(names(MEIS1_chip[, 2:27])) #GSM_names correspond to AML patient samples
sample_names <- gsub(" ", "_", patient_info$title) #sample_names are names given to each element of list that correspond to a patient sample
patient_names <- sample_names[1:26] #patient_names are strings for aml patients
patient_names <- as.list(patient_names) #convert to list so that it may be used in bedgraph_from_txt function

#lapply list of patient names to 
histone_bed <- bedgraph_from_txt(GSM_names, patient_names)
list2env(histone_bed, .GlobalEnv)

#Now stating for CD34+ control
GSM_cd34_names <- as.list(names(MEIS1_chip[, 28:37]))
sample_cd34_names <- sample_names[27:36]
sample_cd34_names <- as.list(sample_cd34_names)

cd34_chip <- bedgraph_from_txt(GSM_cd34_names, sample_cd34_names)
list2env(cd34_chip, .GlobalEnv)

#defining neutrophil control
GSM_neutrophil_names <- as.list(names(MEIS1_chip[, 38:48]))
sample_neutrophil_names <- sample_names[37:47]
sample_neutrophil_names <- as.list(sample_neutrophil_names)

neutrophil_chip <- bedgraph_from_txt(GSM_neutrophil_names, sample_neutrophil_names)
list2env(neutrophil_chip, .GlobalEnv)

#defining t cell control
GSM_tcell_names <- as.list(names(MEIS1_chip[, 49:57]))
sample_tcell_names <- sample_names[48:56]
sample_tcell_names <- as.list(sample_tcell_names)

tcell_chip <- bedgraph_from_txt(GSM_tcell_names, sample_tcell_names)
list2env(tcell_chip, .GlobalEnv)
#determining max and min score for each patient to determine get a grasp on data distribution
for(i in seq_along(histone_bed)){
  print(c(max(score(histone_bed[[i]])), min(score(histone_bed[[i]]))))
}

#### Plotting H3K9ac along E2 enhancer region ####
#setting E2 region 
#old region is 
MEIS1_E2_region <- toGRanges("chr2:66525936-66527140", genome = "hg18")

#setting plot parameters for karyoploteR
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

#can do this to compare normal karyotype to cd34+ healthy controls
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

#CD34+ normal control
#normal donor 1
kpHeatmap(kp, data = CD34_normal_control_1, r0 = 0.7, r1 = 0.73, y = CD34_normal_control_1$score,colors = c("blue", "white", "red"))
#CD34_normal_control 3
kpHeatmap(kp, data = CD34_normal_control_3, r0 = 0.75, r1 = 0.78, y = CD34_normal_control_3$score,colors = c("blue", "white", "red"))
#CD34_normal_control 4
kpHeatmap(kp, data = CD34_normal_control_4, r0 = 0.8, r1 = 0.83, y = CD34_normal_control_4$score,colors = c("blue", "white", "red"))
#CD34_normal_control 5
kpHeatmap(kp, data = CD34_normal_control_5, r0 = 0.85, r1 = 0.88, y = CD34_normal_control_5$score,colors = c("blue", "white", "red"))
#CD34_normal_control 7
kpHeatmap(kp, data = CD34_normal_control_7, r0 = 0.9, r1 = 0.93, y = CD34_normal_control_7$score,colors = c("blue", "white", "red"))


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
  dplyr::select(sample_id,  characteristics_ch1.5, characteristics_ch2.8, characteristics_ch2.9)

#combining both to make dictionary of patients with sample ids for both chip and rna microarrays
patient_ids <- cbind(patient_ids, temp_var)

#changing colnames
colnames(patient_ids) <- c("patient", "rna", "chip", "karyotype", "npm1", "flt3")

#plotting distribution of karyotypes
patient_ids %>%
  dplyr::group_by(karyotype) %>%
  ggplot(aes(x = karyotype)) +
  geom_bar(stat = "count", position = "identity") +
  ggtitle("karyotype distribution of patients with RNA microarray data") +
  theme(plot.title = element_text(hjust = 2)) +
  coord_flip() 

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
#will focus on FLI1 for now
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

#Now to correlate ets factor expression with h3k9ac expression
#use this region "chr2:66525936-66527140" hg18
#first must find the genomic region of chromosome: 66544400 - 66546800 of chr2
chr2_e2 <- chip_annotation %>%
  dplyr::filter(CHROMOSOME %in% "chr2", RANGE_START > 66525936, RANGE_END < 66527140)

#fetching chip microarray data corresponding to this region
ch2_e2_chip <- MEIS1_chip %>%
  dplyr::filter(ID_REF %in% chr2_e2$ID)

#selecting patient samples
#first convert list GSM_names to string
chip_sample_ids <- unlist(GSM_names)

#selecting patients that have chip-seq data
patient_e2_chip <- ch2_e2_chip %>%
  dplyr::select(patient_ids$chip)

#getting mean signal intensity across this region of genome
aml_mean_signal_intensity <- sapply(patient_e2_chip, FUN = mean)

#correlating with FLI1
#have to find what sample ids correspond to which patient
#identifying samples that are present in both rna and chip microarrays
rna_and_chip <- intersect(colnames(ets_microarray), patient_ids)

#extracting fli1 from rna microarray
fli1 <- ets_microarray %>%
  dplyr::filter(symbol %in% "FLI1") %>%
  dplyr::select(patient_ids$rna)
fli1 <- as.numeric(fli1)

#creating a dataframe that contains rna and chip microarray data, for correlation
rna_chip_fli1 <- data.frame(
  chip = aml_mean_signal_intensity,
  rna = fli1
)

#identifying the presence of outliers
#creating a function to plot cooks distance
plot_cooks_distance <- function(cooks_distance){
  plot(cooks_distance, pch="*", cex=2, main="Influential Obs by Cooks distance")  # plot cook's distance
  abline(h = 4*mean(cooks_distance, na.rm=T), col="red")  # add cutoff line
  text(x=1:length(cooks_distance)+1, y=cooks_distance, labels=ifelse(cooks_distance>4*mean(cooks_distance, na.rm=T),names(cooks_distance),""), col="red")  # add labels
}

#fitting lm to plot cook's distance of rna and chip microarray
#this identifies GSM950962 as outlier
chip_rna_lm <- lm(chip ~ rna, data = rna_chip_fli1)
cookd_chip_rna <- cooks.distance(chip_rna_lm)
plot_cooks_distance(cookd_chip_rna)

#removing outlier
rna_chip_fli1 <- rna_chip_fli1[rownames(rna_chip_fli1) != "GSM950958",]
rna_chip_fli1 <- rna_chip_fli1[rownames(rna_chip_fli1) != "GSM950964",]
rna_chip_fli1 <- rna_chip_fli1[rownames(rna_chip_fli1) != "GSM950957",]

#looking to plot relationship of only normal karyotype samples
normal_karyotype <- c("GSM950948", "GSM950957", "GSM950962", "GSM950964", "GSM950965", "GSM950972")

nk <- rna_chip_fli1 %>% rownames_to_column(var = "x")
nk <- nk %>% dplyr::filter(x %in% normal_karyotype)

#plotting correlation
#chip sample is GSM950962
library(ggpubr)
ggscatter(
  rna_chip_fli1, 
  x = "chip", 
  y = "rna", 
  cor.method = "kendall", 
  alternative = "two.sided", 
  xlab = "H3K9ac mean signal intensity", 
  ylab = "FLI1 signal intensity",
  conf.int = TRUE,
  add = "reg.line", 
  cor.coef = TRUE, 
  ggtheme = theme_pubr()
)

ls#
df <- ch2_e2_chip
df <- df[-1]

#creating df for plotting
df <- df %>% dplyr::select(1, 7, 10, 13, 15, 18, 22, 23, 25, 26, 2, 3, 5, 6, 16, 19, 28, 30, 32, 34, 36, 38, 40, 42, 44, 46, 48, 50, 52, 54, 56)

#creating group labels
group_labels <- c(rep("Normal Karyotype", 10), rep("Non-Normal Karyotype", 6), rep("CD34+ normal", 5), rep("neutrophil normal", 5), rep("T cell normal", 5))

#plotting h3k9ac levels across different tissue types
#fist getting mean signal intensity for all samples
mean_signal_intensity <- sapply(df, FUN = mean)

#creating a dataframe for plotting
to_plot <- data.frame(
  group = group_labels, 
  mean_signal_intensity = mean_signal_intensity
)

give.n <- function(x){
  return(c(y = mean(x), label = length(x)))
}

ggplot(to_plot, aes(x = group, y = mean_signal_intensity)) +
  geom_boxplot() +
  stat_summary(fun.data = give.n, geom = "text") +
  ggtitle("mean signal intensity of MEIS1 E2 region") +
  theme_minimal() + theme(plot.title = element_text(hjust = 0.5))

