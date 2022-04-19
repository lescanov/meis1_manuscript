#March 23 2022
#Analyzing gRNA for MEIS1 E2.2 for FLI1 binding
#First will start with the exact gRNA sequence
library(JASPAR2020)
library(TFBSTools)
library(Biostrings)
library(tidyverse)

#Retriving TF binding matrix from JASPAR
FLI1_pfm <- getMatrixByName(JASPAR2020, name = "FLI1")

#Converting pfm to pwm
FLI1_pwm <- toPWM(FLI1_pfm)

#Plotting sequence logo
seqLogo(toICM(FLI1_pfm))

#Testing binding for FLI1 on MEIS1 E2.2 gRNA region
E2.2 <- DNAString("ATACTAGGCGGTATCCCGGA")
#                 CTTTGATACTAGGCGGTATCCCGGAGGGCT
siteseq <- searchSeq(FLI1_pwm, E2.2, seqname = "E2.2")
head(writeGFF3(siteseq))

#I want to discover all transcription factors that can bind region
#will create a list of human tfs outlined in jaspar
opts <- list()
opts[["species"]] <- 9606 #this is the taxonomy ID for humans
human_tfs <- getMatrixSet(JASPAR2020, opts)

#converting human_tfs to pwm
human_tfs_pwm <- toPWM(human_tfs)

#running tf binding analysis on e2.2 grna sequence, on +/- 5bp of locus
#writing a function that will write a csv for numerous scores going down by intervals of .1
FiveUpDownstream <- DNAString("CTTTGATACTAGGCGGTATCCCGGAGGGCT")
score_sequence <- seq(from = 0.8, to = 0, by = -0.1)

for (i in seq_along(score_sequence)){
  score_tf <- searchSeq(human_tfs_pwm, FiveUpDownstream, seqname = "E2.2-5bp-up-downstream", 
                        min.score = score_sequence[i])
  result_tf <- writeGFF3(score_tf)
  result_tf <- result_tf %>%
    mutate(
      gene_symbol = result_tf$attributes
    )
  result_tf$gene_symbol <- gsub(";.*", "", result_tf$gene_symbol)
  result_tf$gene_symbol <- gsub("TF=", "", result_tf$gene_symbol)
  
  write.csv(result_tf, paste(score_sequence[i], ".csv", sep ="" ))
}

