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
#Could use the enrichment tool on the webiste here:
#https://bitbucket.org/CBGR/jaspar_enrichment/src/master/
#Or could just list all transcription factors in JASPAR
opts <- list()
opts[["species"]] <- 9606 #this is the taxonomy ID for humans
human_tfs <- getMatrixSet(JASPAR2020, opts)

#converting human_tfs to pwm
human_tfs_pwm <- toPWM(human_tfs)

#Looking for binding to E2.2 region
siteseq_human_tfs <- searchSeq(human_tfs_pwm, E2.2, seqname = "E2.2")
head(writeGFF3(siteseq_human_tfs))
result <- writeGFF3(siteseq_human_tfs)

write.csv(result, "JASPAR2020_TF_binding_MEIS1_E2_2.csv")

#This result shows that there is MEIS1

#### Running analysis on the area around E2.2 gRNA ####

E2.2_region_around <- DNAString("GTTGGCGGCACGGTTTATTTTATTTTACCTGTTTTCCTCGGAGGGCGCGAAGACTGCCACCCGCGCGGGGACCTGGGATCGACGACTTTGATACTAGGCGGTATCCCGGAGGGCTAAGTCGGCGGAAATCCACTTGACCTTGTAGCGTTAGTCCTTTCTTTTCCTTTCCTTTCCTTTTTCTTTCTTCTCTCTTTCCTATT")

#Running analysis with min score of 0
minscorezero <- searchSeq(human_tfs_pwm, E2.2, seqname = "E2.2", min.score = 0.8)
min_result <- writeGFF3(minscorezero)

region_around_E2_score <- searchSeq(human_tfs_pwm, E2.2_region_around, seqname = "E2.2_around", min.score = 0.8)
proximity_result <- writeGFF3(region_around_E2_score)
proximity_result <- proximity_result %>%
  mutate(
    gene_symbol = proximity_result[["attributes"]]
  )
proximity_result$gene_symbol <- gsub(";.*", "", proximity_result$gene_symbol)
proximity_result$gene_symbol <- gsub("TF=", "", proximity_result$gene_symbol)

write.csv(proximity_result, "100bp_upstream_downstream_E2gRNA.csv")

#Running analysis on +/- 5bp E2.2 gRNA region
FiveUpDownstream <- DNAString("CTTTGATACTAGGCGGTATCCCGGAGGGCT")
score_five <- searchSeq(human_tfs_pwm, FiveUpDownstream, seqname = "E2.2", min.score = 0.75)
result_five <- writeGFF3(score_five)

result_five <- result_five %>%
  mutate(
    gene_symbol =  result_five$attributes
  )
result_five$gene_symbol <- gsub(";.*", "", result_five$gene_symbol)
result_five$gene_symbol <- gsub("TF=", "", result_five$gene_symbol)

write.csv(result_five, "JASPAR_E2_5bp_updownstream.csv")

#### Writing a function that will write a csv for numerous score going down by intervals of .1

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

