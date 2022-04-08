#loading packages for chipseq visualization
library(BiocManager)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(karyoploteR)
library(org.Hs.eg.db)
library(rtracklayer)
library(BSgenome.Hsapiens.UCSC.hg38)

#Building this plot is layering it with multiple functions, starting with plotKaryotype
#The plan for this project is to look into MEIS1 regulation
#Starting with KG-1 cell line as I found ChIP-seq data for FLI1 from this paper:
#https://www.sciencedirect.com/science/article/pii/S0006497121011526?via%3Dihub#sec1
#And KG-1 has high MEIS1 expression, according to CCLE. 

#Importing and normalizing bigwig files
#Bigwig files downloaded from GEO need to have peaks called against IgG
#tutorial for normalizing files is found here: 
#https://compgenomr.github.io/book/peak-calling.html#peak-quality-control

#I'm running into the problem that you cannot normalize from bigwig files. 
#You will need bed files
#Most likely will have to process raw data myself

#importing bigwig files for analysis
wdir <- getwd()
fli1_kg1_bw <- BigWigFile(paste(wdir, "/GSE158794/FLI1_KG1.bw", sep = ""))
h3k27ac_kg1_bw <- BigWigFile(paste(wdir, "/GSE158794/H3K27ac_KG1.bw", sep =""))
fli1_me1_bw <- BigWigFile(paste(wdir, "/GSE158794/FLI1_ME1.bw", sep = ""))
h3k27ac_me1_bw <- BigWigFile(paste(wdir, "/GSE158794/H3K27ac_ME1.bw", sep = ""))

#Creating plot with overlays for KG-1 bed files
pp <- getDefaultPlotParams(plot.type=1)
pp$leftmargin <- 0.15
pp$topmargin <- 15
pp$bottommargin <- 15
pp$ideogramheight <- 5
pp$data1inmargin <- 10
pp$max <- 10

#creating granges object for meis1 e2 enhancer
MEIS1_E2_region <- toGRanges("chr2:66445000-66448000", genome = "hg38")
kp <- plotKaryotype(zoom = MEIS1_E2_region, plot.params = pp, genome = "hg38")
kpAddBaseNumbers(kp, tick.dist = 1000, minor.tick.dist = 200,
                 add.units = TRUE, digits = 6)

#Creating plots of MEIS1 locus
genes_data <- makeGenesDataFromTxDb(TxDb.Hsapiens.UCSC.hg38.knownGene,
                                    karyoplot=kp,
                                    plot.transcripts = TRUE, 
                                    plot.transcripts.structure = TRUE)
#Merging MEIS1 transcripts
genes_data <- addGeneNames(genes_data)
genes_data <- mergeTranscripts(genes_data)

kpPlotGenes(kp, data=genes_data, r0=0, r1=0.1, gene.name.cex = 0.5, 
            avoid.overlapping = TRUE)

#KG-1 cell line
kp <- kpPlotBigWig(kp, data=fli1_kg1_bw, r0=0.15, r1=0.3, ymax = 70, col = "#f2b195")
kpAxis(kp, ymin=0, ymax= 70 , r0=0.15, r1=0.3, cex = 0.5)
kpAddLabels(kp, labels = "FLI1", r0=0.15, r1=0.3, label.margin = 0.05, cex = 0.5)

kp <- kpPlotBigWig(kp, data=h3k27ac_kg1_bw, r0=0.35, r1=0.5, ymax = 70, col = "#f2b195")
kpAxis(kp, ymin=0, ymax= 70, r0=0.35, r1=0.5, cex = 0.5)
kpAddLabels(kp, labels = "H3K27ac", r0=0.35, r1=0.5, label.margin = 0.05, cex = 0.5)

#ME1 cell line
kp <- kpPlotBigWig(kp, data=fli1_me1_bw, r0=0.55, r1=0.70, ymax = 70, col = "#abd7eb")
kpAxis(kp, ymin=0, ymax=70, r0=0.55, r1=0.70, cex = 0.5)
kpAddLabels(kp, labels = "FLI1", r0=0.55, r1=0.70, label.margin = 0.05, cex = 0.5)

kp <- kpPlotBigWig(kp, data=h3k27ac_me1_bw, r0=0.75, r1=0.9, ymax = 70, col = "#abd7eb")
kpAxis(kp, ymin=0, ymax=70, r0=0.75, r1=0.9, cex = 0.5)
kpAddLabels(kp, labels = "H3K27ac", r0=0.75, r1=0.9, label.margin = 0.05, cex = 0.5)


