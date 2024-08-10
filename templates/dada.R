library("ggplot2")
args <-  commandArgs(TRUE)

plot_list = list()

library(dada2); packageVersion("dada2")

# dada parameters
setDadaOpt(DETECT_SINGLETONS = FALSE, USE_KMERS = TRUE, OMEGA_A = 1e-40, BAND_SIZE = 9 )

# Max Espected Error parameter 
max_EE <- 8

# file paths from NF
Samplefq <- file.path(args[1])
SampleKfq <- file.path(args[3]) 
samplemusclfq <- file.path(args[2]) 

# file names for saving plots and tables
Samplename <- basename(args[1])
SampleKname<- basename(args[3])
Samplemuscle<- basename(args[2])
#Quality profile check 

plot1 <-plotQualityProfile(Samplefq)
ggsave(paste0(Samplename,".pdf"),plot = plot1,device = "pdf")
plot2 <-plotQualityProfile(SampleKfq)
ggsave(paste0(SampleKname,".pdf"),plot = plot2,device = "pdf")
plot3 <- plotQualityProfile(samplemusclfq)
ggsave(paste0(Samplemuscle,".pdf"),plot = plot3,device = "pdf")

#Quality filtering          
            
filtSample <- file.path( paste0(Samplename,"_filtered.fastq"))
filtKSample <- file.path( paste0(SampleKname,"_filtered.fastq"))
filtmuscleSample <- file.path(paste0(Samplemuscle,"_filtered.fastq"))

#quality filtering check
out <- filterAndTrim(Samplefq, filtSample,maxN=0,maxEE=max_EE,compress=FALSE, multithread=FALSE)
outk <- filterAndTrim(SampleKfq, filtKSample,maxN=0,maxEE=max_EE,compress=FALSE, multithread=FALSE)
outmuscle <- filterAndTrim(samplemusclfq, filtmuscleSample,maxN=0,maxEE=max_EE,compress=FALSE, multithread=FALSE)

# Quality check po filtrowaniu 

plot4 <- plotQualityProfile(filtSample)
ggsave(paste0(Samplename,"filtered.pdf"),plot = plot4,device = "pdf")
plot5 <- plotQualityProfile(filtKSample)
ggsave(paste0(SampleKname,"filtered.pdf"),plot = plot5,device = "pdf")
plot6 <- plotQualityProfile(filtmuscleSample)
ggsave(paste0(Samplemuscle,"filtered.pdf"),plot = plot6,device = "pdf")


#Dada ASV call control first
errF <- learnErrors(filtKSample, multithread=TRUE)
plot7 <- plotErrors(errF, nominalQ=TRUE)
ggsave(paste0(SampleKname,"errorates.pdf"),plot = plot7,device = "pdf")

dadaFs_k <- dada(filtKSample, err=errF, USE_QUALS = TRUE, multithread=TRUE, verbose = 2,selfConsist =FALSE)
dadaFs_k$clustering
seqtab_k <- makeSequenceTable(list(control = dadaFs_k))  

#DADA ASV call 1st biological sample with errorates from control           
dadaFs <- dada(filtSample, err=errF, USE_QUALS = TRUE, multithread=TRUE, verbose = 2,selfConsist =FALSE)
dadaFs$clustering
seqtab <- makeSequenceTable(list(sample = dadaFs))
#dim(dadaFs)

# DADA ASV call 2nd biological sample with errorates from control
dadaFsmu <- dada(filtmuscleSample, err=errF, USE_QUALS = TRUE, multithread=TRUE, verbose = 2,selfConsist =FALSE)
dadaFsmu$clustering
seqtabmuscle <- makeSequenceTable(list(muscle = dadaFsmu))

# merging ASV tables
merged_seq_tab <- mergeSequenceTables(seqtab_k, seqtab,seqtabmuscle, repeats = "error",orderBy = "abundance")

# writing table
table <-write.csv(merged_seq_tab,paste0(Samplename,'_sequence_table'))



