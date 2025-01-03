library(dada2);# packageVersion("dada2")

path <- "/data/projects/2024/Effenberger-Diabetes/data/20011/Fastq" 
figures_out = "/data/projects/2024/Effenberger-Diabetes/out/figures/"
list.files(path)

# Forward and reverse fastq filenames have format: SAMPLENAME_R1_001.fastq and SAMPLENAME_R2_001.fastq
fnFs <- sort(list.files(path, pattern="_R1_001.fastq", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="_R2_001.fastq", full.names = TRUE))
# Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)

#Inspect read quality profiles
p_forward = plotQualityProfile(fnFs[1:2])
p_reverse = plotQualityProfile(fnRs[1:2])
#png(filename=paste0(figures_out,"read_quality_profiles.png"))
#plot(p)
#dev.off()


#Filter and trim
# Place filtered files in filtered/ subdirectory
filtFs <- file.path(path,"filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(path,"filtered", paste0(sample.names, "_R_filt.fastq.gz"))
names(filtFs) <- sample.names
names(filtRs) <- sample.names


out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(240,160),
                     maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,
                     compress=TRUE, multithread=TRUE) # On Windows set multithread=FALSE
head(out)

#Learn the Error Rates
errF <- learnErrors(filtFs, multithread=TRUE)
errR <- learnErrors(filtRs, multithread=TRUE)
plotErrors(errF, nominalQ=TRUE)

#Sample Inference
dadaFs <- dada(filtFs, err=errF, multithread=TRUE)
