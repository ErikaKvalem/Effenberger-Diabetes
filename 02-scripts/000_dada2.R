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
#plotErrors(errF, nominalQ=TRUE)

#Sample Inference
dadaFs <- dada(filtFs, err=errF, multithread=TRUE)
dadaRs <- dada(filtRs, err=errR, multithread=TRUE)
#inspect 
dadaFs[[1]]

#Merge paired reads 
mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose=TRUE)
# Inspect the merger data.frame from the first sample
head(mergers[[1]])

#Construct sequence table
seqtab <- makeSequenceTable(mergers)
dim(seqtab)

# Inspect distribution of sequence lengths
table(nchar(getSequences(seqtab)))

#Remove chimeras
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
dim(seqtab.nochim)

sum(seqtab.nochim)/sum(seqtab)

#Track reads through the pipeline
getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))
# If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names
head(track)


#Assign taxonomy
taxa <- assignTaxonomy(seqtab.nochim, "/data/projects/2024/Effenberger-Diabetes/data/20011/Fastq/silva/silva_nr99_v138.2_toSpecies_trainset.fa.gz", multithread=TRUE)


#inspect taxonomic assignment
taxa.print <- taxa # Removing sequence rownames for display only
rownames(taxa.print) <- NULL
head(taxa.print)

#Handoff to phyloseq
library(phyloseq); packageVersion("phyloseq")
library(Biostrings); packageVersion("Biostrings")
library(ggplot2); packageVersion("ggplot2")

theme_set(theme_bw())

# simple sample data.frame from the information encoded in the filenames. 
samdf <- read.csv(paste0(path,"/20011_SampleInfo.csv"))
samples.out <- rownames(seqtab.nochim)
# Assuming your dataframe is named 'df'
samdf <- rbind(samdf, data.frame(IMGM.ID = "20002-0053", Sample.Information = "Control", Type = "Control"))
# Assuming your dataframe is named 'df'
samdf <- rbind(samdf, data.frame(IMGM.ID = "20002-0054", Sample.Information = "Control", Type = "Control"))
rownames(samdf) <- samdf$IMGM.ID
rownames(samdf) <- samples.out

#construct a phyloseq object directly from the dada2 outputs.
ps <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows=FALSE), 
               sample_data(samdf), 
               tax_table(taxa))
ps <- prune_samples(sample_names(ps) != "Mock", ps) # Remove mock sample


dna <- Biostrings::DNAStringSet(taxa_names(ps))
names(dna) <- taxa_names(ps)
ps <- merge_phyloseq(ps, dna)
taxa_names(ps) <- paste0("ASV", seq(ntaxa(ps)))
ps

#Visualize alpha-diversity:
plot_richness(ps, x="Type", measures=c("Shannon", "Simpson"), color="Sample.Information")


## Transform data to proportions as appropriate for Bray-Curtis distances
ps.prop <- transform_sample_counts(ps, function(otu) otu/sum(otu))

# Access the OTU table
otu_table_matrix <- as.matrix(ps.prop@otu_table)

# Replace NaN with 0
otu_table_matrix[is.nan(otu_table_matrix)] <- 0

# Update the OTU table in the phyloseq object
ps.prop@otu_table <- otu_table(otu_table_matrix, taxa_are_rows = taxa_are_rows(ps.prop))

#ord.nmds.bray <- ordinate(ps.prop, method="NMDS", distance="bray")

plot_ordination(ps.prop, ord.nmds.bray, color="Sample.Information", title="Bray NMDS")

#Bar plot:
top20 <- names(sort(taxa_sums(ps), decreasing=TRUE))[1:20]
ps.top20 <- transform_sample_counts(ps, function(OTU) OTU/sum(OTU))
ps.top20 <- prune_taxa(top20, ps.top20)
plot_bar(ps.top20, x="Family", fill="Sample.Information")   + facet_wrap(~Type, scales="free_y")
plot_bar(ps.top20, x="Genus", fill="Sample.Information")   + facet_wrap(~Type, scales="free_y")



##########fastq to fasta 
whereitsat <- "/data/projects/2024/Effenberger-Diabetes/data/20011/Fastq/"
inny <- dir(whereitsat, pattern = "fastq.gz$")
outie <- paste0(whereitsat, "/",  gsub("fastq.gz","fa", inny))

for(i in seq(along = inny)) writeFasta(readFastq(whereitsat, inny[i]), outie[i])
