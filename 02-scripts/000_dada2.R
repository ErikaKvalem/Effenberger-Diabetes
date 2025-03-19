library(dada2);# packageVersion("dada2")

library(dplyr)

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

#`%not_in%` <- purrr::negate(`%in%`)
df_taxa <- as.data.frame(taxa, stringsAsFactors = FALSE)
#
df_taxa_f=df_taxa[!is.na(df_taxa$Kingdom),]
df_taxa_f=df_taxa[!is.na(df_taxa$Phylum),]
df_taxa_f=df_taxa[!is.na(df_taxa$Class),]
df_taxa_f=df_taxa[!is.na(df_taxa$Order),]
df_taxa_f=df_taxa[!is.na(df_taxa$Family),]
df_taxa_f=df_taxa[!is.na(df_taxa$Genus),]
df_taxa_f=df_taxa[!is.na(df_taxa$Species),]



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

samdf_filtered <- samdf %>%
  filter(!grepl("Control|Kontrolle", Type))
# Modify IMGM.ID: Replace "-" with "_" and add "s" at the beginning
samdf_filtered$IMGM.ID <- paste0("s", gsub("-", "_", samdf_filtered$IMGM.ID))

# Set the modified IMGM.ID as row names
rownames(samdf_filtered) <- samdf_filtered$IMGM.ID
#
otu_names <- colnames(seqtab.nochim)
#
## Taxa names from df_taxa_f (assuming they are set as rownames)
taxa_names <- rownames(df_taxa_f)
#
## See which names are common
common_taxa <- intersect(otu_names, taxa_names)
#
common_taxa <- intersect(colnames(seqtab.nochim), rownames(df_taxa_f))
#
## Check if any common taxa exist
if (length(common_taxa) == 0) {
  stop("No matching taxa names found between the OTU table and taxonomy table.")
}
#
## Subset the OTU table to keep only columns corresponding to common taxa
seqtab.nochim_filtered <- seqtab.nochim[, common_taxa, drop = FALSE]
#
## Subset the taxonomy table to keep only rows corresponding to common taxa
df_taxa_f_filtered <- df_taxa_f[common_taxa, , drop = FALSE]
#
## Convert the filtered taxonomy data frame to a matrix
tax_mat <- as.matrix(df_taxa_f_filtered)

ps <- phyloseq(
  otu_table(seqtab.nochim_filtered, taxa_are_rows = FALSE),
  sample_data(samdf_filtered),
  tax_table(tax_mat)
)

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


top20 <- names(sort(taxa_sums(ps), decreasing=TRUE))[1:100]


ps.top20 <- transform_sample_counts(ps, function(OTU) OTU/sum(OTU))
ps.top20 <- prune_taxa(top20, ps.top20)
plot_bar(ps.top20, x="Family", fill="Sample.Information")   + facet_wrap(~Type, scales="free_y")
plot_bar(ps.top20, x="Genus", fill="Sample.Information")   + facet_wrap(~Type, scales="free_y")
plot_bar(ps.top20, x="Species", fill="Sample.Information")   + facet_wrap(~Type, scales="free_y")



plot_bar(ps.top20, fill="Genus", x="Sample.Information") 


library(phyloseq)
library(ggplot2)

# Convert phyloseq object to a data frame for ggplot2
ps_melt <- psmelt(ps.top20)



###############################################################
scale_y_continuous(labels = scales::percent)
library(scales)

get_taxa_unique(ps, "Family")
sample_data(ps)$Type<-factor(sample_data(ps)$Type)
sample_data(ps)$Sample.Information<-factor(sample_data(ps)$Sample.Information)

ps


phy2 <- ps %>% aggregate_taxa(level = "Phylum") %>%
  
  microbiome::transform(transform = "compositional")

phy2 %>%  plot_composition()+
  scale_y_continuous(labels = scales::percent) 

phy2 %>% 
  plot_composition() +
  scale_y_continuous(labels = scales::percent) +
  labs(color = "Phylum", fill = "Phylum")

phy2 %>%  plot_composition(group_by = "Type")+ scale_y_continuous(labels = percent)


# Remove taxa where Phylum is NA or "Unknown"

# Remove taxa where any rank has "Unknown" or NA
phy2_clean <- subset_taxa(phy2, apply(tax_table(phy2), 1, function(x) all(!is.na(x) & x != "Unknown")))

#phy2_clean <- subset_taxa(phy2, !is.na(Kingdom) & Kingdom != "Unknown")
#phy2_clean <- subset_taxa(phy2, !is.na(Phylum) & Phylum != "Unknown")
#phy2_clean <- subset_taxa(phy2, !is.na(Class) & Class != "Unknown")
#phy2_clean <- subset_taxa(phy2, !is.na(Order) & Order != "Unknown")
#phy2_clean <- subset_taxa(phy2, !is.na(Family) & Family != "Unknown")


phy2_clean %>%  plot_composition()+
  scale_y_continuous(labels = scales::percent) 

############################### out of ampliseq pipeline 

ps = readRDS("/data/projects/2024/Effenberger-Diabetes/out/nf_core_ampliseq/phyloseq/dada2_phyloseq.rds")
