library(tidyverse)

library(phyloseq)

library(microbiomeMarker)

library(knitr)

library(dplyr)
library(microbiomeMarker)

library(knitr)


phy <- readRDS("/data/projects/2024/Effenberger-Diabetes/out/nf_core_ampliseq_003/phyloseq/dada2_phyloseq.rds")
print(phy)

phy

# Remove NA phylum
phy <- subset_taxa(phy, !is.na(Phylum) & Phylum != "")

# Minimum prevalence filter: keep features in >3 samples
phy <- filter_taxa(phy, function(x) sum(x > 0) > 3, TRUE)

# Abundance threshold filter (e.g., 3%)
phy <- transform_sample_counts(phy, function(x) x / sum(x))
phy <- filter_taxa(phy, function(x) mean(x > 0.03) > 0.05, TRUE) # adjust 0.03 as needed

# Create "Type" column based on pattern in "sample_information"
sample_data(phy)$Type <- ifelse(grepl("PDM", sample_data(phy)$sample_information), "PDM",
                                ifelse(grepl("K", sample_data(phy)$sample_information), "K", "DM"))


# Extract and clean taxonomy
tax <- as.data.frame(tax_table(phy))

# Keep only the 7 standard taxonomic ranks
#tax <- tax[, c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species_exact")]

# Rename 'Species_exact' to 'Species'
#colnames(tax)[colnames(tax) == "Species_exact"] <- "Species"

# Reassign as tax_table and preserve rownames
tax_fixed <- tax_table(as.matrix(tax))
rownames(tax_fixed) <- taxa_names(phy)

# Set it back into phyloseq object
tax_table(phy) <- tax_fixed

# Get sample metadata
meta <- data.frame(sample_data(phy))

# Subset to only PDM and DM
meta_sub <- meta %>% filter(Type %in% c("PDM", "DM"))
abund_sub <- abund_sel[rownames(meta_sub), ]

# Response variable
meta_sub$Type <- factor(meta_sub$Type, levels = c("DM", "PDM"))  # DM = 0, PDM = 1

