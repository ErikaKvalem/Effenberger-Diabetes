library(MicrobiomeAGORA)
library(ggplot2)
source("dat/init.agora.R")
rm(agora)

library(data.table)

# Load your abundance table (ASV IDs as rownames)
otu <- fread("/data/projects/2024/Effenberger-Diabetes/out/nf_core_ampliseq_003/qiime2/abundance_tables/feature-table-copy.tsv")
setnames(otu, 1, "FeatureID")

# Load the feature-to-sequence map
seq.map <- fread("/data/projects/2024/Effenberger-Diabetes/out/nf_core_ampliseq_003/qiime2/representative_sequences/feature-id-to-seq.tsv", header = FALSE)
setnames(seq.map, c("FeatureID", "Representative_Sequence"))

# Merge
otu.seq <- merge(seq.map, otu, by = "FeatureID")

# Drop FeatureID and move Representative_Sequence to the front
otu.seq[, FeatureID := NULL]
setcolorder(otu.seq, c("Representative_Sequence", setdiff(names(otu.seq), "Representative_Sequence")))

# Write out in correct format
fwrite(otu.seq, "/data/scratch/kvalem/projects/2024/Effenberger-Diabetes/02-scripts/tables/otutab.txt", sep = "\t")



mic <- new("Microbiome",
           uniq.table.file = "/data/scratch/kvalem/projects/2024/Effenberger-Diabetes/02-scripts/tables/otutab.txt",
           agora.mapping.file = "dat/16s/agora.mapping.output.m8",
           sample.description.file = "/data/scratch/kvalem/projects/2024/Effenberger-Diabetes/02-scripts/metadata_abundance_copy.txt"
)

mic <- filter.mapping(mic,allowed.identity.deviation=0,method.resolve.multiple = "user") # always choose 1 if asked
mic <- create.AgoraTable(mic)
mic <- filter.samples(mic, min.seqs = 250, max.unclassified = 0.5)
