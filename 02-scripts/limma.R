# Load libraries
library(limma)
library(edgeR) # For voom normalization
library(dplyr)
library(tidyverse)



# Leer el archivo con fila de taxones como rownames
counts_raw <- read.csv(
  "/data/projects/2024/Effenberger-Diabetes/out/nf_core_ampliseq_003/qiime2/barplot/level-6.csv",
  row.names = 1,
  check.names = FALSE
)

cols_to_remove <- c(
  "sample_information", "age", "KHK1", "KHK2", "CA1", "CA2",
  "HbA1C (DCCT/NGSP)1", "HbA1C (DCCT/NGSP)2", "Glukose1", "Glukose2",
  "BMI1", "BMI2", "Pankreatektomie", "id", "Type",
  "CA_diff", "KHK_diff", "BMI_diff", "Glukose_diff", "HbA1C_diff"
)

# Elimina esas columnas (asegura que existen)
counts_raw <- counts_raw[, !(names(counts_raw) %in% cols_to_remove)]

# Transponer para tener taxones como filas y muestras como columnas
counts_t <- t(counts_raw)

keep <- rowSums(counts_t > 20) >= 2
counts_t <- counts_t[keep, ]

# Convertir a matriz numérica asegurando que los valores son números
counts_numeric <- apply(counts_t, 2, function(x) as.numeric(as.character(x)))

# Restaurar rownames (sample IDs)
rownames(counts_numeric) <- rownames(counts_t)
counts_numeric[is.na(counts_numeric)] <- 0

counts_numeric <- t(counts_numeric)



meta<- read.csv(
  "/data/scratch/kvalem/projects/2024/Effenberger-Diabetes/02-scripts/metadata_abundance.csv"
)
rownames(meta) <- meta$id

meta <- meta %>%
  select(
    sample_information, age, KHK1, KHK2, CA1, CA2,
    HbA1C..DCCT.NGSP.1, HbA1C..DCCT.NGSP.2, Glukose1, Glukose2,
    BMI1, BMI2, Pankreatektomie_encoded, id, Type
  )

meta <- meta %>%
  mutate(Type = recode(Type,
                       "pankreopriver Diabetes" = "PDM",
                       "Diabetes mellitus Typ1" = "DM",
                       "Kontrolle" = "K"))

# Define your design matrix
#meta$Type <- factor(meta$Type, levels = c("pankreopriver Diabetes", "Diabetes mellitus Typ1"))
meta <- meta %>%
  filter(Type %in% c("PDM", "K")) %>%
  filter(id %in% rownames(counts_numeric))

meta$Type <- factor(meta$Type, levels = c("K", "PDM"))
meta$Type <- relevel(factor(meta$Type), ref = "K")


design <- model.matrix(~ Type , data = meta)



counts_filtered <- counts_numeric[rownames(counts_numeric) %in% meta$id, ]
counts_filtered <- t(counts_filtered)
dge <- DGEList(counts = counts_filtered)
# Voom transformation
v <- voom(dge, design, plot = TRUE)

fit <- lmFit(v, design)
fit <- eBayes(fit)
results <- topTable(fit, coef = "TypePDM"  , number = Inf, adjust = "fdr")


results_sig_raw <- results[results$P.Value < 0.01, ]
head(results_sig_raw[order(results_sig$P.Value), ])

library(ggplot2)

results$significant <- results$P.Value < 0.01

ggplot(results, aes(x = logFC, y = -log10(P.Value), color = significant)) +
  geom_point(alpha = 0.7) +
  scale_color_manual(values = c("grey", "red")) +
  theme_minimal() +
  labs(title = "Volcano Plot", x = "log2 Fold Change", y = "-log10(p-value)")

subset(results, abs(logFC) > 1 & P.Value < 0.1)


####################################ANCOM-BC Tutorial
#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install("ANCOMBC")

######################################corncob
#remotes::install_github("statdivlab/corncob")

library(corncob)
library(phyloseq)
library(magrittr)

#library(ANCOMBC)

pseq <- readRDS("/data/projects/2024/Effenberger-Diabetes/out/nf_core_ampliseq_003/phyloseq/dada2_phyloseq.rds")

sample_df <- as.data.frame(sample_data(pseq))

# Create new Type column based on sample_information
sample_df$Type <- ifelse(grepl("K", sample_df$sample_information), "K",
                         ifelse(grepl("PDM", sample_df$sample_information), "PDM", "DM"))

sample_data(pseq) <- sample_df



pseq <- pseq %>% 
  phyloseq::subset_samples(Type %in% c("K","DM")) %>%
  phyloseq::tax_glom("Phylum") 


corncob <- bbdml(formula = OTU.1 ~ 1,
                 phi.formula = ~ 1,
                 data = soil)

otu_counts <- as(otu_table(pseq), "matrix")
rowSums(otu_counts["3f3a0eaeea9c0690b6ede1b17b4fd8ce", ])

sum(otu_counts["3f3a0eaeea9c0690b6ede1b17b4fd8ce", ])

single_taxon <- prune_taxa("3f3a0eaeea9c0690b6ede1b17b4fd8ce", pseq)


df <- data.frame(
  count = as.numeric(otu_table(single_taxon)),
  total_reads = sample_sums(single_taxon),
  Type = sample_data(single_taxon)$Type
)

fit <- bbdml(
  formula = count / total_reads ~ Type,
  phi.formula = ~ 1,
  data = df,
  weights = df$total_reads
)
summary(fit)

result <- differentialTest(formula = ~ Type,
                           phi.formula = ~ Type,
                           formula_null = ~ 1,
                           phi.formula_null = ~ 1,
                           test = "Wald",
                           boot = FALSE,
                           data = single_taxon,
                           fdr_cutoff = 1)  # 


corncob <- bbdml(formula = "3f3a0eaeea9c0690b6ede1b17b4fd8ce" ~ 1,
                 phi.formula = ~ 1,
                 data = pseq)
