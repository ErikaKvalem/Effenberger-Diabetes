library(tidyverse)

library(phyloseq)

library(microbiomeMarker)


library(dplyr)
library(microbiomeMarker)

#library(knitr)


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



phy_sub <- subset_samples(phy, Type %in% c("DM", "PDM"))


table(sample_data(phy_sub)$group)

# Extract and clean taxonomy
tax <- as.data.frame(tax_table(phy_sub))

# Keep only the 7 standard taxonomic ranks
tax <- tax[, c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species_exact")]

# Rename 'Species_exact' to 'Species'
colnames(tax)[colnames(tax) == "Species_exact"] <- "Species"

# Reassign as tax_table and preserve rownames
tax_fixed <- tax_table(as.matrix(tax))
rownames(tax_fixed) <- taxa_names(phy_sub)

# Set it back into phyloseq object
tax_table(phy_sub) <- tax_fixed

# Confirm
rank_names(phy_sub)
lef_out <- run_lefse(
  phy_sub,
  group = "Type",
  norm = "CPM", 
  kw_cutoff = 0.05,
  lda_cutoff = 2
)

dat <- marker_table(lef_out) %>%
  data.frame() %>%
  dplyr::select(1:4)

head(dat)

#dat %>% table(align = "c")



plot_ef_bar(lef_out)

dat <- marker_table(lef_out) %>%
  data.frame() %>%
  dplyr::select(feature, enrich_group, ef_lda, pvalue)

dat <- dat %>%
  mutate(signed_lda = ifelse(enrich_group == "DM", -ef_lda, ef_lda))  # Flip for diverging barplot

# Order features for plotting
dat$feature <- factor(dat$feature, levels = dat$feature[order(dat$signed_lda)])


dat$feature_mod <- str_extract(dat$feature, "[^|]+$") 

dat <- dat %>%
  arrange(enrich_group, desc(signed_lda)) %>%
  mutate(feature_mod = factor(feature_mod, levels = rev(feature_mod)))  # reverse for coord_flip





p <- ggplot(dat, aes(x = feature_mod, y = signed_lda, fill = enrich_group)) +
  geom_bar(stat = "identity", width = 0.7) +
  coord_flip() +
  scale_y_continuous(name = "LDA SCORE (log 10)") +
  scale_fill_manual(values = c("DM" = "#F8766D", "PDM" = "#00BFC4")) +
  theme_minimal(base_size = 12) +
  theme(
    axis.title.y = element_blank(),
    axis.text.y = element_text(size = 20),
    legend.title = element_blank()
  )


p<- p + theme(legend.key = element_blank(), 
                         strip.background = element_rect(colour="black", fill="white"))

p
#plot_cladogram(lef_out, color = c("red","blue"), clade_label_level = 2)

#ggsave(plot=p,"/data/scratch/kvalem/projects/2024/Effenberger-Diabetes/02-scripts/figures/v02/barplot_lefse_dm_pdm.svg", height = 10, width = 10)
#ggsave(plot=p,"/data/scratch/kvalem/projects/2024/Effenberger-Diabetes/02-scripts/figures/v02/barplot_lefse_dm_pdm.png", height = 10, width = 10)
###

##################################

tax <- as.data.frame(tax_table(phy))

# Keep only the 7 standard taxonomic ranks
tax <- tax[, c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species_exact")]

# Rename 'Species_exact' to 'Species'
colnames(tax)[colnames(tax) == "Species_exact"] <- "Species"

# Reassign as tax_table and preserve rownames
tax_fixed <- tax_table(as.matrix(tax))
rownames(tax_fixed) <- taxa_names(phy)

# Set it back into phyloseq object
tax_table(phy) <- tax_fixed



# PDM vs DM
phy_dm_pdm <- subset_samples(phy, Type %in% c("DM", "PDM"))
lef_dm_pdm <- run_lefse(phy_dm_pdm, group = "Type", norm = "CPM", kw_cutoff = 0.05, taxa_rank = "Phylum",lda_cutoff = 2)

# DM vs K
phy_dm_k <- subset_samples(phy, Type %in% c("DM", "K"))
lef_dm_k <- run_lefse(phy_dm_k, group = "Type", norm = "CPM", kw_cutoff = 0.05,taxa_rank = "Phylum", lda_cutoff = 2)

# PDM vs K
phy_pdm_k <- subset_samples(phy, Type %in% c("PDM", "K"))
lef_pdm_k <- run_lefse(phy_pdm_k, group = "Type", norm = "CPM", kw_cutoff = 0.05,taxa_rank = "Phylum", lda_cutoff = 2)



get_lefse_df <- function(lef_out, comp_name, group1, group2) {
  df <- marker_table(lef_out) %>%
    data.frame()
  
  # Check if kw_pval exists and has non-NA values
  if ("pvalue" %in% colnames(df) && any(!is.na(df$pvalue))) {
    df$padj <- p.adjust(df$pvalue, method = "fdr")
  } else {
    df$padj <- NA_real_
  }
  
  df <- df %>%
    mutate(
      feature_mod = stringr::str_extract(feature, "[^|]+$"),
      signed_lda = ifelse(enrich_group == group1, -ef_lda, ef_lda),
      comparison = comp_name
    )
}


dat_dm_pdm <- get_lefse_df(lef_dm_pdm, "PDM vs DM", "DM", "PDM")
dat_dm_k    <- get_lefse_df(lef_dm_k, "DM vs K", "DM", "K")
dat_pdm_k   <- get_lefse_df(lef_pdm_k, "PDM vs K", "PDM", "K")

# Combine all
dat_all <- bind_rows(dat_dm_pdm, dat_dm_k, dat_pdm_k)

dat_all <- dat_all %>%
  group_by(comparison) %>%
  arrange(enrich_group, desc(signed_lda), .by_group = TRUE) %>%
  mutate(feature_mod = factor(feature_mod, levels = rev(unique(feature_mod))))

dat_all$enrich_group <- factor(dat_all$enrich_group, levels = c("DM", "PDM", "K"))
dat_all$enrich_group <- factor(dat_all$enrich_group, levels = c("DM", "PDM", "K"))

dat_all <- dat_all %>%
  group_by(comparison) %>%
  arrange(enrich_group, desc(signed_lda), .by_group = TRUE) %>%
  mutate(feature_mod = factor(feature_mod, levels = rev(unique(feature_mod)))) %>%
  ungroup()

dat_all$comparison <- factor(
  dat_all$comparison,
  levels = c("PDM vs DM", "DM vs K", "PDM vs K")
)

dat_all <- dat_all %>%
  mutate(feature_id = paste(comparison, feature_mod, sep = " | "))

dat_all <- dat_all %>%
  group_by(comparison) %>%
  arrange(enrich_group, desc(signed_lda), .by_group = TRUE) %>%
  mutate(feature_id = factor(feature_id, levels = rev(unique(feature_id)))) %>%
  ungroup()

legend_df <- data.frame(
  feature_id = NA,
  signed_lda = NA,
  enrich_group = NA,
  padj_label = "p.adj < 0.05"
)

# Plot
p_all <- ggplot(dat_all, aes(x = feature_id, y = signed_lda, fill = enrich_group)) +
  geom_bar(stat = "identity", width = 0.7) +
  coord_flip() +
  scale_y_continuous(name = "LDA SCORE (log10)                  p.adj < 0.05") +
  scale_fill_manual(values = c("DM" = "#E1812C", "PDM" = "#3A923A", "K" = "#3274A1")) +
  scale_x_discrete(labels = dat_all$feature_mod) +   # 👈 custom axis labels
facet_wrap(~ comparison, scales = "free_y", ncol = 1) +

  theme_minimal(base_size = 12) +
  theme(
    axis.title.y = element_blank(),
    axis.text.y = element_text(size = 12),
    strip.text = element_text(size = 14, face = "bold"),
    legend.title = element_blank()
  )

p_all

#p_all<- p_all + theme(legend.key = element_blank(), 
#              strip.background = element_rect(colour="black", fill="white"))

p_all
write.csv(dat_all, file = "/data/scratch/kvalem/projects/2024/Effenberger-Diabetes/02-scripts/tables/results/lefse_data_all_pvalue_padj.csv", row.names = FALSE)

ggsave(plot=p_all,"/data/scratch/kvalem/projects/2024/Effenberger-Diabetes/02-scripts/figures/v02/barplot_lefse_pall_genus_padj.svg", height = 5, width = 5)
ggsave(plot=p_all,"/data/scratch/kvalem/projects/2024/Effenberger-Diabetes/02-scripts/figures/v02/barplot_lefse_pall_genus_padj.png", height = 5, width = 5)
###

######################### volcano plot 

# Ensure padj is non-zero to avoid -Inf in log scale
dat_all <- dat_all %>%
  mutate(
    neglog10_padj = -log10(padj),
    enrich_group = factor(enrich_group, levels = c("DM", "PDM", "K"))
  )

# Volcano plot
p_all <- ggplot(dat_all, aes(x = signed_lda, y = neglog10_padj, color = enrich_group)) +
  geom_point(size = 3, alpha = 0.8) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "gray50") +
  facet_wrap(~ comparison, scales = "free", ncol = 1) +
  scale_color_manual(values = c("DM" = "#E1812C", "PDM" = "#3A923A", "K" = "#3274A1")) +
  labs(
    x = "LDA Score (signed)",
    y = expression(-log[10](p.adj)),
    color = "Enriched in"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    strip.text = element_text(size = 14, face = "bold"),
    legend.title = element_blank()
  )

p_all
#write.csv(dat_all, file = "/data/scratch/kvalem/projects/2024/Effenberger-Diabetes/02-scripts/tables/v02/dat_all_results.csv", row.names = FALSE)