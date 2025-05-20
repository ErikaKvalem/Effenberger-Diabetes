library(tidyverse)

library(phyloseq)

library(microbiomeMarker)

library(knitr)

library(dplyr)
library(microbiomeMarker)

library(knitr)


library(dplyr)
library(tidyr)
library(stringr)


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
lef_dm_pdm <- run_lefse(phy_dm_pdm, group = "Type", norm = "CPM", kw_cutoff = 0.05, taxa_rank = "Genus",lda_cutoff = 2)

# DM vs K
phy_dm_k <- subset_samples(phy, Type %in% c("DM", "K"))
lef_dm_k <- run_lefse(phy_dm_k, group = "Type", norm = "CPM", kw_cutoff = 0.05,taxa_rank = "Genus", lda_cutoff = 2)

# PDM vs K
phy_pdm_k <- subset_samples(phy, Type %in% c("PDM", "K"))
lef_pdm_k <- run_lefse(phy_pdm_k, group = "Type", norm = "CPM", kw_cutoff = 0.05,taxa_rank = "Genus", lda_cutoff = 2)

get_lefse_df <- function(lef_out, comp_name, group1, group2) {
  df <- marker_table(lef_out) %>%
    data.frame() %>%
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



df <- read_csv("/data/projects/2024/Effenberger-Diabetes/data/PDM merged 3.0_modified.csv")%>%
  clean_names()


name_map_full <- c(
  "probennummer" = "sample_information",
  "verstorben" = "Deceased",
  "pankreatektomie" = "Pancreatectomy",
  "c2" = "C-Peptide (FU)",
  "nikotin" = "Smoking Status",
  "sex" = "Sex",
  "age" = "Age",
  
  # Baseline (BS)
  "hyperlipid_mie_20181" = "Hyperlipidemia (BS)",
  "arterielle_hyperotnie1" = "Arterial Hypertension (BS)",
  "bmi1" = "Body Mass Index (BMI) (BS)",
  "gr_e1" = "Height (BS)",
  "gewicht1" = "Weight (BS)",
  "lipidsenker1" = "Lipid-lowering Medication (BS)",
  "art_von_lipidsenker1" = "Type of Lipid-lowering Medication (BS)",
  "rr_medikation1" = "Blood Pressure Medication (BS)",
  "b_blocker1" = "Beta Blockers (BS)",
  "ace_hemmer1" = "ACE Inhibitors (BS)",
  "diuretika1" = "Diuretics (BS)",
  "insulin1" = "Insulin Therapy (BS)",
  "langzeit_insulin1" = "Long-acting Insulin (BS)",
  "kurzzeit_insulin1" = "Short-acting Insulin (BS)",
  "misch_inuslin1" = "Mixed Insulin (BS)",
  "masld1" = "MASLD (BS)",
  "khk1" = "Coronary Heart Disease (CHD) (BS)",
  "ca1" = "Cancer (unspecified) (BS)",
  
  # Follow-up (FU)
  "hyperlipid_mie_20182" = "Hyperlipidemia (FU)",
  "arterielle_hyperotnie2" = "Arterial Hypertension (FU)",
  "bmi2" = "Body Mass Index (BMI) (FU)",
  "gr_e2" = "Height (FU)",
  "gewicht2" = "Weight (FU)",
  "lipidsenker2" = "Lipid-lowering Medication (FU)",
  "art_von_lipidsenker2" = "Type of Lipid-lowering Medication (FU)",
  "rr_medikation2" = "Blood Pressure Medication (FU)",
  "b_blocker2" = "Beta Blockers (FU)",
  "ace_hemmer2" = "ACE Inhibitors (FU)",
  "diuretika2" = "Diuretics (FU)",
  "insulin2" = "Insulin Therapy (FU)",
  "langzeit_insulin2" = "Long-acting Insulin (FU)",
  "kurzzeit_insulin2" = "Short-acting Insulin (FU)",
  "misch_inuslin2" = "Mixed Insulin (FU)",
  "masld2" = "MASLD (FU)",
  "khk2" = "Coronary Heart Disease (CHD) (FU)",
  "ca2" = "Cancer (unspecified) (FU)",
  
  # Lab values (BS)
  "leukozyten1" = "Leukocytes (BS)",
  "h_moglobin1" = "Hemoglobin (BS)",
  "h_matokrit1" = "Hematocrit (BS)",
  "thrombozyten1" = "Platelets (BS)",
  "harnstoff1" = "Urea (BS)",
  "creatinin_enzym_idms_1" = "Creatinine (IDMS) (BS)",
  "glomerul_re_filtrationsrate1" = "Glomerular Filtration Rate (GFR) (BS)",
  "bilirubin_gesamt1" = "Total Bilirubin (BS)",
  "natrium1" = "Sodium (BS)",
  "got_asat_1" = "AST (BS)",
  "gpt_alat_1" = "ALT (BS)",
  "gamma_gt1" = "GGT (BS)",
  "alkalische_phosphatase1" = "Alkaline Phosphatase (BS)",
  "lactat_dehydrogenase_ldh_1" = "LDH (BS)",
  "c_reaktives_prot_crp_1" = "CRP (BS)",
  "quicktest_pt_1" = "Prothrombin Time (BS)",
  "inr_pt_1" = "INR (BS)",
  "part_thrombopl_zeit_a_ptt_1" = "aPTT (BS)",
  "fibrinogen_funkt_n_clauss1" = "Fibrinogen (Clauss) (BS)",
  "albumin1" = "Albumin (BS)",
  "glukose1" = "Glucose (BS)",
  "hb_a1c_dcct_ngsp_1" = "HbA1c (DCCT/NGSP) (BS)",
  "hb_a1c_ifcc_1" = "HbA1c (IFCC) (BS)",
  "cholesterin1" = "Total Cholesterol (BS)",
  "non_hdl_cholesterin1" = "Non-HDL Cholesterol (BS)",
  "triglyceride1" = "Triglycerides (BS)",
  "hdl_cholesterin1" = "HDL Cholesterol (BS)",
  "ldl_cholesterin1" = "LDL Cholesterol (BS)",
  "eisen1" = "Iron (BS)",
  "ferritin1" = "Ferritin (BS)",
  "transferrin1" = "Transferrin (BS)",
  "transferrins_ttigung1" = "Transferrin Saturation (BS)",
  "troponin_t_hoch_sens_1" = "Troponin T (High Sensitivity) (BS)",
  "nt_pro_bnp1" = "NT-proBNP (BS)",
  
  # Lab values (FU)
  "leukozyten2" = "Leukocytes (FU)",
  "h_moglobin2" = "Hemoglobin (FU)",
  "h_matokrit2" = "Hematocrit (FU)",
  "thrombozyten2" = "Platelets (FU)",
  "harnstoff2" = "Urea (FU)",
  "creatinin_enzym_idms_2" = "Creatinine (IDMS) (FU)",
  "glomerul_re_filtrationsrate2" = "Glomerular Filtration Rate (GFR) (FU)",
  "bilirubin_gesamt2" = "Total Bilirubin (FU)",
  "natrium2" = "Sodium (FU)",
  "got_asat_2" = "AST (FU)",
  "gpt_alat_2" = "ALT (FU)",
  "gamma_gt2" = "GGT (FU)",
  "alkalische_phosphatase2" = "Alkaline Phosphatase (FU)",
  "lactat_dehydrogenase_ldh_2" = "LDH (FU)",
  "c_reaktives_prot_crp_2" = "CRP (FU)",
  "quicktest_pt_2" = "Prothrombin Time (FU)",
  "inr_pt_2" = "INR (FU)",
  "part_thrombopl_zeit_a_ptt_2" = "aPTT (FU)",
  "fibrinogen_funkt_n_clauss2" = "Fibrinogen (Clauss) (FU)",
  "albumin2" = "Albumin (FU)",
  "glukose2" = "Glucose (FU)",
  "hb_a1c_dcct_ngsp_2" = "HbA1c (DCCT/NGSP) (FU)",
  "hb_a1c_ifcc_2" = "HbA1c (IFCC) (FU)",
  "cholesterin2" = "Total Cholesterol (FU)",
  "non_hdl_cholesterin2" = "Non-HDL Cholesterol (FU)",
  "triglyceride2" = "Triglycerides (FU)",
  "hdl_cholesterin2" = "HDL Cholesterol (FU)",
  "ldl_cholesterin2" = "LDL Cholesterol (FU)",
  "eisen2" = "Iron (FU)",
  "ferritin2" = "Ferritin (FU)",
  "transferrin2" = "Transferrin (FU)",
  "transferrins_ttigung2" = "Transferrin Saturation (FU)",
  "troponin_t_hoch_sens_2" = "Troponin T (High Sensitivity) (FU)",
  "nt_pro_bnp2" = "NT-proBNP (FU)"
)


df <- df %>%
  dplyr::rename_with(.cols = all_of(names(name_map_full)), .fn = ~ name_map_full[.x])

sample_data(phy) <- sample_data(phy)[, "sample_information", drop = FALSE]

# 1. Extract sample data from phy as data.frame
sam <- as(sample_data(phy), "data.frame")

# 2. Join with df based on sample_information
sam_new <- dplyr::left_join(sam, df, by = "sample_information")

# 3. Set rownames to match phyloseq sample names
rownames(sam_new) <- rownames(sam)  # Ensures sample_names match exactly

# 4. Convert and assign back to phyloseq
sample_data(phy) <- phyloseq::sample_data(sam_new)


############################################## 



# Add taxonomy columns to dat_all without losing other columns
dat_all <- dat_all %>%
  mutate(
    feature_clean = str_replace_all(feature, "s__$", "s__")  # optional cleanup
  ) %>%
  separate(
    col = feature_clean,
    into = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"),
    sep = "\\|",
    fill = "right",
    remove = FALSE
  ) %>%
  mutate(across(Kingdom:Species, ~str_remove(., "^[a-z]__")))  # strip "k__", "p__", etc.

# Sample metadata (clinical features)
clinical_df <- data.frame(sample_data(phy), check.names = FALSE)

# Microbial abundance table
abund <- as.data.frame(otu_table(phy))
if (taxa_are_rows(phy)) {
  abund <- t(abund)
}
all(rownames(clinical_df) == rownames(abund))  # should be TRUE

shared_samples <- intersect(rownames(clinical_df), rownames(abund))
clinical_df <- clinical_df[shared_samples, , drop = FALSE]
abund <- abund[shared_samples, , drop = FALSE]

# Extract taxonomy
tax <- as.data.frame(tax_table(phy))
tax$tax_id <- rownames(tax)

# Build full feature label like "g__Subdoligranulum", "f__Ruminococcaceae", etc.
tax_long <- tax %>%
  mutate(across(everything(), as.character)) %>%
  pivot_longer(cols = Genus, names_to = "rank", values_to = "name") %>%
  filter(!is.na(name), name != "") %>%
  mutate(feature_mod = paste0(tolower(substr(rank, 1, 1)), "__", name)) %>%
  dplyr::select(tax_id, feature_mod)



abund <- as.data.frame(otu_table(phy))
if (taxa_are_rows(phy)) {
  abund <- t(abund)
}
abund$tax_id <- colnames(otu_table(phy))
abund_tidy <- as.data.frame(t(otu_table(phy)))
abund_tidy$sample_id <- rownames(abund_tidy)


abund_long <- abund_tidy %>%
  pivot_longer(-sample_id, names_to = "tax_id", values_to = "abundance")

abund_long <- abund_long %>%
  left_join(tax_long, by = "tax_id") %>%
  filter(!is.na(feature_mod))

abund_grouped <- abund_long %>%
  group_by(sample_id, feature_mod) %>%
  summarise(abundance = sum(abundance), .groups = "drop")

abund_wide <- abund_grouped %>%
  pivot_wider(names_from = feature_mod, values_from = abundance)

abund_matrix <- as.data.frame(abund_wide)
rownames(abund_matrix) <- abund_matrix$sample_id
abund_matrix$sample_id <- NULL

# Get list of 40 features to correlate
features_to_use <- unique(dat_all$feature_mod)
features_to_use <- paste0("g__", features_to_use)


# Make sure they exist in the abundance matrix
features_to_use <- intersect(features_to_use, colnames(abund_matrix))

# Subset abundance data
abund_sel <- abund_matrix[, features_to_use, drop = FALSE]

clinical_df <- data.frame(sample_data(phy), check.names = FALSE)



clinical_vars <- c(
  "NT-proBNP (BS)", "NT-proBNP (FU)",
  "Ferritin (BS)", "Ferritin (FU)",
  "Platelets (BS)", "Platelets (FU)",
  "LDH (BS)", "LDH (FU)",
  "LDL Cholesterol (BS)", "LDL Cholesterol (FU)","Alkaline Phosphatase (BS)","Alkaline Phosphatase (FU)","AST (BS)", "AST (FU)",
  "ALT (BS)", "ALT (FU)", "Short-acting Insulin (BS)", "Short-acting Insulin (FU)", "Troponin T (High Sensitivity) (BS)","Troponin T (High Sensitivity) (FU)", 
  "Urea (BS)", "Urea (FU)", "Weight (BS)", "Weight (FU)", "Iron (BS)", "Iron (FU)", "Body Mass Index (BMI) (BS)", "Body Mass Index (BMI) (FU)", "Lipid-lowering Medication (BS)", "Lipid-lowering Medication (FU)","Hyperlipidemia (BS)", "Hyperlipidemia (FU)", 
  "Creatinine (IDMS) (BS)", "Creatinine (IDMS) (FU)", "Prothrombin Time (BS)", "Prothrombin Time (FU)", "Non-HDL Cholesterol (BS)", "Non-HDL Cholesterol (FU)"
)

bs_vars <- clinical_vars[grepl("\\(BS\\)$", clinical_vars)]
fu_vars <- gsub("\\(BS\\)$", "(FU)", bs_vars)
valid_pairs <- bs_vars[fu_vars %in% colnames(clinical_df)]
fu_vars <- gsub("\\(BS\\)$", "(FU)", valid_pairs)

# Compute deltas safely
for (i in seq_along(valid_pairs)) {
  delta_name <- gsub(" \\(BS\\)$", " (Delta)", valid_pairs[i])
  clinical_df[[delta_name]] <- as.numeric(clinical_df[[fu_vars[i]]]) - as.numeric(clinical_df[[valid_pairs[i]]])
}

delta_vars <- grep(" \\(Delta\\)$", colnames(clinical_df), value = TRUE)

shared_samples <- intersect(rownames(clinical_df), rownames(abund_sel))
clinical_df <- clinical_df[shared_samples, ]
abund_sel <- abund_sel[shared_samples, ]



# Filter to available ones
clinical_vars <- clinical_vars[clinical_vars %in% colnames(clinical_df)]

# Compute correlations
cor_results <- expand.grid(
  microbe = colnames(abund_sel),
  clinical_var = delta_vars,
  stringsAsFactors = FALSE
) %>%
  mutate(
    cor = purrr::map2_dbl(microbe, clinical_var, ~{
      x <- abund_sel[[.x]]
      y <- clinical_df[[.y]]
      if (length(unique(na.omit(x))) < 2 || length(unique(na.omit(y))) < 2) return(NA)
      cor.test(x, y, method = "pearson")$estimate
    }),
    p_value = purrr::map2_dbl(microbe, clinical_var, ~{
      x <- abund_sel[[.x]]
      y <- clinical_df[[.y]]
      if (length(unique(na.omit(x))) < 2 || length(unique(na.omit(y))) < 2) return(NA)
      cor.test(x, y, method = "pearson")$p.value
    }),
    p_adj = p.adjust(p_value, method = "fdr")
  )

cor_results <- cor_results %>%
  mutate(signif = case_when(
    p_value < 0.001 ~ "***",
    p_value < 0.01 ~ "**",
    p_value < 0.05 ~ "*",
    TRUE ~ ""
  ))

ggplot(cor_results, aes(x = clinical_var, y = microbe, fill = cor)) +
  geom_tile(color = "white") +
  geom_text(aes(label = signif), size = 4, color = "black") +
  scale_fill_gradient2(
    low = "blue", mid = "white", high = "red", midpoint = 0,
    name = "Spearman R"
  ) +
  labs(
    x = "",
    y = "",
    title = "Correlation Heatmap with Significance"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1),
    panel.grid = element_blank()
  )

##################################

clinical_df$Type <- ifelse(grepl("PDM", clinical_df$sample_information), "PDM",
                                ifelse(grepl("K", clinical_df$sample_information), "K", "DM"))

library(dplyr)
library(tidyr)
library(purrr)

types <- unique(clinical_df$Type)

cor_results_all <- purrr::map_dfr(types, function(type_level) {
  # Subset to group
  clinical_sub <- clinical_df[clinical_df$Type == type_level, ]
  abund_sub <- abund_sel[rownames(clinical_sub), ]
  
  # Compute correlations
  df <- expand.grid(
    microbe = colnames(abund_sub),
    clinical_var = delta_vars,
    stringsAsFactors = FALSE
  ) %>%
    mutate(
      cor = map2_dbl(microbe, clinical_var, ~{
        x <- abund_sub[[.x]]
        y <- clinical_sub[[.y]]
        if (length(unique(na.omit(x))) < 2 || length(unique(na.omit(y))) < 2) return(NA)
        cor.test(x, y, method = "spearman")$estimate
      }),
      p_value = map2_dbl(microbe, clinical_var, ~{
        x <- abund_sub[[.x]]
        y <- clinical_sub[[.y]]
        if (length(unique(na.omit(x))) < 2 || length(unique(na.omit(y))) < 2) return(NA)
        cor.test(x, y, method = "pearson")$p.value
      }),
      p_adj = p.adjust(p_value, method = "fdr"),
      Type = type_level
    )
})

cor_results_all <- cor_results_all %>%
  mutate(signif = case_when(
    p_value < 0.001 ~ "***",
    p_value < 0.01 ~ "**",
    p_value < 0.05 ~ "*",
    TRUE ~ ""
  ))

cor_results_all$Type <- factor(cor_results_all$Type, levels = c("K", "DM", "PDM"))

cor_results_all <- cor_results_all %>%
  mutate(clinical_var_clean = gsub(" \\(Delta\\)$", "", clinical_var))

#cor_results_all <- cor_results_all %>%
#  filter(Type != "K")

ggplot(cor_results_all, aes(x = clinical_var_clean, y = microbe, fill = cor)) +
  geom_tile(color = "white") +
  geom_text(aes(label = signif), size = 3) +
  facet_wrap(~Type,ncol = 1) +
  scale_fill_gradient2(
    low = "blue", high = "red", mid = "white", midpoint = 0,
    name = "Spearman R"
  ) +
  labs(
    title = "Correlation Heatmaps by Type",
    x = "Clinical Variable (Delta)",
    y = "Microbial Feature"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid = element_blank()
  )


p_all <- ggplot(cor_results_all, aes(x = clinical_var_clean, y = microbe)) +
  geom_point(aes(color = cor, size = -log10(p_value)), alpha = 0.8) +
  facet_wrap(~Type, ncol = 1) +
  scale_color_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +
  labs(
    x = "Clinical Variable (Δ)",
    y = "Microbial Feature",
    title = "Spearman Correlation Dot Plot by Type",
    size = "-log10(p)",
    color = "Spearman R"
  ) +
  theme_minimal(base_size = 11) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


cor_top <- cor_results_all %>%
  filter(!is.na(cor)) %>%
  arrange(desc(abs(cor))) %>%
  group_by(Type) %>%
  slice_head(n = 5) %>%
  ungroup()

ggplot(cor_top, aes(x = reorder_within(microbe, abs(cor), Type), y = cor, fill = clinical_var_clean)) +
  geom_col(position = position_dodge()) +
  facet_wrap(~Type, scales = "free_y", ncol = 1) +
  coord_flip() +
  scale_x_reordered() +
  labs(title = "Top 10 Microbial Correlations by Type",
       x = "Microbial Feature", y = "Spearman R") +
  theme_minimal()


library(ggplot2)

cor_top10 <- cor_results_all %>%
  filter(!is.na(cor)) %>%
  mutate(cor_abs = abs(cor)) %>%
  arrange(desc(cor_abs)) %>%
  slice_head(n = 500)

ggplot(cor_results_all, aes(x = reorder(clinical_var_clean, cor), y = reorder(microbe, cor))) +
  geom_point(aes(size = -log10(p_value), color = cor)) +
  facet_wrap(~Type, ncol = 1) +
  scale_color_gradient2(
    low = "blue", high = "red", mid = "white", midpoint = 0,
    name = "Spearman R"
  ) +
  scale_size_continuous(name = "-log10(p)") +
  labs(
    x = "Clinical Variable (Δ)",
    y = "Microbial Feature",
    title = "Top 10 Strongest Correlations by Type"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid = element_blank()
  )


#ggsave(plot=p_all,"/data/scratch/kvalem/projects/2024/Effenberger-Diabetes/02-scripts/figures/v02/dotplot_correlation_micro_clinical.svg", height = 5, width = 5)
#ggsave(plot=p_all,"/data/scratch/kvalem/projects/2024/Effenberger-Diabetes/02-scripts/figures/v02/dotplot_correlation_micro_clinical.png", height = 5, width = 5)

############################################# linear model 

sample_data(phy)$Type <- ifelse(grepl("PDM", sample_data(phy)$sample_information), "PDM",
                                ifelse(grepl("K", sample_data(phy)$sample_information), "K", "DM"))


# Ensure the response is a binary factor
meta_sub$Type <- factor(meta_sub$Type, levels = c("DM", "PDM"))
#
data_model <- cbind(abund_sub, Type = meta_sub$Type)

model <- glm(Type ~ ., data = data_model, family = binomial)
summary(model)

# Predict probabilities
pred_probs <- predict(model, type = "response")

# Convert to class predictions
pred_class <- ifelse(pred_probs > 0.5, "PDM", "DM")

# Confusion matrix
table(True = meta_sub$Type, Predicted = pred_class)

library(pROC)
roc_obj <- roc(meta_sub$Type, pred_probs)
plot(roc_obj, main = paste("AUC:", round(auc(roc_obj), 2)))

meta_sub$Type <- factor(meta_sub$Type, levels = c("DM", "PDM"))  # PDM = 1, DM = 0
data_model <- cbind(abund_sub, Type = meta_sub$Type)

# Fit full logistic regression
model <- glm(Type ~ ., data = data_model, family = binomial)
summary(model)

library(broom)

model_summary <- tidy(model) %>%
  filter(term != "(Intercept)") %>%
  mutate(odds_ratio = exp(estimate)) %>%
  arrange(p.value)

# View top 10 features
top_features <- model_summary$term[1:10]
print(top_features)

top_features <- gsub("`", "", top_features)  # remove backticks
top_features <- trimws(top_features)         # remove whitespace

top_features <- intersect(top_features, colnames(abund_sub))  # keep only matching ones
abund_top <- abund_sub[, top_features, drop = FALSE]




abund_top <- abund_sub[, top_features, drop = FALSE]
data_model_top <- cbind(abund_top, Type = meta_sub$Type)

model_top <- glm(Type ~ ., data = data_model_top, family = binomial)
summary(model_top)

tidy(model_top) %>%
  mutate(odds_ratio = exp(estimate)) %>%
  dplyr::select(term, estimate, odds_ratio, std.error, p.value)

glimpse(tidy(model_top))

library(glmnet)

x <- as.matrix(abund_sub)
y <- as.factor(meta_sub$Type)

cv_fit <- cv.glmnet(x, y, family = "binomial", alpha = 1)  # LASSO
plot(cv_fit)

best_lambda <- cv_fit$lambda.min
lasso_model <- glmnet(x, y, family = "binomial", alpha = 1, lambda = best_lambda)

# Extract non-zero features
lasso_coefs <- coef(lasso_model)
selected_features <- rownames(lasso_coefs)[which(lasso_coefs != 0)]
print(selected_features)

selected_features <- rownames(lasso_coefs)[which(lasso_coefs != 0)]
print(selected_features)

plot(cv_fit)



library(ggplot2)
library(broom)

model_plot_data <- tidy(model_top) %>%
  filter(term != "(Intercept)") %>%
  mutate(
    term = gsub("`", "", term),
    odds_ratio = exp(estimate)
  )

ggplot(model_plot_data, aes(x = reorder(term, estimate), y = estimate)) +
  geom_col(fill = "steelblue") +
  geom_errorbar(aes(ymin = estimate - std.error, ymax = estimate + std.error), width = 0.2) +
  coord_flip() +
  labs(
    title = "Logistic Regression Coefficients",
    x = "Microbial Feature",
    y = "Log-Odds (Coefficient)"
  ) +
  theme_minimal()

library(pROC)
library(PRROC)

# ROC curve (already computed)
roc_obj <- roc(meta_sub$Type, pred_probs)
plot(roc_obj, main = paste("AUC:", round(auc(roc_obj), 3)))

# Precision-Recall
fg <- pred_probs[meta_sub$Type == "PDM"]
bg <- pred_probs[meta_sub$Type == "DM"]
pr <- pr.curve(scores.class0 = fg, scores.class1 = bg, curve = TRUE)
plot(pr, main = "Precision-Recall Curve")


# Pearson correlation matrix of top features
cor_mat <- cor(abund_sub[, top_features])
heatmap(cor_mat, main = "Feature Correlation Heatmap")

# VIF - needs all numeric and no perfect multicollinearity
library(car)
vif_model <- glm(Type ~ ., data = data_model_top, family = binomial)
vif(vif_model)

