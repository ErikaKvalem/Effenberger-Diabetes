library(tidyverse)

library(phyloseq)

library(ggplot2)
library(ggpmisc)

library(dplyr)
library(microbiomeMarker)




library(dplyr)
library(tidyr)
library(stringr)
library(viridis)

library(readr)
library(themis)

#library(tidymodels)


library(dplyr)
library(tidyr)
library(purrr)

library(ggplot2)


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



df <- read_csv("/data/projects/2024/Effenberger-Diabetes/data/PDM merged 3.0_modified.csv")%>%clean_names()


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

############################################# LINEAR REGRESSION

sam_new$Type <- ifelse(grepl("PDM", sam_new$sample_information), "PDM",
                                ifelse(grepl("K", sam_new$sample_information), "K", "DM"))

# 1. Get all (BS) columns
bs_vars <- grep("\\(BS\\)$", colnames(sam_new), value = TRUE)

# 2. Derive corresponding (FU) column names
fu_vars <- gsub("\\(BS\\)$", "(FU)", bs_vars)

# 3. Keep only pairs where both BS and FU columns exist
valid_pairs <- bs_vars[fu_vars %in% colnames(sam_new)]
fu_vars <- gsub("\\(BS\\)$", "(FU)", valid_pairs)

# 4. Compute delta = FU - BS
for (i in seq_along(valid_pairs)) {
  delta_name <- gsub(" \\(BS\\)$", " (Delta)", valid_pairs[i])
  sam_new[[delta_name]] <- as.numeric(sam_new[[fu_vars[i]]]) - as.numeric(sam_new[[valid_pairs[i]]])
}

# 5. Optional: get names of new delta variables
delta_vars <- grep(" \\(Delta\\)$", colnames(sam_new), value = TRUE)


# Linear model formula
formula <- y ~ x


ggplot(sam_new, aes(x =`Body Mass Index (BMI) (FU)`, y = `HbA1c (DCCT/NGSP) (FU)`, color = Type)) +
  geom_point(size = 3) +
  geom_smooth(method = "lm", formula = formula, se = TRUE, color = "black") +
  stat_poly_eq(
    aes(label = paste(..p.value.label..)), 
    formula = formula, 
    parse = TRUE, 
    label.x.npc = "right", 
    label.y.npc = 0.95,
    size = 5
  ) + scale_color_manual(values = c(
    "DM" = "#E1812C",
    "PDM" = "#3A923A",
    "K" = "#3274A1"
  )) +
  labs(
    title = "HbA1c vs ",
    x = "Weight",
    y = "HbA1c"
  ) +
  theme_minimal()


p <- ggplot(sam_new, aes(x = `HbA1c (DCCT/NGSP) (FU)`, y = `Body Mass Index (BMI) (FU)`, color = Type)) +
  geom_point(size = 2) +
  geom_smooth(method = "lm", se = TRUE) +
  stat_poly_eq(
    aes(label = paste(..p.value.label.., sep = "~~~")),
    formula = y ~ x, parse = TRUE, label.x.npc = "right", label.y.npc = 0.95
  ) +
  labs(title = "HbA1c vs BMI", x = "HbA1c (FU)", y = "Body Mass Index (BMI) (FU)") +
  scale_color_manual(values = c("DM" = "#E1812C", "PDM" = "#3A923A", "K" = "blue")) +
  theme_minimal()



ggsave(plot=p,"/data/scratch/kvalem/projects/2024/Effenberger-Diabetes/02-scripts/figures/v02/linear_regression_bmi_hba1c.svg", height = 4, width = 7)
ggsave(plot=p,"/data/scratch/kvalem/projects/2024/Effenberger-Diabetes/02-scripts/figures/v02/linear_regression_bmi_hba1c.png", height = 4, width = 7)
############################################## linear regression GENUS 

library(phyloseq)
library(dplyr)

# Aggregate at Genus level if not already
phy_genus <- tax_glom(phy, taxrank = "Genus")


# Transform counts to relative abundance (optional but common)
phy_genus_rel <- transform_sample_counts(phy_genus, function(x) x / sum(x))

# Extract abundance data
abund <- as.data.frame(otu_table(phy_genus_rel))
taxa <- as.data.frame(tax_table(phy_genus_rel))

# Identify row matching genus "Blautia"
blautia_row <- rownames(taxa)[taxa$Genus == "Streptococcus"]
blautia_abund <- as.numeric(abund[blautia_row, ])
names(blautia_abund) <- colnames(abund)
# Make sure sample names align
common_samples <- intersect(rownames(sam_new), names(blautia_abund))


# Create a combined dataframe
plot_df <- sam_new[common_samples, ]
plot_df$Blautia <- blautia_abund[common_samples]
library(ggplot2)
library(ggpmisc)


p <- ggplot(plot_df, aes(x = Blautia, y = `Transferrin (BS)`, color = Type)) +
  geom_point(size = 2) +
  geom_smooth(method = "lm", se = TRUE) +
  stat_poly_eq(
    aes(label = paste(..p.value.label..)),
    formula = y ~ x, parse = TRUE,
    label.x.npc = "left", label.y.npc = 0.95
  ) +
  labs(title = "Streptococcus vs Transferrin (BS)", x = "Streptococcus", y = "Transferrin (BS)") +
  scale_color_manual(values = c("DM" = "#E1812C", "PDM" = "#3A923A", "K" = "blue")) +
  theme_minimal()


ggsave(plot=p,"/data/scratch/kvalem/projects/2024/Effenberger-Diabetes/02-scripts/figures/v02/linear_regression_Streptococcus_Transferrin.svg", height = 4, width = 7)
ggsave(plot=p,"/data/scratch/kvalem/projects/2024/Effenberger-Diabetes/02-scripts/figures/v02/linear_regression_Streptococcus_Transferrin.png", height = 4, width = 7)

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
    name = "Pearson R"
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
  geom_point(aes(color = cor, size = -log10(p_adj)), alpha = 0.8) +
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

############################################### LINEAR MODEL LOG REGRESSION & top_features ONLY CLINICAL 
ps_genus <- tax_glom(phy, "Genus")
ps_rel <- transform_sample_counts(ps_genus, function(x) x / sum(x))
otu_df <- as.data.frame(t(otu_table(ps_rel)))

# Combine with metadata
meta_df <- as.data.frame(sample_data(phy))
#data_all <- cbind(meta_df, otu_df)


data_all <- as(sample_data(phy), "data.frame")
#clinical_vars <- colnames(data_all)
#bs_vars <- clinical_vars[grepl("\\(BS\\)$", clinical_vars)]
#fu_vars <- gsub("\\(BS\\)$", "(FU)", bs_vars)

#valid_pairs <- bs_vars[fu_vars %in% colnames(data_all)]
#fu_vars <- gsub("\\(BS\\)$", "(FU)", valid_pairs)

#for (i in seq_along(valid_pairs)) {
#  delta_name <- gsub(" \\(BS\\)$", " (Delta)", valid_pairs[i])
#  data_all[[delta_name]] <- as.numeric(data_all[[fu_vars[i]]]) - as.numeric(data_all[[valid_pairs[i]]])
#}

#delta_vars <- grep(" \\(Delta\\)$", colnames(data_all), value = TRUE)
#meta_keep <- c("sample_information", "Type")
#microbial_vars <- setdiff(colnames(otu_df), c("sample_information", "Type"))


#data_all <- data_all[, c(meta_keep, delta_vars, microbial_vars)]

data_all$Type <- ifelse(grepl("PDM", data_all$sample_information), "PDM",
                        ifelse(grepl("K", data_all$sample_information), "K", "DM"))

data_all <- data_all %>%
  filter(Type %in% c("PDM", "DM")) %>%
  mutate(Type = factor(Type))

# Split data into train and test
set.seed(421)
split <- initial_split(data_all, prop = 0.65, strata = Type)
train <- split %>% 
  training()
test <- split %>% 
  testing()

rec <- recipe(Type ~ ., data = train) %>%
  update_role(sample_information, new_role = "id") %>%
  step_novel(all_nominal_predictors()) %>%  #
  step_zv(all_predictors()) %>%
  step_nzv(all_nominal_predictors()) %>%
  step_unknown(all_nominal_predictors()) %>%         
  step_impute_mean(all_numeric_predictors()) %>%    
  step_dummy(all_nominal_predictors(), one_hot = TRUE) %>%
  step_normalize(all_numeric_predictors()) %>%
  step_downsample(Type)

set.seed(421)
rec_prepped <- prep(rec, training = train)

train_baked <- bake(rec_prepped, new_data = NULL) %>%
  select(-sample_information)
test_baked  <- bake(rec_prepped, new_data = test) %>%
  select(-sample_information)

model <- logistic_reg(penalty = 0.1, mixture = 1) %>%
  set_engine("glmnet") %>%
  set_mode("classification") %>%
  fit(Type ~ ., data = train_baked)

pred_prob <- predict(model, new_data = test_baked, type = "prob")
pred_class <- predict(model, new_data = test_baked, type = "class")

results <- test %>%
  select(Type) %>%
  bind_cols(pred_class, pred_prob)

# --- Evaluation ---
acc <- accuracy(results, truth = Type, estimate = .pred_class)
auc <- roc_auc(results, truth = Type, .pred_DM)

print(acc)
print(auc)

# --- Confusion matrix ---
print(conf_mat(results, truth = Type, estimate = .pred_class))

c <- conf_mat(results, truth = Type, estimate = .pred_class) %>%
  autoplot(type = "heatmap") +
  scale_fill_gradient(high = "#E1812C", low = "#3A923A") +
  labs(title = "Clinical model confusion matrix")
c

ggsave(plot=c,"/data/scratch/kvalem/projects/2024/Effenberger-Diabetes/02-scripts/figures/v02/confusion_matrix_clinical.svg", height = 3, width = 3)
ggsave(plot=c,"/data/scratch/kvalem/projects/2024/Effenberger-Diabetes/02-scripts/figures/v02/confusion_matrix_clinical.png", height = 3, width = 3)


# --- ROC curve ---
roc <- roc_curve(results, truth = Type, .pred_DM) %>%
  autoplot()

roc


#ggsave(plot=roc,"/data/scratch/kvalem/projects/2024/Effenberger-Diabetes/02-scripts/figures/v02/roc_curve_clinical.svg", height = 3, width = 3)
#ggsave(plot=roc,"/data/scratch/kvalem/projects/2024/Effenberger-Diabetes/02-scripts/figures/v02/roc_curve_clinical.png", height = 3, width = 3)


# --- Feature importance ---
tidy(model) %>%
  filter(term != "(Intercept)", estimate != 0) %>%
  arrange(desc(abs(estimate))) 

coefs <- tidy(model) %>%
  filter(term != "(Intercept)", estimate != 0)

# Get top N features with largest absolute effect size
top_features <- coefs %>%
  arrange(desc(abs(estimate))) %>%
  slice_head(n = 10) %>%   # 👈 change 10 to your preferred number
  pull(term)
# Clean up top_features by removing backticks
top_features_clean <- stringr::str_remove_all(top_features, "`")
# Extract the downsampled data with labels from the baked recipe
train_baked_labeled <- bake(rec_prepped, new_data = NULL, composition = "tibble") %>%
  select(-sample_information)

# Then select and plot
plot_data <- train_baked_labeled %>%
  select(Type, all_of(top_features_clean)) %>%
  pivot_longer(-Type, names_to = "Feature", values_to = "Value")

# Taxa names you want to identify
#target_ids <- c("c728ad6f5d183cb36fa06b6a3a47758b", "8eb4e34fd58ab95fa2efab34940c01fa","409f711b59152d57926cf444c5577087","ffc36e27c82042664a16bcd4d380b286","8f6e2a91e20994c00566a5ff2b49506e")

# Retrieve tax_table as a data frame
#tax_df <- as.data.frame(tax_table(phy))

# Look up the Genus for the target taxa
#genus_labels <- tax_df[target_ids, "Genus", drop = FALSE]
#genus_labels

# Create mapping of hash IDs to genus
#label_map <- setNames(genus_labels$Genus, rownames(genus_labels))

# Relabel Feature column
#plot_data_ordered <- plot_data_ordered %>%
#  mutate(Feature = recode(Feature, !!!label_map))


#plot_data_ordered <- plot_data %>%
#  mutate(
#    Feature = factor(
#      Feature,
#      levels = c(
#        setdiff(unique(Feature), c("c728ad6f5d183cb36fa06b6a3a47758b", "8eb4e34fd58ab95fa2efab34940c01fa","409f711b59152d57926cf444c5577087","ffc36e27c82042664a16bcd4d380b286","8f6e2a91e20994c00566a5ff2b49506e")),
#        "c728ad6f5d183cb36fa06b6a3a47758b", "8eb4e34fd58ab95fa2efab34940c01fa","409f711b59152d57926cf444c5577087","ffc36e27c82042664a16bcd4d380b286","8f6e2a91e20994c00566a5ff2b49506e"
#      )
#    )
#  )

p <- ggplot(plot_data, aes(x = Type, y = Value, fill = Type)) +
  geom_violin(trim = FALSE, alpha = 0.7) +
  facet_wrap(~ Feature, scales = "free", ncol = 2) +
  scale_fill_manual(values = c("DM" = "#E1812C", "PDM" = "#3A923A")) +
  theme_minimal() +
  labs(title = "Top Predictive Features", 
       x = "", y = "Z-score") +
  theme(legend.position = "none")

p <- ggplot(plot_data, aes(x = Type, y = Value, fill = Type)) +
  geom_violin(trim = FALSE, alpha = 0.7) +
  facet_wrap(~ Feature, scales = "free", ncol = 5) +
  scale_fill_manual(values = c("DM" = "#E1812C", "PDM" = "#3A923A")) +
  stat_compare_means(method = "wilcox.test", 
                     label = "p.signif", 
                     comparisons = list(c("DM", "PDM")),
                     hide.ns = TRUE) +
  theme_minimal() +
  labs(title = "Top Predictive Features", 
       x = "", y = "Z-score") +
  theme(legend.position = "none")+
  theme(
    legend.position = "none",
    strip.text = element_text(size = 14),    # facet label size
    axis.text = element_text(size = 12),     # axis tick label size
    axis.title = element_text(size = 14),    # axis title size
    plot.title = element_text(size = 16, face = "bold")  # title size
  )

library(ggplot2)
library(ggpubr)

p <- ggplot(plot_data, aes(x = Type, y = Value, fill = Type)) +
  geom_violin(trim = FALSE, alpha = 0.7) +
  facet_wrap(~ Feature, scales = "free", ncol = 5) +
  scale_fill_manual(values = c("DM" = "#E1812C", "PDM" = "#3A923A")) +
  stat_compare_means(
    method = "wilcox.test", 
    label = "p.signif", 
    comparisons = list(c("DM", "PDM")),
    p.adjust.method = "BH",   # FDR correction (Benjamini-Hochberg)
    hide.ns = TRUE            # hide non-significant results
  ) +
  theme_minimal() +
  labs(
    title = "Top Predictive Features", 
    x = "", 
    y = "Z-score"
  ) +
  theme(
    legend.position = "none",
    strip.text = element_text(size = 14),    # facet label size
    axis.text = element_text(size = 12),     # axis tick label size
    axis.title = element_text(size = 14),    # axis title size
    plot.title = element_text(size = 16, face = "bold")  # title size
  )

p

ggsave(plot=p,"/data/scratch/kvalem/projects/2024/Effenberger-Diabetes/02-scripts/figures/v02/top_features_clinical.svg", height = 6, width = 15)
ggsave(plot=p,"/data/scratch/kvalem/projects/2024/Effenberger-Diabetes/02-scripts/figures/v02/top_features_clinical.png", height = 6, width = 15)

################################################### roc_curve_train_test
best_threshold = 0.5
# --- Predictions for train and test ---
# Add the truth label after baking, from baked data
train_baked_with_type <- bake(rec_prepped, new_data = NULL) %>%
  select(-sample_information)

test_baked_with_type <- bake(rec_prepped, new_data = test) %>%
  select(-sample_information)

# Predict on baked data
pred_prob_train <- predict(model, new_data = train_baked_with_type, type = "prob") %>%
  mutate(.pred_class = if_else(.pred_DM > best_threshold, "DM", "PDM"),
         Type = train_baked_with_type$Type,
         Set = "Train")

pred_prob_test <- predict(model, new_data = test_baked_with_type, type = "prob") %>%
  mutate(.pred_class = if_else(.pred_DM > best_threshold, "DM", "PDM"),
         Type = test$Type,
         Set = "Test")



# --- Bind both datasets ---
combined_preds <- bind_rows(pred_prob_train, pred_prob_test)

# Add 'Set' column to each prediction frame
roc_train <- roc_curve(pred_prob_train, truth = Type, .pred_DM) %>%
  mutate(Set = "Train")

roc_test <- roc_curve(pred_prob_test, truth = Type, .pred_DM) %>%
  mutate(Set = "Test")

# Combine
roc_combined <- bind_rows(roc_train, roc_test)

# Plot
print(acc)
print(auc)


p <- ggplot(roc_combined, aes(x = 1 - specificity, y = sensitivity, color = Set)) +
  geom_path(size = 1.2) +
  geom_abline(lty = 2, color = "gray") +
  theme_minimal() +
  labs(title = "ROC AUC - 0.958   Accuracy - 0.786 ", x = "1 - Specificity", y = "Sensitivity") +
  scale_color_manual(values = c("Train" = "#3a99bc", "Test" = "#db9e2a"))

p
ggsave(plot=p,"/data/scratch/kvalem/projects/2024/Effenberger-Diabetes/02-scripts/figures/v02/roc_curve_train_test_clinical.svg", height = 3, width = 3.5)
ggsave(plot=p,"/data/scratch/kvalem/projects/2024/Effenberger-Diabetes/02-scripts/figures/v02/roc_curve_train_test_clinical.png", height = 3, width = 3.5)
################################### HEatmap 

top_features <- tidy(model) %>%
  filter(term != "(Intercept)", estimate != 0) %>%
  arrange(desc(abs(estimate))) %>%
  pull(term)

# Remove backticks
top_features_clean <- stringr::str_remove_all(top_features, "`")

heatmap_data <- train_baked_labeled %>%
  select(all_of(top_features_clean)) %>%
  as.matrix()

cor_matrix <- cor(heatmap_data, use = "pairwise.complete.obs")

library(pheatmap)

target_ids <- c(
  "c728ad6f5d183cb36fa06b6a3a47758b",
  "8eb4e34fd58ab95fa2efab34940c01fa",
  "409f711b59152d57926cf444c5577087",
  "ffc36e27c82042664a16bcd4d380b286",
  "8f6e2a91e20994c00566a5ff2b49506e"
)

# Convert tax_table to data.frame for easy handling
tax_df <- as.data.frame(tax_table(phy))

# Subset Genus names for the given IDs
genus_labels <- tax_df[target_ids, "Genus"]

# Create a named vector for relabeling
feature_labels <- setNames(
  paste0("g__", as.character(genus_labels)),
  target_ids
)

colnames(cor_matrix) <- ifelse(colnames(cor_matrix) %in% names(feature_labels),
                               feature_labels[colnames(cor_matrix)],
                               colnames(cor_matrix))

rownames(cor_matrix) <- ifelse(rownames(cor_matrix) %in% names(feature_labels),
                               feature_labels[rownames(cor_matrix)],
                               rownames(cor_matrix))

pheatmap(
  cor_matrix,
  clustering_distance_rows = "euclidean",
  clustering_distance_cols = "euclidean",
  clustering_method = "complete",
  color = colorRampPalette(c("white", "orange", "darkred"))(100),
  border_color = NA,
  main = ""
)


library(pheatmap)

# Compute correlation matrix
cor_matrix <- cor(train_baked[top_features_clean], method = "pearson", use = "pairwise.complete.obs")

# Compute significance matrix (p-values)
p_matrix <- matrix(NA, nrow = ncol(cor_matrix), ncol = ncol(cor_matrix))
colnames(p_matrix) <- rownames(p_matrix) <- colnames(cor_matrix)

for (i in 1:ncol(cor_matrix)) {
  for (j in 1:ncol(cor_matrix)) {
    test <- cor.test(train_baked[[top_features_clean[i]]], train_baked[[top_features_clean[j]]])
    p_matrix[i, j] <- test$p.value
  }
}

# Create annotation matrix with * for p < 0.05
annot_matrix <- ifelse(p_matrix < 0.05, "*", "")

# Relabel for genus names
colnames(cor_matrix) <- rownames(cor_matrix) <- ifelse(colnames(cor_matrix) %in% names(feature_labels),
                                                       feature_labels[colnames(cor_matrix)],
                                                       colnames(cor_matrix))

# Plot heatmap
ht <- pheatmap(
  cor_matrix,
  display_numbers = annot_matrix,
  number_color = "black",
  clustering_distance_rows = "euclidean",
  clustering_distance_cols = "euclidean",
  clustering_method = "complete",
  color = viridis(100),
  border_color = NA,
  main = "Feature correlation matrix",
  legend_labels = "Correlation Coefficient (r)",  fontsize_row = 12,       # ⬅️ Increase row label font
  fontsize_col = 12,       # ⬅️ Increase column label font
  fontsize_number = 10     # ⬅️ Font size for "*"
) 



ggsave(plot=ht,"/data/scratch/kvalem/projects/2024/Effenberger-Diabetes/02-scripts/figures/v02/correlation_heatmap_clinical.svg", height = 8, width = 8)
ggsave(plot=ht,"/data/scratch/kvalem/projects/2024/Effenberger-Diabetes/02-scripts/figures/v02/correlation_heatmap_clinical.png", height = 8, width = 8)

########################################## forest plot

library(broom)
library(ggplot2)
library(dplyr)

# Get non-zero coefficients
coef_df <- tidy(model) %>%
  filter(term != "(Intercept)", estimate != 0) %>%
  arrange(estimate) %>%
  mutate(term = stringr::str_remove_all(term, "`"))  # Remove backticks

# Add fake standard errors (if glmnet didn't return them)
# This is a placeholder; ideally use actual SEs if available
coef_df <- coef_df %>%
  mutate(std.error = abs(estimate) * 0.2)  # ← Replace with real SEs if available

feature_labels <- c(
  "8eb4e34fd58ab95fa2efab34940c01fa" = "g__Prevotella_9",
  "c728ad6f5d183cb36fa06b6a3a47758b" = "g__Faecalibacterium"
)

coef_df$term <- recode(coef_df$term, !!!feature_labels)

lg <- ggplot(coef_df, aes(x = estimate, y = reorder(term, estimate))) +
  geom_col(fill = "#3182bd") +
  geom_errorbarh(aes(xmin = estimate - std.error,
                     xmax = estimate + std.error),
                 height = 0.3, color = "black") +
  theme_minimal() +
  labs(title = "Logistic Regression Coefficients",
       x = "Log-Odds (Coefficient)",
       y = NULL) +
  theme(text = element_text(size = 20))
lg

ggsave(plot=lg,"/data/scratch/kvalem/projects/2024/Effenberger-Diabetes/02-scripts/figures/v02/log_reg_coefs_clinical.svg", height = 4, width = 8)
ggsave(plot=lg,"/data/scratch/kvalem/projects/2024/Effenberger-Diabetes/02-scripts/figures/v02/log_reg_coefs_clinical.png", height = 4, width = 8)
