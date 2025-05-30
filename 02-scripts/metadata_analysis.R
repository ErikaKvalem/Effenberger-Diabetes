library(tidyverse)
library(readxl)
library(janitor)
library(GGally)
library(ggpubr)
library(dplyr)
library(stringr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(broom.mixed)
library(lmerTest)



df <- read_csv("/data/projects/2024/Effenberger-Diabetes/data/PDM merged 3.0_modified.csv")%>%
  clean_names()  # makes column names consistent and easier to handle
glimpse(df)
summary(df)

df <- df %>%
  mutate(Type = case_when(
    str_detect(probennummer, "PDM") ~ "PDM",
    str_detect(probennummer, "K") ~ "K",
    TRUE ~ "DM"
  ))


df <- df %>%
  mutate(bmi_change = bmi2 - bmi1,
         glucose_change = glukose2 - glukose1,
         hba1c_change = hb_a1c_dcct_ngsp_2 - hb_a1c_dcct_ngsp_1)

ggpaired(df, cond1 = "glukose1", cond2 = "glukose2", id = "probennummer", line.color = "gray") +
  labs(title = "BMI before and after")


df_bmi <- df %>%
  filter(!is.na(bmi1) & !is.na(bmi2))  # Keep only complete cases for BMI

type_colors <- c("DM" = "#E1812C", "PDM" = "#3A923A", "K" = "#3274A1")

# Plot
ggpaired(df_bmi,
         cond1 = "bmi1",
         cond2 = "bmi2",
         id = "probennummer",
         color = "Type",  # Add this line
         line.color = "gray",
         line.size = 0.4,
         palette = "jco") +
  stat_compare_means(paired = TRUE, method = "wilcox.test", label = "p.signif") +
  labs(title = "BMI before and after surgery")

table(df$verstorben)

table(df$pankreatektomie)

table(df$nikotin)

table(df$sex)

summary(df$age)

hist(df$age, main="Age Distribution", xlab="Age")

df$BMI_change <- df$bmi2 - df$bmi1
summary(df$BMI_change)
hist(df$BMI_change, main="Change in BMI", xlab="ΔBMI")

table(df$hyperlipidamie_2018)
table(df$hyperlipidamie_20181)

table(df$lipidsenker1)
table(df$lipidsenker2)


colSums(is.na(df))

barplot(table(df$verstorben),
        main = "verstorben Distribution",
        xlab = "verstorben",
        ylab = "Count",
        col = "skyblue")


cols_ja_nein <- names(df)[
  sapply(df, function(col) {
    is.character(col) && all(c("ja", "nein") %in% unique(col))
  })
]

################################################ EPIDEMIOLOGICAL 
library(pheatmap)
df_epi <- read_csv("/data/projects/2024/Effenberger-Diabetes/data/epidemiological_variables.csv")%>%
  clean_names() 

df_epi <- df_epi %>%
  mutate(Type = case_when(
    grepl("^PDM", probennummer) ~ "PDM",
    grepl("^K", probennummer) ~ "K",
    TRUE ~ "DM"
  ))

# Define columns for df_subset
subset_cols <- c("probennummer", "verstorben", "pankreatektomie", 
                 "c2", "nikotin", "sex","age")

# df_subset with only those columns
df_subset <- df_epi[, subset_cols]
df_subset$c2[df_subset$c2 == 20] <- 1

df_subset$sex[df_subset$sex == "m"] <- 1
df_subset$sex[df_subset$sex == "f"] <- 0

# Optionally, convert to numeric
df_subset$sex <- as.numeric(df_subset$sex)

df_subset_heat <- df_subset %>%
  column_to_rownames("probennummer") 


annotation_row <- data.frame(age = df_subset$age)
rownames(annotation_row) <- df_subset$probennummer

subset_cols <- c("probennummer", "verstorben", "pankreatektomie", 
                 "c2", "nikotin", "sex")

# df_subset with only those columns
binary_data <- df_subset[, subset_cols]
names(binary_data)[names(binary_data) == "c2"] <- "C2 (20g/d)"
names(binary_data)[names(binary_data) == "verstorben"] <- "Deceased"
names(binary_data)[names(binary_data) == "nikotin"] <- "Nicotine"
names(binary_data)[names(binary_data) == "sex"] <- "Sex"
names(binary_data)[names(binary_data) == "pankreatektomie"] <- "Pancreatectomy"


binary_data <- binary_data %>%
  column_to_rownames("probennummer") 

binary_colors <- colorRampPalette(c("purple", "orange"))(100)

# 4. Define viridis color scale for age annotation
annotation_colors <- list(
  age = viridis(100)
)

custom_order <- c("K1",  "K2",  "K3", "K4", "K5",   "K6",  "K7",  "K8",   "K9", "K10",
                  "DM1", "DM2",  "DM3", "DM4", "DM5",  "DM6", "DM7", "DM8", "DM9", "DM10", 
                  "DM11", "DM12", "DM13", "DM14", "DM15", "DM16", "DM17", "DM18", "DM19", "DM20", "DM21",
                  "PDM1", "PDM2", "PDM4", "PDM5", "PDM6", "PDM7", "PDM8", "PDM9", "PDM10", 
                  "PDM11", "PDM12", "PDM13", "PDM14", "PDM15", "PDM16", "PDM17", "PDM18")

# Ensure custom_order matches rownames
custom_order <- custom_order[custom_order %in% rownames(binary_matrix)]

# Reorder binary_matrix and annotation_row
binary_matrix <- binary_matrix[custom_order, , drop = FALSE]
annotation_row <- annotation_row[custom_order, , drop = FALSE]

rownames(binary_data) <- as.character(rownames(binary_data))
rownames(annotation_row) <- as.character(rownames(annotation_row))
binary_matrix <- as.matrix(binary_data)


#svg("/data/scratch/kvalem/projects/2024/Effenberger-Diabetes/02-scripts/figures/v02/heatmap_binary.svg", height = 10, width = 5)


pheatmap(
  mat = binary_matrix,
  color = binary_colors,
  annotation_row = annotation_row,
  annotation_colors = annotation_colors,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  main = ""
)
dev.off()

cols_to_exclude <- c("verstorben", "pankreatektomie", "c2", "nikotin", "sex", "age")
df_longitudinal <- df_epi[, !(names(df_epi) %in% cols_to_exclude)]
df_longitudinal$insulin1[df_longitudinal$insulin1 == "nein"] <- 0
df_longitudinal$insulin2[df_longitudinal$insulin2 == "nein"] <- 0


# 1. Identify boolean-like columns (2 unique values only)
bool_like_cols <- names(df_longitudinal)[sapply(df_longitudinal, function(col) {
  length(unique(na.omit(col))) == 2
})]



# 2. Subset and reshape data
df_bool <- dplyr::select(df_longitudinal, dplyr::all_of(bool_like_cols))
names(df_bool) <- sub("(.*)(.)$", "\\1_\\2", names(df_bool))

rownames(df_bool) <- as.character(df_longitudinal$probennummer)

df_longitudinal <- df_longitudinal %>%
  mutate(Type = case_when(
    grepl("^PDM", probennummer) ~ "PDM",
    grepl("^K", probennummer) ~ "K",
    TRUE ~ "DM"
  ))

df_longitudinal <- df_longitudinal[, !(names(df_longitudinal) %in%
                                         c("art_von_lipidsenker1", "b_blocker1", "ace_hemmer1", "diuretika1",
                                           "art_von_lipidsenker2", "b_blocker2", "ace_hemmer2", "diuretika2"))]


df_bool$probennummer <- df_longitudinal$probennummer

df_bool <- df_bool %>%
  left_join(df_longitudinal %>% dplyr::select(probennummer, Type), by = "probennummer")


# Assuming your boolean dataframe is named df_bool

df_bool <- df_bool %>%
  mutate(across(everything(), ~ as.numeric(as.character(.))))

df_bool_long <- df_bool %>%
  pivot_longer(cols = -c(probennummer, Type), names_to = "Variable", values_to = "Value") %>%
  mutate(Value = as.character(Value))

df_bool_long <- df_bool_long %>%
  filter(Type != "K")


df_bool_counts <- df_bool_long %>%
  group_by(Variable, Value, Type) %>%
  summarise(Count = n(), .groups = "drop")

# 3. Create stacked barplot

p <- ggplot(df_bool_counts, aes(x = Variable, y = Count, fill = Value)) +
  geom_bar(stat = "identity", position = "stack") +
  facet_wrap(~ Type) +
  labs(
    title = "",
    x = "", y = "Number patients", fill = ""
  ) +
  scale_fill_manual(values = c("0" = "purple", "1" = "orange")) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
p <- p + theme(
  panel.background = element_rect(fill = "white", color = NA),
  plot.background = element_rect(fill = "white", color = NA)
)
#ggsave(plot = p,"/data/scratch/kvalem/projects/2024/Effenberger-Diabetes/02-scripts/figures/v02/boolean_count.svg", height = 10, width = 10)
#ggsave(plot = p,"/data/scratch/kvalem/projects/2024/Effenberger-Diabetes/02-scripts/figures/v02/boolean_count.png", height = 10, width = 10)


################################################ LAB BLOOD TEST

df_lab <- read_csv("/data/projects/2024/Effenberger-Diabetes/data/laboratory_variables.csv")%>%
  clean_names() 

df_lab <- df_lab %>%
  mutate(Type = case_when(
    grepl("^PDM", probennummer) ~ "PDM",
    grepl("^K", probennummer) ~ "K",
    TRUE ~ "DM"
  ))



######################################################################### ALL 

df_all <- read_csv("/data/projects/2024/Effenberger-Diabetes/data/PDM merged 3.0_modified.csv")%>%
  clean_names() 

df_lab <- read_csv("/data/projects/2024/Effenberger-Diabetes/data/laboratory_variables.csv")%>%
  clean_names() 

var_names_lab <- colnames(df_lab)
base_vars_lab <- unique(str_remove(var_names_lab, "[12]$"))

df_epi <- read_csv("/data/projects/2024/Effenberger-Diabetes/data/epidemiological_variables.csv")%>%
  clean_names() 

var_names_epi <- colnames(df_epi)
base_vars_epi <- unique(str_remove(var_names_epi, "[12]$"))

#df_all <- df_all %>%
#  dplyr::select(where(~ !all(.x %in% c(0, 1), na.rm = TRUE)))



df_all <- df_all %>%
  mutate(Type = case_when(
    grepl("^PDM", probennummer) ~ "PDM",
    grepl("^K", probennummer) ~ "K",
    TRUE ~ "DM"
  ))


#Get all unique base variable names
var_names <- colnames(df_all)
base_vars <- unique(str_remove(var_names, "[12]$"))

base_vars <- setdiff(base_vars, c("probennummer", "Type","verstorben","pankreatektomie"))


var_names <- unlist(lapply(base_vars, function(v) paste0(v, c("1", "2"))))
var_names <- var_names[var_names %in% colnames(df_all)]  # keep existing

df_long <- df_all %>%
  mutate(ID = probennummer) %>%
  pivot_longer(
    cols = all_of(var_names),
    names_to = c(".value", "time"),
    names_pattern = "(.*)([12])"
  ) %>%
  mutate(time = as.integer(time))


# Step 4: Fit GLMMs for each variable and store results
results <- list()

for (var in base_vars) {
  # Check for existence of variable in long format
  if (!all(c(var, "time", "probennummer") %in% colnames(df_long))) next
  
  df_sub <- df_long %>%
    dplyr::select(probennummer, time,age, sex,Type, !!sym(var)) %>%
    dplyr::filter(!is.na(!!sym(var)), !is.na(Type))
  
  # Skip if not enough data
  if (nrow(df_sub) < 10) next
  
  # Fit model depending on type
  model <- NULL
  result <- NULL
  
  # Numeric -> LMM
  if (is.numeric(df_sub[[var]])) {
    model <- tryCatch(
      lmer(as.formula(paste(var, "~ time + Type + age + sex + (1 | probennummer)")), data = df_sub),
      error = function(e) NULL
    )
  }
  # Binary -> GLMM
  else if (all(df_sub[[var]] %in% c(0, 1))) {
    model <- tryCatch(
      glmer(as.formula(paste(var, "~ time + Type + age + sex + (1 | probennummer)")), data = df_sub, family = binomial),
      error = function(e) NULL
    )
  }
  
  if (!is.null(model)) {
    est <- broom.mixed::tidy(model, effects = "fixed", conf.int = TRUE) %>%
      dplyr::filter(term == "time") %>%
      dplyr::mutate(
        variable = var,
        p.value = 2 * (1 - pnorm(abs(estimate / std.error)))  # add this line
      ) %>%
      dplyr::select(variable, estimate, std.error, conf.low, conf.high, p.value)
    
    results[[var]] <- est
  }
}

# Combine all results
forest_df <- bind_rows(results)
forest_df <- forest_df %>%
  dplyr::mutate(
    sig = case_when(
      p.value < 0.001 ~ "***",
      p.value < 0.01  ~ "**",
      p.value < 0.05  ~ "*",
      TRUE            ~ ""
    )
  )



forest_df <- forest_df %>%
  dplyr::mutate(
    source = case_when(
      variable %in% base_vars_epi ~ "epi",
      variable %in% base_vars_lab ~ "lab",
      TRUE ~ "other"
    )
  )

name_map <- c(
  "age" = "Age",
  "hyperlipid_mie_2018" = "Hyperlipidemia (2018)",
  "arterielle_hyperotnie" = "Arterial Hypertension",
  "bmi" = "Body Mass Index (BMI)",
  "gr_e" = "Height",
  "gewicht" = "Weight",
  "lipidsenker" = "Lipid-lowering Medication",
  "rr_medikation" = "Blood Pressure Medication",
  "langzeit_insulin" = "Long-acting Insulin",
  "kurzzeit_insulin" = "Short-acting Insulin",
  "misch_inuslin" = "Mixed Insulin",
  "masld" = "MASLD",
  "khk" = "Coronary Heart Disease (CHD)",
  "ca" = "Cancer (unspecified)",
  "leukozyten" = "Leukocytes",
  "h_moglobin" = "Hemoglobin",
  "h_matokrit" = "Hematocrit",
  "thrombozyten" = "Platelets",
  "harnstoff" = "Urea",
  "creatinin_enzym_idms_" = "Creatinine (IDMS)",
  "glomerul_re_filtrationsrate" = "Glomerular Filtration Rate (GFR)",
  "bilirubin_gesamt" = "Total Bilirubin",
  "natrium" = "Sodium",
  "got_asat_" = "AST",
  "gpt_alat_" = "ALT",
  "gamma_gt" = "GGT",
  "alkalische_phosphatase" = "Alkaline Phosphatase",
  "lactat_dehydrogenase_ldh_" = "LDH",
  "c_reaktives_prot_crp_" = "CRP",
  "quicktest_pt_" = "Prothrombin Time",
  "inr_pt_" = "INR",
  "part_thrombopl_zeit_a_ptt_" = "aPTT",
  "fibrinogen_funkt_n_clauss" = "Fibrinogen (Clauss)",
  "albumin" = "Albumin",
  "glukose" = "Glucose",
  "hb_a1c_dcct_ngsp_" = "HbA1c (DCCT/NGSP)",
  "hb_a1c_ifcc_" = "HbA1c (IFCC)",
  "cholesterin" = "Total Cholesterol",
  "non_hdl_cholesterin" = "Non-HDL Cholesterol",
  "triglyceride" = "Triglycerides",
  "hdl_cholesterin" = "HDL Cholesterol",
  "ldl_cholesterin" = "LDL Cholesterol",
  "eisen" = "Iron",
  "ferritin" = "Ferritin",
  "transferrin" = "Transferrin",
  "transferrins_ttigung" = "Transferrin Saturation",
  "troponin_t_hoch_sens_" = "Troponin T (High Sensitivity)",
  "nt_pro_bnp" = "NT-proBNP"
)

forest_df <- forest_df %>%
  mutate(variable = recode(variable, !!!name_map))

label_colors <- setNames(
  ifelse(forest_df$source == "longitudinal", "red",
         ifelse(forest_df$source == "lab", "blue", "black")),
  forest_df$variable
)


# Step 5: Plot forest plot
p <- ggplot(forest_df, aes(x = estimate, y = reorder(variable, estimate))) +
  geom_point() +
  geom_errorbarh(aes(xmin = conf.low, xmax = conf.high), height = 0.2) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray") +
  labs(
    x = "Effect of Time (Follow up vs Baseline)",
    y = "Variable",
    title = "Forest Plot of Longitudinal GLMMs"
  ) +
  theme_minimal() 
 
p <- ggplot(forest_df, aes(x = estimate, y = reorder(variable, estimate))) +
  geom_point() +
  geom_errorbarh(aes(xmin = conf.low, xmax = conf.high), height = 0.2) +
  geom_text(aes(label = sig), hjust = -0.5, size = 4) +  # Add stars
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray") +
  labs(
    x = "Effect of Time (Follow up vs Baseline)",
    y = "",
    title = ""
  ) +
  theme_minimal()

p <- p + theme(
  panel.background = element_rect(fill = "white", color = NA),
  plot.background = element_rect(fill = "white", color = NA)
)
p

#ggsave(plot = p,"/data/scratch/kvalem/projects/2024/Effenberger-Diabetes/02-scripts/figures/v02/forest_plot_numerical.svg", height = 10, width = 5)
#ggsave(plot = p,"/data/scratch/kvalem/projects/2024/Effenberger-Diabetes/02-scripts/figures/v02/forest_plot_numerical.png", height = 10, width = 5)

########################################################################################
run_models_by_group <- function(df_group, group_label, base_vars) {
  results <- list()
  
  for (var in base_vars) {
    if (!all(c(var, "time", "probennummer", "age", "sex") %in% colnames(df_group))) next
    
    df_sub <- df_group %>%
      dplyr::select(probennummer, time, age, sex, !!sym(var)) %>%
      dplyr::filter(!is.na(!!sym(var)), !is.na(age), !is.na(sex))
    
    if (nrow(df_sub) < 10) next
    
    model <- NULL
    
    # Numeric variable
    if (is.numeric(df_sub[[var]])) {
      model <- tryCatch(
        lmer(as.formula(paste(var, "~ time + age + sex + (1 | probennummer)")), data = df_sub),
        error = function(e) NULL
      )
    }
    # Binary variable
    else if (all(df_sub[[var]] %in% c(0, 1), na.rm = TRUE)) {
      model <- tryCatch(
        glmer(as.formula(paste(var, "~ time + age + sex + (1 | probennummer)")), data = df_sub, family = binomial),
        error = function(e) NULL
      )
    }
    
    if (!is.null(model)) {
      est <- broom.mixed::tidy(model, effects = "fixed", conf.int = TRUE) %>%
        dplyr::filter(term == "time") %>%
        dplyr::mutate(
          variable = var,
          group = group_label,
          p.value = 2 * (1 - pnorm(abs(estimate / std.error)))
        ) %>%
        dplyr::select(variable, estimate, std.error, conf.low, conf.high, p.value, group)
      
      # Convert to odds ratio if GLMM
      if ("glmerMod" %in% class(model)) {
        est <- est %>%
          dplyr::mutate(
            estimate = exp(estimate),
            conf.low = exp(conf.low),
            conf.high = exp(conf.high)
          )
      }
      
      results[[var]] <- est
    }
  }
  
  bind_rows(results)
}

# 2. Prepare subsets and base_vars
df_pdm <- df_long %>% filter(Type == "PDM")
df_dm  <- df_long %>% filter(Type == "DM")

# Base variable names (you should define this earlier or extract from df_long)
all_vars <- colnames(df_long)
base_vars <- setdiff(all_vars, c("probennummer", "time", "Type", "sex", "age"))

# 3. Run the function for both groups
forest_pdm <- run_models_by_group(df_pdm, "PDM", base_vars)
forest_dm  <- run_models_by_group(df_dm,  "DM",  base_vars)

# 4. Combine results
forest_df <- bind_rows(forest_pdm, forest_dm)

# 5. Add significance stars
forest_df <- forest_df %>%
  mutate(
    sig = case_when(
      p.value < 0.001 ~ "***",
      p.value < 0.01  ~ "**",
      p.value < 0.05  ~ "*",
      TRUE            ~ ""
    )
  )
name_map <- c(
  "age" = "Age",
  "hyperlipid_mie_2018" = "Hyperlipidemia (2018)",
  "arterielle_hyperotnie" = "Arterial Hypertension",
  "bmi" = "Body Mass Index (BMI)",
  "gr_e" = "Height",
  "gewicht" = "Weight",
  "lipidsenker" = "Lipid-lowering Medication",
  "rr_medikation" = "Blood Pressure Medication",
  "langzeit_insulin" = "Long-acting Insulin",
  "kurzzeit_insulin" = "Short-acting Insulin",
  "misch_inuslin" = "Mixed Insulin",
  "masld" = "MASLD",
  "khk" = "Coronary Heart Disease (CHD)",
  "ca" = "Cancer (unspecified)",
  "leukozyten" = "Leukocytes",
  "h_moglobin" = "Hemoglobin",
  "h_matokrit" = "Hematocrit",
  "thrombozyten" = "Platelets",
  "harnstoff" = "Urea",
  "creatinin_enzym_idms_" = "Creatinine (IDMS)",
  "glomerul_re_filtrationsrate" = "Glomerular Filtration Rate (GFR)",
  "bilirubin_gesamt" = "Total Bilirubin",
  "natrium" = "Sodium",
  "got_asat_" = "AST",
  "gpt_alat_" = "ALT",
  "gamma_gt" = "GGT",
  "alkalische_phosphatase" = "Alkaline Phosphatase",
  "lactat_dehydrogenase_ldh_" = "LDH",
  "c_reaktives_prot_crp_" = "CRP",
  "quicktest_pt_" = "Prothrombin Time",
  "inr_pt_" = "INR",
  "part_thrombopl_zeit_a_ptt_" = "aPTT",
  "fibrinogen_funkt_n_clauss" = "Fibrinogen (Clauss)",
  "albumin" = "Albumin",
  "glukose" = "Glucose",
  "hb_a1c_dcct_ngsp_" = "HbA1c (DCCT/NGSP)",
  "hb_a1c_ifcc_" = "HbA1c (IFCC)",
  "cholesterin" = "Total Cholesterol",
  "non_hdl_cholesterin" = "Non-HDL Cholesterol",
  "triglyceride" = "Triglycerides",
  "hdl_cholesterin" = "HDL Cholesterol",
  "ldl_cholesterin" = "LDL Cholesterol",
  "eisen" = "Iron",
  "ferritin" = "Ferritin",
  "transferrin" = "Transferrin",
  "transferrins_ttigung" = "Transferrin Saturation",
  "troponin_t_hoch_sens_" = "Troponin T (High Sensitivity)",
  "nt_pro_bnp" = "NT-proBNP"
)

forest_df <- forest_df %>%
  mutate(variable = recode(variable, !!!name_map))


forest_df <- forest_df %>%
  group_by(group) %>%
  mutate(p.adj = p.adjust(p.value, method = "fdr")) %>%
  ungroup() %>%
  mutate(
    sig = case_when(
      p.adj < 0.001 ~ "***",
      p.adj < 0.01  ~ "**",
      p.adj < 0.05  ~ "*",
      TRUE          ~ ""
    )
  )


forest_df_sig <- forest_df %>%
  filter(p.adj < 0.05) %>%
  mutate(variable = recode(variable, !!!name_map))



# 6. Plot
p <- ggplot(forest_df_sig, aes(x = estimate, y = reorder(variable, estimate), color = group)) +
  geom_point(position = position_dodge(width = 0.6)) +
  geom_errorbarh(aes(xmin = conf.low, xmax = conf.high), height = 0.2, position = position_dodge(width = 0.6)) +
  geom_text(aes(label = sig), hjust = -0.5, size = 4, position = position_dodge(width = 0.6)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray") +
  scale_color_manual(
    values = c("DM" = "#E1812C", "PDM" = "#3A923A")
  ) +
  labs(
    x = "Mean effect (Follow up vs Baseline)",
    y = "",
    title = "",
    color = ""
  ) +
  theme_minimal() + 
  theme(
    axis.text.y = element_text(size = 12),  # 👈 Increase this number as needed
    axis.text.x = element_text(size = 10),  # optional: adjust x-axis too
    legend.text = element_text(size = 11)   # optional: adjust legend text
  )

p <- p + theme(
  panel.background = element_rect(fill = "white", color = NA),
  plot.background = element_rect(fill = "white", color = NA)
)

p <- p + 
  annotate("text", x = Inf, y = -Inf, hjust = 1.05, vjust = -1.5,
           label = "* p < 0.05   ** p < 0.01   *** p < 0.001",
           size = 4, color = "black", fontface = "italic")



p
#ggsave(plot = p,"/data/scratch/kvalem/projects/2024/Effenberger-Diabetes/02-scripts/figures/v02/forest_plot_numerical_groups_sig.svg", height = 5, width = 7)
#ggsave(plot = p,"/data/scratch/kvalem/projects/2024/Effenberger-Diabetes/02-scripts/figures/v02/forest_plot_numerical_groups_sig.png", height = 5, width = 7)

#write.csv(forest_df_sig, file = "/data/scratch/kvalem/projects/2024/Effenberger-Diabetes/02-scripts/tables/v02/forest_df_sig.csv", row.names = FALSE)
