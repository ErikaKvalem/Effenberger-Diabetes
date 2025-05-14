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

ggpaired(df, cond1 = "glukose1", cond2 = "glukose2", id = "pat_id", line.color = "gray") +
  labs(title = "BMI before and after")


df_bmi <- df %>%
  filter(!is.na(bmi1) & !is.na(bmi2))  # Keep only complete cases for BMI

type_colors <- c("DM" = "#E1812C", "PDM" = "#3A923A", "K" = "#3274A1")

# Plot
ggpaired(df_bmi,
         cond1 = "bmi1",
         cond2 = "bmi2",
         id = "pat_id",
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
rownames(binary_data) <- as.character(rownames(binary_data))
rownames(annotation_row) <- as.character(rownames(annotation_row))
binary_matrix <- as.matrix(binary_data)


svg("/data/scratch/kvalem/projects/2024/Effenberger-Diabetes/02-scripts/figures/v02/heatmap_binary.svg", height = 10, width = 5)


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

library(dplyr)
library(tidyr)
library(ggplot2)

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

df_bool$probennummer <- rownames(df_bool)

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
ggsave(plot = p,"/data/scratch/kvalem/projects/2024/Effenberger-Diabetes/02-scripts/figures/v02/boolean_count.svg", height = 10, width = 10)
ggsave(plot = p,"/data/scratch/kvalem/projects/2024/Effenberger-Diabetes/02-scripts/figures/v02/boolean_count.png", height = 10, width = 10)


################################################ LAB BLOOD TEST

df_lab <- read_csv("/data/projects/2024/Effenberger-Diabetes/data/laboratory_variables.csv")%>%
  clean_names() 


