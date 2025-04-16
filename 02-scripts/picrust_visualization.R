# install.packages("devtools")
#devtools::install_github("cafferychen777/ggpicrust2")


# If you want to analyze the abundance of KEGG pathways instead of KO within the pathway, please set `ko_to_kegg` to TRUE.
# KEGG pathways typically have more descriptive explanations.

library(readr)
library(ggpicrust2)
library(tibble)
library(tidyverse)
library(ggprism)
library(patchwork)

# Load necessary data: abundance data and metadata

abundance_file <- "/data/scratch/kvalem/projects/2024/Effenberger-Diabetes/02-scripts/picrust2-2.6.1/picrust2_out_pipeline/KO_metagenome_out/filtered_abundance.tsv"
abundance_df <- read_delim(
  abundance_file,
  delim = "\t",
  escape_double = FALSE,
  trim_ws = TRUE
)

metacyc_abundance <- "/data/scratch/kvalem/projects/2024/Effenberger-Diabetes/02-scripts/picrust2-2.6.1/picrust2_out_pipeline/EC_metagenome_out/metacyc_abundance_filtered.tsv"
metacyc_abundance_df <- read_delim(
  metacyc_abundance,
  delim = "\t",
  escape_double = FALSE,
  trim_ws = TRUE
)

ec_name_map <- read.delim("/data/scratch/kvalem/projects/2024/Effenberger-Diabetes/02-scripts/picrust2-2.6.1/picrust2/default_files/description_mapfiles/ec_name.txt.gz", header = FALSE, sep = "\t", col.names = c("feature", "description"))
metadata <- read_delim(
  "/data/scratch/kvalem/projects/2024/Effenberger-Diabetes/02-scripts/tables/original/20011_SampleInfo.csv",
  delim = ",",
  escape_double = FALSE,
  trim_ws = TRUE
)


kegg_abundance <- ko2kegg_abundance("/data/scratch/kvalem/projects/2024/Effenberger-Diabetes/02-scripts/picrust2-2.6.1/picrust2_out_pipeline/KO_metagenome_out/pred_metagenome_unstrat.tsv")



############################################################################################# KEGG barplot
metadata_filtered <- metadata[metadata$Type %in% c("Diabetes mellitus Typ1", "pankreopriver Diabetes"), ]

# If sample IDs are in 'SampleID' column
sample_ids <- metadata_filtered$`IMGM ID`

# Filter columns in kegg_abundance to keep only matching samples
# Keep the first column (pathway/KO ID), and matching samples
kegg_abundance_filtered <- kegg_abundance[, c(TRUE, colnames(kegg_abundance)[-1] %in% sample_ids)]
kegg_abundance_filtered <- kegg_abundance_filtered[, colnames(kegg_abundance_filtered) != "s20011_0001"]


daa_results_df <- pathway_daa(abundance = kegg_abundance_filtered, metadata = metadata_filtered, group = "Type", daa_method = "DESeq2", select = NULL, reference = "Kontrolle")
daa_sub_method_results_df <- daa_results_df[daa_results_df$method == "DESeq2" , ]

daa_annotated_sub_method_results_df <- pathway_annotation(pathway = "KO", daa_results_df = daa_sub_method_results_df, ko_to_kegg = TRUE)
filtered_df <- daa_annotated_sub_method_results_df[grepl("Metabolism;", daa_annotated_sub_method_results_df$pathway_class), ]

p <- pathway_errorbar(abundance = kegg_abundance_filtered, daa_results_df = filtered_df, Group = metadata_filtered$Type, p_values_threshold = 0.01, order = "pathway_class", select = NULL, ko_to_kegg = TRUE, p_value_bar = TRUE, colors = c("#3685BC","#D93C3E"), x_lab = "pathway_name")

############################################################################################# MetaCyc Pathway and EC




metadata_filtered <- metadata[metadata$Type %in% c("Diabetes mellitus Typ1", "pankreopriver Diabetes"), ]

# If sample IDs are in 'SampleID' column
sample_ids <- metadata_filtered$`IMGM ID`


metacyc_abundance <-  read_delim("/data/scratch/kvalem/projects/2024/Effenberger-Diabetes/02-scripts/picrust2-2.6.1/picrust2_out_pipeline/EC_metagenome_out/pred_metagenome_unstrat.tsv", delim = "\t",
                                 escape_double = FALSE,
                                 trim_ws = TRUE)
metacyc_abundance_filtered <- metacyc_abundance[, c(TRUE, colnames(metacyc_abundance)[-1] %in% sample_ids)]
metacyc_abundance_filtered <- metacyc_abundance_filtered[, colnames(metacyc_abundance_filtered) != "s20011_0001"]
names(metacyc_abundance_filtered)[names(metacyc_abundance_filtered) == "function"] <- "pathway"

metacyc_daa_results_df <- pathway_daa(abundance = metacyc_abundance_filtered %>% column_to_rownames("pathway"), metadata = metadata_filtered, group = "Type", daa_method = "DESeq2")
metacyc_daa_annotated_results_df <- pathway_annotation(pathway = "MetaCyc", daa_results_df = metacyc_daa_results_df, ko_to_kegg = FALSE)
daa_sub_method_metacyc_daa_annotated_results_df <- metacyc_daa_annotated_results_df[metacyc_daa_annotated_results_df$method == "DESeq2" , ]


filtered_daa_results_annotated <- metacyc_daa_annotated_results_df %>%
  left_join(ec_name_map, by = "feature")

filtered_daa_results_annotated <- filtered_daa_results_annotated %>%
  mutate(description = description.y) %>%
  select(-description.x, -description.y)


metacyc_abundance_filtered_named <- metacyc_abundance_filtered %>%
  column_to_rownames("pathway")  # or "pathway" or "description", depending on column name

top_features <- filtered_daa_results_annotated %>%
  arrange(p_adjust) %>%
  slice_head(n = 30) %>%
  pull(feature)

pathway_errorbar_patched <- ggpicrust2::pathway_errorbar

pathway_errorbar_patched <- function (abundance, daa_results_df, Group, ko_to_kegg = FALSE, 
                                      p_values_threshold = 0.05, order = "group", select = NULL, 
                                      p_value_bar = TRUE, colors = NULL, x_lab = NULL) 
{
  if (!is.matrix(abundance) && !is.data.frame(abundance)) {
    stop("'abundance' must be a matrix or data frame")
  }
  if (!is.data.frame(daa_results_df)) {
    stop("'daa_results_df' must be a data frame")
  }
  required_cols <- c("feature", "method", "group1", "group2", "p_adjust")
  missing_cols <- setdiff(required_cols, colnames(daa_results_df))
  if (length(missing_cols) > 0) {
    stop("Missing required columns in daa_results_df: ", paste(missing_cols, collapse = ", "))
  }
  if (length(Group) != ncol(abundance)) {
    stop("Length of Group must match number of columns in abundance matrix")
  }
  
  # ⛔ method check REMOVED
  
  if (nlevels(factor(daa_results_df$group1)) != 1 || nlevels(factor(daa_results_df$group2)) != 1) {
    stop("The 'group1' or 'group2' column in the 'daa_results_df' data frame contains more than one group. Please filter each to contain only one group.")
  }
  
  if (is.null(x_lab)) {
    x_lab <- ifelse(ko_to_kegg, "pathway_name", "description")
  }
  
  if (!(x_lab %in% colnames(daa_results_df))) {
    stop("The 'x_lab' you defined does not exist as a column in the 'daa_results_df' data frame.")
  }
  
  daa_results_df <- daa_results_df[!is.na(daa_results_df[, x_lab]), ]
  daa_results_filtered_df <- daa_results_df[daa_results_df$p_adjust < p_values_threshold, ]
  
  if (!is.null(select)) {
    daa_results_filtered_sub_df <- daa_results_filtered_df[daa_results_filtered_df$feature %in% select, ]
  } else {
    daa_results_filtered_sub_df <- daa_results_filtered_df
  }
  
  if (nrow(daa_results_filtered_sub_df) > 30) {
    stop("Too many significant features (>30). Use 'select' to reduce the number of features.")
  }
  if (nrow(daa_results_filtered_sub_df) == 0) {
    stop("No significant features found.")
  }
  
  relative_abundance_mat <- apply(t(as.matrix(abundance)), 1, function(x) x / sum(x))
  sub_relative_abundance_mat <- relative_abundance_mat[rownames(relative_abundance_mat) %in% daa_results_filtered_sub_df$feature, ]
  error_bar_matrix <- cbind(sample = colnames(sub_relative_abundance_mat), group = Group, t(sub_relative_abundance_mat))
  error_bar_df <- as.data.frame(error_bar_matrix)
  error_bar_df$group <- factor(Group, levels = levels(as.factor(Group)))
  error_bar_pivot_longer_df <- tidyr::pivot_longer(error_bar_df, -c(sample, group))
  error_bar_pivot_longer_tibble <- dplyr::mutate(error_bar_pivot_longer_df,
                                                 group = as.factor(group),
                                                 sample = factor(sample),
                                                 name = factor(name),
                                                 value = as.numeric(value))
  
  error_bar_summary <- error_bar_pivot_longer_tibble %>%
    dplyr::group_by(name, group) %>%
    dplyr::summarise(mean = mean(value), sd = stats::sd(value), .groups = "drop")
  
  if (order == "group") {
    daa_results_filtered_sub_df$pro <- 1
    for (i in levels(error_bar_summary$name)) {
      sub <- error_bar_summary[error_bar_summary$name == i, ]
      top_group <- as.vector(sub[sub$mean == max(sub$mean), ]$group)
      daa_results_filtered_sub_df[daa_results_filtered_sub_df$feature == i, ]$pro <- top_group
    }
    order_idx <- order(daa_results_filtered_sub_df$pro, daa_results_filtered_sub_df$p_adjust)
  } else if (order == "p_values") {
    order_idx <- order(daa_results_filtered_sub_df$p_adjust)
  } else {
    order_idx <- order(daa_results_filtered_sub_df$feature)
  }
  
  daa_results_filtered_sub_df <- daa_results_filtered_sub_df[order_idx, ]
  error_bar_summary_ordered <- do.call(rbind, lapply(daa_results_filtered_sub_df$feature, function(i) {
    error_bar_summary[error_bar_summary$name == i, ]
  }))
  
  error_bar_summary_ordered[, x_lab] <- rep(daa_results_filtered_sub_df[, x_lab], each = length(unique(error_bar_summary$group)))
  error_bar_summary_ordered$name <- factor(error_bar_summary_ordered$name, levels = rev(daa_results_filtered_sub_df$feature))
  
  if (is.null(colors)) {
    colors <- c("#d93c3e", "#3685bc", "#6faa3e", "#e8a825", "#c973e6", "#ee6b3d", "#2db0a7", "#f25292")[1:nlevels(as.factor(Group))]
  }
  
  error_bar_summary_ordered$name <- factor(
    daa_results_filtered_sub_df[[x_lab]][match(as.character(error_bar_summary_ordered$name),
                                               as.character(daa_results_filtered_sub_df$feature))],
    levels = rev(daa_results_filtered_sub_df[[x_lab]])
  )
  
  # Debug print to confirm!
  print("✅ Label mapping applied:")
  print(head(error_bar_summary_ordered$name))
  
  bar_plot <- ggplot2::ggplot(error_bar_summary_ordered,
                              ggplot2::aes(mean, name, fill = group)) +
    ggplot2::geom_errorbar(ggplot2::aes(xmax = mean + sd, xmin = 0),
                           position = ggplot2::position_dodge(width = 0.8),
                           width = 0.5, linewidth = 0.5, color = "black") +
    ggplot2::geom_bar(stat = "identity", position = ggplot2::position_dodge(width = 0.8), width = 0.8) +
    ggplot2::scale_fill_manual(values = colors) +
    ggplot2::labs(x = "Relative Abundance", y = NULL) +
    ggplot2::theme_minimal() +
    ggplot2::theme(legend.position = "top",
                   axis.text.y = ggplot2::element_text(size = 10, color = "black"),
                   axis.title.x = ggplot2::element_text(size = 10, color = "black"))

  
  
  
  
  return(bar_plot)
}





p2 <- pathway_errorbar_patched(
  abundance = metacyc_abundance_filtered_named,
  daa_results_df = filtered_daa_results_annotated,
  Group = metadata_filtered$Type,
  ko_to_kegg = FALSE,
  p_values_threshold = 0.05,
  order = "group",
  select = top_features,
  p_value_bar = TRUE,
  colors = NULL,
  x_lab = "description"
)


################################################################################