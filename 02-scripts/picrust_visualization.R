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

library(KEGGREST)
library(dplyr)
library(purrr)
library(tibble)
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

filtered_df <- daa_annotated_sub_method_results_df %>%
  as.data.frame() %>%
  filter(p_adjust < 0.05) %>%
  slice_head(n = 50)

filtered_df <- daa_annotated_sub_method_results_df[
  grepl("Organismal Systems;|Metabolism; Nucleotide metabolism", daa_annotated_sub_method_results_df$pathway_class),
]

################################################################################FUNCTION PATHWAY ERRORBAR
my_pathway_errorbar <- function(abundance, daa_results_df, Group, 
                                p_values_threshold = 0.05, 
                                order = "group", select = NULL, 
                                colors = NULL, x_lab = NULL) {
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(patchwork)
  
  if (!is.matrix(abundance) && !is.data.frame(abundance)) {
    stop("'abundance' must be a matrix or data frame")
  }
  if (!is.data.frame(daa_results_df)) {
    stop("'daa_results_df' must be a data frame")
  }
  if (length(Group) != ncol(abundance)) {
    stop("Length of Group must match number of columns in abundance matrix")
  }
  
  if (is.null(x_lab)) x_lab <- "description"
  if (!(x_lab %in% colnames(daa_results_df))) {
    stop("The column specified in 'x_lab' does not exist in daa_results_df.")
  }
  
  if (is.null(colors)) {
    colors <- c("#d93c3e", "#3685bc", "#6faa3e")[1:nlevels(as.factor(Group))]
  }
  
  # Filter significant features
  daa_filtered <- daa_results_df[daa_results_df$p_adjust < p_values_threshold, ]
  daa_filtered <- daa_filtered[!is.na(daa_filtered[[x_lab]]), ]
  
  if (!is.null(select)) {
    daa_filtered <- daa_filtered[daa_filtered$feature %in% select, ]
  }
  if (nrow(daa_filtered) == 0) stop("No significant features to plot.")
  
  # Filter abundance matrix
  abundance_sub <- abundance[rownames(abundance) %in% daa_filtered$feature, , drop = FALSE]
  
  # Compute relative abundance per sample
  rel_abund <- apply(abundance_sub, 2, function(x) x / sum(x))
  
  # Transpose and convert to long format
  rel_abund_df <- as.data.frame(t(rel_abund))
  rel_abund_df$sample <- colnames(abundance)
  rel_abund_df$group <- Group
  
  long_df <- rel_abund_df %>%
    pivot_longer(
      cols = -c(sample, group),
      names_to = "feature",
      values_to = "value"
    )
  
  # Summary stats
  summary_df <- long_df %>%
    group_by(feature, group) %>%
    summarise(mean = mean(value), sd = sd(value), .groups = "drop")
  
  # Reorder features and match p-values
  feature_order <- rev(unique(summary_df$feature))
  summary_df$feature <- factor(summary_df$feature, levels = feature_order)
  daa_filtered <- daa_filtered[match(levels(summary_df$feature), daa_filtered$feature), ]
  
  # Labels for y-axis
  summary_df$label <- rep(daa_filtered[[x_lab]], each = length(unique(summary_df$group)))
  #daa_filtered$p_label <- substr(formatC(daa_filtered$p_adjust, format = "e", digits = 2), 1, 6)
  
  daa_filtered$p_label <- formatC(daa_filtered$p_adjust, format = "e", digits = 2)
  
  daa_filtered$feature <- factor(daa_filtered$feature, levels = feature_order)
  
  # Bar Plot
  bar_plot <- ggplot(summary_df, aes(x = mean, y = feature, fill = group)) +
    geom_errorbar(aes(xmin = 0, xmax = mean + sd), 
                  position = position_dodge(width = 0.8),
                  width = 0.4, linewidth = 0.5, color = "black") +
    geom_bar(stat = "identity",
             position = position_dodge(width = 0.8),
             width = 0.7) +
    scale_fill_manual(values = colors) +
    scale_x_continuous(expand = c(0, 0)) +
    scale_y_discrete(labels = rev(daa_filtered[[x_lab]])) +
    labs(x = "Relative Abundance", y = NULL) +
    theme_bw(base_size = 16) +
    theme(
      axis.text = element_text(size = 14, color = "black"),
      axis.text.x = element_text(angle = 45, hjust = 1),
      axis.title.x = element_text(size = 16),
      axis.text.y = element_text(size = 14, color = "black"),
      legend.position = "top",
      legend.title = element_blank(),
      legend.text = element_text(size = 12),
      panel.grid.major.y = element_blank(),
      panel.grid.major.x = element_blank()
    )
  
  # P-value text plot
  pval_plot <- ggplot(daa_filtered, aes(x = 1, y = feature)) +
    geom_text(aes(label = p_label), size = 4.2, hjust = 0) +
    scale_y_discrete(labels = NULL) +
    labs(x = "adj. p-value", y = NULL) +
    xlim(1, 1.5) +  # Add extra space to the right
    theme_void(base_size = 16) +
    theme(
      axis.title.x = element_text(size = 14, hjust = 0),
      plot.margin = margin(5, 20, 5, 5)  # extra right margin
    )
  
  # Combine with patchwork
  combined_plot <- bar_plot + pval_plot + plot_layout(ncol = 2, widths = c(3, 0.8))
  return(combined_plot)
}




################################################################################

filtered_df %>% 
  filter(p_adjust < 0.05) %>% 
  pull(pathway_name) %>% 
  is.na() %>% 
  all()


p <- my_pathway_errorbar(
  abundance = kegg_abundance_filtered,
  daa_results_df = filtered_df,
  Group = metadata_filtered$Type,
  x_lab = "pathway_name"
)



p <- my_pathway_errorbar(abundance = kegg_abundance_filtered, daa_results_df = filtered_df, Group = metadata_filtered$Type, p_values_threshold = 0.01, order = "pathway_class", select = NULL,  colors = c("#3685BC","#D93C3E"), x_lab = "pathway_name")

p <- pathway_errorbar(abundance = kegg_abundance_filtered, daa_results_df = filtered_df, Group = metadata_filtered$Type, p_values_threshold = 0.01, order = "pathway_class", select = NULL, ko_to_kegg = TRUE, p_value_bar = TRUE, colors = c("#3685BC","#D93C3E"), x_lab = "pathway_name")



ggsave("/data/scratch/kvalem/projects/2024/Effenberger-Diabetes/02-scripts/figures/v02/pathway_plot.png", plot = p, width = 20, height = 10, dpi = 300)

# Save as SVG
ggsave("/data/scratch/kvalem/projects/2024/Effenberger-Diabetes/02-scripts/figures/v02/pathway_plot.svg", plot = p, width = 8, height = 10, device = "svg")
############################################################################################# MetaCyc Pathway and EC BARPLOT




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
  ko_to_kegg = TRUE,
  p_values_threshold = 0.05,
  order = "group",
  select = top_features,
  p_value_bar = TRUE,
  colors = NULL,
  x_lab = "description"
)

p2


# STEP 1: Get top EC numbers from your annotated DAA results
ec_features <- filtered_daa_results_annotated %>%
  filter(!is.na(feature)) %>%
  select(feature, description) %>%
  distinct() %>%
  filter(grepl("^\\d+\\.\\d+\\.\\d+\\.\\d+$", feature))  # Keep valid EC format

# STEP 2: Create a lookup function to query KEGG for each EC number
get_kegg_info <- function(ec_number) {
  kegg_id <- paste0("ec:", ec_number)
  tryCatch({
    kegg_entry <- keggGet(kegg_id)[[1]]
    tibble(
      EC = ec_number,
      KEGG_Name = kegg_entry$NAME %>% paste(collapse = "; "),
      Pathways = kegg_entry$PATHWAY %>% paste(collapse = "; ")
    )
  }, error = function(e) {
    tibble(EC = ec_number, KEGG_Name = NA, Pathways = NA)
  })
}

# STEP 3: Apply the lookup to all EC numbers
ec_annotations <- map_dfr(ec_features$feature, get_kegg_info)



################################################################################ MetaCyc Pathway and EC HEATMAP

feature_with_p_0.05 <- filtered_daa_results_annotated %>% filter(p_adjust < 0.05)

feature_with_p_0.05 <- filtered_daa_results_annotated %>%
  as.data.frame() %>%
  mutate(across(where(~ inherits(.x, "Rle")), as.vector)) %>%  # convert Rle columns
  filter(p_adjust < 0.05) %>%
  arrange(p_adjust) %>%
  slice_head(n = 30)


# 2. Subset the abundance matrix to include only significant pathways
sig_abundance <- metacyc_abundance_filtered_named[rownames(metacyc_abundance_filtered_named) %in% feature_with_p_0.05$feature, ]

# 3. Replace rownames with descriptions
# Ensure the rownames are matched correctly to the descriptions
# Replace rownames with unique descriptions
descriptions <- feature_with_p_0.05$description[match(rownames(sig_abundance), feature_with_p_0.05$feature)]
rownames(sig_abundance) <- make.unique(descriptions)
custom_colors <- c( "#FFA555","#6ABC6A")
p3 <- my_pathway_heatmap(abundance = sig_abundance, metadata = metadata_filtered, group = "Type",show_legend = TRUE,colors = custom_colors)
p3
################################################################################HEATMAP patched
my_pathway_heatmap <- function (abundance, metadata, group, colors = NULL, font_size = 12, 
          show_row_names = TRUE, show_legend = TRUE, custom_theme = NULL) 
{
  if (!is.matrix(abundance) && !is.data.frame(abundance)) {
    stop("abundance must be a data frame or matrix")
  }
  abundance <- as.matrix(abundance)
  if (ncol(abundance) < 2) {
    stop("At least two samples are required for creating a heatmap")
  }
  if (nrow(abundance) < 1) {
    stop("At least one pathway is required")
  }
  if (is.null(colnames(abundance))) {
    colnames(abundance) <- paste0("Sample", seq_len(ncol(abundance)))
  }
  if (is.null(rownames(abundance))) {
    rownames(abundance) <- paste0("Pathway", seq_len(nrow(abundance)))
  }
  if (!is.data.frame(metadata)) {
    stop("metadata must be a data frame")
  }
  if (!is.character(group) || length(group) != 1) {
    stop("group must be a single character string")
  }
  if (!is.null(colors) && !is.character(colors)) {
    stop("colors must be NULL or a character vector of color codes")
  }
  group_levels <- unique(metadata[[group]])
  if (length(group_levels) < 2) {
    stop("At least two groups are required for comparison")
  }
  if (!group %in% colnames(metadata)) {
    stop(paste("group:", group, "must be a column in metadata"))
  }
  sample_name_col <- colnames(metadata)[sapply(colnames(metadata), 
                                               function(x) all(colnames(abundance) %in% metadata[[x]]))]
  metadata$sample_name <- metadata %>% select(all_of(c(sample_name_col))) %>% 
    pull()
  if (!all(colnames(abundance) %in% metadata$sample_name)) {
    stop("Samples in abundance and metadata must match")
  }
  #z_abundance <- t(apply(abundance, 1, scale))
  z_abundance <- t(apply(abundance, 1, function(x) {
    z <- scale(x)
    z[z > 1] <- 1
    z[z < -1] <- -1
    return(z)
  }))
  
  colnames(z_abundance) <- colnames(abundance)
  z_df <- as.data.frame(z_abundance)
  metadata <- metadata %>% as.data.frame()
  ordered_metadata <- metadata[order(metadata[, group]), ]
  ordered_sample_names <- ordered_metadata$sample_name
  order <- ordered_metadata$sample_name
  ordered_group_levels <- ordered_metadata %>% select(all_of(c(group))) %>% 
    pull()
  long_df <- z_df %>% tibble::rownames_to_column() %>% tidyr::pivot_longer(cols = -rowname, 
                                                                           names_to = "Sample", values_to = "Value") %>% left_join(metadata %>% 
                                                                                                                                     select(all_of(c("sample_name", group))), by = c(Sample = "sample_name"))
  long_df$Sample <- factor(long_df$Sample, levels = order)
  breaks <- range(long_df$Value, na.rm = TRUE)
  if (is.null(colors)) {
    colors <- c("#d93c3e", "#3685bc", "#6faa3e", "#e8a825", 
                "#c973e6", "#ee6b3d", "#2db0a7", "#f25292")
  }
  p <- ggplot2::ggplot(data = long_df, mapping = ggplot2::aes(x = Sample, 
                                                              y = rowname, fill = Value)) + ggplot2::geom_tile() + 
    ggplot2::scale_fill_gradient2(low = "#2166ac", mid = "white", 
                                  high = "#b2182b", midpoint = 0,limits = c(-1, 1), name = "Z Score") + ggplot2::labs(x = NULL, 
                                                                                  y = NULL) + ggplot2::scale_y_discrete(expand = c(0, 0), 
                                                                                                                        position = "left") + ggplot2::scale_x_discrete(expand = c(0, 
                                                                                                                                                                                  0)) + ggplot2::theme(axis.text.x = ggplot2::element_blank(), 
                                                                                                                                                                                                       axis.text.y = ggplot2::element_text(size = font_size, 
                                                                                                                                                                                                                                           color = "black"), axis.ticks = ggplot2::element_blank(), 
                                                                                                                                                                                                       axis.text = ggplot2::element_text(color = "black", size = 10, 
                                                                                                                                                                                                                                         face = "bold"), panel.spacing = unit(0, "lines"), 
                                                                                                                                                                                                  legend.title = ggplot2::element_text(size = 12, color = "black", 
                                                                                                                                                                                                                                            face = "bold"), legend.text = ggplot2::element_text(size = 12, 
                                                                                                                                                                                                                                                                                                color = "black", face = "bold"), panel.background = ggplot2::element_blank(), 
                                                                                                                                                                                                       legend.margin = ggplot2::margin(l = 0, unit = "cm"), 
                                                                                                                                                                                                       strip.text = element_text(size = 12, face = "bold")) + 
    ggplot2::guides(fill = ggplot2::guide_colorbar(direction = "vertical", 
                                                   reverse = F, barwidth = unit(0.6, "cm"), barheight = unit(9, 
                                                                                                             "cm"), title = "Z Score", title.position = "top", 
                                                   title.hjust = -1, ticks = TRUE, label = TRUE)) + 
    ggh4x::facet_nested(cols = vars(!!sym(group)), space = "free", 
                        scale = "free", switch = "x", strip = ggh4x::strip_nested(background_x = ggh4x::elem_list_rect(fill = colors)))
  if (!show_row_names) {
    p <- p + theme(axis.text.y = element_blank())
  }
  if (!show_legend) {
    p <- p + theme(legend.position = "none")
  }
  if (!is.null(custom_theme)) {
    p <- p + custom_theme
  }
  cat("The Sample Names in order from left to right are:\n")
  cat(ordered_sample_names, sep = ", ")
  cat("\n")
  cat("The Group Levels in order from left to right are:\n")
  cat(ordered_group_levels, sep = ", ")
  cat("\n")
  return(p)
}


################################################################################


names(metadata_filtered)[names(metadata_filtered) == "IMGM ID"] <- "sample_name"

p4 <- pathway_pca(abundance = sig_abundance, metadata = metadata_filtered, group = "Type")
p4
