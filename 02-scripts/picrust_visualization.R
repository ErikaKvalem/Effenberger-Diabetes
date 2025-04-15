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


# Ensure it's a data.frame
#metadata <- as.data.frame(metadata)

# Filter metadata for the two groups
#metadata_filtered <- metadata[metadata$Type %in% c("pankreopriver Diabetes", "Diabetes mellitus Typ1"), ]

#write.table(metadata_filtered, "/data/scratch/kvalem/projects/2024/Effenberger-Diabetes/02-scripts/tables/original/metadata_filtered.tsv", sep = "\t", quote = FALSE, row.names = FALSE)

metadata2 <- read_delim(
  "/data/scratch/kvalem/projects/2024/Effenberger-Diabetes/02-scripts/tables/original/20011_SampleInfo.csv",
  delim = ",",
  escape_double = FALSE,
  trim_ws = TRUE
)


# Get matching sample IDs
#selected_samples <- metadata_filtered$`IMGM ID`

#abundance_filtered <- abundance_df[, c("function", selected_samples)]
#write.table(abundance_filtered, "/data/scratch/kvalem/projects/2024/Effenberger-Diabetes/02-scripts/picrust2-2.6.1/picrust2_out_pipeline/KO_metagenome_out/filtered_abundance.tsv", sep = "\t", quote = FALSE, row.names = FALSE)

# Filter abundance data: keep only 'function' column + selected samples
#metacyc_abundance_df_filtered <- metacyc_abundance_df[, c("function", selected_samples)]
#write.table(metacyc_abundance_df_filtered, "/data/scratch/kvalem/projects/2024/Effenberger-Diabetes/02-scripts/picrust2-2.6.1/picrust2_out_pipeline/EC_metagenome_out/metacyc_abundance_filtered.tsv", sep = "\t", quote = FALSE, row.names = FALSE)





# If you want to analyze the abundance of KEGG pathways instead of KO within the pathway, please set `ko_to_kegg` to TRUE.
# KEGG pathways typically have more descriptive explanations.



# Run ggpicrust2 with input file path
results_file_input <- ggpicrust2(file = abundance_file,
                                 metadata = metadata,
                                 group = "Type", # For example dataset, group = "Environment"
                                 pathway = "KO",
                                 daa_method = "DESeq2",
                                 ko_to_kegg = TRUE,
                                 order = "pathway_class",
                                 p_values_bar = TRUE,
                                 x_lab = "pathway_name")

# Run ggpicrust2 with imported data.frame
abundance_data <- read_delim(abundance_file, delim = "\t", col_names = TRUE, trim_ws = TRUE)

# Run ggpicrust2 with input data
results_data_input <- ggpicrust2(data = abundance_data,
                                 metadata = metadata,
                                 group = "Type", # For example dataset, group = "Environment"
                                 pathway = "KO",
                                 daa_method = "LinDA",
                                 ko_to_kegg = TRUE,
                                 order = "pathway_class",
                                 p_values_bar = TRUE,
                                 x_lab = "pathway_name")

# Access the plot and results dataframe for the first DA method
example_plot <- results_file_input[[1]]$plot
example_results <- results_file_input[[1]]$results


# Analyze the EC or MetaCyc pathway

results_file_input <- ggpicrust2(
  data = metacyc_abundance_df,
  metadata = metadata,
  group = "Type",
  pathway = "MetaCyc",
  daa_method = "ALDEx2",
  ko_to_kegg = FALSE,
  order = "group",
  p_values_bar = FALSE,  # Disable plotting for now
  x_lab = "description"
)
results_file_input[[1]]$plot
results_file_input[[1]]$results




#############################################################################################
metadata <- read_delim(
  "/data/scratch/kvalem/projects/2024/Effenberger-Diabetes/02-scripts/tables/original/20011_SampleInfo.csv",
  delim = ",",
  escape_double = FALSE,
  trim_ws = TRUE
)



kegg_abundance <- ko2kegg_abundance("/data/scratch/kvalem/projects/2024/Effenberger-Diabetes/02-scripts/picrust2-2.6.1/picrust2_out_pipeline/KO_metagenome_out/pred_metagenome_unstrat.tsv")

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


# Workflow for MetaCyc Pathway and EC

metacyc_abundance <-  read_delim("/data/scratch/kvalem/projects/2024/Effenberger-Diabetes/02-scripts/picrust2-2.6.1/picrust2_out_pipeline/EC_metagenome_out/pred_metagenome_unstrat.tsv", delim = "\t",
                                 escape_double = FALSE,
                                 trim_ws = TRUE)
metacyc_abundance_filtered <- metacyc_abundance[, c(TRUE, colnames(metacyc_abundance)[-1] %in% sample_ids)]
metacyc_abundance_filtered <- metacyc_abundance_filtered[, colnames(metacyc_abundance_filtered) != "s20011_0001"]

metacyc_daa_results_df <- pathway_daa(abundance = metacyc_abundance_filtered %>% column_to_rownames("function"), metadata = metadata_filtered, group = "Type", daa_method = "DESeq2")
metacyc_daa_annotated_results_df <- pathway_annotation(pathway = "MetaCyc", daa_results_df = metacyc_daa_results_df, ko_to_kegg = FALSE)
daa_sub_method_metacyc_daa_annotated_results_df<- metacyc_daa_annotated_results_df[metacyc_daa_annotated_results_df$method == "DESeq2" , ]


#pathway_errorbar(abundance = metacyc_abundance_filtered %>% column_to_rownames("function"), daa_results_df = daa_sub_method_metacyc_daa_annotated_results_df, Group = metadata_filtered$Type, ko_to_kegg = FALSE, p_values_threshold = 0.05, order = "group", select = NULL, p_value_bar = TRUE, colors = NULL, x_lab = "description")


feature_with_p_0.05 <- metacyc_daa_results_df %>% filter(p_adjust < 0.05)
pathway_heatmap(abundance = metacyc_daa_results_df %>% filter("function" %in% feature_with_p_0.05$feature) %>% column_to_rownames("function"), metadata = metadata_filtered, group = "Type")
