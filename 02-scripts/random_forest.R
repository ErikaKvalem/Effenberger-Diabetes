library(randomForest)
library(readr)
library(dplyr)
library(caret)
library(randomForest)
library(pROC)
library(tidyverse)
library(corrplot)
library(randomForest)
library(pROC)
library(caret)
library(dplyr)
library(devtools)
install_github("jokergoo/ComplexHeatmap")
library(ComplexHeatmap)
library(ComplexHeatmap)
library(circlize)


df <- read_csv("metadata_abundance.csv")
df <- df %>%
  filter(Type %in% c("pankreopriver Diabetes", "Diabetes mellitus Typ1"))
df <- df %>%
  filter(Type %in% c("pankreopriver Diabetes", "Diabetes mellitus Typ1")) %>%
  mutate(Type = recode(Type,
                       "pankreopriver Diabetes" = "PDM",
                       "Diabetes mellitus Typ1" = "DM"))

# Convert Type to factor for classification
df$Type <- as.factor(df$Type)
# Remove columns starting with "..."
df <- df[, !grepl("^\\.\\.\\.", names(df))]
df <- df %>% select(-sample_information)  
df <- df %>% select(-CA_diff)  
df <- df %>% select(-KHK_diff)  
#df <- df %>% select(-`Pat ID`)  
# Example: Convert all "increase"/"decrease" to numeric


#df$CA1 <- ifelse(df$CA1 == "ja", 1,
#                 ifelse(df$CA1 == "nein", 0, NA))
#df$CA2 <- ifelse(df$CA2 == "ja", 1,
#                 ifelse(df$CA2 == "nein", 0, NA))


colnames(df) <- make.names(colnames(df))


set.seed(123)
train_index <- createDataPartition(df$Type, p = 0.7, list = FALSE)
train_data <- df[train_index, ]
test_data  <- df[-train_index, ]

colnames(train_data) <- make.names(colnames(train_data))
colnames(test_data)  <- make.names(colnames(test_data))  # Also update test data


col_na_fraction <- colMeans(is.na(train_data))
sort(col_na_fraction, decreasing = TRUE)
# Drop columns where more than 50% of values are NA
train_data <- train_data[, col_na_fraction < 0.5]
test_data  <- test_data[, colnames(train_data)]  # match columns


preProc <- preProcess(train_data, method = "medianImpute")
train_data <- predict(preProc, train_data)
test_data <- predict(preProc, test_data)

train_data <- train_data %>% select(-id)
test_data  <- test_data %>% select(-id)

tune_grid <- expand.grid(mtry = c(2, 5, 10, 20, 50))
fit_control <- trainControl(method = "LOOCV", number = 5)
rf_model <- train(Type ~ ., data = train_data, method = "rf",
                  trControl = fit_control,
                  tuneGrid = tune_grid, importance = TRUE)

predict<- predict(rf_model, test_data)
confusionMatrix(predict, test_data$Type)

####################################### 
importance_df <- varImp(rf_model)$importance

# If it has multiple columns (e.g., one per class), create an overall score
if (ncol(importance_df) > 1) {
  importance_df$Overall <- rowSums(importance_df)
}

# Now sort by importance
importance_df <- importance_df[order(-importance_df$Overall), , drop = FALSE]

# Assign cleaned row names as a new Feature column
importance_df$Feature <- make.names(rownames(importance_df))

# Preview
head(importance_df, 20)


top_features <- importance_df$Feature[1:10]
top_features_clean <- make.names(top_features)


# Subset training and test data to include only top features + label
train_data_selected <- train_data[, c("Type", top_features)]
test_data_selected  <- test_data[, c("Type", top_features)]

fit_control <- trainControl(method = "LOOCV", number = 5)


rf_selected <- train(Type ~ ., data = train_data_selected, method = "rf",
                     trControl = fit_control, importance = TRUE)

# Evaluate
print(rf_selected)
predict_selected <- predict(rf_selected, test_data_selected)
confusionMatrix(predict_selected, test_data_selected$Type)

top10_df <- importance_df[1:10, c("Feature", "Overall")]

################# PLOT 



# Create long-format data for plotting
df$CA1 <- ifelse(df$CA1 == "ja", 1,
                 ifelse(df$CA1 == "nein", 0, NA))
df$CA2 <- ifelse(df$CA2 == "ja", 1,
                 ifelse(df$CA2 == "nein", 0, NA))
top_features

plot_data <- df %>%
  select(all_of(c("Type", top_features))) %>%
  pivot_longer(-Type, names_to = "Feature", values_to = "Abundance")



################### PDM vs DM 
microbial_features <- c("Lachnospiraceae.UCG.009", "Parabacteroides", "Oscillospira", 
                        "Family.XIII.UCG.001", "Pankreatektomie_encoded", 
                        "Faecalibacterium", "Limosilactobacillus")

clinical_features <- c("CA2", "age", "CA1")

# Add feature_type column
plot_data <- plot_data %>%
  mutate(feature_type = case_when(
    Feature %in% microbial_features ~ "Microbial taxa",
    Feature %in% clinical_features ~ "Clinical metadata",
    TRUE ~ "Other"
  ))

p <- ggplot(plot_data, aes(x = Type, y = Abundance, fill = feature_type)) +
  geom_violin(trim = TRUE, alpha = 0.7) +
  geom_boxplot(width = 0.1, outlier.shape = NA, color = "black") +
  facet_wrap(~ Feature, scales = "free_y", ncol = 2) +
#  scale_fill_manual(values = c("DM" = "#E4A165", "PDM" = "#75B375")) +
  theme_minimal(base_size = 12) +
  labs(title = "Top Features",
       y = "Feature vallue", x = "Diabetes Type") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none")


p
ggsave("/data/scratch/kvalem/projects/2024/Effenberger-Diabetes/02-scripts/figures/v02/top10_features.png", plot = p, width = 10, height = 10, dpi = 300)
ggsave("/data/scratch/kvalem/projects/2024/Effenberger-Diabetes/02-scripts/figures/v02/top10_features.svg", plot = p, width = 10, height = 10, dpi = 300)

ggplot(plot_data, aes(x = Type, y = Abundance, fill = Type)) 

##########################################
library(caret)
library(ggplot2)

# Generate confusion matrix

cm <- confusionMatrix(predict_selected, test_data_selected$Type)

# Convert to data frame for ggplot
cm_df <- as.data.frame(cm$table)
colnames(cm_df) <- c("Reference", "Prediction", "Freq")


p_CMatrix <- ggplot(cm_df, aes(x = Reference, y = Prediction, fill = Freq)) +
  geom_tile(color = "white") +
  geom_text(aes(label = Freq), size = 5) +
  scale_fill_gradient(low = "lightblue", high = "steelblue") +
  labs(title = "Confusion   \n Matrix Accuracy : 0.8182  \n Sensitivity : 0.6000 Specificity : 1.0000   ", x = "Actual", y = "Predicted") +
  theme_minimal()

ggsave("/data/scratch/kvalem/projects/2024/Effenberger-Diabetes/02-scripts/figures/v02/confusion_matrix.png", plot = p_CMatrix, width = 3.5, height = 3, dpi = 300)
ggsave("/data/scratch/kvalem/projects/2024/Effenberger-Diabetes/02-scripts/figures/v02/confusion_matrix.svg", plot = p_CMatrix, width = 3.5, height = 3, dpi = 300)
#################################################

library(pROC)

# Predict probabilities from random forest model
pred_probs <- predict(rf_selected, test_data_selected, type = "prob")

# Ensure correct order of levels (e.g., "PDM" = 0, "DM" = 1)
test_data_selected$Type <- factor(test_data_selected$Type, levels = c("PDM", "DM"))

# Compute ROC curve for the positive class ("DM")
roc_obj <- roc(test_data_selected$Type, pred_probs[,"DM"])

library(pROC)
library(ggplot2)

# Compute ROC
roc_obj <- roc(test_data_selected$Type, pred_probs[,"DM"])

# Use ggroc() to create a ggplot object
plot_auc <- ggroc(roc_obj, color = "darkorange", size = 1.2) +
  ggtitle("ROC Curve - Accuracy ") +
  theme_minimal() 
  annotate("text", x = 0.6, y = 0.4, label = paste("AUC =", round(auc(roc_obj), 3)), size = 2)
  
plot_auc

# Save to PNG and SVG
ggsave("/data/scratch/kvalem/projects/2024/Effenberger-Diabetes/02-scripts/figures/v02/roc_auc.png",
       plot = plot_auc, width = 2.5, height = 2, dpi = 300)

ggsave("/data/scratch/kvalem/projects/2024/Effenberger-Diabetes/02-scripts/figures/v02/roc_auc.svg",
       plot = plot_auc, width = 2.5, height = 2, dpi = 300)

###############################

library(ggplot2)
library(caret)

# Train confusion matrix
train_pred <- predict(rf_selected, train_data_selected)
cm_train <- as.data.frame(confusionMatrix(train_pred, train_data_selected$Type)$table)

ggplot(cm_train, aes(Prediction, Reference, fill = Freq)) +
  geom_tile() +
  geom_text(aes(label = Freq), color = "white", size = 4) +
  scale_fill_gradient(low = "lightblue", high = "steelblue") +
  labs(title = "Confusion Matrix - Train", x = "Predicted", y = "Actual") +
  theme_minimal()
##############################

library(pROC)

# Probabilities
train_probs <- predict(rf_selected, train_data_selected, type = "prob")
test_probs  <- predict(rf_selected, test_data_selected, type = "prob")

# ROC curves
roc_train <- roc(train_data_selected$Type, train_probs[, "DM"])
roc_test  <- roc(test_data_selected$Type, test_probs[, "DM"])

# Plot overlay
plot(roc_train, col = "blue", lwd = 2, main = "ROC Curves")
plot(roc_test, col = "red", lwd = 2, add = TRUE)
legend("bottomright", legend = c("Train", "Test"),
       col = c("blue", "red"), lwd = 2)

varImpPlot(rf_selected$finalModel, type = 1, main = "Top Variables (Train)")


#################### venn diagram pdm vs dm vs k 

pdm_k_microbial = c( "X.Ruminococcus..gauvreauii.group", "Faecalibacterium","CAG.56",  "Lachnospiraceae.NK4A136.group",    "Lachnospiraceae.NC2004.group", "Anaerostipes" )
pdm_k_clinical =c("HbA1C..DCCT.NGSP.1","HbA1C..DCCT.NGSP.2","Glukose1" ,"Glukose2")

dm_k_microbial =  c(  "Bilophila"  , "X.Eubacterium..nodatum.group", "Tyzzerella", "Anaerostipes", "Monoglobus", "UCG.002", "CAG.56")
dm_k_clinical =  c(  "HbA1C..DCCT.NGSP.1","HbA1C..DCCT.NGSP.2","Glukose1" ,"Glukose2" )

pdm_dm_microbial = c("Lachnospiraceae.UCG.009", "Parabacteroides", "Oscillospira", 
                        "Family.XIII.UCG.001", 
                        "Faecalibacterium", "Limosilactobacillus")

pdm_dm_clinical = c("CA2", "age", "CA1","Pankreatektomie_encoded")



shared_microbial <- intersect(pdm_k_microbial, dm_k_microbial)
unique_pdm_dm_microbial <- setdiff(pdm_dm_microbial, union(pdm_k_microbial, dm_k_microbial))


library(VennDiagram)

venn.plot <- draw.triple.venn(
  area1 = length(pdm_k_microbial),
  area2 = length(dm_k_microbial),
  area3 = length(pdm_dm_microbial),
  n12 = length(intersect(pdm_k_microbial, dm_k_microbial)),
  n23 = length(intersect(dm_k_microbial, pdm_dm_microbial)),
  n13 = length(intersect(pdm_k_microbial, pdm_dm_microbial)),
  n123 = length(Reduce(intersect, list(pdm_k_microbial, dm_k_microbial, pdm_dm_microbial))),
  category = c("PDM vs K", "DM vs K", "PDM vs DM"),
  fill = c("skyblue", "orange", "pink")
)

ggsave("/data/scratch/kvalem/projects/2024/Effenberger-Diabetes/02-scripts/figures/v02/venn_microbial.png",
       plot = venn.plot , width = 2.5, height = 2, dpi = 300)

ggsave("/data/scratch/kvalem/projects/2024/Effenberger-Diabetes/02-scripts/figures/v02/venn_microbial.svg",
       plot = venn.plot , width = 2.5, height = 2, dpi = 300)


shared_microbial <- intersect(pdm_k_clinical, dm_k_clinical)
unique_pdm_dm_microbial <- setdiff(pdm_dm_clinical, union(pdm_k_clinical, dm_k_clinical))


library(VennDiagram)

venn.plot <- draw.triple.venn(
  area1 = length(pdm_k_clinical),
  area2 = length(dm_k_clinical),
  area3 = length(pdm_dm_clinical),
  n12 = length(intersect(pdm_k_clinical, dm_k_clinical)),
  n23 = length(intersect(dm_k_clinical, pdm_dm_clinical)),
  n13 = length(intersect(pdm_k_clinical, pdm_dm_clinical)),
  n123 = length(Reduce(intersect, list(pdm_k_clinical, dm_k_clinical, pdm_dm_clinical))),
  category = c("PDM vs K", "DM vs K", "PDM vs DM"),
  fill = c("skyblue", "orange", "pink")
)

ggsave("/data/scratch/kvalem/projects/2024/Effenberger-Diabetes/02-scripts/figures/v02/venn_clinical.png",
       plot = venn.plot , width = 2.5, height = 2, dpi = 300)

ggsave("/data/scratch/kvalem/projects/2024/Effenberger-Diabetes/02-scripts/figures/v02/venn_clinical.svg",
       plot = venn.plot , width = 2.5, height = 2, dpi = 300)

###########################


select_features <-unique(c ( "X.Ruminococcus..gauvreauii.group", "Faecalibacterium","CAG.56",  "Lachnospiraceae.NK4A136.group",    "Lachnospiraceae.NC2004.group", "Anaerostipes","Bilophila"  , "X.Eubacterium..nodatum.group", "Tyzzerella", "Anaerostipes", "Monoglobus", "UCG.002", "CAG.56","Lachnospiraceae.UCG.009", "Parabacteroides", "Oscillospira", 
    "Family.XIII.UCG.001",  
    "Faecalibacterium", "Limosilactobacillus"))
df_matrix <- as.matrix(df[, select_features])
rownames(df_matrix) <- df$id
df_matrix_norm <- scale(df_matrix)


Heatmap(df_matrix_norm)


###########################

df <- read_csv("metadata_abundance.csv")
#df <- df %>%
#  filter(Type %in% c("pankreopriver Diabetes", "Diabetes mellitus Typ1"))
#df <- df %>%
#  filter(Type %in% c("pankreopriver Diabetes", "Diabetes mellitus Typ1")) %>%
df <- df %>% mutate(Type = recode(Type,
                      "pankreopriver Diabetes" = "PDM",
                       "Diabetes mellitus Typ1" = "DM", "Kontrolle"="K"))

# Convert Type to factor for classification
df$Type <- as.factor(df$Type)
# Remove columns starting with "..."
df <- df[, !grepl("^\\.\\.\\.", names(df))]
df <- df %>% select(-sample_information)  
df <- df %>% select(-CA_diff)  
df <- df %>% select(-KHK_diff)  

colnames(df) <- make.names(colnames(df))


# Load necessary libraries


ordered_indices <- order(df$Type)
df$Type <- factor(df$Type, levels = c("K", "PDM", "DM"))
ordered_indices <- order(df$Type) 
select_features <-unique(c ( "X.Ruminococcus..gauvreauii.group", "Faecalibacterium","CAG.56",  "Lachnospiraceae.NK4A136.group",    "Lachnospiraceae.NC2004.group", "Anaerostipes","Bilophila"  , "X.Eubacterium..nodatum.group", "Tyzzerella", "Anaerostipes", "Monoglobus", "UCG.002", "CAG.56","Lachnospiraceae.UCG.009", "Parabacteroides", "Oscillospira", 
                             "Family.XIII.UCG.001",  
                             "Faecalibacterium", "Limosilactobacillus"))
select_features_clinical = c("CA2", "age", "CA1","Pankreatektomie_encoded", "HbA1C..DCCT.NGSP.1","HbA1C..DCCT.NGSP.2","Glukose1" ,"Glukose2")


df_matrix <- as.matrix(df[, select_features])
rownames(df_matrix) <- df$id
df_matrix_norm <- scale(df_matrix)

df_matrix <- log10(as.matrix(df[, select_features]) + 1e-6)


df_matrix_ordered <- df_matrix_norm[ordered_indices, ]

# Update Type accordingly
ordered_type <- df$Type[ordered_indices]

type_colors <- c("DM" = "#E1812C", "PDM" = "#3A923A", "K"="#3274A1")
row_annot <- HeatmapAnnotation(
  Type = ordered_type,
  col = list(Type = type_colors),
  which = "row"
)
df_matrix_transposed <- t(df_matrix_ordered)
col_annot <- HeatmapAnnotation(
  Type = ordered_type,
  col = list(Type = type_colors),
  which = "column"
)
h <- Heatmap(
  df_matrix_transposed,
  name = "Z-score",
  show_row_names = TRUE,
  show_column_names = FALSE,
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  top_annotation = col_annot,
  row_names_gp = gpar(fontsize = 12),  # optional: adjust label size
  heatmap_legend_param = list(direction = "vertical")  # optional: vertical legend
)
h

png("/data/scratch/kvalem/projects/2024/Effenberger-Diabetes/02-scripts/figures/v02/heatmap.png", 
    width = 2000, height = 2000, res = 300)
draw(h)
dev.off()

png("/data/scratch/kvalem/projects/2024/Effenberger-Diabetes/02-scripts/figures/v02/heatmap.svg", 
    width = 2000, height = 2000, res = 300)
draw(h)
dev.off()

################################################################# CLINICAL
# I need to remove Pankreatektomie_encoded because it is encoded, i will plot it in another way
select_features_clinical = c( "age","CA2", "CA1", "HbA1C..DCCT.NGSP.1","HbA1C..DCCT.NGSP.2","Glukose1" ,"Glukose2")


df_matrix <- as.matrix(df[, select_features_clinical])
rownames(df_matrix) <- df$id
df_matrix_norm <- scale(df_matrix)

df_matrix <- log10(as.matrix(df[, select_features_clinical]) + 1e-6)


df_matrix_ordered <- df_matrix_norm[ordered_indices, ]

# Update Type accordingly
ordered_type <- df$Type[ordered_indices]

type_colors <- c("DM" = "#E1812C", "PDM" = "#3A923A", "K"="#3274A1")
row_annot <- HeatmapAnnotation(
  Type = ordered_type,
  col = list(Type = type_colors),
  which = "row"
)
df_matrix_transposed <- t(df_matrix_ordered)
col_annot <- HeatmapAnnotation(
  Type = ordered_type,
  col = list(Type = type_colors),
  which = "column"
)
h <- Heatmap(
  df_matrix_transposed,
  name = "Z-score",
  show_row_names = TRUE,
  show_column_names = FALSE,
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  top_annotation = col_annot,
  row_names_gp = gpar(fontsize = 12),  # optional: adjust label size
  heatmap_legend_param = list(direction = "vertical")  # optional: vertical legend
)
h

png("/data/scratch/kvalem/projects/2024/Effenberger-Diabetes/02-scripts/figures/v02/heatmap_clinical.png", 
    width = 2000, height = 2000, res = 300)
draw(h)
dev.off()

png("/data/scratch/kvalem/projects/2024/Effenberger-Diabetes/02-scripts/figures/v02/heatmap_clinical.svg", 
    width = 2000, height = 2000, res = 300)
draw(h)
dev.off()
############################## correlation 

# Extract and log-transform microbial data
microbe_data <- log10(df[, select_features] + 1e-6)

# Extract clinical data
clinical_data <- df[, select_features_clinical]

# Make sure row counts match
stopifnot(nrow(microbe_data) == nrow(clinical_data))
# Compute Spearman correlation
cor_matrix <- cor(microbe_data, clinical_data, method = "spearman", use = "pairwise.complete.obs")
#col_fun <- colorRamp2(c(-0.5, 0, 0.5), c("blue", "white", "red"))
cor_matrix <- cor(microbe_data, clinical_data, method = "spearman", use = "pairwise.complete.obs")

Heatmap(cor_matrix,
        name = "Spearman\nCorrelation",
 #       col = col_fun,
        column_names_rot = 45,
        cluster_rows = FALSE,
        cluster_columns = FALSE,
        row_names_gp = gpar(fontsize = 10),
        column_names_gp = gpar(fontsize = 10))

my_palette <- colorRampPalette(c("blue", "white", "red"))(200)

# Optional: clip correlation values to [-0.5, 0.5] for better contrast
cor_matrix_clipped <- cor_matrix
cor_matrix_clipped[cor_matrix_clipped > 0.5] <- 0.5
cor_matrix_clipped[cor_matrix_clipped < -0.5] <- -0.5

# Define your color scale and breakpoints
my_palette <- colorRampPalette(c("blue", "white", "red"))(200)
my_breaks <- seq(-0.5, 0.5, length.out = 201)  # match palette length

# Plot without text, with custom legend
corrplot(cor_matrix_clipped,
         method = "color",
         type = "lower",
         col = my_palette,
         tl.col = "black",
         tl.srt = 45,
         addCoef.col = NULL,     
         cl.breaks = c(-0.5, -0.25, 0, 0.25, 0.5),  
         cl.length = 5,          # make sure it matches number of breaks
         is.corr = FALSE,        # to respect your custom scale
         mar = c(0,0,1,0),)



par(mar = c(1, 1, 1, 7))  # bottom, left, top, right margins
svg("/data/scratch/kvalem/projects/2024/Effenberger-Diabetes/02-scripts/figures/v02/correlation.svg", width = 8, height = 8)
# Plot the correlation matrix
corrplot(cor_matrix_clipped,
         method = "color",
         type = "lower",
         col = my_palette,
         tl.col = "black",
         tl.srt = 45,
         addCoef.col = NULL,     
         cl.breaks = c(-0.5, -0.25, 0, 0.25, 0.5),  
         cl.length = 5,          # make sure it matches number of breaks
         is.corr = FALSE,        # to respect your custom scale
         mar = c(0,0,1,0),)


# Add legend label on the right side (side = 4)
mtext("Spearman correlation", side = 2, line = 0.5, cex = 1)
dev.off()

