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
# Example: Convert all "increase"/"decrease" to numeric
df <- df %>%
  mutate(across(where(is.character), ~ ifelse(. == "increase", 1,
                                              ifelse(. == "decrease", -1, NA))))

colnames(df) <- make.names(colnames(df))


set.seed(123)
train_index <- createDataPartition(df$Type, p = 0.7, list = FALSE)
train_data <- df[train_index, ]
test_data  <- df[-train_index, ]

col_na_fraction <- colMeans(is.na(train_data))
sort(col_na_fraction, decreasing = TRUE)
# Drop columns where more than 50% of values are NA
train_data <- train_data[, col_na_fraction < 0.5]
test_data  <- test_data[, colnames(train_data)]  # match columns


preProc <- preProcess(train_data, method = "medianImpute")
train_data <- predict(preProc, train_data)
test_data <- predict(preProc, test_data)

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
plot_data <- df %>%
  select(all_of(c("Type", top_features))) %>%
  pivot_longer(-Type, names_to = "Feature", values_to = "Abundance")

p <- ggplot(plot_data, aes(x = Type, y = Abundance, fill = Type)) +
  geom_violin(trim = FALSE, alpha = 0.7) +
  geom_boxplot(width = 0.1, outlier.shape = NA, color = "black") +
  facet_wrap(~ Feature, scales = "free_y", ncol = 2) +
  theme_minimal(base_size = 12) +
  labs(title = "Top Features by Class",
       y = "Abundance", x = "Diabetes Type") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none")

#ggsave("/data/scratch/kvalem/projects/2024/Effenberger-Diabetes/02-scripts/figures/v02/top10_features.png", plot = p, width = 10, height = 10, dpi = 300)
#ggsave("/data/scratch/kvalem/projects/2024/Effenberger-Diabetes/02-scripts/figures/v02/top10_features.svg", plot = p, width = 10, height = 10, dpi = 300)



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


