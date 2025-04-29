library(randomForest)
library(readr)
library(dplyr)
library(caret)
library(randomForest)



df <- read_csv("metadata_abundance.csv")
df <- df %>%
  filter(Type %in% c("pankreopriver Diabetes", "Diabetes mellitus Typ1"))
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
rf_model <- train(Type ~ ., data = train_data, method = "rf",
                  trControl = fit_control,
                  tuneGrid = tune_grid, importance = TRUE)



####################################### 



importance_df <- varImp(rf_model)$importance

# If it has multiple columns (e.g., one per class), create an overall score
if (ncol(importance_df) > 1) {
  importance_df$Overall <- rowSums(importance_df)
}

# Now sort by importance
importance_df <- importance_df[order(-importance_df$Overall), , drop = FALSE]

# Add feature names
importance_df$Feature <- rownames(importance_df)

# Preview
head(importance_df, 20)


top_features <- importance_df$Feature[1:20]

# Subset training and test data to include only top features + label
train_data_selected <- train_data[, c("Type", top_features)]
test_data_selected  <- test_data[, c("Type", top_features)]


# Match feature names safely using exact matching
selected_cols <- match(top_features, names(train_data))

# Add Type column manually
selected_cols <- c(match("Type", names(train_data)), selected_cols)


