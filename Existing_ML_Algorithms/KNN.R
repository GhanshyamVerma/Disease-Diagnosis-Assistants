# Required Programming Language
# R version 3.6.2 or above

# Required packages
required_packages <- c("caret", "dplyr", "e1071", "ggplot2", "class", "pROC")

# Check if required r libraries are already installed or not
new_packages <- required_packages[!(required_packages %in% installed.packages()[,"Package"])]

# If required r packages or libraries are not already installed, Install them.
if(length(new_packages)) install.packages(new_packages)

# load packages
library(caret) 
library(dplyr) 
library(e1071) 
library(ggplot2)
library(class) 
library(randomForest) 
library(pROC) 
library(kernlab) 
library(xgboost)

# Provide seed values for data partition and ml model building

data_partition_seed <- 7
ml_model_seed <- 1234


# Define an input path
input_path <- "/Datasets/Gene_Expression/"

# Define an output path
output_path <- "/All_ML_Models_Results/"

# Read Gene Expression Datasets
read_gene_expression_data <- function(filename,input_path) {
  path_filename <- paste0(input_path, filename)
  return(read.csv(file = path_filename, header = TRUE, sep=",", check.names = FALSE))
}

# Split data into train, validation, and test
split_data <- function(Gene_Exp_Data, data_partition_seed) {
  total_data <- Gene_Exp_Data
  Data_target <- total_data %>% filter(Time > 0)
  
  # Split into train and test
  set.seed(data_partition_seed)
  index_test <- createDataPartition(y = Data_target$True_Class_Label,
                                    p = 0.50, list = FALSE)
  test_data <- Data_target[index_test, ]
  train_data <- Data_target[-index_test, ]
  
  dim(train_data)
  dim(test_data)
  
  # Split test data into holdout test and validation
  set.seed(data_partition_seed)
  index_test <- createDataPartition(y = test_data$True_Class_Label,
                                    p = 0.50, list = FALSE)
  hold_out_test <- test_data[index_test, ]
  valid_data <- test_data[-index_test, ]
  dim(hold_out_test)
  dim(valid_data)
  
  # Full train data for final model building
  full_train_data <- rbind(train_data,valid_data)
  
  # train test all time points
  g_train_data <- total_data %>% filter(Super_Subject_ID %in% train_data$Super_Subject_ID)
  g_valid_data <- total_data %>% filter(Super_Subject_ID %in% valid_data$Super_Subject_ID)
  g_test_data <- total_data %>% filter(Super_Subject_ID %in% hold_out_test$Super_Subject_ID)
  g_full_train_data <- total_data %>% filter(Super_Subject_ID %in% full_train_data$Super_Subject_ID)
  
  # Display the data
  g_test_data[c(1:6),c(1:7)] # show first 6 rows
  
  return(list(train_data = g_train_data, full_train_data = g_full_train_data, hold_out_test = hold_out_test, valid_data = g_valid_data))
}

# Train KNN model
train_knn_model <- function(train_data, Labels_train_data, ml_model_seed) {
  set.seed(ml_model_seed)
  metric <- "Accuracy"
  grid <- expand.grid(k = c(1:30))
  trained_model <- train(x= train_data,
                         y = Labels_train_data,
                         method = "knn",
                         metric = metric,
                         tuneGrid = grid)
  return(trained_model)
}



# Final Training function for KNN using best parameters
final_knn_training_function <- function(train_data, Labels_train_data, ml_model_seed, final_k) {
  set.seed(ml_model_seed)
  metric <- "Accuracy"
  grid <- expand.grid(k = final_k)
  trained_model <- train(x= train_data,
                         y = Labels_train_data,
                         method = "knn",
                         metric = metric,
                         tuneGrid = grid)
  return(trained_model)
}


# Write confusion matrix results to TXT
write_confusion_to_txt <- function(predictions, actual_labels, model_name, dataset_name, output_path) {
  # Ensure that both predictions and actual_labels are factors
  predictions <- factor(predictions)
  actual_labels <- factor(actual_labels)
  
  dataset_name <- substr(dataset_name, 1, nchar(dataset_name) - 4)
  
  # Create the filename
  filename <- paste0(output_path, model_name, "_", dataset_name, "_confusion_matrix.txt")
  
  # Open a connection for writing
  con <- file(filename, "w")
  
  if(length(levels(predictions)) == 1 && length(levels(actual_labels)) == 1) {
    # Calculate accuracy for a single class
    accuracy <- sum(predictions == actual_labels) / length(actual_labels)
    write("Only one class present in both predictions and actual labels. Calculated accuracy for single class.", con)
    write("Accuracy:", con)
    write(accuracy, con)
    
    print("Only one class present in both predictions and actual labels.")
    print(paste0("Calculated accuracy for single class: ", accuracy))
    print(accuracy)
  } else {
    # Ensure both have the same levels
    levels_union <- union(levels(predictions), levels(actual_labels))
    predictions <- factor(predictions, levels = levels_union)
    actual_labels <- factor(actual_labels, levels = levels_union)
    
    # Compute the confusion matrix
    conf_matrix_result <- confusionMatrix(predictions, actual_labels)
    
    # Extract the matrix and statistics
    conf_matrix_table <- conf_matrix_result$table
    overall_stats <- conf_matrix_result$overall
    class_stats <- conf_matrix_result$byClass
    
    # Print and write confusion matrix
    print(paste0("Print confusion matrix for hold-out test set for ", model_name, " ", "for ", dataset_name, ":"))
    print(conf_matrix_result)
    
    write("Confusion Matrix:", con)
    write.table(conf_matrix_table, con, sep = "\t")
    
    write("\nOverall Statistics:", con)
    write.table(as.data.frame(overall_stats), con, sep = "\t", row.names = TRUE)
    
    write("\nClass Statistics:", con)
    write.table(as.data.frame(t(class_stats)), con, sep = "\t", row.names = TRUE)
  }
  
  # Close the connection
  close(con)
}



# Main function
main_function <- function() {
  # Read datasets
  dataset_names <- c("Gene_Expression_Dataset_1_GSE73072.csv", "Gene_Expression_Dataset_2_GSE68310.csv", "Gene_Expression_Dataset_3_GSE90732.csv", "Gene_Expression_Dataset_4_GSE61754.csv")
  for(dataset_name in dataset_names) {
    Gene_Exp_Data <- read_gene_expression_data(dataset_name, input_path)
    
    # Split data
    splits <- split_data(Gene_Exp_Data, data_partition_seed)
    
    # Model building using training data
    trained_knn_model <- train_knn_model(as.matrix(splits$train_data[,-c(1:6)]),
                                         as.factor(splits$train_data$True_Class_Label), ml_model_seed)

    # Print KNN training results
    print("KNN training results:")
    print(trained_knn_model)

    # Performing validation and hyper parameter selection using validation data
    validation_knn_model <- train_knn_model(as.matrix(splits$valid_data[,-c(1:6)]),
                                            as.factor(splits$valid_data$True_Class_Label), ml_model_seed)


    # Print KNN validation results
    print("KNN validation results:")
    print(validation_knn_model)

    # Selecting final model parameters
    final_k <- validation_knn_model$finalModel$tuneValue[1]

    # Print final parameter values
    print("Final value of k:")
    print(final_k)

    # Final model building using KNN
    final_knn_trained_model <- final_knn_training_function(as.matrix(splits$full_train_data[,-c(1:6)]),
                                                           as.factor(splits$full_train_data$True_Class_Label),
                                                           ml_model_seed, final_k)


    # Test KNN
    knn_predictions <- predict(final_knn_trained_model, newdata = as.matrix(splits$hold_out_test[,-c(1:6)]))

    # Write results to a text file
    write_confusion_to_txt(knn_predictions, as.factor(splits$hold_out_test$True_Class_Label), "KNN", dataset_name, output_path)
  
  }
}

main_function()
