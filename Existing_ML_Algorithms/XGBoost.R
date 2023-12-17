# Required Programming Language
# R version 3.6.2 or above

# Required packages
required_packages <- c("caret", "dplyr", "e1071", "ggplot2", "class", 
                       "randomForest", "pROC", "kernlab", "xgboost")

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
input_path <- "./Datasets/" 

# Define an output path
output_path <- "./All_ML_Models_Results/" 

# Read Gene Expression Datasets
read_gene_expression_data <- function(filename,input_path) {
  path_filename <- paste0(input_path, filename)
  print(paste0("Reading CSV file called ", path_filename))
  return(read.csv(file = path_filename, header = TRUE, sep=",", check.names = FALSE))
}


# Split data into train, validation, and test
split_data <- function(Gene_Exp_Data, data_partition_seed, dataset_name) {
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
  
  
  # Extract dataset name only by removing .csv from the dataset file name
  dataset_name <- substr(dataset_name, 1, nchar(dataset_name) - 4)
  
  # Partitioning hold_out_test data further into Testset 1a, Testset 1b or  Testset 2a, Testset 2b only for Gene_Expression_Dataset_1 and Gene_Expression_Dataset_2
  if(dataset_name == "Gene_Expression_Dataset_1_GSE73072" | dataset_name == "Gene_Expression_Dataset_2_GSE68310") {
    set.seed(data_partition_seed)
    # Dividing test data set further into holdout test set (25%) and validation set (25%) using createDataPartition function of caret package
    index_test1 <- createDataPartition(y = hold_out_test$True_Class_Label,
                                       p = 0.50, list = FALSE)
    hold_out_test_a <- hold_out_test[index_test1, ]
    hold_out_test_b <- hold_out_test[-index_test1, ]
    
    # Testset all time points
    holdout_test_a <- total_data %>% filter(Super_Subject_ID %in% hold_out_test_a$Super_Subject_ID)
    holdout_test_b <- total_data %>% filter(Super_Subject_ID %in% hold_out_test_b$Super_Subject_ID)
    
    return(list(train_data = g_train_data, full_train_data = g_full_train_data, hold_out_test = hold_out_test, valid_data = g_valid_data, hold_out_test_a = hold_out_test_a, hold_out_test_b = hold_out_test_b, holdout_test = g_test_data, holdout_test_a = holdout_test_a, holdout_test_b = holdout_test_b))
    
  } else {
    return(list(train_data = g_train_data, full_train_data = g_full_train_data, hold_out_test = hold_out_test, valid_data = g_valid_data, holdout_test = g_test_data))
  }
  
}


# Train XGBoost model 
train_XGBoost_model <- function(train_data, Labels_train_data, ml_model_seed, max_depth_vals = 6, eta_vals = 0.3) {
  set.seed(ml_model_seed)
  Labels_train_data <- recode(Labels_train_data,'RVI'=1, 'Not RVI'=0)
  
  # training procedure
  trained_model <- xgboost(data = train_data, 
                           label = Labels_train_data, 
                           max.depth = max_depth_vals, 
                           eta = eta_vals, 
                           nrounds = 100, 
                           eval_metric = "error",
                           objective = "binary:logistic",
                           verbose = 1)
  
  # Print trained model
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
  print(paste0("Writing results to ", output_path))
  
  if(length(levels(predictions)) == 1 && length(levels(actual_labels)) == 1) {
    # Calculate accuracy for a single class
    accuracy <- sum(predictions == actual_labels) / length(actual_labels)
    write("Only one class present in both predictions and actual labels. Calculated accuracy for single class.", con)
    write(print(paste0("Accuracy: ", accuracy*100,"%")), con)
    
    print("Only one class present in both predictions and actual labels.")
    print(paste0("Calculated accuracy for single class. Accuracy: ", accuracy*100,"%"))
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


##################################################################################################
####################################### Main Program  ############################################
##################################################################################################

# Read datasets
dataset_names <- c("Gene_Expression_Dataset_1_GSE73072.csv", "Gene_Expression_Dataset_2_GSE68310.csv", "Gene_Expression_Dataset_3_GSE90732.csv", "Gene_Expression_Dataset_4_GSE61754.csv")

for(dataset_name in dataset_names) {
  Gene_Exp_Data <- read_gene_expression_data(dataset_name, input_path)
  
  # Split data
  splits <- split_data(Gene_Exp_Data, data_partition_seed, dataset_name)
  
  
  # Code for XGBoost model building, validation and evaluation
  # Model building using training data
  print("###########  Starting XGBoost learning using training data ########### ")
  trained_XGBoost_model <- train_XGBoost_model(as.matrix(splits$train_data[,-c(1:6)]),
                                               as.factor(splits$train_data$True_Class_Label), ml_model_seed)
  
  # Print XGBoost training results
  print("XGBoost training results:")
  print(trained_XGBoost_model)
  
  # best performance on train data (error)
  print(min(trained_XGBoost_model$evaluation_log$train_error))
  
  # best performance on train data (Accuracy)
  print(1-min(trained_XGBoost_model$evaluation_log$train_error))
  
  # Performing validation and hyper parameter selection using validation data
  validation_XGBoost_model <- train_XGBoost_model(as.matrix(splits$valid_data[,-c(1:6)]),
                                                  as.factor(splits$valid_data$True_Class_Label), ml_model_seed)
  
  # Print XGBoost validation results
  print("########### XGBoost validation results ###########")
  print(validation_XGBoost_model)
  
  # Selecting final model parameters
  final_max_depth <- validation_XGBoost_model$params$max_depth
  final_eta <- validation_XGBoost_model$params$eta
  
  # best performance on valid (error)
  print(min(validation_XGBoost_model$evaluation_log$train_error))
  
  # best performance on valid (Accuracy)
  print(1-min(validation_XGBoost_model$evaluation_log$train_error))
  
  
  # Print final parameter values
  print("Final value of parameter max depth of XGBoost:")
  print(final_max_depth)
  print("Final value of parameter eta of XGBoost:")
  print(final_eta)
  
  
  # Final model building using XGBoost
  final_XGBoost_trained_model <- train_XGBoost_model(as.matrix(splits$full_train_data[,-c(1:6)]),
                                                     as.factor(splits$full_train_data$True_Class_Label),
                                                     ml_model_seed, final_max_depth, final_eta)
  
  # Test XGBoost model
  print("########### Starting prediction for holdout testset ###########")
  XGBoost_predictions <- predict(final_XGBoost_trained_model, newdata = as.matrix(splits$hold_out_test[,-c(1:6)]))
  
  for (i in 1:length(XGBoost_predictions)) {
    if(XGBoost_predictions[i] < 0.5){
      XGBoost_predictions[i] <- 'Not RVI'
    }else{
      XGBoost_predictions[i] <- 'RVI'
    }
    
  }
  
  # Write results to a text file for full holdout testset
  write_confusion_to_txt(as.factor(XGBoost_predictions), as.factor(splits$hold_out_test$True_Class_Label), "XGBoost", dataset_name, output_path)
  
  # Result for Testset 1a or  Testset 2a only for Gene_Expression_Dataset_1 and Gene_Expression_Dataset_2
  if(dataset_name == "Gene_Expression_Dataset_1_GSE73072.csv" | dataset_name == "Gene_Expression_Dataset_2_GSE68310.csv") {
    
    # Test XGBoost model
    XGBoost_predictions <- predict(final_XGBoost_trained_model, newdata = as.matrix(splits$hold_out_test_a[,-c(1:6)]))
    
    for (i in 1:length(XGBoost_predictions)) {
      if(XGBoost_predictions[i] < 0.5){
        XGBoost_predictions[i] <- 'Not RVI'
      }else{
        XGBoost_predictions[i] <- 'RVI'
      }
    }
    
    data_name <- substr(dataset_name, 1, nchar(dataset_name) - 4)
    # Create the filename
    result_filename <- paste0(data_name, "_holdout_testset_a_results.txt")
    
    # Write results to a text file
    write_confusion_to_txt(as.factor(XGBoost_predictions), as.factor(splits$hold_out_test_a$True_Class_Label), "XGBoost", result_filename, output_path)
    
  }
  
  # Result for  Testset 1b or Testset 2b only for Gene_Expression_Dataset_1 and Gene_Expression_Dataset_2
  if(dataset_name == "Gene_Expression_Dataset_1_GSE73072.csv" | dataset_name == "Gene_Expression_Dataset_2_GSE68310.csv") {
    
    # Test XGBoost model
    XGBoost_predictions <- predict(final_XGBoost_trained_model, newdata = as.matrix(splits$hold_out_test_b[,-c(1:6)]))
    
    for (i in 1:length(XGBoost_predictions)) {
      if(XGBoost_predictions[i] < 0.5){
        XGBoost_predictions[i] <- 'Not RVI'
      }else{
        XGBoost_predictions[i] <- 'RVI'
      }
    }
    
    data_name <- substr(dataset_name, 1, nchar(dataset_name) - 4)
    # Create the filename
    result_filename <- paste0(data_name, "_holdout_testset_b_results.txt")
    
    # Write results to a text file
    write_confusion_to_txt(as.factor(XGBoost_predictions), as.factor(splits$hold_out_test_b$True_Class_Label), "XGBoost", result_filename, output_path)
    
  } # End if condition
  
}

