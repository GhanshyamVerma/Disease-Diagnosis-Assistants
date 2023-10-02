# Required Programming Language
# R version 3.6.2 or above

# Required packages
# caret 6.0.93
# dplyr 1.1.0
# e1071 1.7.11
# xgboost 1.7.3.1

# load packages
library(caret) # for various machine learning functions
library(dplyr) # for data manipulation
library(e1071) # for various functions like confusion matrix
library(xgboost) # for building XGboost model


############################### Input from user ##########################################

# read Gene expression data collected at two required time-points for all those subjects for which disease diagnosis need to be performed
# required two time-points (t_1 = 0 hours i.e. healthy state and at time-point t_D i.e. diseased state or the time-point at which disease diagnosis has been requested)
# if you want to use another another Gene Expression Dataset, just read another Gene Expression Dataset by giving path of that dataset in below statement.
Gene_Exp_Data <- read.csv(file = "Gene_Expression_Dataset_4_GSE61754.csv", header = TRUE, sep=",", check.names = FALSE)

# Display the data
Gene_Exp_Data[c(1:7),c(1:7)] # show first 7 rows


# go inside the train test split
total_data <- Gene_Exp_Data

Data_target <- total_data %>% filter(Time > 0) 

## First dividing data into training and test set (single time point - target time point)
set.seed(7)
# Dividing data set into train (50%) and test (50%) using createDataPartition function of caret package
index_test <- createDataPartition(y = Data_target$True_Class_Label,
                                  p = 0.50, list = FALSE)
test_data <- Data_target[index_test, ]
train_data <- Data_target[-index_test, ]


dim(train_data)
dim(test_data)

set.seed(7)
# Dividing test data set further into holdout test set (25%) and validation set (25%) using createDataPartition function of caret package
index_test <- createDataPartition(y = test_data$True_Class_Label,
                                  p = 0.50, list = FALSE)
hold_out_test <- test_data[index_test, ]
valid_data <- test_data[-index_test, ]

dim(hold_out_test)
dim(valid_data)

# full train data for final model building
full_train_data <- rbind(train_data,valid_data)

# train test all time points
g_train_data <- total_data %>% filter(Super_Subject_ID %in% train_data$Super_Subject_ID)
g_valid_data <- total_data %>% filter(Super_Subject_ID %in% valid_data$Super_Subject_ID)
g_test_data <- total_data %>% filter(Super_Subject_ID %in% hold_out_test$Super_Subject_ID)
g_full_train_data <- total_data %>% filter(Super_Subject_ID %in% full_train_data$Super_Subject_ID)



# Display the dimensions of training and hold-out test set(rows columns)
(dim(g_train_data))
(dim(g_valid_data))
(dim(g_test_data))
(dim(g_full_train_data)) # for final model building

# Display the data
g_test_data[c(1:6),c(1:7)] # show first 6 rows


# Train and test data
g_train_data_all_genes <- g_train_data
g_valid_data_all_genes <- g_valid_data
g_test_data_all_genes <- g_test_data
g_full_train_data_all_genes <- g_full_train_data

# Converting Label vector into factor as per the requirement of train function of caret package
Labels_g_train_data_all_genes <- as.factor(g_train_data$True_Class_Label)
Labels_g_valid_data_all_genes <- as.factor(g_valid_data$True_Class_Label)
Labels_g_test_data_all_genes <- as.factor(g_test_data$True_Class_Label)
Labels_g_full_train_data_all_genes <- as.factor(g_full_train_data_all_genes$True_Class_Label)

# Converting training and test data into matrix as per the requirement of train function of caret package
g_train_data_all_genes <- as.matrix(g_train_data_all_genes[,-c(1:6)])
g_valid_data_all_genes <- as.matrix(g_valid_data_all_genes[,-c(1:6)])
g_test_data_all_genes <- as.matrix(g_test_data_all_genes[,-c(1:6)])
g_full_train_data_all_genes <- as.matrix(g_full_train_data_all_genes[,-c(1:6)])


Labels_g_train_data_all_genes <- recode(Labels_g_train_data_all_genes,'RVI'=1, 'Not RVI'=0)

# training procedure

set.seed(1234)
train_full_Feature <- xgboost(data = g_train_data_all_genes, 
                              label = Labels_g_train_data_all_genes, 
                              max.depth = 6, 
                              eta = 0.3, 
                              nrounds = 100, 
                              eval_metric = "error",
                              objective = "binary:logistic",
                              verbose = 1)





# Print trained model
print(train_full_Feature)

# best performance on train (error)
print(min(train_full_Feature$evaluation_log$train_error))

# best performance on train (Accuracy)
print(1-min(train_full_Feature$evaluation_log$train_error))

##################################### Validation data prediction ####################################

Labels_g_valid_data_all_genes <- recode(Labels_g_valid_data_all_genes,'RVI'=1, 'Not RVI'=0)

# training procedure

set.seed(1234)
Valid_full_Feature <- xgboost(data = g_valid_data_all_genes, 
                              label = Labels_g_valid_data_all_genes, 
                              max.depth = 6, 
                              eta = 0.3, 
                              nrounds = 100, 
                              eval_metric = "error",
                              objective = "binary:logistic",
                              verbose = 1)



# Print valid model
print(Valid_full_Feature)

# best performance on valid (error)
print(min(Valid_full_Feature$evaluation_log$train_error))

# best performance on valid (Accuracy)
print(1-min(Valid_full_Feature$evaluation_log$train_error))


Valid_full_Feature$params


###################### Final model building using full train data ####################################

Labels_g_full_train_data_all_genes <- recode(Labels_g_full_train_data_all_genes,'RVI'=1, 'Not RVI'=0)

# training procedure

set.seed(1234)
final_trained_model <- xgboost(data = g_full_train_data_all_genes, 
                               label = Labels_g_full_train_data_all_genes, 
                               max.depth = Valid_full_Feature$params$max_depth, 
                               eta = Valid_full_Feature$params$eta, 
                               nrounds = 100, 
                               eval_metric = Valid_full_Feature$params$eval_metric,
                               objective = Valid_full_Feature$params$objective,
                               verbose = 1)



# Print valid model
print(final_trained_model)

# best performance on valid (error)
print(min(final_trained_model$evaluation_log$train_error))

# best performance on valid (Accuracy)
print(1-min(final_trained_model$evaluation_log$train_error))


##################################### test data prediction Dataset 4 ####################################


### Test set predicition
# Converting Label vector into factor as per the requirement of train function of caret package
Labels_g_test_data_time_t <- as.factor(hold_out_test$True_Class_Label)


# Converting training and test data into matrix as per the requirement of train function of caret package
g_test_data_time_t_all_genes <- as.matrix(hold_out_test[,-c(1:6)])

# Predicting Test Data
# Passing test data without labels (without fist column which contains labels)
(testPrediction1 <- predict(final_trained_model, newdata = g_test_data_time_t_all_genes))

for (i in 1:length(testPrediction1)) {
  if(testPrediction1[i] < 0.5){
    testPrediction1[i] <- 'Not RVI'
  }else{
    testPrediction1[i] <- 'RVI'
  }
  
}

# Show predicted labels
print(testPrediction1)
# Show Test Data
print(Labels_g_test_data_time_t)


### Performance measure

# Display confusion matrix Dataset 4
print(confusionMatrix(as.factor(testPrediction1), Labels_g_test_data_time_t))


