# Required Programming Language
# R version 3.6.2 or above

# Required packages
caret 6.0.93
dplyr 1.1.0
e1071 1.7.11
ggplot2 3.4.1
class 7.3.20
randomForest  4.7.1.1
pROC  1.18.0
kernlab 0.9.32

# load packages
library(caret) # for various machine learning functions
library(dplyr) # for data manupulation
library(e1071) # for various functions like confusion matrix
library(ggplot2) # for plots
library(class) # # For various classification functions
library(randomForest) # For random forest algorithm
library(pROC) # For plotting the ROC curves
library(kernlab) # For kernals

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

# training procedure

#### Model building using all the Genes
set.seed(1234)

# Assigning values to the parameter C
grid <- expand.grid(C = c(2^-5, 2^-3, 2^-1, 1, 2^1, 2^3, 2^5,
                          2^7, 2^9, 2^11, 2^13, 2^15))

# training LSVM classifier
train_full_Feature <- train(x= g_train_data_all_genes, # Training Data
                            y = Labels_g_train_data_all_genes,  # Class labels of training data
                            method = "svmLinear", # Train using LSVM
                            tuneGrid = grid) # Passing training control parameters


# Print trained model
print(train_full_Feature)

# standard deviation for 12023 genes
print((sd(train_full_Feature$resample$Accuracy)))


##################################### Validation data prediction ####################################

set.seed(1234)

# Assigning values to the parameter C
grid <- expand.grid(C = c(2^-5, 2^-3, 2^-1, 1, 2^1, 2^3, 2^5,
                          2^7, 2^9, 2^11, 2^13, 2^15))

# training LSVM classifier
Valid_full_Feature <- train(x= g_valid_data_all_genes, # Training Data
                            y = Labels_g_valid_data_all_genes,  # Class labels of training data
                            method = "svmLinear", # Train using LSVM
                            tuneGrid = grid) # Passing training control parameters
# Print trained model
print(Valid_full_Feature)

# standard deviation for 12023 genes
print((sd(Valid_full_Feature$resample$Accuracy)))



###################### Final model building using full train data ####################################

set.seed(1234)

# Assigning values to the parameter C
grid <- expand.grid(C = Valid_full_Feature$finalModel@param$C)

# training LSVM classifier
final_trained_model <- train(x= g_full_train_data_all_genes, # Training Data
                             y = Labels_g_full_train_data_all_genes,  # Class labels of training data
                             method = "svmLinear", # Train using LSVM
                             tuneGrid = grid) # Passing training control parameters
# Print trained model
print(final_trained_model)

# standard deviation for 12023 genes
print((sd(final_trained_model$resample$Accuracy)))
##################################### test data prediction ####################################


### Test set predicition
# Converting Label vector into factor as per the requirement of train function of caret package
Labels_g_test_data_time_t <- as.factor(hold_out_test$True_Class_Label)


# Converting training and test data into matrix as per the requirement of train function of caret package
g_test_data_time_t_all_genes <- as.matrix(hold_out_test[,-c(1:6)])

# Predicting Test Data
# Passing test data without labels (without fist column which contains labels)
(testPrediction1 <- predict(final_trained_model, newdata = g_test_data_time_t_all_genes))


# Show Test Data
print(Labels_g_test_data_time_t)


### Performance measure

# Display confusion matrix
print(confusionMatrix(testPrediction1, Labels_g_test_data_time_t))
