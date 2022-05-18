################ A Personalized Disease Diagnosis Assistant ######################

# load packages
library(caret) # for various machine learning functions
library(dplyr) # for data manupulation
library(e1071) # for various functions like confusion matrix
library(ggplot2) # for plots


############################### Input from user ##########################################

# read CTD knowledge graph that contains gene disease links 
KG_ctd_gene_disease <- read.csv(file = "/home/ghaver/Ghanshyam/all_data/KG/Gene_Diseases_Knowledge_Graph_with_Curated_Links_v3.csv", header = TRUE, sep=",")

# read Gene expression data collected at two required time-points for all those subjects for which disease diagnosis need to be performed
# required two time-points (t_1 = 0 hours i.e. healthy state and at time-point t_D i.e. diseased state or the time-point at which disease diagnosis has been requested)
Gene_expression_data <- read.csv(file = "/home/ghaver/Ghanshyam/all_data/Gene_Expression/Gene_Expression_Dataset_3_GSE90732.csv", header = TRUE, sep=",", check.names = FALSE)

Gene_expression_data[1:10,1:10]
dim(Gene_expression_data)

# go inside the train test split
total_data <- Gene_expression_data

# Collecting samples other than 0 time stamp (samples at day 2 or 3)
Data_target <- total_data %>% filter(Time > 0) 


## First dividing data into training and test set (single time point - target time point)
set.seed(7)
# Dividing data set into train (50%) and test (50%) using createDataPartition function of caret package
index_test <- createDataPartition(y = Data_target$True_Class_Label,
                                  p = 0.50, list = FALSE)
g_test_data <- Data_target[index_test, ]
g_train_data <- Data_target[-index_test, ]


dim(g_train_data)
dim(g_test_data)

set.seed(7)
# Dividing test data set further into holdout test set (25%) and validation set (25%) using createDataPartition function of caret package
index_test <- createDataPartition(y = g_test_data$True_Class_Label,
                                  p = 0.50, list = FALSE)
g_hold_out_test <- g_test_data[index_test, ]
g_valid_data <- g_test_data[-index_test, ]

dim(g_hold_out_test)
dim(g_valid_data)

# train test all time points
g_total_train_data <- total_data %>% filter(Super_Subject_ID %in% g_train_data$Super_Subject_ID)
g_total_valid_data <- total_data %>% filter(Super_Subject_ID %in% g_valid_data$Super_Subject_ID)
g_total_holdout_test_data <- total_data %>% filter(Super_Subject_ID %in% g_hold_out_test$Super_Subject_ID)

# enter the number of top "m" most likely diseases you want to see with computed probabilites for each patient 
m = 5

# enter values of MDEGs and LDEGs you want to consider or enter a range
# if you don't want to perform grid search and want to perform the experiment on single P and Q value 
# then assign the required value to variable P_start, and P_end and assign 0 to the P_step. Similarly the desired value (value of Q) to variable Q_start, and Q_end and 0 to the Q_step as shown below
# P_start <- 200
# P_end <- 200
# P_step <- 0
# Q_start <- 200
# Q_end <- 200
# Q_step <- 0


###################################    Training code   ##############################################################

# Grid search: alternatively to perform grid search, enter appropriate start, end, and step size for P and Q
P_start <- 25
P_end <- 500
P_step <- 25
Q_start <- 25
Q_end <- 500
Q_step <- 25

#################################################################################################

#################################################################################################

# extract all unique diseases from KG to assign computed disease weight
Data_Unique_Disease <- KG_ctd_gene_disease[!duplicated(KG_ctd_gene_disease[,c('Disease_ID')]), c(2,3)]

# reorder the columns
Data_Unique_Disease <- Data_Unique_Disease[ , c(2,1)]

# add a new column named Disease_Weight
Data_Unique_Disease <- Data_Unique_Disease %>% mutate(Disease_Weight = 0)

# temporary variable of Data_Unique_Disease
Data_Unique_Disease_initial_weights <- Data_Unique_Disease

# gene expression values start index in gene expression data
s_index <- 7

# computing length for grid search
PS <- length(seq(P_start,P_end,P_step))
QS <- length(seq(Q_start,Q_end,Q_step))

# creating data frame to compute and store accuracy at different values of P and Q and, and different value of top n diseases
Accuracy_matrix <- data.frame("P" = 1:(PS*QS), "Q" = 1:(PS*QS), "Acc_Top_1_Dis" = 1:(PS*QS), "Acc_Top_2_Dis" = 1:(PS*QS), "Acc_Top_3_Dis" = 1:(PS*QS), "Acc_Top_4_Dis" = 1:(PS*QS), "Acc_Top_5_Dis" = 1:(PS*QS), "Acc_Top_10_Dis" = 1:(PS*QS))

Accuracy_matrix <- Accuracy_matrix %>% mutate(Avg_Acc = 0)
Accuracy_matrix <- Accuracy_matrix %>% mutate(Time = 0)
# initializing the increamentor for Accuracy_matrix
Acc_index <- 1


#### Code for assigning weights to the diseases in KG based on changes observed in genes

for(p in seq(P_start,P_end,P_step)){  # loop of p is for top genes / P MDEGs
  
  for(q in seq(Q_start,Q_end,Q_step)){  # loop of q is for bottom genes / Q LDEGs
    
    # total number of subjects
    s <- dim(g_total_train_data)[1]/2
    
    # creating dataframe to store vlaues of predicted class labels
    
    predicted_info <- g_total_train_data[ , c(1:s_index-1)]
    predicted_info <- predicted_info %>% mutate(predicted_label_top_1 = 1:(2*s))
    predicted_info <- predicted_info %>% mutate(predicted_label_top_2 = 1:(2*s))
    predicted_info <- predicted_info %>% mutate(predicted_label_top_3 = 1:(2*s))
    predicted_info <- predicted_info %>% mutate(predicted_label_top_4 = 1:(2*s))
    predicted_info <- predicted_info %>% mutate(predicted_label_top_5 = 1:(2*s))
    predicted_info <- predicted_info %>% mutate(predicted_label_top_10 = 1:(2*s))
    
    Gene_Data_All_ti_prediction <- g_total_train_data[!g_total_train_data$Time == 0, c(1:(s_index-1))]
    Gene_Data_All_ti_prediction <- Gene_Data_All_ti_prediction %>% mutate(predicted_label_top_1 = 1:s)
    Gene_Data_All_ti_prediction <- Gene_Data_All_ti_prediction %>% mutate(predicted_label_top_2 = 1:s)
    Gene_Data_All_ti_prediction <- Gene_Data_All_ti_prediction %>% mutate(predicted_label_top_3 = 1:s)
    Gene_Data_All_ti_prediction <- Gene_Data_All_ti_prediction %>% mutate(predicted_label_top_4 = 1:s)
    Gene_Data_All_ti_prediction <- Gene_Data_All_ti_prediction %>% mutate(predicted_label_top_5 = 1:s)
    Gene_Data_All_ti_prediction <- Gene_Data_All_ti_prediction %>% mutate(predicted_label_top_10 = 1:s)
    
    All_Sub_temp_prediction <- data.frame("Top_1"=1:s, "Top_2"=1:s,"Top_3"=1:s, "Top_4"=1:s,"Top_5"=1:s, "Top_10"=1:s)
    
    # loop of l for number of subjects
    
    for(l in 1:s){ 
      
      # extracting data of lth subject 
      # making super_subject_id serial wise
      g_total_train_data$Super_Subject_ID <- Gene_expression_data$Super_Subject_ID[1:(2*s)]
      Gene_expression_data_sub_l <- g_total_train_data %>% filter(Super_Subject_ID == l)
      Gene_expression_data_sub_l <- Gene_expression_data_sub_l[ , -c(1:(s_index-1))]
      
      cat("\n")
      print("################################# New training data subject's computation start from here ##################################")
      cat("\n")
      print("Train Data: Subject/Patient id is:")
      print(l)
      cat("\n")
      
      # for each subject, again initialize the Data_Unique_Disease variable dataframe
      Data_Unique_Disease <- Data_Unique_Disease_initial_weights
      
      # compute changes in gene expression values (Target sample (TD) - Reference sample (T1)
      Gene_Transition_Matrix <- Gene_expression_data_sub_l[2, ] - Gene_expression_data_sub_l[1, ]
      
      # extracting top P MDEGs
      Gene_Transition_Matrix_top_p_Genes <- Gene_Transition_Matrix[ , order(-abs(Gene_Transition_Matrix[ , ]))]
      Gene_Transition_Matrix_top_p_Genes <- Gene_Transition_Matrix_top_p_Genes[ , c(1:p)]
      
      print("Top 5 Most Differencially Expressed Genes (MDEGs):")
      print(Gene_Transition_Matrix_top_p_Genes[ , 1:5])
      cat("\n")
      
      for(j in 1:p){ # loop of j for number of genes
        
        # extract all diseases' ids from KG that are associated with top P genes/ MDEGs of the subject
        Disease_IDs <- KG_ctd_gene_disease[KG_ctd_gene_disease$GeneSymbol == names(Gene_Transition_Matrix_top_p_Genes)[j], "Disease_ID" ]
        
        for(k in 1:length(Disease_IDs)){ # loop for every disease id
          
          # computing the weight/score for each disease for lth subject
          Data_Unique_Disease[Data_Unique_Disease$Disease_ID ==  Disease_IDs[k], "Disease_Weight"] <-   Data_Unique_Disease[Data_Unique_Disease$Disease_ID ==  Disease_IDs[k], "Disease_Weight"] + abs(Gene_Transition_Matrix_top_p_Genes[,j])
          
        }  # end for loop k
      } # end for loop j
      
      #### Code for down-weighting the diseases based on the associated Q LDEGs to them
      
      # extracting Q LDEGs
      Gene_Transition_Matrix_bottom_q_Genes <- Gene_Transition_Matrix[ , order(abs(Gene_Transition_Matrix[ , ]))]
      Gene_Transition_Matrix_bottom_q_Genes <- Gene_Transition_Matrix_bottom_q_Genes[ , c(1:q)]
      
      print("5 Least Differencially Expressed Genes (LDEGs):")
      print(Gene_Transition_Matrix_bottom_q_Genes[ , 1:5])
      cat("\n")
      
      for(j in 1:q){ # loop of j for number of bottom genes
        
        # extract all diseases' ids from KG that are associated with bottom Q genes / LDEGs of the subject
        Disease_IDs <- KG_ctd_gene_disease[KG_ctd_gene_disease$GeneSymbol == names(Gene_Transition_Matrix_bottom_q_Genes)[j], "Disease_ID" ]
        
        for(k in 1:length(Disease_IDs)){ # loop for every disease id
          
          # down-weighting the diseases based on the associated Q LDEGs to them
          Data_Unique_Disease[Data_Unique_Disease$Disease_ID ==  Disease_IDs[k], "Disease_Weight"] <-  Data_Unique_Disease[Data_Unique_Disease$Disease_ID ==  Disease_IDs[k], "Disease_Weight"] - (abs(Gene_Transition_Matrix_top_p_Genes[,1]) - abs(Gene_Transition_Matrix_bottom_q_Genes[,j]))
          
        }  # end for loop k
      } # end for loop j
      
      
      # create file name to write data into csv file
      #file_name <- paste("Disease_Weight_Train_Sub_",l,"_p_",p,"_q_",q,".csv", collapse = "",sep="")
      
      # write data into csv file
      #write.csv(Data_Unique_Disease[order(-Data_Unique_Disease$Disease_Weight), ], file = file_name, row.names = FALSE)
      
      print("Value of P for train data (number of MDEGs) is :")
      print(p)
      cat("\n")
      
      print("Value of Q for train data (number of LDEGs) is :")
      print(q)
      cat("\n")
      
      print("This train subject at this time point has following True Class Label:")
      print(g_total_train_data[ g_total_train_data$Super_Subject_ID == l , ]$True_Class_Label[2])
      cat("\n")
      
      for(i in 1:6){ # loop for how many top disease you want to look for Acc calc (current loop is for top 1 to 5 and top 10 predicted diseases )
        
        if(i<6){
          Top_Disease_Names <- Data_Unique_Disease[order(-Data_Unique_Disease$Disease_Weight), ][1:i, "DiseaseName"]
        }else{
          Top_Disease_Names <- Data_Unique_Disease[order(-Data_Unique_Disease$Disease_Weight), ][1:10, "DiseaseName"]
        }
        
        if(any(Top_Disease_Names == "Respiratory Viral Infection")){
          predicted_info[predicted_info$Super_Subject_ID == l , (i+s_index-1)][2] <- "RVI"
          Gene_Data_All_ti_prediction[Gene_Data_All_ti_prediction$Super_Subject_ID == l , (i+s_index-1)][1] <- "RVI"
          All_Sub_temp_prediction[l,i] <- "RVI"
        }else{
          predicted_info[predicted_info$Super_Subject_ID == l , (i+s_index-1)][2] <- "Not RVI"
          Gene_Data_All_ti_prediction[Gene_Data_All_ti_prediction$Super_Subject_ID == l , (i+s_index-1)][1] <- "Not RVI"
          All_Sub_temp_prediction[l,i] <- "Not RVI"
        }
        
        if(i<6){
          print(paste("Predicted label using top ", i, "disease is:"))
          print(All_Sub_temp_prediction[l,i])
          cat("\n")
          
        }else{
          print("Predicted label using top 10 disease is:")
          print(All_Sub_temp_prediction[l,i])
          cat("\n")
          
        }
        
      }
      
      # creating dataframe to store computed disease probabilites
      Data_Unique_Disease_with_Probability <- Data_Unique_Disease[order(-Data_Unique_Disease$Disease_Weight), ][1:m, ]
      
      # function to compute probability of predicted disease using generalized softmax
      generalized_softmax <- function(Dis_weight){
        S_Prob <- Dis_weight
        e_Dis_weight <- exp(Dis_weight)
        for(i in 1:length(Dis_weight)){
          S_Prob[i] <- (exp(Dis_weight[i])/sum(e_Dis_weight))*100
        }
        return(S_Prob)
      }
      
      # calling probability function
      Disease_Prob <- generalized_softmax(Data_Unique_Disease_with_Probability$Disease_Weight)
      Data_Unique_Disease_with_Probability <- Data_Unique_Disease_with_Probability %>% mutate(Disease_Probability = Disease_Prob)
      
      print("Top 5 most likely diseases predicted for this subject are:")
      print(Data_Unique_Disease_with_Probability)
      cat("\n")
      
      cat("\n")
      print("################################## This subject's predictions ends here #####################################")
      cat("\n")
    } # end for loop 
    
    
    
    # create file name to write data into csv file
    file_name_1 <- paste("predicted_info_train_Sub_",l,"_p_",p,"_q_",q,".csv", collapse = "",sep="")
    write.csv(Gene_Data_All_ti_prediction, file = file_name_1, row.names = FALSE)
    
    cat("\n")
    print("Accuracy calculated using all subjects considering different values of top n diseases:")
    cat("\n")
    
    ## computing accuracy using all subjects of test data 
    # computing accuraacy using confusion matrix if both positive and negative subjects are there in the data
    if(any(Gene_Data_All_ti_prediction$True_Class_Label == "Not RVI")){
      for(i in 1:6){
        confusion_mat <- confusionMatrix( as.factor(All_Sub_temp_prediction[,i]), as.factor(Gene_Data_All_ti_prediction$True_Class_Label), positive = "RVI")
        
        if(i<6){
          cat("\n")
          print(paste("Accuracy using top", i, "disease is:"))
          print(confusion_mat)
          cat("\n")
          
        }else{
          cat("\n")
          print("Accuracy using top 10 diseases is:")
          print(confusion_mat)
          cat("\n")
          
        }
        
        Accuracy_matrix[Acc_index, (2+i)] <- confusion_mat$overall[1]
      }
    }else{ #computing accuraacy using stardard method (hit rate) if only positive subjects are there in the data
      for(i in 1:6){
        hit <- 0
        for(k in 1:dim(All_Sub_temp_prediction)[1]){
          if(Gene_Data_All_ti_prediction[k,"True_Class_Label"] == All_Sub_temp_prediction[k,i]){
            hit <- hit + 1
          }
        }
        Acc <- hit/dim(All_Sub_temp_prediction)[1]
        cat("\n")
        print(paste("Accuracy using top", i, "disease is:"))
        print(Acc)
        if(i<6){
          cat("\n")
          print(paste("Accuracy using top", i, "disease is:"))
          print(Acc)
          cat("\n")
          
        }else{
          cat("\n")
          print("Accuracy using top 10 diseases is:")
          print(Acc)
          cat("\n")
          
        }
        Accuracy_matrix[Acc_index, (2+i)] <- Acc
        
      }
      
    }
    
    # assign current value of P and Q 
    Accuracy_matrix$P[Acc_index] <- p
    Accuracy_matrix$Q[Acc_index] <- q
    
    cat("\n")
    print("Accuracy calculated using all subjects on different values of top n diseases and varing P and Q:")
    cat("\n")
    print("Accuracy Matrix (Grid Search):")
    print(Accuracy_matrix[1:Acc_index, ])
    cat("\n")
    
    
    # adding time to the Accuracy_matrix
    Accuracy_matrix$Time[Acc_index] <- each_itteration_time[3]
    
    # adding average acc to the Accuracy_matrix
    Accuracy_matrix$Avg_Acc <- rowMeans(Accuracy_matrix[,c(3:8)])
    
    
    
    
    
    # create file name of accuracy matrix to write data into csv file
    file_name_1 <- paste("Accuracy_matrix_Train_Sub",l,"_p_",p,"_q_",q,".csv", collapse = "",sep="")
    write.csv(Accuracy_matrix, file = file_name_1, row.names = FALSE)
    
    # increamentor 
    Acc_index <- Acc_index +1
    
    
    
    
    print("################################## Current itteration of the loop for Q ends here ##########################################")
    cat("\n")
    
  } # ending loop q 
  
  print("################################## Current itteration of the loop for P ends here ##########################################")
  cat("\n")
  
} # ending loop p

print("Accuracy Matrix (Training):")
print(Accuracy_matrix)



#######################################################################################################################################

##############################################    Validation code   ###################################################################

#######################################################################################################################################

# creating data frame to compute and store accuracy at different values of P and Q and, and different value of top n diseases
Accuracy_matrix <- data.frame("P" = 1:(PS*QS), "Q" = 1:(PS*QS), "Acc_Top_1_Dis" = 1:(PS*QS), "Acc_Top_2_Dis" = 1:(PS*QS), "Acc_Top_3_Dis" = 1:(PS*QS), "Acc_Top_4_Dis" = 1:(PS*QS), "Acc_Top_5_Dis" = 1:(PS*QS), "Acc_Top_10_Dis" = 1:(PS*QS))

Accuracy_matrix <- Accuracy_matrix %>% mutate(Avg_Acc = 0)
Accuracy_matrix <- Accuracy_matrix %>% mutate(Time = 0)
# initializing the increamentor for Accuracy_matrix
Acc_index <- 1

#### Code for assigning weights to the diseases in KG based on changes observed in genes

for(p in seq(P_start,P_end,P_step)){  # loop of p is for top genes / P MDEGs
  
  for(q in seq(Q_start,Q_end,Q_step)){  # loop of q is for bottom genes / Q LDEGs
    
    # total number of subjects
    s <- dim(g_total_valid_data)[1]/2
    
    # creating dataframe to store vlaues of predicted class labels
    
    predicted_info <- g_total_valid_data[ , c(1:s_index-1)]
    predicted_info <- predicted_info %>% mutate(predicted_label_top_1 = 1:(2*s))
    predicted_info <- predicted_info %>% mutate(predicted_label_top_2 = 1:(2*s))
    predicted_info <- predicted_info %>% mutate(predicted_label_top_3 = 1:(2*s))
    predicted_info <- predicted_info %>% mutate(predicted_label_top_4 = 1:(2*s))
    predicted_info <- predicted_info %>% mutate(predicted_label_top_5 = 1:(2*s))
    predicted_info <- predicted_info %>% mutate(predicted_label_top_10 = 1:(2*s))
    
    Gene_Data_All_ti_prediction <- g_total_valid_data[!g_total_valid_data$Time == 0, c(1:(s_index-1))]
    Gene_Data_All_ti_prediction <- Gene_Data_All_ti_prediction %>% mutate(predicted_label_top_1 = 1:s)
    Gene_Data_All_ti_prediction <- Gene_Data_All_ti_prediction %>% mutate(predicted_label_top_2 = 1:s)
    Gene_Data_All_ti_prediction <- Gene_Data_All_ti_prediction %>% mutate(predicted_label_top_3 = 1:s)
    Gene_Data_All_ti_prediction <- Gene_Data_All_ti_prediction %>% mutate(predicted_label_top_4 = 1:s)
    Gene_Data_All_ti_prediction <- Gene_Data_All_ti_prediction %>% mutate(predicted_label_top_5 = 1:s)
    Gene_Data_All_ti_prediction <- Gene_Data_All_ti_prediction %>% mutate(predicted_label_top_10 = 1:s)
    
    All_Sub_temp_prediction <- data.frame("Top_1"=1:s, "Top_2"=1:s,"Top_3"=1:s, "Top_4"=1:s,"Top_5"=1:s, "Top_10"=1:s)
    
    # loop of l for number of subjects
    
    for(l in 1:s){ 
      
      # extracting data of lth subject 
      # making super_subject_id serial wise
      g_total_valid_data$Super_Subject_ID <- Gene_expression_data$Super_Subject_ID[1:(2*s)]
      Gene_expression_data_sub_l <- g_total_valid_data %>% filter(Super_Subject_ID == l)
      Gene_expression_data_sub_l <- Gene_expression_data_sub_l[ , -c(1:(s_index-1))]
      
      cat("\n")
      print("################################# New training data subject's computation start from here ##################################")
      cat("\n")
      print("Train Data: Subject/Patient id is:")
      print(l)
      cat("\n")
      
      # for each subject, again initialize the Data_Unique_Disease variable dataframe
      Data_Unique_Disease <- Data_Unique_Disease_initial_weights
      
      # compute changes in gene expression values (Target sample (TD) - Reference sample (T1)
      Gene_Transition_Matrix <- Gene_expression_data_sub_l[2, ] - Gene_expression_data_sub_l[1, ]
      
      # extracting top P MDEGs
      Gene_Transition_Matrix_top_p_Genes <- Gene_Transition_Matrix[ , order(-abs(Gene_Transition_Matrix[ , ]))]
      Gene_Transition_Matrix_top_p_Genes <- Gene_Transition_Matrix_top_p_Genes[ , c(1:p)]
      
      print("Top 5 Most Differencially Expressed Genes (MDEGs):")
      print(Gene_Transition_Matrix_top_p_Genes[ , 1:5])
      cat("\n")
      
      for(j in 1:p){ # loop of j for number of genes
        
        # extract all diseases' ids from KG that are associated with top P genes/ MDEGs of the subject
        Disease_IDs <- KG_ctd_gene_disease[KG_ctd_gene_disease$GeneSymbol == names(Gene_Transition_Matrix_top_p_Genes)[j], "Disease_ID" ]
        
        for(k in 1:length(Disease_IDs)){ # loop for every disease id
          
          # computing the weight/score for each disease for lth subject
          Data_Unique_Disease[Data_Unique_Disease$Disease_ID ==  Disease_IDs[k], "Disease_Weight"] <-   Data_Unique_Disease[Data_Unique_Disease$Disease_ID ==  Disease_IDs[k], "Disease_Weight"] + abs(Gene_Transition_Matrix_top_p_Genes[,j])
          
        }  # end for loop k
      } # end for loop j
      
      #### Code for down-weighting the diseases based on the associated Q LDEGs to them
      
      # extracting Q LDEGs
      Gene_Transition_Matrix_bottom_q_Genes <- Gene_Transition_Matrix[ , order(abs(Gene_Transition_Matrix[ , ]))]
      Gene_Transition_Matrix_bottom_q_Genes <- Gene_Transition_Matrix_bottom_q_Genes[ , c(1:q)]
      
      print("5 Least Differencially Expressed Genes (LDEGs):")
      print(Gene_Transition_Matrix_bottom_q_Genes[ , 1:5])
      cat("\n")
      
      for(j in 1:q){ # loop of j for number of bottom genes
        
        # extract all diseases' ids from KG that are associated with bottom Q genes / LDEGs of the subject
        Disease_IDs <- KG_ctd_gene_disease[KG_ctd_gene_disease$GeneSymbol == names(Gene_Transition_Matrix_bottom_q_Genes)[j], "Disease_ID" ]
        
        for(k in 1:length(Disease_IDs)){ # loop for every disease id
          
          # down-weighting the diseases based on the associated Q LDEGs to them
          Data_Unique_Disease[Data_Unique_Disease$Disease_ID ==  Disease_IDs[k], "Disease_Weight"] <-  Data_Unique_Disease[Data_Unique_Disease$Disease_ID ==  Disease_IDs[k], "Disease_Weight"] - (abs(Gene_Transition_Matrix_top_p_Genes[,1]) - abs(Gene_Transition_Matrix_bottom_q_Genes[,j]))
          
        }  # end for loop k
      } # end for loop j
      
      
      # create file name to write data into csv file
      #file_name <- paste("Disease_Weight_Train_Sub_",l,"_p_",p,"_q_",q,".csv", collapse = "",sep="")
      
      # write data into csv file
      #write.csv(Data_Unique_Disease[order(-Data_Unique_Disease$Disease_Weight), ], file = file_name, row.names = FALSE)
      
      print("Value of P for validation data (number of MDEGs) is :")
      print(p)
      cat("\n")
      
      print("Value of Q for validation data (number of LDEGs) is :")
      print(q)
      cat("\n")
      
      print("This validation subject at this time point has following True Class Label:")
      print(g_total_valid_data[ g_total_valid_data$Super_Subject_ID == l , ]$True_Class_Label[2])
      cat("\n")
      
      for(i in 1:6){ # loop for how many top disease you want to look for Acc calc (current loop is for top 1 to 5 and top 10 predicted diseases )
        
        if(i<6){
          Top_Disease_Names <- Data_Unique_Disease[order(-Data_Unique_Disease$Disease_Weight), ][1:i, "DiseaseName"]
        }else{
          Top_Disease_Names <- Data_Unique_Disease[order(-Data_Unique_Disease$Disease_Weight), ][1:10, "DiseaseName"]
        }
        
        if(any(Top_Disease_Names == "Respiratory Viral Infection")){
          predicted_info[predicted_info$Super_Subject_ID == l , (i+s_index-1)][2] <- "RVI"
          Gene_Data_All_ti_prediction[Gene_Data_All_ti_prediction$Super_Subject_ID == l , (i+s_index-1)][1] <- "RVI"
          All_Sub_temp_prediction[l,i] <- "RVI"
        }else{
          predicted_info[predicted_info$Super_Subject_ID == l , (i+s_index-1)][2] <- "Not RVI"
          Gene_Data_All_ti_prediction[Gene_Data_All_ti_prediction$Super_Subject_ID == l , (i+s_index-1)][1] <- "Not RVI"
          All_Sub_temp_prediction[l,i] <- "Not RVI"
        }
        
        if(i<6){
          print(paste("Predicted label using top ", i, "disease is:"))
          print(All_Sub_temp_prediction[l,i])
          cat("\n")
          
        }else{
          print("Predicted label using top 10 disease is:")
          print(All_Sub_temp_prediction[l,i])
          cat("\n")
          
        }
        
      }
      
      # creating dataframe to store computed disease probabilites
      Data_Unique_Disease_with_Probability <- Data_Unique_Disease[order(-Data_Unique_Disease$Disease_Weight), ][1:m, ]
      
      # function to compute probability of predicted disease using generalized softmax
      generalized_softmax <- function(Dis_weight){
        S_Prob <- Dis_weight
        e_Dis_weight <- exp(Dis_weight)
        for(i in 1:length(Dis_weight)){
          S_Prob[i] <- (exp(Dis_weight[i])/sum(e_Dis_weight))*100
        }
        return(S_Prob)
      }
      
      # calling probability function
      Disease_Prob <- generalized_softmax(Data_Unique_Disease_with_Probability$Disease_Weight)
      Data_Unique_Disease_with_Probability <- Data_Unique_Disease_with_Probability %>% mutate(Disease_Probability = Disease_Prob)
      
      print("Top 5 most likely diseases predicted for this subject are:")
      print(Data_Unique_Disease_with_Probability)
      cat("\n")
      
      cat("\n")
      print("################################## This subject's predictions ends here #####################################")
      cat("\n")
    } # end for loop l
    
    # create file name to write data into csv file
    file_name_1 <- paste("predicted_info_train_Sub_",l,"_p_",p,"_q_",q,".csv", collapse = "",sep="")
    write.csv(Gene_Data_All_ti_prediction, file = file_name_1, row.names = FALSE)
    
    cat("\n")
    print("Accuracy calculated using all subjects considering different values of top n diseases:")
    cat("\n")
    
    ## computing accuracy using all subjects of test data 
    # computing accuraacy using confusion matrix if both positive and negative subjects are there in the data
    if(any(Gene_Data_All_ti_prediction$True_Class_Label == "Not RVI")){
      for(i in 1:6){
        confusion_mat <- confusionMatrix( as.factor(All_Sub_temp_prediction[,i]), as.factor(Gene_Data_All_ti_prediction$True_Class_Label), positive = "RVI")
        
        if(i<6){
          cat("\n")
          print(paste("Accuracy using top", i, "disease is:"))
          print(confusion_mat)
          cat("\n")
          
        }else{
          cat("\n")
          print("Accuracy using top 10 diseases is:")
          print(confusion_mat)
          cat("\n")
          
        }
        
        Accuracy_matrix[Acc_index, (2+i)] <- confusion_mat$overall[1]
      }
    }else{ #computing accuraacy using stardard method (hit rate) if only positive subjects are there in the data
      for(i in 1:6){
        hit <- 0
        for(k in 1:dim(All_Sub_temp_prediction)[1]){
          if(Gene_Data_All_ti_prediction[k,"True_Class_Label"] == All_Sub_temp_prediction[k,i]){
            hit <- hit + 1
          }
        }
        Acc <- hit/dim(All_Sub_temp_prediction)[1]
        cat("\n")
        print(paste("Accuracy using top", i, "disease is:"))
        print(Acc)
        if(i<6){
          cat("\n")
          print(paste("Accuracy using top", i, "disease is:"))
          print(Acc)
          cat("\n")
          
        }else{
          cat("\n")
          print("Accuracy using top 10 diseases is:")
          print(Acc)
          cat("\n")
          
        }
        Accuracy_matrix[Acc_index, (2+i)] <- Acc
        
      }
      
    }
    
    # assign current value of P and Q 
    Accuracy_matrix$P[Acc_index] <- p
    Accuracy_matrix$Q[Acc_index] <- q
    
    cat("\n")
    print("Accuracy calculated using all subjects on different values of top n diseases and varing P and Q:")
    cat("\n")
    print("Accuracy Matrix (Grid Search):")
    print(Accuracy_matrix[1:Acc_index, ])
    cat("\n")
    
    # adding time to the Accuracy_matrix
    Accuracy_matrix$Time[Acc_index] <- each_itteration_time[3]
    
    # adding average acc to the Accuracy_matrix
    Accuracy_matrix$Avg_Acc <- rowMeans(Accuracy_matrix[,c(3:8)])
    
    # create file name of accuracy matrix to write data into csv file
    file_name_1 <- paste("Accuracy_matrix_Valid_Sub",l,"_p_",p,"_q_",q,".csv", collapse = "",sep="")
    write.csv(Accuracy_matrix, file = file_name_1, row.names = FALSE)
    
    # increamentor 
    Acc_index <- Acc_index +1
    
    print("################################## Current itteration of the loop for Q ends here ##########################################")
    cat("\n")
    
  } # ending loop q 
  
  print("################################## Current itteration of the loop for P ends here ##########################################")
  cat("\n")
  
} # ending loop p

print("Accuracy Matrix (Training):")
print(Accuracy_matrix)



#######################################################################################################################################

####################################################    Test code   ###################################################################

#######################################################################################################################################

# Find Best parameter
P <- Accuracy_matrix[Accuracy_matrix$Avg_Acc == max(Accuracy_matrix$Avg_Acc), "P"][1]
Q <- Accuracy_matrix[Accuracy_matrix$Avg_Acc == max(Accuracy_matrix$Avg_Acc), "Q"][1]

# Grid search: alternatively to perform grid search, enter appropriate start, end, and step size for P and Q
P_start <- P
P_end <- P
P_step <- 0
Q_start <- Q
Q_end <- Q
Q_step <- 0


p <- P
q <- Q
#################################################################################################

#################################################################################################

# extract all unique diseases from KG to assign computed disease weight
Data_Unique_Disease <- KG_ctd_gene_disease[!duplicated(KG_ctd_gene_disease[,c('Disease_ID')]), c(2,3)]

# reorder the columns
Data_Unique_Disease <- Data_Unique_Disease[ , c(2,1)]

# add a new column named Disease_Weight
Data_Unique_Disease <- Data_Unique_Disease %>% mutate(Disease_Weight = 0)

# temporary variable of Data_Unique_Disease
Data_Unique_Disease_initial_weights <- Data_Unique_Disease

# gene expression values start index in gene expression data
s_index <- 7

# computing length for grid search
PS <- length(seq(P_start,P_end,P_step))
QS <- length(seq(Q_start,Q_end,Q_step))

# creating data frame to compute and store accuracy at different values of P and Q and, and different value of top n diseases
Accuracy_matrix <- data.frame("P" = 1:(PS*QS), "Q" = 1:(PS*QS), "Acc_Top_1_Dis" = 1:(PS*QS), "Acc_Top_2_Dis" = 1:(PS*QS), "Acc_Top_3_Dis" = 1:(PS*QS), "Acc_Top_4_Dis" = 1:(PS*QS), "Acc_Top_5_Dis" = 1:(PS*QS), "Acc_Top_10_Dis" = 1:(PS*QS))

# initializing the increamentor for Accuracy_matrix
Acc_index <- 1

#### Code for assigning weights to the diseases in KG based on changes observed in genes

# total number of subjects
s <- dim(g_total_holdout_test_data)[1]/2

# creating dataframe to store vlaues of predicted class labels

predicted_info <- g_total_holdout_test_data[ , c(1:s_index-1)]
predicted_info <- predicted_info %>% mutate(predicted_label_top_1 = 1:(2*s))
predicted_info <- predicted_info %>% mutate(predicted_label_top_2 = 1:(2*s))
predicted_info <- predicted_info %>% mutate(predicted_label_top_3 = 1:(2*s))
predicted_info <- predicted_info %>% mutate(predicted_label_top_4 = 1:(2*s))
predicted_info <- predicted_info %>% mutate(predicted_label_top_5 = 1:(2*s))
predicted_info <- predicted_info %>% mutate(predicted_label_top_10 = 1:(2*s))

Gene_Data_All_ti_prediction <- g_total_holdout_test_data[!g_total_holdout_test_data$Time == 0, c(1:(s_index-1))]
Gene_Data_All_ti_prediction <- Gene_Data_All_ti_prediction %>% mutate(predicted_label_top_1 = 1:s)
Gene_Data_All_ti_prediction <- Gene_Data_All_ti_prediction %>% mutate(predicted_label_top_2 = 1:s)
Gene_Data_All_ti_prediction <- Gene_Data_All_ti_prediction %>% mutate(predicted_label_top_3 = 1:s)
Gene_Data_All_ti_prediction <- Gene_Data_All_ti_prediction %>% mutate(predicted_label_top_4 = 1:s)
Gene_Data_All_ti_prediction <- Gene_Data_All_ti_prediction %>% mutate(predicted_label_top_5 = 1:s)
Gene_Data_All_ti_prediction <- Gene_Data_All_ti_prediction %>% mutate(predicted_label_top_10 = 1:s)

All_Sub_temp_prediction <- data.frame("Top_1"=1:s, "Top_2"=1:s,"Top_3"=1:s, "Top_4"=1:s,"Top_5"=1:s, "Top_10"=1:s)

# loop of l for number of subjects

for(l in 1:s){ 
  
  # extracting data of lth subject 
  # making super_subject_id serial wise
  g_total_holdout_test_data$Super_Subject_ID <- Gene_expression_data$Super_Subject_ID[1:(2*s)]
  Gene_expression_data_sub_l <- g_total_holdout_test_data %>% filter(Super_Subject_ID == l)
  Gene_expression_data_sub_l <- Gene_expression_data_sub_l[ , -c(1:(s_index-1))]
  
  cat("\n")
  print("################################# New subject's computation start from here ##################################")
  cat("\n")
  print("Subject/Patient id is:")
  print(l)
  cat("\n")
  
  # for each subject, again initialize the Data_Unique_Disease variable dataframe
  Data_Unique_Disease <- Data_Unique_Disease_initial_weights
  
  # compute changes in gene expression values (Target sample (TD) - Reference sample (T1)
  Gene_Transition_Matrix <- Gene_expression_data_sub_l[2, ] - Gene_expression_data_sub_l[1, ]
  
  # extracting top P MDEGs
  Gene_Transition_Matrix_top_p_Genes <- Gene_Transition_Matrix[ , order(-abs(Gene_Transition_Matrix[ , ]))]
  Gene_Transition_Matrix_top_p_Genes <- Gene_Transition_Matrix_top_p_Genes[ , c(1:p)]
  
  print("Top 5 Most Differencially Expressed Genes (MDEGs):")
  print(Gene_Transition_Matrix_top_p_Genes[ , 1:5])
  cat("\n")
  
  for(j in 1:p){ # loop of j for number of genes
    
    # extract all diseases' ids from KG that are associated with top P genes/ MDEGs of the subject
    Disease_IDs <- KG_ctd_gene_disease[KG_ctd_gene_disease$GeneSymbol == names(Gene_Transition_Matrix_top_p_Genes)[j], "Disease_ID" ]
    
    for(k in 1:length(Disease_IDs)){ # loop for every disease id
      
      # computing the weight/score for each disease for lth subject
      Data_Unique_Disease[Data_Unique_Disease$Disease_ID ==  Disease_IDs[k], "Disease_Weight"] <-   Data_Unique_Disease[Data_Unique_Disease$Disease_ID ==  Disease_IDs[k], "Disease_Weight"] + abs(Gene_Transition_Matrix_top_p_Genes[,j])
      
    }  # end for loop k
  } # end for loop j
  
  #### Code for down-weighting the diseases based on the associated Q LDEGs to them
  
  # extracting Q LDEGs
  Gene_Transition_Matrix_bottom_q_Genes <- Gene_Transition_Matrix[ , order(abs(Gene_Transition_Matrix[ , ]))]
  Gene_Transition_Matrix_bottom_q_Genes <- Gene_Transition_Matrix_bottom_q_Genes[ , c(1:q)]
  
  print("5 Least Differencially Expressed Genes (LDEGs):")
  print(Gene_Transition_Matrix_bottom_q_Genes[ , 1:5])
  cat("\n")
  
  for(j in 1:q){ # loop of j for number of bottom genes
    
    # extract all diseases' ids from KG that are associated with bottom Q genes / LDEGs of the subject
    Disease_IDs <- KG_ctd_gene_disease[KG_ctd_gene_disease$GeneSymbol == names(Gene_Transition_Matrix_bottom_q_Genes)[j], "Disease_ID" ]
    
    for(k in 1:length(Disease_IDs)){ # loop for every disease id
      
      # down-weighting the diseases based on the associated Q LDEGs to them
      Data_Unique_Disease[Data_Unique_Disease$Disease_ID ==  Disease_IDs[k], "Disease_Weight"] <-  Data_Unique_Disease[Data_Unique_Disease$Disease_ID ==  Disease_IDs[k], "Disease_Weight"] - (abs(Gene_Transition_Matrix_top_p_Genes[,1]) - abs(Gene_Transition_Matrix_bottom_q_Genes[,j]))
      
    }  # end for loop k
  } # end for loop j
  
  
  # create file name to write data into csv file
  #file_name <- paste("Disease_Weight_Test_Sub_",l,"_p_",p,"_q_",q,".csv", collapse = "",sep="")
  
  # write data into csv file
  #write.csv(Data_Unique_Disease[order(-Data_Unique_Disease$Disease_Weight), ], file = file_name, row.names = FALSE)
  
  print("Value of P for test (number of MDEGs) is :")
  print(p)
  cat("\n")
  
  print("Value of Q for test (number of LDEGs) is :")
  print(q)
  cat("\n")
  
  print("This test subject at this time point has following True Class Label:")
  print(g_total_holdout_test_data[ g_total_holdout_test_data$Super_Subject_ID == l , ]$True_Class_Label[2])
  cat("\n")
  
  for(i in 1:6){ # loop for how many top disease you want to look for Acc calc (current loop is for top 1 to 5 and top 10 predicted diseases )
    
    if(i<6){
      Top_Disease_Names <- Data_Unique_Disease[order(-Data_Unique_Disease$Disease_Weight), ][1:i, "DiseaseName"]
    }else{
      Top_Disease_Names <- Data_Unique_Disease[order(-Data_Unique_Disease$Disease_Weight), ][1:10, "DiseaseName"]
    }
    
    if(any(Top_Disease_Names == "Respiratory Viral Infection")){
      predicted_info[predicted_info$Super_Subject_ID == l , (i+s_index-1)][2] <- "RVI"
      Gene_Data_All_ti_prediction[Gene_Data_All_ti_prediction$Super_Subject_ID == l , (i+s_index-1)][1] <- "RVI"
      All_Sub_temp_prediction[l,i] <- "RVI"
    }else{
      predicted_info[predicted_info$Super_Subject_ID == l , (i+s_index-1)][2] <- "Not RVI"
      Gene_Data_All_ti_prediction[Gene_Data_All_ti_prediction$Super_Subject_ID == l , (i+s_index-1)][1] <- "Not RVI"
      All_Sub_temp_prediction[l,i] <- "Not RVI"
    }
    
    if(i<6){
      print(paste("Predicted test data label using top ", i, "disease is:"))
      print(All_Sub_temp_prediction[l,i])
      cat("\n")
      
    }else{
      print("Predicted test data label using top 10 disease is:")
      print(All_Sub_temp_prediction[l,i])
      cat("\n")
      
    }
    
  }
  
  # creating dataframe to store computed disease probabilites
  Data_Unique_Disease_with_Probability <- Data_Unique_Disease[order(-Data_Unique_Disease$Disease_Weight), ][1:m, ]
  
  # function to compute probability of predicted disease using generalized softmax
  generalized_softmax <- function(Dis_weight){
    S_Prob <- Dis_weight
    e_Dis_weight <- exp(Dis_weight)
    for(i in 1:length(Dis_weight)){
      S_Prob[i] <- (exp(Dis_weight[i])/sum(e_Dis_weight))*100
    }
    return(S_Prob)
  }
  
  # calling probability function
  Disease_Prob <- generalized_softmax(Data_Unique_Disease_with_Probability$Disease_Weight)
  Data_Unique_Disease_with_Probability <- Data_Unique_Disease_with_Probability %>% mutate(Disease_Probability = Disease_Prob)
  
  print("Top 5 most likely diseases predicted for this test subject are:")
  print(Data_Unique_Disease_with_Probability)
  cat("\n")
  
  cat("\n")
  print("################################## This test subject's predictions ends here #####################################")
  cat("\n")
} # end for loop l

# create file name to write data into csv file
file_name_1 <- paste("predicted_info_holdout_test_Sub_",l,"_p_",p,"_q_",q,".csv", collapse = "",sep="")
write.csv(Gene_Data_All_ti_prediction, file = file_name_1, row.names = FALSE)

cat("\n")
print("Accuracy calculated using all test subjects considering different values of top n diseases:")
cat("\n")

## computing accuracy using all subjects of test data 
# computing accuraacy using confusion matrix if both positive and negative subjects are there in the data
if(any(Gene_Data_All_ti_prediction$True_Class_Label == "Not RVI")){
  for(i in 1:6){
    confusion_mat <- confusionMatrix( as.factor(All_Sub_temp_prediction[,i]), as.factor(Gene_Data_All_ti_prediction$True_Class_Label), positive = "RVI")
    
    if(i<6){
      cat("\n")
      print(paste("Accuracy using top", i, "disease is:"))
      print(confusion_mat)
      cat("\n")
      
    }else{
      cat("\n")
      print("Accuracy using top 10 diseases is:")
      print(confusion_mat)
      cat("\n")
      
    }
    
    Accuracy_matrix[Acc_index, (2+i)] <- confusion_mat$overall[1]
  }
}else{ #computing accuraacy using stardard method (hit rate) if only positive subjects are there in the data
  for(i in 1:6){
    hit <- 0
    for(k in 1:dim(All_Sub_temp_prediction)[1]){
      if(Gene_Data_All_ti_prediction[k,"True_Class_Label"] == All_Sub_temp_prediction[k,i]){
        hit <- hit + 1
      }
    }
    Acc <- hit/dim(All_Sub_temp_prediction)[1]
    cat("\n")
    print(paste("Accuracy using top", i, "disease is:"))
    print(Acc)
    if(i<6){
      cat("\n")
      print(paste("Accuracy using top", i, "disease is:"))
      print(Acc)
      cat("\n")
      
    }else{
      cat("\n")
      print("Accuracy using top 10 diseases is:")
      print(Acc)
      cat("\n")
      
    }
    Accuracy_matrix[Acc_index, (2+i)] <- Acc
    
  }
  
}

# assign current value of P and Q 
Accuracy_matrix$P[Acc_index] <- p
Accuracy_matrix$Q[Acc_index] <- q

cat("\n")
print("Accuracy calculated using all test subjects on different values of top n diseases and varing P and Q:")
cat("\n")
print("Accuracy Matrix (Grid Search):")
print(Accuracy_matrix[1:Acc_index, ])
cat("\n")

# adding average acc to the Accuracy_matrix
Accuracy_matrix <- Accuracy_matrix %>% mutate(Avg_Acc = 0)
Accuracy_matrix$Avg_Acc <- rowMeans(Accuracy_matrix[,c(3:8)])

# create file name of accuracy matrix to write data into csv file
file_name_1 <- paste("Accuracy_matrix_holdout_test_Sub",l,"_p_",p,"_q_",q,".csv", collapse = "",sep="")
write.csv(Accuracy_matrix, file = file_name_1, row.names = FALSE)

# increamentor 
Acc_index <- Acc_index +1

print("################################## Current itteration of the loop for Q ends here ##########################################")
cat("\n")

#  } # ending loop q 

print("################################## Current itteration of the loop for P ends here ##########################################")
cat("\n")

#} # ending loop p

print("Accuracy Matrix (Test Data):")
print(Accuracy_matrix)
