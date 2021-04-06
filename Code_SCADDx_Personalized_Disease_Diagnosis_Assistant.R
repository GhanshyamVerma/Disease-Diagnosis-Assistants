################ SCADDx: A Personalized Disease Diagnosis Assistant ######################

# load packages
library(caret) # for various machine learning functions
library(dplyr) # for data manupulation
library(e1071) # for various functions like confusion matrix
library(ggplot2) # for plots


############################# input from user #################################################

# read CTD knowlege graph that contains gene disease links
KG_ctd_gene_disease <- read.csv(file = "KG_ctd_gene_disease.csv", header = TRUE, sep=",")

# read Gene expression data collected at two required time-points for all those subjects for which disease diagnosis need to be performed
# required two time-points (0 hours i.e. healthy state and at time-point T i.e. diseased state or the time-point at which disease diagnosis has been requested)
Gene_expression_data_full <- read.csv(file = "Gene_expr_data_all_Sub_all_features.csv", header = TRUE, sep=",", check.names = FALSE)

# gene start index Gene_expression_data_full dataframe
s_index <- 10

# read data that contains 3 columns (1. unique disease id of all the diseases in KG, 2. unique disease name of all the diseases in KG and 3. Disease weights initialised with 0) 
Data_Unique_Disease <- read.csv(file = "Data_Unique_Disease.csv", header = TRUE, sep=",", check.names = FALSE)

# how many top gene you want for assigning weights to the diseases in KG 

# rang of top P genes and bottom Q genes 
P_start <- 10
P_end <- 500 
P_step <- 10
Q_start <- 10
Q_end <- 500
Q_step <- 10

#################################################################################################

# removing duplicaate genes if any
Gene_expression_data_full_no_duplicate <- Gene_expression_data_full[, !duplicated(colnames(Gene_expression_data_full), fromLast = TRUE)]

# total number of subjects
s <- Gene_expression_data_full_no_duplicate[dim(Gene_expression_data_full_no_duplicate)[1] , "Super_Subject_ID"]

# gene end index for the dataframe
e_index <- dim(Gene_expression_data_full_no_duplicate)[2]

# performing 0-1 normalization
for(i in 1:s){  # loop of i for subject wise normalisation 
  Gene_expression_data_full_no_duplicate[Gene_expression_data_full_no_duplicate$Super_Subject_ID == i, s_index:e_index] <- ((Gene_expression_data_full_no_duplicate[Gene_expression_data_full_no_duplicate$Super_Subject_ID == i, s_index:e_index] - min(Gene_expression_data_full_no_duplicate[Gene_expression_data_full_no_duplicate$Super_Subject_ID == i, s_index:e_index])) / ( max(Gene_expression_data_full_no_duplicate[Gene_expression_data_full_no_duplicate$Super_Subject_ID == i, s_index:e_index]) - min(Gene_expression_data_full_no_duplicate[Gene_expression_data_full_no_duplicate$Super_Subject_ID == i, s_index:e_index])))
}

# extract initial columns that do not contains genes columns 
Gene_expression_data_inicial_columns <- Gene_expression_data_full_no_duplicate[ , c(1:(s_index-1))]


# find common genes in KG and gene expression data
common_gene_KG_and_exp_data <- intersect(KG_ctd_gene_disease$GeneSymbol, names(Gene_expression_data_full_no_duplicate[ ,-c(1:(s_index-1))]))


# taking only genes 
Gene_expression_data_only_genes <- Gene_expression_data_full_no_duplicate[ , -c(1:(s_index-1))]

# keep only common genes
Gene_expression_data_full_no_duplicate <- Gene_expression_data_only_genes[ , colnames(Gene_expression_data_only_genes) %in% common_gene_KG_and_exp_data]


# combine first initial columns again which has information except genes like label, title, time point, subject id etc.
Gene_expression_data_full_no_duplicate <- cbind(Gene_expression_data_inicial_columns, Gene_expression_data_full_no_duplicate)

# temporary dataframe
Data_Unique_Disease_initial_weights <- Data_Unique_Disease


#### Code for assigning weights to the diseases in KG based on changes observed in genes

PS <- length(seq(P_start,P_end,P_step))
QS <- length(seq(Q_start,Q_end,Q_step))

Accuracy_matrix <- data.frame("P" = 1:(PS*QS), "Q" = 1:(PS*QS), "Acc_Top_1_Dis" = 1:(PS*QS), "Acc_Top_2_Dis" = 1:(PS*QS), "Acc_Top_3_Dis" = 1:(PS*QS), "Acc_Top_4_Dis" = 1:(PS*QS), "Acc_Top_5_Dis" = 1:(PS*QS), "Acc_Top_10_Dis" = 1:(PS*QS), "Max_Score" = 1:(PS*QS), "Min_Score" = 1:(PS*QS))
Acc_index <- 1


# loop starts
for(p in seq(P_start,P_end,P_step)){  # loop of p is for top genes 
  
  for(q in seq(Q_start,Q_end,Q_step)){  # loop of q is for bottom genes 
    # for all P and Q values
    # Data_Unique_Disease_all_times <- Data_Unique_Disease_initial_weights
    
    G_Max_Score <- 0
    G_Min_Score <- 0
    
    
    # total number of subjects
    s <- Gene_expression_data_full_no_duplicate[dim(Gene_expression_data_full_no_duplicate)[1] , "Super_Subject_ID"]
    
    predicted_info <- Gene_expression_data_full_no_duplicate[ , c(1:(s_index-1))]
    predicted_info <- predicted_info %>% mutate(predicted_label_top_1 = 1:(2*s))
    predicted_info <- predicted_info %>% mutate(predicted_label_top_2 = 1:(2*s))
    predicted_info <- predicted_info %>% mutate(predicted_label_top_3 = 1:(2*s))
    predicted_info <- predicted_info %>% mutate(predicted_label_top_4 = 1:(2*s))
    predicted_info <- predicted_info %>% mutate(predicted_label_top_5 = 1:(2*s))
    predicted_info <- predicted_info %>% mutate(predicted_label_top_10 = 1:(2*s))
    
    Gene_Data_All_ti_prediction <- Gene_expression_data_full_no_duplicate[Gene_expression_data_full_no_duplicate$Time_Point_Adjusted == 72, c(1:(s_index-1))]
    
    Gene_Data_All_ti_prediction <- Gene_Data_All_ti_prediction %>% mutate(predicted_label_top_1 = 1:s)
    Gene_Data_All_ti_prediction <- Gene_Data_All_ti_prediction %>% mutate(predicted_label_top_2 = 1:s)
    Gene_Data_All_ti_prediction <- Gene_Data_All_ti_prediction %>% mutate(predicted_label_top_3 = 1:s)
    Gene_Data_All_ti_prediction <- Gene_Data_All_ti_prediction %>% mutate(predicted_label_top_4 = 1:s)
    Gene_Data_All_ti_prediction <- Gene_Data_All_ti_prediction %>% mutate(predicted_label_top_5 = 1:s)
    Gene_Data_All_ti_prediction <- Gene_Data_All_ti_prediction %>% mutate(predicted_label_top_10 = 1:s)
    
    All_Sub_temp_prediction <- data.frame("Top_1"=1:s, "Top_2"=1:s,"Top_3"=1:s, "Top_4"=1:s,"Top_5"=1:s, "Top_10"=1:s)
    
    for(l in 1:s){ # loop of l for number of subjects
      
      Gene_expression_data_sub_l <- Gene_expression_data_full_no_duplicate %>% filter(Super_Subject_ID == l)
      Gene_expression_data_sub_l <- Gene_expression_data_sub_l[ , -c(1:(s_index-1))]
      
      print("########################## New subject computation start from here ############################")
      print("Subject id is:")
      print(l)
      
      # for each time point initialize again
      Data_Unique_Disease <- Data_Unique_Disease_initial_weights
      
      # Ti - T1 (infected - healthy)
      Gene_Transition_Matrix <- Gene_expression_data_sub_l[2, ] - Gene_expression_data_sub_l[1, ]
      
      # print("show value of Gene_Transition_Matrix:")
      # print(Gene_Transition_Matrix[, 1:20])
      
      Gene_Transition_Matrix_top_p_Genes <- Gene_Transition_Matrix[ , order(-abs(Gene_Transition_Matrix[ , ]))]
      Gene_Transition_Matrix_top_p_Genes <- Gene_Transition_Matrix_top_p_Genes[ , c(1:p)]
      
      
      print("show values of top 5 Gene_Transition_Matrix_top_p_Genes:")
      print(Gene_Transition_Matrix_top_p_Genes[ , 1:5])
      
      # code for assigning reward based on top P genes
      for(j in 1:p){ # loop of j for number of genes
        Disease_IDs <- KG_ctd_gene_disease[KG_ctd_gene_disease$GeneSymbol == names(Gene_Transition_Matrix_top_p_Genes)[j], "Disease_ID" ]
        for(k in 1:length(Disease_IDs)){ # loop for every disease id
          Data_Unique_Disease[Data_Unique_Disease$Disease_ID ==  Disease_IDs[k], "Disease_Weight"] <-   Data_Unique_Disease[Data_Unique_Disease$Disease_ID ==  Disease_IDs[k], "Disease_Weight"] + abs(Gene_Transition_Matrix_top_p_Genes[,j])
        }  # end for loop k
      } # end for loop j
     
      # code for penalty based on bottom Q genes
      Gene_Transition_Matrix_bottom_q_Genes <- Gene_Transition_Matrix[ , order(abs(Gene_Transition_Matrix[ , ]))]
      Gene_Transition_Matrix_bottom_q_Genes <- Gene_Transition_Matrix_bottom_q_Genes[ , c(1:q)]
      
      print("show values of Gene_Transition_Matrix_bottom_q_Genes:")
      print(Gene_Transition_Matrix_bottom_q_Genes[ , 1:5])
      
      for(j in 1:q){ # loop of j for number of bottom genes
        Disease_IDs <- KG_ctd_gene_disease[KG_ctd_gene_disease$GeneSymbol == names(Gene_Transition_Matrix_bottom_q_Genes)[j], "Disease_ID" ]
        for(k in 1:length(Disease_IDs)){ # loop for every disease id
          Data_Unique_Disease[Data_Unique_Disease$Disease_ID ==  Disease_IDs[k], "Disease_Weight"] <-  Data_Unique_Disease[Data_Unique_Disease$Disease_ID ==  Disease_IDs[k], "Disease_Weight"] - (abs(Gene_Transition_Matrix_top_p_Genes[,1]) - abs(Gene_Transition_Matrix_bottom_q_Genes[,j]))
        }  # end for loop k
      } # end for loop j
      
      Max_Score <- Data_Unique_Disease[order(-Data_Unique_Disease$Disease_Weight),  "Disease_Weight"][1]
      Min_Score <- Data_Unique_Disease[order(Data_Unique_Disease$Disease_Weight),  "Disease_Weight"][1]
      
      # create file name to write data into csv file
      file_name <- paste("Disease_Weight_Sub_",l,"_p_",p,"_q_",q,".csv", collapse = "",sep="")
      
      # write data into csv file
      write.csv(Data_Unique_Disease[order(-Data_Unique_Disease$Disease_Weight), ], file = file_name, row.names = FALSE)
      
      print("Subject id is:")
      print(l)
      
      print("Value of p is :")
      print(p)
      
      print("Value of q is :")
      print(q)
      
      print("This subject at this time point has following label:")
      print(Gene_expression_data_full_no_duplicate[ Gene_expression_data_full_no_duplicate$Super_Subject_ID == l , ]$Label[2])
      
      for(i in 1:6){ # loop for how many top disease you want to look for Acc calc
        
        if(i<6){
          Top_Disease_Names <- Data_Unique_Disease[order(-Data_Unique_Disease$Disease_Weight), ][1:i, "Disease_Name"]
        }else{
          Top_Disease_Names <- Data_Unique_Disease[order(-Data_Unique_Disease$Disease_Weight), ][1:10, "Disease_Name"]
        }
        
        if(any(Top_Disease_Names == "Respiratory Viral Infection")){
          predicted_info[predicted_info$Super_Subject_ID == l , (i+s_index-1)][2] <- "RVI"
          Gene_Data_All_ti_prediction[Gene_Data_All_ti_prediction$Super_Subject_ID == l , (i+s_index-1)][1] <- "RVI"
          All_Sub_temp_prediction[l,i] <- "RVI"
        }else{
          predicted_info[predicted_info$Super_Subject_ID == l , (i+s_index-1)][2] <- "NOT RVI"
          Gene_Data_All_ti_prediction[Gene_Data_All_ti_prediction$Super_Subject_ID == l , (i+s_index-1)][1] <- "NOT RVI"
          All_Sub_temp_prediction[l,i] <- "NOT RVI"
        }
        print(paste("Predicted label using top ", i, "disease is:"))
        print(All_Sub_temp_prediction[l,i])
        
        
      }
      
      
      print("Top 10 Diseases are:")
      print(Data_Unique_Disease[order(-Data_Unique_Disease$Disease_Weight), ][1:10, ])
      print("Top 5 Genes are:")
      print(Gene_Transition_Matrix_top_p_Genes[ ,1:5])
      
      print("############################################################################")
      
      if(G_Min_Score > Min_Score){
        G_Min_Score <- Min_Score
      }
      if(G_Max_Score < Max_Score){
        G_Max_Score <- Max_Score
      }
      
    } # end for loop l
    
    # create file name to write data into csv file
    file_name_1 <- paste("predicted_info_60hr_Sub_",l,"_p_",p,"_q_",q,".csv", collapse = "",sep="")
    write.csv(Gene_Data_All_ti_prediction, file = file_name_1, row.names = FALSE)
    
    if(any(Gene_Data_All_ti_prediction$Label == "Not RVI")){
      for(i in 1:6){
        confusion_mat <- confusionMatrix( as.factor(All_Sub_temp_prediction[,i]), as.factor(Gene_Data_All_ti_prediction$Label), positive = "RVI")
        print(paste("Accuracy using top", i, "disease is:"))
        print(confusion_mat)
    
        Accuracy_matrix[Acc_index, (2+i)] <- confusion_mat$overall[1]
      }
  
    }
    else{
      for(i in 1:6){
        hit <- 0
          for(k in 1:dim(All_Sub_temp_prediction)[1]){
            if(Gene_Data_All_ti_prediction[k,"Label"] == All_Sub_temp_prediction[k,i]){
              hit <- hit + 1
            }
          }
        Acc <- hit/dim(All_Sub_temp_prediction)[1]
        print(paste("Accuracy using top", i, "disease is:"))
        print(Acc)
    
        Accuracy_matrix[Acc_index, (2+i)] <- Acc
      }
    }
    
    Accuracy_matrix$P[Acc_index] <- p
    Accuracy_matrix$Q[Acc_index] <- q
    Accuracy_matrix$Max_Score[Acc_index] <- G_Max_Score
    Accuracy_matrix$Min_Score[Acc_index] <- G_Min_Score
    
    Acc_index <- Acc_index +1
    
    # create file name to write data into csv file
    file_name_1 <- paste("Accuracy_matrix_60hr_Sub",l,"_p_",p,"_q_",q,".csv", collapse = "",sep="")
    write.csv(Accuracy_matrix, file = file_name_1, row.names = FALSE)
    
    print("Accuracy Matrix:")
    print(Accuracy_matrix[1:Acc_index, ])
    
  } # ending loop q 
} # ending loop p
