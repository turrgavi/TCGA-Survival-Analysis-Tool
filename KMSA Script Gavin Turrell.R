#Gavin Turrell

#Import all required packages
library("TCGAbiolinks")
library("stringr")
library("dplyr")
library("survival")
library("ggplot2")
library("survminer")
library(cowplot)

#Set cutoff of survival analysis in days
survival_cutoff <-1825

#Set directory where expression data is located
setwd("/Users/gavinturrell/Documents")

#Load expression data
print("Loading in expression data....")
raw_expression_data <- as.matrix(read.csv("TCGABLCAexp.csv", header = FALSE))

#Open file prompt to select gene comparison list
print("Select comparison gene list as CSV format")
filename <- file.choose()
desired_gene <- as.matrix(read.csv(filename, header = FALSE))
name <- basename(filename)

#Remove duplicate genes in pathway list
desired_gene <- unique(desired_gene)
pathway_name <- str_sub(name, end=-5)

#Create a matrix to hold all transformed data
test_buffer = matrix(data = NA, (nrow=length(desired_gene[,1])), ncol=length(raw_expression_data[,1]))

#Match genes from selected list to genes found in expression data
test_buffer <- merge(desired_gene, raw_expression_data, by = "V1")
test_buffer <- rbind(raw_expression_data[1,], test_buffer)

#Removes all non-Primary Tumor Samples by excluding all barcodes that do NOT end in "-01"
run_time_ptf = length(test_buffer[1,])
p = 1
t = 1
test_buffer <- na.omit(test_buffer)

for (i in 1:run_time_ptf){
  p=p+1
  
  if (str_sub(test_buffer[1,p], -4) != "-01A" & 
      str_sub(test_buffer[1,p], -4) != "-01B" &
      str_sub(test_buffer[1,p], -4) != "-01C")
  {
    test_buffer <- test_buffer[,-p]
  }

  if (p == length(test_buffer[1,])){
    break
  }
}

#Create a matrix for holding upper and lower quartile calculations for each gene in the test_buffer matrix
quart <- matrix(data = NA, ncol= 3, nrow = length(test_buffer[,1]))

#Name columns by gene name, LQ for lower quartile, UQ for upper quartile
quart[,1] <- test_buffer[,1]
quart[1,2] <- "LQ"
quart[1,3] <- "UQ"

#For loop for calculating the UQ and LQ for each gene present in test_buffer and adding it to the appropriate cell in the "quart" matrix
run_time_quart = length(test_buffer[,1])
p = 1
q = 1
o = 1
for (i in 1:run_time_quart){
  p <- p+1
  
  if (p>length(test_buffer[,1])){
    break()
  }
  
  LQ <- round(quantile(as.numeric(test_buffer[p,2:length(test_buffer[2,])]), 0.25, na.rm = TRUE), 2)
  UQ <- round(quantile(as.numeric(test_buffer[p,2:length(test_buffer[2,])]), 0.75, na.rm = TRUE), 2)
  
  quart[p,2] <- LQ
  quart[p,3] <- UQ
  
}

#For loop to cycle through all each gene's raw expression data in test_buffer and label as UQ or LQ according to the quart matrix values and gene names
#Each gene will have unique UQ and LQ values 
#Some cells will not fall into UQ or LQ category and will not be labelled
patient_countLQ = 0
patient_countUQ = 0
v = 1
sample_count = 1
row_count = 2
for (sample_count in 1:length(test_buffer[,1])){
  for (q in 1:length(test_buffer[1,])){
    v=v+1
    if(row_count > length(test_buffer[,1])){
      break()
    }

    if (v > length(test_buffer[1,])){
      v = 2
    }
    
    if (round(as.numeric(test_buffer[row_count,v]),4) < as.numeric(quart[row_count,2]) | round(as.numeric(test_buffer[row_count,v]),4) == 0){
      test_buffer[row_count,v] <- "LQ"
      patient_countLQ = patient_countLQ+1
    }
  
    if (round(as.numeric(test_buffer[row_count,v]),4) > round(as.numeric(quart[row_count,3]),4) &
       test_buffer[row_count,v] != "LQ" &
       test_buffer[row_count,v] !=  test_buffer[row_count,1]){
      
       test_buffer[row_count,v] <- "UQ"
       patient_countUQ = patient_countUQ+1
    }
    
    if (v == length(test_buffer[1,])){
      row_count = row_count+1
    }
    }
  cat("Assigning Quartiles... % Complete: ", round((sample_count-1)/length(test_buffer[,1])*100,2), "%\n\n")
}

#Create a matrix to hold P-values, expression quartiles and patient count for upper and lower quartiles
pval <- matrix(data=NA, ncol=5, nrow=length(test_buffer[,1])) 
pval[,1] <- test_buffer[,1]
pval[1,2] <- "P-value"
pval[1,3] <- "Decreased Survival Expression Quartile"
pval[1,4] <- "UQ Total"
pval[1,5] <- "LQ Total"
gene_inc <-1

#This is the main KMSA loop that runs for every matched gene
for (gene_inc in 1:length(test_buffer[,1])){
  
  #Create a matrix called "dataf" for loading in only UQ and LQ labelled data from test_buffer
  dataf <- matrix(data = NA, ncol = 3, nrow=length(test_buffer[1,]))
  
  #All test_buffer data is transfered to "dataf"
  dataf[,1] <- t(test_buffer[1,])
  dataf[,2] <- t(test_buffer[gene_inc+1,])
  
  #All non-LQ or non-UQ values are removed from the list and added to a new object called dataf_sub
  dataf_sub <- subset(dataf, dataf[,2] == "LQ" | dataf[,2] == "UQ")
  
  #Sys.sleep(3)
  #PLEASE NOTE THERE IS AN ISSUE WITH LOADING IN THIS DATA SEQUENTIALLY, MAY NEED TO SLEEP() PRIOR TO THIS
  clinical_data_all <- read.csv("TCGA-BLCA.GDC_phenotype.csv", header = TRUE)
  
  #Save clinical data as matrix for easier handling into new data sets
  clinical_data_all <- as.matrix(clinical_data_all)
  
  #Remove NAs from clinical dataset
  clinical_data_all <- subset(clinical_data_all, !is.na(clinical_data_all[,2]))
  
  #Create new dataset for survival analysis formatting
  survival_data_all <- matrix(data = NA, nrow=nrow(clinical_data_all), ncol= 4)
  
  #Set first column of survival_data_all to be list of submitter_id from clinical_data_all
  survival_data_all[,1] <- clinical_data_all[,1]
  survival_data_all[,2] <- as.numeric(clinical_data_all[,"days_to_death.demographic"])
  survival_data_all[,3] <- as.numeric(clinical_data_all[,"days_to_last_follow_up.diagnoses"])
  survival_data_all[,4] <- clinical_data_all[,"vital_status.demographic"]
  
  #For loop to replace days_to_death values with days_to_last_follow_up values if days_to_death is NA and days_to_last_follow_up is NOT NA
  i = 1
  for (i in 1:length(survival_data_all[,3])){
    if (!is.na(survival_data_all[i,3]) & is.na(survival_data_all[i,2])){
      survival_data_all[i,2] <- survival_data_all[i,3]
    }
  }
  
  #For loop for assigning a status to each pateint based on their days_to_death of days_to_last_follow_up
  i = 1
  for (i in 1:length(survival_data_all[,3])){
    
    if (survival_data_all[i,4] == "Alive"){
      survival_data_all[i,3] <- 0
    }
    
    if (survival_data_all[i,4] == "Dead"){
      survival_data_all[i,3] <- 1
    }
  }
  
  i = 1
  for (i in 1:length(survival_data_all[,1])){
    print(survival_data_all[i,2])
    if (as.numeric(survival_data_all[i,2]) > survival_cutoff){
      survival_data_all[i,3] <- 0
      survival_data_all[i,2] <- survival_cutoff
    }
  }
  
  #create a dataframe merging survival_data_all and dataf_sub by the IDs
  #This creates a dataframe containing all patients with a UQ or LQ label alongside their time and status and quartile label 
  merged_dataf <- merge(survival_data_all, dataf_sub, by = "V1")
  merged_dataf <- subset(merged_dataf, as.numeric(as.character(merged_dataf[,2])) <= survival_cutoff)
  
  #This is the actual KMSA function which assigns the ouput to the "fit" variable
  fit <- survfit(Surv(as.numeric(as.vector(merged_dataf$V2.x)),as.numeric(merged_dataf$V3.x))~merged_dataf$V2.y, data = merged_dataf)
  
  #Perform a log rank test called diff_fit using  time, status and quartile label as a factor using rho=0 as advised by xenabrowser
  try({ #TRY function used to continue loop following error
    diff_fit = survdiff(Surv(as.numeric(as.vector(merged_dataf$V2.x)),as.numeric(merged_dataf$V3.x))~merged_dataf$V2.y, rho=0)
    
    i = 1
    LQ_total <- sum(merged_dataf[,5] == "LQ")
    UQ_total <- sum(merged_dataf[,5] == "UQ")
    LQ_deaths = 0
    UQ_deaths = 0
    for (i in 1:length(merged_dataf[,1])){
      if (merged_dataf[i,4] == "Dead" & merged_dataf[i,5] == "LQ"){
        LQ_deaths = LQ_deaths+1
      }
      
      if (merged_dataf[i,4] == "Dead" & merged_dataf[i,5] == "UQ"){
        UQ_deaths = UQ_deaths+1
      }
      
    }
    
    if (LQ_total-LQ_deaths < UQ_total-UQ_deaths){
      pval[gene_inc+1,3] <- "Lower"
    }
    
    if (LQ_total-LQ_deaths > UQ_total-UQ_deaths){
      pval[gene_inc+1,3] <- "Upper"
    }
    
    pval[gene_inc+1,4] <- UQ_total
    pval[gene_inc+1,5] <- LQ_total
    #return only the p-value from the Kaplan Meier analysis
    pval[gene_inc+1,2] <- round(pchisq(diff_fit$chisq, length(diff_fit$n)-1, lower.tail = FALSE),5)
    gene_inc <- gene_inc+1
    
    cat("Gene: ", test_buffer[gene_inc,1], " analyzed\n")
    cat("Percentage Complete: ", round((gene_inc)/(length(test_buffer[,1]))*100,2), "%\n\n")
    
    if (gene_inc == length(test_buffer[,1])){
      break()
    }
  }
  ) #End of Try function
  
}


#pval <- subset(pval, pval[,2] < 0.05 | pval[,2] == "P-value" )
write.table(pval, file = paste0(pathway_name, "_Pvalues.csv"), sep = ",", col.names = FALSE, row.names = FALSE)
cat("All Genes analyzed and written to file [List Name]_Pvalues.csv")

ggsurvplot(fit, 
           conf.int = F, 
           pval = TRUE, 
           risk.table = FALSE, 
           ggtheme = theme_bw(), 
           xlim=c(0,1825), 
           break.time.by=365, 
           font.legend = c(9),
           legend = c("top"),
           legend.labs = c("Lower", "Upper"),
           legend.title = c("Quartile: "),
           font.x = c(12), 
           font.y = c(12),
           font.tickslab = c(9))
