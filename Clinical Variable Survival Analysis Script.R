#Gavin Turrell

#Import all required packages
#library("TCGAbiolinks")
library("stringr")
library("dplyr")
library("survival")
library("ggplot2")
library("survminer")
#library(cowplot)

#Set cutoff of survival analysis in days
survival_cutoff <-1825

#Read config file
config_file <- read.csv("Clinical Configuration File.csv", header = FALSE)
#What variables do you want to analyse? [TRUE/FALSE]
# check_bmi = config_file[1, 2]
# check_smoking_status = config_file[2, 2]
# check_age_at_diagnosis = config_file[3,2]
# check_stage = config_file[4,2]
# check_stage_break = config_file[5,2]

#Only have one variable true at a time, aside from check_stage and check_stage_break
check_bmi = FALSE
check_smoking_status = FALSE
check_age_at_diagnosis = TRUE
check_stage = FALSE
check_stage_break = TRUE #if true, check_stage must also be true

#set the working directory for files to be saved and loaded from
#setwd("/Users/gavinturrell/Tresors/Gavin's tresor/University of Queensland/Research/Dr Janin Chandra/R/TCGA-CESC R/Post Analysis Results")
setwd("C:/Users/Gavin Turrell/My Tresors/Gavin's tresor/University of Queensland/Research/Dr Janin Chandra/R/TCGA-CESC R/Post Analysis Results")


clinical_data_all <- read.csv("TCGA-CESC_clinical.csv", header = TRUE)

#create vectors of clinical variables
age_years <- clinical_data_all[,"age_at_diagnosis"]/365
bmi <- clinical_data_all[,"bmi"]
smoking_status <- clinical_data_all[,"pack_years_smoked"]


#Define quartiles for clinical variables
LQ_age <- quantile(age_years, 0.25, na.rm = TRUE)
UQ_age <- quantile(age_years, 0.75, na.rm = TRUE)
# LQ_bmi <- quantile(bmi, 0.25, na.rm = TRUE)
# UQ_bmi <- quantile(bmi, 0.75, na.rm = TRUE)
UW_bmi <- 18.5
OW_bmi <- 25
Obese_bmi <- 30
LQ_smoking_status <- quantile(smoking_status, 0.25, na.rm = TRUE)
UQ_smoking_status <- quantile(smoking_status, 0.75, na.rm = TRUE)

#Create matrix for final analysis to read from
survival_data = matrix(data = NA, nrow=length(clinical_data_all[,1]), ncol=5)


colnames(survival_data) <- c("Status", "Event", "Quartile", "death", "fu")

if (check_bmi == TRUE){
  colnames(survival_data) <- c("Status", "Event", "BMI", "death", "fu")
}
  



#1 = Alive, 2 = dead
survival_data[,1] <- clinical_data_all$vital_status
survival_data[,4] <- clinical_data_all$days_to_death
survival_data[,5] <- clinical_data_all$days_to_last_follow_up

#Loop for assigning quartile labels to patients
i = 1
LQ_count = 0
UQ_count = 0

if (check_age_at_diagnosis == TRUE){
  for (i in 1:length(clinical_data_all[,1])){
    
    if (!is.na(clinical_data_all[i,"age_at_diagnosis"]) && as.numeric(clinical_data_all[i,"age_at_diagnosis"])/365 < as.numeric(LQ_age)){
      survival_data[i,"Quartile"] <- "LQ"
      LQ_count = LQ_count+1
    }
    
    if (!is.na(clinical_data_all[i,"age_at_diagnosis"]) && as.numeric(clinical_data_all[i,"age_at_diagnosis"])/365 > as.numeric(UQ_age)){
      survival_data[i,"Quartile"] <- "UQ"
      UQ_count = UQ_count+1
    }
    
    print(i)
    i = i+1
  }
}

if (check_bmi == TRUE){

  for (i in 1:length(clinical_data_all[,1])){
    
    if (!is.na(clinical_data_all[i,"bmi"]) && as.numeric(clinical_data_all[i,"bmi"]) < as.numeric(UW_bmi)){
      survival_data[i,"BMI"] <- "UW"
      
    }
    
    if (!is.na(clinical_data_all[i,"bmi"]) && as.numeric(clinical_data_all[i,"bmi"]) > as.numeric(Obese_bmi)){
      survival_data[i,"BMI"] <- "O"
     
    }
    
    if (!is.na(clinical_data_all[i,"bmi"]) && 
        as.numeric(clinical_data_all[i,"bmi"]) < as.numeric(Obese_bmi) &&
        as.numeric(clinical_data_all[i,"bmi"]) > as.numeric(OW_bmi)){
      survival_data[i,"BMI"] <- "OW"
    }
    
    if (!is.na(clinical_data_all[i,"bmi"]) && 
        as.numeric(clinical_data_all[i,"bmi"]) < as.numeric(OW_bmi) &&
        as.numeric(clinical_data_all[i,"bmi"]) > as.numeric(UW_bmi)){
      survival_data[i,"BMI"] <- "NW"
    }
    print(i)
    i = i+1
  }
}

if (check_smoking_status == TRUE){
  
  for (i in 1:length(clinical_data_all[,1])){
    
    if (!is.na(clinical_data_all[i,"pack_years_smoked"]) && as.numeric(clinical_data_all[i,"pack_years_smoked"]) < as.numeric(LQ_smoking_status)){
      survival_data[i,"Quartile"] <- "LQ"
      LQ_count = LQ_count+1
    }
    
    if (!is.na(clinical_data_all[i,"pack_years_smoked"]) && as.numeric(clinical_data_all[i,"pack_years_smoked"]) > as.numeric(UQ_smoking_status)){
      survival_data[i,"Quartile"] <- "UQ"
      UQ_count = UQ_count+1
    }
    
    print(i)
    i = i+1
  }
}

if (check_stage == TRUE){
  survival_data[,"Quartile"] <- as.character(clinical_data_all[,"figo_stage"])
  survival_data[,"Quartile"] <- gsub("[A-C]", "", survival_data[,"Quartile"])
  survival_data[,"Quartile"] <- gsub("[1-2]", "", survival_data[,"Quartile"])
  
  if (check_stage_break == TRUE){
    i = 2 
    for (i in 1:length(survival_data[,1])){
      if (!is.na(survival_data[i,"Quartile"]) && as.character(survival_data[i,"Quartile"]) != "Stage I"){
        survival_data[i,"Quartile"] <- "Stage II+"
      }
      i = i+1
      print(i)
    }
  }
}
#For loop for assigning event time (days)
i = 2

for (i in 1:length(survival_data[,1])){
  
  if (is.na(survival_data[i,"death"])){
    survival_data[i,"Event"] <- survival_data[i,"fu"]
  }
  
  if (!is.na(survival_data[i, "death"]) && is.na(survival_data[i,"fu"])){
    survival_data[i,"Event"] <- survival_data[i,"death"]
  }
  
  if (!is.na(survival_data[i, "death"]) && !is.na(survival_data[i,"fu"])){
    
    if (as.numeric(survival_data[i,"death"]) > as.numeric(survival_data[i,"fu"])){
      survival_data[i,"Event"] <- survival_data[i,"death"]
    }
    
    if (as.numeric(survival_data[i,"death"]) <= as.numeric(survival_data[i,"fu"])){
      survival_data[i,"Event"] <- survival_data[i,"fu"]
    }
  }
  i = i+1
  print(i)
}

i = 2
for (i in 1:length(survival_data[,1])){
  if (as.numeric(survival_data[i,"Event"]) > survival_cutoff){
    survival_data[i,"Status"] <- 1
    survival_data[i,"Event"] <- survival_cutoff
  }
  
}

#Remove all non-quartile assigned patients and create cutoff
if (check_bmi == FALSE){
  survival_data <- subset(survival_data, !is.na(survival_data[,"Quartile"]))
  #survival_data <- subset(survival_data, survival_data[,"Event"] <= survival_cutoff)
}

if (check_bmi == TRUE){
  survival_data <- subset(survival_data, !is.na(survival_data[,"BMI"]))
}


#convert to dataframe as required by the Surv package
survival_data <- as.data.frame(survival_data)

if (check_stage == TRUE){
  colnames(survival_data) <- c("Status", "Event", "FIGO_Stage", "death", "fu")
  fit <- survfit(Surv(as.numeric(as.vector(survival_data$Event)),as.numeric(survival_data$Status))~survival_data$FIGO_Stage, data = survival_data)
}

#Perform SA analysis
if (check_stage != TRUE & check_bmi != TRUE){
  fit <- survfit(Surv(as.numeric(as.vector(survival_data$Event)),as.numeric(survival_data$Status))~survival_data$Quartile, data = survival_data)
}


if (check_bmi == TRUE){
  fit <- survfit(Surv(as.numeric(as.vector(survival_data$Event)),as.numeric(survival_data$Status))~survival_data$BMI, data = survival_data)
}

#Plot the results of the analysis (very basic plot for now but incudes necessary information)

if (check_stage == TRUE & check_stage_break == FALSE){
  ggsurvplot(fit, 
             conf.int = F, 
             pval = TRUE, 
             risk.table = FALSE, 
             ggtheme = theme_bw(), 
             xlim=c(0,1825), 
             break.time.by=365, 
             font.legend = c(9),
             legend = c("top"),
             legend.labs = c("I", "II", "III", "IV"),
             legend.title = c("FIGO Stage:"),
             font.x = c(12), 
             font.y = c(12),
             font.tickslab = c(9))
}

if (check_bmi == TRUE){
  ggsurvplot(fit, 
             conf.int = F, 
             pval = TRUE, 
             risk.table = FALSE, 
             ggtheme = theme_bw(), 
             xlim=c(0,1825), 
             break.time.by=365, 
             font.legend = c(9),
             legend = c("top"),
             legend.labs = c("Normal", "Obese", "Overweight", "Underweight"),
             legend.title = c("BMI: ", "bold"),
             font.x = c(12), 
             font.y = c(12),
             font.tickslab = c(9))
}

if (check_bmi == FALSE & check_stage == FALSE){
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
}