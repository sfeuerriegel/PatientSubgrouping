
rm(list = ls())
library(dplyr)
library(tidyr)
library(foreign)
library(psych)
library(forcats)
library(caret)
#First, we exclude all patients only completing less than 10 of 52 weekly responses 
#(76 out of 928) and additionally excluded all patients whose pain intensity levels had a standard deviation of 
#0 (5 out of 852), leaving us with a total of 847 patients.



lbp.data <- as.data.frame(read.dta(file = "../../../Data/casey_master_11.dta"))

lbp.intensity <- lbp.data[, paste0("nrs", 1:52)]

lbp.data$NA_NUM <- rowSums(is.na(lbp.intensity))

lbp.data$SD <- sapply(1:nrow(lbp.data),function(k){sd(lbp.intensity[k,],na.rm = TRUE)})

lbp.data_filtered <- lbp.data %>% 
  filter (NA_NUM <= 42) %>%
  filter (SD != 0) %>%
  droplevels()
row.names(lbp.data_filtered)<-seq(nrow(lbp.data_filtered))

nrow(lbp.data[lbp.data$NA_NUM>42,])
nrow(lbp.data[lbp.data$SD==0 & lbp.data$NA_NUM<=42,])
nrow(lbp.data_filtered)
rm(lbp.data,lbp.intensity)

#Next, we recode the variables to be potentially used as risk factors in subsequent analysis

lbp.data_filtered <- lbp.data_filtered %>%
  mutate(
    dlva0 = fct_recode(
      dlva0 ,
      "0-2 weeks" = "0-2" ,
      "2-4 weeks" = "2-4" ,
      "1-3 months" = "1_månde-3_måneder" ,
      "More than 3 months" = "mere_end_3_måneder"
    )
  ) %>%
  mutate(bsex0 = fct_recode(bsex0 ,
                            "Female" = "female" ,
                            "Male" = "male")) %>%
  mutate(
    nootherdisease = fct_recode(
      nootherdisease ,
      "No other disease" = "ja" ,
      "Other chronic disease" = "nej"
    )
  ) %>%
  mutate(
    bryg0 = fct_recode(
      bryg0 ,
      "Never smoked" = "ikke_ryger" ,
      "Ex-smoker" = "eks_ryger" ,
      "Smoker" = "ryger"
    )
  ) %>%
  mutate(
    tlda0 = fct_recode(
      tlda0 ,
      "Less than 30 days" = "30_dage_eller_færre" ,
      "More than 30 days" = "over_30_dage"
    )
  )%>%
  mutate(
    domin_bp = fct_recode(
      domin_bp ,
      "Back pain dominating" = "nej" ,
      "Back pain not dominating" = "ja"
    )
  )%>%
  mutate(
    bfbe0 = fct_recode(
      bfbe0 ,
      "Sitting" = "siddende" ,
      "Sitting and walking" = "siddende_gående" ,
      "Light physical work" = "let_fysisk_arbejde" ,
      "Heavy physical work" = "tungt_fysisk_arbejde"
    )
  ) %>%
  mutate(
    budd0 = fct_recode(
      budd0 ,
      "No vocational training" = "Ingen_erhvervsudd" ,
      "Vocational training" = "Faglig_udd" ,
      "Short Higher Education" = "Kort_videregående_udd" ,
      "Middle Higher Education" = "Mellemlang_videregående_udd" ,
      "Long Higher Education" = "Lang_videregående_udd"
    )
  ) %>% 
  mutate(
    barb0 = fct_collapse(
      barb0 ,
      "Self-employed" = "selvstændig" ,
      "Full-time" = "fuldtidsansat" ,
      "Part-time" = "deltidsansat" ,
      "Student" = "studerende" , 
      "Unemployed" = "arbejdsledig" , 
      "Early retirement/retired" = "efterløn_pension" , 
      "Other" = c("hjemmegående" , "medarbejdende")
    )
  ) %>%
  mutate(tlep0 = fct_recode(
    tlep0 ,
    "0" = "0",
    "1-3" = "1-3",
    "More than  3" = "flere_end_3"
  ))

factors <- c("bsex0" , 
             "vasb0_dic" , 
             "dlva0" , 
             "tlep0" , 
             "bfbe0" , 
             "tlda0" ,
             "nootherdisease" , 
             "bryg0" , 
             "pain_dis" , 
             "domin_bp")
cont <- c("vasl0" , 
          "age" , 
          "bhoej0" , 
          "bmi" , 
          "htil0" , 
          "rmprop")

temp_df<-cbind(lbp.data_filtered[,cont],na.tree.replace(lbp.data_filtered[,factors]))
#Imputation of NA values for continuous variable 
temp_df[,cont] <- sapply(cont,function(k){temp_df[,k] <- 
  ifelse(is.na(temp_df[,k]) , 
         mean(temp_df[,k],na.rm = TRUE) , 
         temp_df[,k])})
risk_factors <- model.matrix(formula(paste("~", paste(c(factors, cont), collapse = "+"))), data = temp_df)

rm(temp_df,cont,factors)


N <- 425
set.seed(1)

train_idx <- sample(1:nrow(lbp.data_filtered), N, replace = FALSE)
risk_factors<-risk_factors[,-1]
PreProcess<-preProcess(risk_factors[train_idx,],c("center"))
risk_factors[train_idx,] <-predict(PreProcess,risk_factors[train_idx,])
risk_factors[-train_idx,] <-predict(PreProcess,risk_factors[-train_idx,])
lbp.duration <- lbp.data_filtered[,paste0("lbpdays",1:52)]
lbp.intensity <- lbp.data_filtered[,paste0("nrs",1:52)]
lbp.limit <- lbp.data_filtered[,paste0("actlimi",1:52)]
T <- 52
Y_dur <- lbp.duration 
Y_pain <- lbp.intensity 
Y_limit <- lbp.limit 
Y_dur[is.na(Y_dur)] <-8
Y_pain[is.na(Y_pain)] <-11
Y_limit[is.na(Y_limit)] <-8
Y_dur<-matrix(as.integer(unlist(Y_dur)),ncol=T,byrow=FALSE)
Y_pain<-matrix(as.integer(unlist(Y_pain)),ncol=T,byrow=FALSE)
Y_limit<-matrix(as.integer(unlist(Y_limit)),ncol=T,byrow=FALSE)

rm(lbp.intensity,lbp.limit,lbp.duration)
save(list = ls(), file = "../Preprocessed.MHMM.Data.R")
rm(list=ls())





