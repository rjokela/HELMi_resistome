
# Housekeeping #####
library(tidyverse)
library(vegan)
library(skimr)
library(coda.base)

setwd("/path/to/folder/")

# Source code modified from https://github.com/kdyson/R_Scripts/blob/master/AICc_PERMANOVA.R
source("../folder2/r-codes/AIC_PERMANOVA.R")

meta <- read_csv("metadata/metadata.csv")
colnames(meta)
colnames(meta) <- gsub(" ", "_", colnames(meta))
colnames(meta) <- gsub("-", "", colnames(meta))

arg <- read_csv("arg_files/normalised_simple_bla.csv")
arg <- arg %>%
  column_to_rownames("simplified bla")
targ <- as.data.frame(t(arg))
rtarg <- targ/rowSums(targ)

skim(meta[meta$Sample_type == "B4",])


# B4 ####
meta_b4 <- meta_df[meta_df$Sample_type == "B4",]
meta_b4 <- meta_b4 %>%
  select_if(~ sum(!is.na(.)) > 0.9 * 475)
meta_b4 <- na.omit(meta_b4)

summarytools::dfSummary(meta_b4)

targ_b4 <- targ[grepl("B4", rownames(targ)),]
targ_b4 <- targ_b4[,colSums(targ_b4) > 0]

targ_b4 <- targ_b4[rownames(targ_b4) %in% meta_b4$Sample_ID,]

identical(rownames(targ_b4), meta_df$Sample_ID) #?
all(rownames(targ_b4) == meta_b4$Sample_ID) # OK

rtarg_b4 <- rtarg[grepl("B4", rownames(rtarg)),]
rtarg_b4 <- rtarg_b4[,colSums(rtarg_b4) > 0]

dist_b4 <- rtarg_b4 %>%
  vegdist(method = "bray")

ado_b4 <- adonis2(dist_b4 ~ MGI_Plate + inf_DeliveryMode + IAP_substance + 
                    Maternal_ab + Paternal_ab + birth_season + 
                    m_PreviousDeliveries + inf_AgeHome + 
                    inf_Gestational_term + inf_HasHospitalization_past3m + 
                    Formula_1stmonth + Formula + ExclusiveBF + 
                    genotype +
                    m_Education + 
                    inf_cummulative_NbHospi +
                    inf_HospitalLocal + env_FurOrFeathers + 
                    m_AgeDelivery + Penicillins + 
                    inf_SkinContactAge + m_BMI + 
                    ATB_prev6m_corr + familly_TimeSpendNonParents3m + 
                    SUM + familly_NbCaregivers3m_Cat + 
                    breastfeeding + 
                    inf_TravelAbroad3m + inf_Sex + 
                    chao1 + Cluster,
                  meta_b4,
                  permutations = 99) 

AICc.PERMANOVA2(ado_b4, meta_b4)$small
# Small, use AICc
curr_AIC <- AICc.PERMANOVA2(ado_b4, meta_b4)$AICc


mod_vars <- c("MGI_Plate", "inf_DeliveryMode", "IAP_substance",
              "Maternal_ab", "Paternal_ab", "birth_season",
              "m_PreviousDeliveries", "inf_AgeHome",
              "inf_Gestational_term",
              "inf_HasHospitalization_past3m",
              "Formula_1stmonth", "Formula", "ExclusiveBF",
              "genotype",
              "m_Education", "inf_cummulative_NbHospi",
              "inf_HospitalLocal", "env_FurOrFeathers",
              "m_AgeDelivery", "Penicillins",
              "inf_SkinContactAge", "m_BMI",
              "ATB_prev6m_corr", "familly_TimeSpendNonParents3m",
              "SUM", "familly_NbCaregivers3m_Cat",
              "breastfeeding", "inf_TravelAbroad3m", 
              "inf_Sex",
              "chao1", "Cluster") 
curr_model <- ado_b4


model_list <- list()
AIC_list <- list()

### model selection#####
# DIY drop 1
# loop without the technical variable
for (i in 2:length(mod_vars)) {
  curr_vars <- paste(mod_vars[-i], collapse = "+")
  tmp_model <- adonis2(as.formula(paste("dist_b4", "~", curr_vars)),
                       meta_b4, 
                       permutations = 99
  )
  tmp_AIC <- AICc.PERMANOVA2(tmp_model, meta_b4)$AICc
  
  # save models to list
  model_list[[i]] <- tmp_model
  AIC_list[[i]] <- tmp_AIC
  new_AIC <- Reduce(min, AIC_list[2:length(AIC_list)])
  # if new model has AIC of 2 of more smaller, update model and variable list (and use new AIC in next round)
}


AIC_list
indx <- which(AIC_list[2:length(AIC_list)] == min(unlist(AIC_list[2:length(AIC_list)]))) 
indx

model_list[2:length(model_list)][[indx]]

# base the selection on lowest AIC 

if (curr_AIC - new_AIC >= 2) {
  curr_model <- model_list[2:length(model_list)][[indx]]
  curr_AIC <- new_AIC
}

# next round
mod_vars2 <- c("MGI_Plate", "inf_DeliveryMode", #"IAP_substance",
               "Maternal_ab", "Paternal_ab", "birth_season",
               "m_PreviousDeliveries", "inf_AgeHome",
               "inf_Gestational_term",
               "inf_HasHospitalization_past3m",
               "Formula_1stmonth", "Formula", "ExclusiveBF",
               "genotype",
               "m_Education", "inf_cummulative_NbHospi",
               "inf_HospitalLocal", "env_FurOrFeathers",
               "m_AgeDelivery", "Penicillins",
               "inf_SkinContactAge", "m_BMI",
               "ATB_prev6m_corr", "familly_TimeSpendNonParents3m",
               "SUM", "familly_NbCaregivers3m_Cat",
               "breastfeeding", "inf_TravelAbroad3m", 
               "inf_Sex",
               "chao1", "Cluster") 

model_list2 <- list()
AIC_list2 <- list()

### model selection#####
# DIY drop 1
# loop without the technical variable
for (i in 2:length(mod_vars2)) {
  curr_vars <- paste(mod_vars2[-i], collapse = "+")
  tmp_model <- adonis2(as.formula(paste("dist_b4", "~", curr_vars)),
                       meta_b4, 
                       permutations = 99
  )
  tmp_AIC <- AICc.PERMANOVA2(tmp_model, meta_b4)$AICc
  
  # save models to list
  model_list2[[i]] <- tmp_model
  AIC_list2[[i]] <- tmp_AIC
  new_AIC <- Reduce(min, AIC_list2[2:length(AIC_list2)])
  # if new model has AIC of 2 of more smaller, update model and variable list (and use new AIC in next round)
}


AIC_list2
indx <- which(AIC_list2[2:length(AIC_list2)] == min(unlist(AIC_list2[2:length(AIC_list2)]))) 
indx

model_list2[2:length(model_list2)][[indx]]

# clear lowest


if (curr_AIC - new_AIC >= 2) {
  curr_model <- model_list2[2:length(model_list2)][[indx]]
  curr_AIC <- new_AIC
}

# next round
mod_vars3 <- c("MGI_Plate", "inf_DeliveryMode", #"IAP_substance",
               "Maternal_ab", "Paternal_ab", #"birth_season",
               "m_PreviousDeliveries", "inf_AgeHome",
               "inf_Gestational_term",
               "inf_HasHospitalization_past3m",
               "Formula_1stmonth", "Formula", "ExclusiveBF",
               "genotype",
               "m_Education", "inf_cummulative_NbHospi",
               "inf_HospitalLocal", "env_FurOrFeathers",
               "m_AgeDelivery", "Penicillins",
               "inf_SkinContactAge", "m_BMI",
               "ATB_prev6m_corr", "familly_TimeSpendNonParents3m",
               "SUM", "familly_NbCaregivers3m_Cat",
               "breastfeeding", "inf_TravelAbroad3m", 
               "inf_Sex",
               "chao1", "Cluster") 

model_list3 <- list()
AIC_list3 <- list()

### model selection#####
# DIY drop 1
# loop without the technical variable
for (i in 2:length(mod_vars3)) {
  curr_vars <- paste(mod_vars3[-i], collapse = "+")
  tmp_model <- adonis2(as.formula(paste("dist_b4", "~", curr_vars)),
                       meta_b4, 
                       permutations = 99
  )
  tmp_AIC <- AICc.PERMANOVA2(tmp_model, meta_b4)$AICc
  
  # save models to list
  model_list3[[i]] <- tmp_model
  AIC_list3[[i]] <- tmp_AIC
  new_AIC <- Reduce(min, AIC_list3[2:length(AIC_list3)])
  # if new model has AIC of 2 of more smaller, update model and variable list (and use new AIC in next round)
}


AIC_list3
indx <- which(AIC_list3[2:length(AIC_list3)] == min(unlist(AIC_list3[2:length(AIC_list3)]))) 
indx

model_list3[2:length(model_list3)][[indx]]

# 

if (curr_AIC - new_AIC >= 2) {
  curr_model <- model_list3[2:length(model_list3)][[indx]]
  curr_AIC <- new_AIC
}

# next round
mod_vars4 <- c("MGI_Plate", "inf_DeliveryMode", #"IAP_substance",
               "Maternal_ab", #"Paternal_ab",
               # "birth_season",
               "m_PreviousDeliveries", "inf_AgeHome",
               "inf_Gestational_term",
               "inf_HasHospitalization_past3m",
               "Formula_1stmonth", "Formula", "ExclusiveBF",
               "genotype",
               "m_Education", "inf_cummulative_NbHospi",
               "inf_HospitalLocal", "env_FurOrFeathers",
               "m_AgeDelivery", "Penicillins",
               "inf_SkinContactAge", "m_BMI",
               "ATB_prev6m_corr", "familly_TimeSpendNonParents3m",
               "SUM", "familly_NbCaregivers3m_Cat",
               "breastfeeding", "inf_TravelAbroad3m", 
               "inf_Sex",
               "chao1", "Cluster") 

model_list4 <- list()
AIC_list4 <- list()

### model selection#####
# DIY drop 1
# loop without the technical variable
for (i in 2:length(mod_vars4)) {
  curr_vars <- paste(mod_vars4[-i], collapse = "+")
  tmp_model <- adonis2(as.formula(paste("dist_b4", "~", curr_vars)),
                       meta_b4, 
                       permutations = 99
  )
  tmp_AIC <- AICc.PERMANOVA2(tmp_model, meta_b4)$AICc
  
  # save models to list
  model_list4[[i]] <- tmp_model
  AIC_list4[[i]] <- tmp_AIC
  new_AIC <- Reduce(min, AIC_list4[2:length(AIC_list4)])
  # if new model has AIC of 2 of more smaller, update model and variable list (and use new AIC in next round)
}


AIC_list4
indx <- which(AIC_list4[2:length(AIC_list4)] == min(unlist(AIC_list4[2:length(AIC_list4)]))) 
indx

model_list4[2:length(model_list4)][[indx]]

# selection

if (curr_AIC - new_AIC >= 2) {
  curr_model <- model_list4[2:length(model_list4)][[indx]]
  curr_AIC <- new_AIC
}


# next round
mod_vars5 <- c("MGI_Plate", "inf_DeliveryMode", #"IAP_substance",
               "Maternal_ab", #"Paternal_ab",# "birth_season",
               "m_PreviousDeliveries", "inf_AgeHome",
               # "inf_Gestational_term",
               "inf_HasHospitalization_past3m",
               "Formula_1stmonth", "Formula", "ExclusiveBF",
               "genotype",
               "m_Education", "inf_cummulative_NbHospi",
               "inf_HospitalLocal", "env_FurOrFeathers",
               "m_AgeDelivery", "Penicillins",
               "inf_SkinContactAge", "m_BMI",
               "ATB_prev6m_corr", "familly_TimeSpendNonParents3m",
               "SUM", "familly_NbCaregivers3m_Cat",
               "breastfeeding", "inf_TravelAbroad3m", 
               "inf_Sex",
               "chao1", "Cluster") 

model_list5 <- list()
AIC_list5 <- list()

### model selection#####
# DIY drop 1
# loop without the technical variable
for (i in 2:length(mod_vars5)) {
  curr_vars <- paste(mod_vars5[-i], collapse = "+")
  tmp_model <- adonis2(as.formula(paste("dist_b4", "~", curr_vars)),
                       meta_b4, 
                       permutations = 99
  )
  tmp_AIC <- AICc.PERMANOVA2(tmp_model, meta_b4)$AICc
  
  # save models to list
  model_list5[[i]] <- tmp_model
  AIC_list5[[i]] <- tmp_AIC
  new_AIC <- Reduce(min, AIC_list5[2:length(AIC_list5)])
  # if new model has AIC of 2 of more smaller, update model and variable list (and use new AIC in next round)
}


AIC_list5
indx <- which(AIC_list5[2:length(AIC_list5)] == min(unlist(AIC_list5[2:length(AIC_list5)]))) 
indx

model_list5[2:length(model_list5)][[indx]]

# select


if (curr_AIC - new_AIC >= 2) {
  curr_model <- model_list5[2:length(model_list5)][[indx]]
  curr_AIC <- new_AIC
}


# next round
mod_vars6 <- c("MGI_Plate", "inf_DeliveryMode", #"IAP_substance",
               "Maternal_ab", #"Paternal_ab",# "birth_season",
               "m_PreviousDeliveries", "inf_AgeHome",
               # "inf_Gestational_term",
               "inf_HasHospitalization_past3m",
               "Formula_1stmonth", "Formula", "ExclusiveBF",
               # "genotype",
               "m_Education",
               "inf_cummulative_NbHospi",
               "inf_HospitalLocal", "env_FurOrFeathers",
               "m_AgeDelivery", "Penicillins",
               "inf_SkinContactAge", "m_BMI",
               "ATB_prev6m_corr", "familly_TimeSpendNonParents3m",
               "SUM", "familly_NbCaregivers3m_Cat",
               "breastfeeding", "inf_TravelAbroad3m", 
               "inf_Sex",
               "chao1", "Cluster") 


model_list6 <- list()
AIC_list6 <- list()

### model selection#####
# DIY drop 1
# loop without the technical variable
for (i in 2:length(mod_vars6)) {
  curr_vars <- paste(mod_vars6[-i], collapse = "+")
  tmp_model <- adonis2(as.formula(paste("dist_b4", "~", curr_vars)),
                       meta_b4, 
                       permutations = 99
  )
  tmp_AIC <- AICc.PERMANOVA2(tmp_model, meta_b4)$AICc
  
  # save models to list
  model_list6[[i]] <- tmp_model
  AIC_list6[[i]] <- tmp_AIC
  new_AIC <- Reduce(min, AIC_list6[2:length(AIC_list6)])
  # if new model has AIC of 2 of more smaller, update model and variable list (and use new AIC in next round)
}


AIC_list6
indx <- which(AIC_list6[2:length(AIC_list6)] == min(unlist(AIC_list6[2:length(AIC_list6)]))) 
indx

model_list6[2:length(model_list6)][[indx]]

# select the one selected by drop1

if (curr_AIC - new_AIC >= 2) {
  curr_model <- model_list6[2:length(model_list6)][[indx]]
  curr_AIC <- new_AIC
}

# next round
mod_vars7 <- c("MGI_Plate", "inf_DeliveryMode", #"IAP_substance",
               "Maternal_ab", #"Paternal_ab",# "birth_season",
               "m_PreviousDeliveries", "inf_AgeHome",
               # "inf_Gestational_term",
               "inf_HasHospitalization_past3m",
               "Formula_1stmonth", "Formula", "ExclusiveBF",
               # "genotype",
               # "m_Education",
               "inf_cummulative_NbHospi",
               "inf_HospitalLocal", "env_FurOrFeathers",
               "m_AgeDelivery", "Penicillins",
               "inf_SkinContactAge", "m_BMI",
               "ATB_prev6m_corr", "familly_TimeSpendNonParents3m",
               "SUM", "familly_NbCaregivers3m_Cat",
               "breastfeeding", "inf_TravelAbroad3m", 
               "inf_Sex",
               "chao1", "Cluster") 

model_list7 <- list()
AIC_list7 <- list()

### model selection#####
# DIY drop 1
# loop without the technical variable
for (i in 2:length(mod_vars7)) {
  curr_vars <- paste(mod_vars7[-i], collapse = "+")
  tmp_model <- adonis2(as.formula(paste("dist_b4", "~", curr_vars)),
                       meta_b4, 
                       permutations = 99
  )
  tmp_AIC <- AICc.PERMANOVA2(tmp_model, meta_b4)$AICc
  
  # save models to list
  model_list7[[i]] <- tmp_model
  AIC_list7[[i]] <- tmp_AIC
  new_AIC <- Reduce(min, AIC_list7[2:length(AIC_list7)])
  # if new model has AIC of 2 of more smaller, update model and variable list (and use new AIC in next round)
}


AIC_list7
indx <- which(AIC_list7[2:length(AIC_list7)] == min(unlist(AIC_list7[2:length(AIC_list7)]))) 
indx

model_list7[2:length(model_list7)][[indx]]

# 

if (curr_AIC - new_AIC >= 2) {
  curr_model <- model_list7[2:length(model_list7)][[indx]]
  curr_AIC <- new_AIC
}

# next round
mod_vars8 <- c("MGI_Plate", "inf_DeliveryMode", #"IAP_substance",
               # "Maternal_ab", #"Paternal_ab",# "birth_season",
               "m_PreviousDeliveries", "inf_AgeHome",
               # "inf_Gestational_term",
               "inf_HasHospitalization_past3m",
               "Formula_1stmonth", "Formula", "ExclusiveBF",
               # "genotype",
               # "m_Education",
               "inf_cummulative_NbHospi",
               "inf_HospitalLocal", "env_FurOrFeathers",
               "m_AgeDelivery", "Penicillins",
               "inf_SkinContactAge", "m_BMI",
               "ATB_prev6m_corr", "familly_TimeSpendNonParents3m",
               "SUM", "familly_NbCaregivers3m_Cat",
               "breastfeeding", "inf_TravelAbroad3m", 
               "inf_Sex",
               "chao1", "Cluster") 

model_list8 <- list()
AIC_list8 <- list()

### model selection#####
# DIY drop 1
# loop without the technical variable
for (i in 2:length(mod_vars8)) {
  curr_vars <- paste(mod_vars8[-i], collapse = "+")
  tmp_model <- adonis2(as.formula(paste("dist_b4", "~", curr_vars)),
                       meta_b4, 
                       permutations = 99
  )
  tmp_AIC <- AICc.PERMANOVA2(tmp_model, meta_b4)$AICc
  
  # save models to list
  model_list8[[i]] <- tmp_model
  AIC_list8[[i]] <- tmp_AIC
  new_AIC <- Reduce(min, AIC_list8[2:length(AIC_list8)])
  # if new model has AIC of 2 of more smaller, update model and variable list (and use new AIC in next round)
}


AIC_list8
indx <- which(AIC_list8[2:length(AIC_list8)] == min(unlist(AIC_list8[2:length(AIC_list8)]))) 
indx

model_list8[2:length(model_list8)][[indx]]

# Select, be more strict with AIC 1.5

if (curr_AIC - new_AIC >= 1.5) {
  curr_model <- model_list8[2:length(model_list8)][[indx]]
  curr_AIC <- new_AIC
}

# next round
mod_vars9 <- c("MGI_Plate", "inf_DeliveryMode", #"IAP_substance",
               # "Maternal_ab", #"Paternal_ab",# "birth_season",
               "m_PreviousDeliveries", "inf_AgeHome",
               # "inf_Gestational_term",
               "inf_HasHospitalization_past3m",
               "Formula_1stmonth", "Formula", "ExclusiveBF",
               # "genotype",
               # "m_Education",
               "inf_cummulative_NbHospi",
               "inf_HospitalLocal", "env_FurOrFeathers",
               "m_AgeDelivery", "Penicillins",
               "inf_SkinContactAge", "m_BMI",
               "ATB_prev6m_corr", "familly_TimeSpendNonParents3m",
               "SUM", "familly_NbCaregivers3m_Cat",
               "breastfeeding", #"inf_TravelAbroad3m", 
               "inf_Sex",
               "chao1", "Cluster")

model_list9 <- list()
AIC_list9 <- list()

### model selection#####
# DIY drop 1
# loop without the technical variable
for (i in 2:length(mod_vars9)) {
  curr_vars <- paste(mod_vars9[-i], collapse = "+")
  tmp_model <- adonis2(as.formula(paste("dist_b4", "~", curr_vars)),
                       meta_b4, 
                       permutations = 99
  )
  tmp_AIC <- AICc.PERMANOVA2(tmp_model, meta_b4)$AICc
  
  # save models to list
  model_list9[[i]] <- tmp_model
  AIC_list9[[i]] <- tmp_AIC
  new_AIC <- Reduce(min, AIC_list9[2:length(AIC_list9)])
  # if new model has AIC of 2 of more smaller, update model and variable list (and use new AIC in next round)
}


AIC_list9
indx <- which(AIC_list9[2:length(AIC_list9)] == min(unlist(AIC_list9[2:length(AIC_list9)]))) 
indx

model_list9[2:length(model_list9)][[indx]]

# Select, be more strict with AIC

if (curr_AIC - new_AIC >= 1.5) {
  curr_model <- model_list9[2:length(model_list9)][[indx]]
  curr_AIC <- new_AIC
}

# next round
mod_vars10 <- c("MGI_Plate", "inf_DeliveryMode", #"IAP_substance",
                # "Maternal_ab", #"Paternal_ab",# "birth_season",
                "m_PreviousDeliveries", "inf_AgeHome",
                # "inf_Gestational_term",
                "inf_HasHospitalization_past3m",
                "Formula_1stmonth", "Formula", "ExclusiveBF",
                # "genotype",
                # "m_Education",
                # "inf_cummulative_NbHospi",
                "inf_HospitalLocal", "env_FurOrFeathers",
                "m_AgeDelivery", "Penicillins",
                "inf_SkinContactAge", "m_BMI",
                "ATB_prev6m_corr", "familly_TimeSpendNonParents3m",
                "SUM", "familly_NbCaregivers3m_Cat",
                "breastfeeding", #"inf_TravelAbroad3m", 
                "inf_Sex",
                "chao1", "Cluster")

model_list10 <- list()
AIC_list10 <- list()

### model selection#####
# DIY drop 1
# loop without the technical variable
for (i in 2:length(mod_vars10)) {
  curr_vars <- paste(mod_vars10[-i], collapse = "+")
  tmp_model <- adonis2(as.formula(paste("dist_b4", "~", curr_vars)),
                       meta_b4, 
                       permutations = 99
  )
  tmp_AIC <- AICc.PERMANOVA2(tmp_model, meta_b4)$AICc
  
  # save models to list
  model_list10[[i]] <- tmp_model
  AIC_list10[[i]] <- tmp_AIC
  new_AIC <- Reduce(min, AIC_list10[2:length(AIC_list10)])
  # if new model has AIC of 2 of more smaller, update model and variable list (and use new AIC in next round)
}


AIC_list10
indx <- which(AIC_list10[2:length(AIC_list10)] == min(unlist(AIC_list10[2:length(AIC_list10)]))) 
indx

model_list10[2:length(model_list10)][[indx]]

# Select, be more strict with AIC


if (curr_AIC - new_AIC >= 1.5) {
  curr_model <- model_list10[2:length(model_list10)][[indx]]
  curr_AIC <- new_AIC
}

# next round
mod_vars11 <- c("MGI_Plate", "inf_DeliveryMode", #"IAP_substance",
                # "Maternal_ab", #"Paternal_ab",# "birth_season",
                "m_PreviousDeliveries", "inf_AgeHome",
                # "inf_Gestational_term",
                "inf_HasHospitalization_past3m",
                "Formula_1stmonth", "Formula", "ExclusiveBF",
                # "genotype",
                # "m_Education",
                # "inf_cummulative_NbHospi",
                "inf_HospitalLocal", "env_FurOrFeathers",
                "m_AgeDelivery", "Penicillins",
                "inf_SkinContactAge", "m_BMI",
                "ATB_prev6m_corr", "familly_TimeSpendNonParents3m",
                "SUM", "familly_NbCaregivers3m_Cat",
                "breastfeeding", #"inf_TravelAbroad3m", 
                # "inf_Sex",
                "chao1", "Cluster")

model_list11 <- list()
AIC_list11 <- list()

### model selection#####
# DIY drop 1
# loop without the technical variable
for (i in 2:length(mod_vars11)) {
  curr_vars <- paste(mod_vars11[-i], collapse = "+")
  tmp_model <- adonis2(as.formula(paste("dist_b4", "~", curr_vars)),
                       meta_b4, 
                       permutations = 99
  )
  tmp_AIC <- AICc.PERMANOVA2(tmp_model, meta_b4)$AICc
  
  # save models to list
  model_list11[[i]] <- tmp_model
  AIC_list11[[i]] <- tmp_AIC
  new_AIC <- Reduce(min, AIC_list11[2:length(AIC_list11)])
  # if new model has AIC of 2 of more smaller, update model and variable list (and use new AIC in next round)
}


AIC_list11
indx <- which(AIC_list11[2:length(AIC_list11)] == min(unlist(AIC_list11[2:length(AIC_list11)]))) 
indx

model_list11[2:length(model_list11)][[indx]]

# Done


if (curr_AIC - new_AIC >= 1.5) {
  curr_model <- model_list11[2:length(model_list11)][[indx]]
  curr_AIC <- new_AIC
}

final_ado_b4 <- adonis2(dist_b4 ~ MGI_Plate + inf_DeliveryMode + #IAP_substance + 
                          # Maternal_ab + Paternal_ab + birth_season + 
                          m_PreviousDeliveries + inf_AgeHome + 
                          # inf_Gestational_term + 
                          inf_HasHospitalization_past3m + 
                          Formula_1stmonth +
                          Formula + ExclusiveBF + 
                          # m_Education + 
                          # inf_cummulative_NbHospi +
                          inf_HospitalLocal + env_FurOrFeathers + 
                          m_AgeDelivery + Penicillins + 
                          inf_SkinContactAge + m_BMI + 
                          ATB_prev6m_corr +
                          familly_TimeSpendNonParents3m + 
                          SUM + familly_NbCaregivers3m_Cat +
                          breastfeeding +
                          # inf_TravelAbroad3m + 
                          inf_Sex + 
                          chao1 + Cluster,
                        meta_b4,
                        permutations = 999)



# B5 ######
skim(meta[meta$Sample_type == "B5",])

# B5 ####
meta_b5 <- meta_df[meta_df$Sample_type == "B5",]
meta_b5 <- meta_b5 %>%
  select_if(~ sum(!is.na(.)) > 0.9 * 475)
meta_b5 <- na.omit(meta_b5)

targ_b5 <- targ[grepl("B5", rownames(targ)),]
targ_b5 <- targ_b5[,colSums(targ_b5) > 0]

targ_b5 <- targ_b5[rownames(targ_b5) %in% meta_b5$Sample_ID,]

identical(rownames(targ_b5), meta_df$Sample_ID) #?
all(rownames(targ_b5) == meta_b5$Sample_ID) # OK

dist_b5 <- rtarg_b5 %>%
  vegdist(method = "bray")

ado_b5 <- adonis2(dist_b5 ~ MGI_Plate + inf_DeliveryMode + Solids_ever +
                    Maternal_ab + Paternal_ab + IAP_substance + 
                    birth_season + Solids_quantity +
                    inf_Gestational_term + env_FurOrFeathers + 
                    inf_AgeHome + m_PreviousDeliveries + 
                    inf_TravelAbroad3m + breastfeeding + 
                    ExclusiveBF + genotype +
                    m_Education + familly_NbCaregivers3m_Cat + 
                    Solids + Formula_1stmonth + Formula + 
                    ATB_prev6m_corr + inf_SkinContactAge + 
                    Penicillins + m_AgeDelivery + m_BMI + 
                    inf_IsAway3m + inf_HospitalLocal + 
                    familly_TimeSpendNonParents3m + 
                    inf_HasHospitalization_past3m + 
                    inf_cummulative_NbHospi +
                    inf_Sex + SUM + 
                    chao1 + Cluster,
                  meta_b5,
                  permutations = 99) 

AICc.PERMANOVA2(ado_b5, meta_b5)$small
# Small, use AICc
curr_AIC <- AICc.PERMANOVA2(ado_b5, meta_b5)$AICc


mod_vars <- c("MGI_Plate", "inf_DeliveryMode", "Solids_ever",
              "Maternal_ab", "Paternal_ab", "IAP_substance",
              "birth_season", "Solids_quantity",
              "inf_Gestational_term", "env_FurOrFeathers",
              "inf_AgeHome", "m_PreviousDeliveries",
              "inf_TravelAbroad3m", "breastfeeding",
              "ExclusiveBF", "genotype",
              "m_Education", "familly_NbCaregivers3m_Cat",
              "Solids", "Formula_1stmonth", "Formula",
              "ATB_prev6m_corr", "inf_SkinContactAge", 
              "Penicillins", "m_AgeDelivery", "m_BMI", 
              "inf_IsAway3m", "inf_HospitalLocal", 
              "familly_TimeSpendNonParents3m", 
              "inf_HasHospitalization_past3m", 
              "inf_cummulative_NbHospi", 
              "inf_Sex", "SUM", 
              "chao1", "Cluster") 

curr_model <- ado_b5


model_list <- list()
AIC_list <- list()

### model selection#####
# DIY drop 1
# loop without the technical variable
for (i in 2:length(mod_vars)) {
  curr_vars <- paste(mod_vars[-i], collapse = "+")
  tmp_model <- adonis2(as.formula(paste("dist_b5", "~", curr_vars)),
                       meta_b5, 
                       permutations = 99
  )
  tmp_AIC <- AICc.PERMANOVA2(tmp_model, meta_b5)$AICc
  
  # save models to list
  model_list[[i]] <- tmp_model
  AIC_list[[i]] <- tmp_AIC
  new_AIC <- Reduce(min, AIC_list[2:length(AIC_list)])
  # if new model has AIC of 2 of more smaller, update model and variable list (and use new AIC in next round)
}


AIC_list
indx <- which(AIC_list[2:length(AIC_list)] == min(unlist(AIC_list[2:length(AIC_list)]))) 
indx

model_list[2:length(model_list)][[indx]]

# selection


if (curr_AIC - new_AIC >= 2) {
  curr_model <- model_list[2:length(model_list)][[indx]]
  curr_AIC <- new_AIC
}

# next round
mod_vars2 <- c("MGI_Plate", "inf_DeliveryMode", "Solids_ever",
               "Maternal_ab", "Paternal_ab", "IAP_substance",
               # "birth_season", 
               "Solids_quantity",
               "inf_Gestational_term",
               "env_FurOrFeathers",
               "inf_AgeHome", "m_PreviousDeliveries",
               "inf_TravelAbroad3m", "breastfeeding",
               "ExclusiveBF", "genotype",
               "m_Education", "familly_NbCaregivers3m_Cat",
               "Solids", "Formula_1stmonth", "Formula",
               "ATB_prev6m_corr", "inf_SkinContactAge", 
               "Penicillins", "m_AgeDelivery", "m_BMI", 
               "inf_IsAway3m", "inf_HospitalLocal", 
               "familly_TimeSpendNonParents3m", 
               "inf_HasHospitalization_past3m", 
               "inf_cummulative_NbHospi", 
               "inf_Sex", "SUM", 
               "chao1", "Cluster") 

model_list2 <- list()
AIC_list2 <- list()

### model selection#####
# DIY drop 1
# loop without the technical variable
for (i in 2:length(mod_vars2)) {
  curr_vars <- paste(mod_vars2[-i], collapse = "+")
  tmp_model <- adonis2(as.formula(paste("dist_b5", "~", curr_vars)),
                       meta_b5, 
                       permutations = 99
  )
  tmp_AIC <- AICc.PERMANOVA2(tmp_model, meta_b5)$AICc
  
  # save models to list
  model_list2[[i]] <- tmp_model
  AIC_list2[[i]] <- tmp_AIC
  new_AIC <- Reduce(min, AIC_list2[2:length(AIC_list2)])
  # if new model has AIC of 2 of more smaller, update model and variable list (and use new AIC in next round)
}


AIC_list2
indx <- which(AIC_list2[2:length(AIC_list2)] == min(unlist(AIC_list2[2:length(AIC_list2)]))) 
indx

model_list2[2:length(model_list2)][[indx]]

# selection

if (curr_AIC - new_AIC >= 2) {
  curr_model <- model_list2[2:length(model_list2)][[indx]]
  curr_AIC <- new_AIC
}


# next round
mod_vars3 <- c("MGI_Plate", "inf_DeliveryMode", "Solids_ever",
               "Maternal_ab", "Paternal_ab", #"IAP_substance",
               # "birth_season",
               "Solids_quantity",
               "inf_Gestational_term",
               "env_FurOrFeathers",
               "inf_AgeHome", "m_PreviousDeliveries",
               "inf_TravelAbroad3m", "breastfeeding",
               "ExclusiveBF", "genotype",
               "m_Education", "familly_NbCaregivers3m_Cat",
               "Solids", "Formula_1stmonth", "Formula",
               "ATB_prev6m_corr", "inf_SkinContactAge", 
               "Penicillins", "m_AgeDelivery", "m_BMI", 
               "inf_IsAway3m", "inf_HospitalLocal", 
               "familly_TimeSpendNonParents3m", 
               "inf_HasHospitalization_past3m", 
               "inf_cummulative_NbHospi", 
               "inf_Sex", "SUM", 
               "chao1", "Cluster") 

model_list3 <- list()
AIC_list3 <- list()

### model selection#####
# DIY drop 1
# loop without the technical variable
for (i in 2:length(mod_vars3)) {
  curr_vars <- paste(mod_vars3[-i], collapse = "+")
  tmp_model <- adonis2(as.formula(paste("dist_b5", "~", curr_vars)),
                       meta_b5, 
                       permutations = 99
  )
  tmp_AIC <- AICc.PERMANOVA2(tmp_model, meta_b5)$AICc
  
  # save models to list
  model_list3[[i]] <- tmp_model
  AIC_list3[[i]] <- tmp_AIC
  new_AIC <- Reduce(min, AIC_list3[2:length(AIC_list3)])
  # if new model has AIC of 2 of more smaller, update model and variable list (and use new AIC in next round)
}


AIC_list3
indx <- which(AIC_list3[2:length(AIC_list3)] == min(unlist(AIC_list3[2:length(AIC_list3)]))) 
indx

model_list3[2:length(model_list3)][[indx]]

# go with the selection


if (curr_AIC - new_AIC >= 2) {
  curr_model <- model_list3[2:length(model_list3)][[indx]]
  curr_AIC <- new_AIC
}

# next round
mod_vars4 <- c("MGI_Plate", "inf_DeliveryMode", "Solids_ever",
               "Maternal_ab", "Paternal_ab", #"IAP_substance",
               # "birth_season",
               "Solids_quantity",
               "inf_Gestational_term",
               "env_FurOrFeathers",
               "inf_AgeHome", "m_PreviousDeliveries",
               "inf_TravelAbroad3m", "breastfeeding",
               "ExclusiveBF", "genotype",
               "m_Education", #"familly_NbCaregivers3m_Cat",
               "Solids", "Formula_1stmonth", "Formula",
               "ATB_prev6m_corr", "inf_SkinContactAge", 
               "Penicillins", "m_AgeDelivery", "m_BMI", 
               "inf_IsAway3m", "inf_HospitalLocal", 
               "familly_TimeSpendNonParents3m", 
               "inf_HasHospitalization_past3m", 
               "inf_cummulative_NbHospi", 
               "inf_Sex", "SUM", 
               "chao1", "Cluster") 

model_list4 <- list()
AIC_list4 <- list()

### model selection#####
# DIY drop 1
# loop without the technical variable
for (i in 2:length(mod_vars4)) {
  curr_vars <- paste(mod_vars4[-i], collapse = "+")
  tmp_model <- adonis2(as.formula(paste("dist_b5", "~", curr_vars)),
                       meta_b5, 
                       permutations = 99
  )
  tmp_AIC <- AICc.PERMANOVA2(tmp_model, meta_b5)$AICc
  
  # save models to list
  model_list4[[i]] <- tmp_model
  AIC_list4[[i]] <- tmp_AIC
  new_AIC <- Reduce(min, AIC_list4[2:length(AIC_list4)])
  # if new model has AIC of 2 of more smaller, update model and variable list (and use new AIC in next round)
}


AIC_list4
indx <- which(AIC_list4[2:length(AIC_list4)] == min(unlist(AIC_list4[2:length(AIC_list4)]))) 
indx

model_list4[2:length(model_list4)][[indx]]

# go with the selection

if (curr_AIC - new_AIC >= 2) {
  curr_model <- model_list4[2:length(model_list4)][[indx]]
  curr_AIC <- new_AIC
}


# next round
mod_vars5 <- c("MGI_Plate", "inf_DeliveryMode", "Solids_ever",
               "Maternal_ab", "Paternal_ab", #"IAP_substance",
               # "birth_season",
               "Solids_quantity",
               "inf_Gestational_term",
               "env_FurOrFeathers",
               "inf_AgeHome", "m_PreviousDeliveries",
               "inf_TravelAbroad3m", "breastfeeding",
               "ExclusiveBF", "genotype",
               # "m_Education", #"familly_NbCaregivers3m_Cat",
               "Solids", "Formula_1stmonth", "Formula",
               "ATB_prev6m_corr", "inf_SkinContactAge", 
               "Penicillins", "m_AgeDelivery", "m_BMI", 
               "inf_IsAway3m", "inf_HospitalLocal", 
               "familly_TimeSpendNonParents3m", 
               "inf_HasHospitalization_past3m", 
               "inf_cummulative_NbHospi", 
               "inf_Sex", "SUM", 
               "chao1", "Cluster") 

model_list5 <- list()
AIC_list5 <- list()

### model selection#####
# DIY drop 1
# loop without the technical variable
for (i in 2:length(mod_vars5)) {
  curr_vars <- paste(mod_vars5[-i], collapse = "+")
  tmp_model <- adonis2(as.formula(paste("dist_b5", "~", curr_vars)),
                       meta_b5, 
                       permutations = 99
  )
  tmp_AIC <- AICc.PERMANOVA2(tmp_model, meta_b5)$AICc
  
  # save models to list
  model_list5[[i]] <- tmp_model
  AIC_list5[[i]] <- tmp_AIC
  new_AIC <- Reduce(min, AIC_list5[2:length(AIC_list5)])
  # if new model has AIC of 2 of more smaller, update model and variable list (and use new AIC in next round)
}


AIC_list5
indx <- which(AIC_list5[2:length(AIC_list5)] == min(unlist(AIC_list5[2:length(AIC_list5)]))) 
indx

model_list5[2:length(model_list5)][[indx]]

#  selection

if (curr_AIC - new_AIC >= 2) {
  curr_model <- model_list5[2:length(model_list5)][[indx]]
  curr_AIC <- new_AIC
}

# next round
mod_vars6 <- c("MGI_Plate", "inf_DeliveryMode", "Solids_ever",
               "Maternal_ab", "Paternal_ab", #"IAP_substance",
               # "birth_season",
               "Solids_quantity",
               # "inf_Gestational_term",
               "env_FurOrFeathers",
               "inf_AgeHome", "m_PreviousDeliveries",
               "inf_TravelAbroad3m", "breastfeeding",
               "ExclusiveBF", "genotype",
               # "m_Education", #"familly_NbCaregivers3m_Cat",
               "Solids", "Formula_1stmonth", "Formula",
               "ATB_prev6m_corr", "inf_SkinContactAge", 
               "Penicillins", "m_AgeDelivery", "m_BMI", 
               "inf_IsAway3m", "inf_HospitalLocal", 
               "familly_TimeSpendNonParents3m", 
               "inf_HasHospitalization_past3m", 
               "inf_cummulative_NbHospi", 
               "inf_Sex", "SUM", 
               "chao1", "Cluster") 


model_list6 <- list()
AIC_list6 <- list()

### model selection#####
# DIY drop 1
# loop without the technical variable
for (i in 2:length(mod_vars6)) {
  curr_vars <- paste(mod_vars6[-i], collapse = "+")
  tmp_model <- adonis2(as.formula(paste("dist_b5", "~", curr_vars)),
                       meta_b5, 
                       permutations = 99
  )
  tmp_AIC <- AICc.PERMANOVA2(tmp_model, meta_b5)$AICc
  
  # save models to list
  model_list6[[i]] <- tmp_model
  AIC_list6[[i]] <- tmp_AIC
  new_AIC <- Reduce(min, AIC_list6[2:length(AIC_list6)])
  # if new model has AIC of 2 of more smaller, update model and variable list (and use new AIC in next round)
}


AIC_list6
indx <- which(AIC_list6[2:length(AIC_list6)] == min(unlist(AIC_list6[2:length(AIC_list6)]))) 
indx

model_list6[2:length(model_list6)][[indx]]

#  selection


if (curr_AIC - new_AIC >= 2) {
  curr_model <- model_list6[2:length(model_list6)][[indx]]
  curr_AIC <- new_AIC
}


# next round
mod_vars7 <- c("MGI_Plate", "inf_DeliveryMode", "Solids_ever",
               "Maternal_ab", "Paternal_ab", #"IAP_substance",
               # "birth_season",
               "Solids_quantity",
               # "inf_Gestational_term",
               "env_FurOrFeathers",
               "inf_AgeHome", "m_PreviousDeliveries",
               "inf_TravelAbroad3m", "breastfeeding",
               "ExclusiveBF", #"genotype",
               # "m_Education", #"familly_NbCaregivers3m_Cat",
               "Solids", "Formula_1stmonth", "Formula",
               "ATB_prev6m_corr", "inf_SkinContactAge", 
               "Penicillins", "m_AgeDelivery", "m_BMI", 
               "inf_IsAway3m", "inf_HospitalLocal", 
               "familly_TimeSpendNonParents3m", 
               "inf_HasHospitalization_past3m", 
               "inf_cummulative_NbHospi", 
               "inf_Sex", "SUM", 
               "chao1", "Cluster") 

model_list7 <- list()
AIC_list7 <- list()

### model selection#####
# DIY drop 1
# loop without the technical variable
for (i in 2:length(mod_vars7)) {
  curr_vars <- paste(mod_vars7[-i], collapse = "+")
  tmp_model <- adonis2(as.formula(paste("dist_b5", "~", curr_vars)),
                       meta_b5, 
                       permutations = 99
  )
  tmp_AIC <- AICc.PERMANOVA2(tmp_model, meta_b5)$AICc
  
  # save models to list
  model_list7[[i]] <- tmp_model
  AIC_list7[[i]] <- tmp_AIC
  new_AIC <- Reduce(min, AIC_list7[2:length(AIC_list7)])
  # if new model has AIC of 2 of more smaller, update model and variable list (and use new AIC in next round)
}


AIC_list7
indx <- which(AIC_list7[2:length(AIC_list7)] == min(unlist(AIC_list7[2:length(AIC_list7)]))) 
indx

model_list7[2:length(model_list7)][[indx]]

#  selection


if (curr_AIC - new_AIC >= 2) {
  curr_model <- model_list7[2:length(model_list7)][[indx]]
  curr_AIC <- new_AIC
}

# next round
mod_vars8 <- c("MGI_Plate", "inf_DeliveryMode", "Solids_ever",
               "Maternal_ab", "Paternal_ab", #"IAP_substance",
               # "birth_season",
               # "Solids_quantity",
               # "inf_Gestational_term",
               "env_FurOrFeathers",
               "inf_AgeHome", "m_PreviousDeliveries",
               "inf_TravelAbroad3m", "breastfeeding",
               "ExclusiveBF", # "genotype",
               # "m_Education", #"familly_NbCaregivers3m_Cat",
               "Solids", "Formula_1stmonth", "Formula",
               "ATB_prev6m_corr", "inf_SkinContactAge", 
               "Penicillins", "m_AgeDelivery", "m_BMI", 
               "inf_IsAway3m", "inf_HospitalLocal", 
               "familly_TimeSpendNonParents3m", 
               "inf_HasHospitalization_past3m", 
               "inf_cummulative_NbHospi", 
               "inf_Sex", "SUM", 
               "chao1", "Cluster") 

model_list8 <- list()
AIC_list8 <- list()

### model selection#####
# DIY drop 1
# loop without the technical variable
for (i in 2:length(mod_vars8)) {
  curr_vars <- paste(mod_vars8[-i], collapse = "+")
  tmp_model <- adonis2(as.formula(paste("dist_b5", "~", curr_vars)),
                       meta_b5, 
                       permutations = 99
  )
  tmp_AIC <- AICc.PERMANOVA2(tmp_model, meta_b5)$AICc
  
  # save models to list
  model_list8[[i]] <- tmp_model
  AIC_list8[[i]] <- tmp_AIC
  new_AIC <- Reduce(min, AIC_list8[2:length(AIC_list8)])
  # if new model has AIC of 2 of more smaller, update model and variable list (and use new AIC in next round)
}


AIC_list8
indx <- which(AIC_list8[2:length(AIC_list8)] == min(unlist(AIC_list8[2:length(AIC_list8)]))) 
indx

model_list8[2:length(model_list8)][[indx]]

# selection


if (curr_AIC - new_AIC >= 2) {
  curr_model <- model_list8[2:length(model_list8)][[indx]]
  curr_AIC <- new_AIC
}

# next round
mod_vars9 <- c("MGI_Plate", "inf_DeliveryMode", "Solids_ever",
               "Maternal_ab", #"Paternal_ab", #"IAP_substance",
               # "birth_season",
               # "Solids_quantity",
               # "inf_Gestational_term",
               "env_FurOrFeathers",
               "inf_AgeHome", "m_PreviousDeliveries",
               "inf_TravelAbroad3m", "breastfeeding",
               "ExclusiveBF", #"genotype",
               # "m_Education", #"familly_NbCaregivers3m_Cat",
               "Solids", "Formula_1stmonth", "Formula",
               "ATB_prev6m_corr", "inf_SkinContactAge", 
               "Penicillins", "m_AgeDelivery", "m_BMI", 
               "inf_IsAway3m", "inf_HospitalLocal", 
               "familly_TimeSpendNonParents3m", 
               "inf_HasHospitalization_past3m", 
               "inf_cummulative_NbHospi", 
               "inf_Sex", "SUM", 
               "chao1", "Cluster") 

model_list9 <- list()
AIC_list9 <- list()

### model selection#####
# DIY drop 1
# loop without the technical variable
for (i in 2:length(mod_vars9)) {
  curr_vars <- paste(mod_vars9[-i], collapse = "+")
  tmp_model <- adonis2(as.formula(paste("dist_b5", "~", curr_vars)),
                       meta_b5, 
                       permutations = 99
  )
  tmp_AIC <- AICc.PERMANOVA2(tmp_model, meta_b5)$AICc
  
  # save models to list
  model_list9[[i]] <- tmp_model
  AIC_list9[[i]] <- tmp_AIC
  new_AIC <- Reduce(min, AIC_list9[2:length(AIC_list9)])
  # if new model has AIC of 2 of more smaller, update model and variable list (and use new AIC in next round)
}


AIC_list9
indx <- which(AIC_list9[2:length(AIC_list9)] == min(unlist(AIC_list9[2:length(AIC_list9)]))) 
indx

model_list9[2:length(model_list9)][[indx]]

# selection

if (curr_AIC - new_AIC >= 2) {
  curr_model <- model_list9[2:length(model_list9)][[indx]]
  curr_AIC <- new_AIC
}

# next round
mod_vars10 <- c("MGI_Plate", "inf_DeliveryMode", "Solids_ever",
                # "Maternal_ab", #"Paternal_ab", #"IAP_substance",
                # "birth_season",
                # "Solids_quantity",
                # "inf_Gestational_term",
                "env_FurOrFeathers",
                "inf_AgeHome", "m_PreviousDeliveries",
                "inf_TravelAbroad3m", "breastfeeding",
                "ExclusiveBF", #"genotype",
                # "m_Education", #"familly_NbCaregivers3m_Cat",
                "Solids", "Formula_1stmonth", "Formula",
                "ATB_prev6m_corr", "inf_SkinContactAge", 
                "Penicillins", "m_AgeDelivery", "m_BMI", 
                "inf_IsAway3m", "inf_HospitalLocal", 
                "familly_TimeSpendNonParents3m", 
                "inf_HasHospitalization_past3m", 
                "inf_cummulative_NbHospi", 
                "inf_Sex", "SUM", 
                "chao1", "Cluster") 

model_list10 <- list()
AIC_list10 <- list()

### model selection#####
# DIY drop 1
# loop without the technical variable
for (i in 2:length(mod_vars10)) {
  curr_vars <- paste(mod_vars10[-i], collapse = "+")
  tmp_model <- adonis2(as.formula(paste("dist_b5", "~", curr_vars)),
                       meta_b5, 
                       permutations = 99
  )
  tmp_AIC <- AICc.PERMANOVA2(tmp_model, meta_b5)$AICc
  
  # save models to list
  model_list10[[i]] <- tmp_model
  AIC_list10[[i]] <- tmp_AIC
  new_AIC <- Reduce(min, AIC_list10[2:length(AIC_list10)])
  # if new model has AIC of 2 of more smaller, update model and variable list (and use new AIC in next round)
}


AIC_list10
indx <- which(AIC_list10[2:length(AIC_list10)] == min(unlist(AIC_list10[2:length(AIC_list10)]))) 
indx

model_list10[2:length(model_list10)][[indx]]

# use more strict AIC for selection

if (curr_AIC - new_AIC >= 1.5) {
  curr_model <- model_list10[2:length(model_list10)][[indx]]
  curr_AIC <- new_AIC
}

# next round
mod_vars11 <- c("MGI_Plate", "inf_DeliveryMode", "Solids_ever",
                # "Maternal_ab", #"Paternal_ab", #"IAP_substance",
                # "birth_season",
                # "Solids_quantity",
                # "inf_Gestational_term",
                "env_FurOrFeathers",
                "inf_AgeHome", "m_PreviousDeliveries",
                "inf_TravelAbroad3m", "breastfeeding",
                "ExclusiveBF", #"genotype",
                # "m_Education", #"familly_NbCaregivers3m_Cat",
                "Solids", "Formula_1stmonth", "Formula",
                "ATB_prev6m_corr", "inf_SkinContactAge", 
                "Penicillins", "m_AgeDelivery", "m_BMI", 
                # "inf_IsAway3m", 
                "inf_HospitalLocal", 
                "familly_TimeSpendNonParents3m", 
                "inf_HasHospitalization_past3m", 
                "inf_cummulative_NbHospi", 
                "inf_Sex", "SUM", 
                "chao1", "Cluster") 

model_list11 <- list()
AIC_list11 <- list()

### model selection#####
# DIY drop 1
# loop without the technical variable
for (i in 2:length(mod_vars11)) {
  curr_vars <- paste(mod_vars11[-i], collapse = "+")
  tmp_model <- adonis2(as.formula(paste("dist_b5", "~", curr_vars)),
                       meta_b5, 
                       permutations = 99
  )
  tmp_AIC <- AICc.PERMANOVA2(tmp_model, meta_b5)$AICc
  
  # save models to list
  model_list11[[i]] <- tmp_model
  AIC_list11[[i]] <- tmp_AIC
  new_AIC <- Reduce(min, AIC_list11[2:length(AIC_list11)])
  # if new model has AIC of 2 of more smaller, update model and variable list (and use new AIC in next round)
}


AIC_list11
indx <- which(AIC_list11[2:length(AIC_list11)] == min(unlist(AIC_list11[2:length(AIC_list11)]))) 
indx

model_list11[2:length(model_list11)][[indx]]

# use more strict AIC for selection


if (curr_AIC - new_AIC >= 1.5) {
  curr_model <- model_list11[2:length(model_list11)][[indx]]
  curr_AIC <- new_AIC
}

# next round
mod_vars12 <- c("MGI_Plate", "inf_DeliveryMode", "Solids_ever",
                # "Maternal_ab", #"Paternal_ab", #"IAP_substance",
                # "birth_season",
                # "Solids_quantity",
                # "inf_Gestational_term",
                "env_FurOrFeathers",
                "inf_AgeHome", "m_PreviousDeliveries",
                "inf_TravelAbroad3m", "breastfeeding",
                "ExclusiveBF", #"genotype",
                # "m_Education", #"familly_NbCaregivers3m_Cat",
                "Solids", "Formula_1stmonth", "Formula",
                # "ATB_prev6m_corr",
                "inf_SkinContactAge", 
                "Penicillins", "m_AgeDelivery", "m_BMI", 
                # "inf_IsAway3m", 
                "inf_HospitalLocal", 
                "familly_TimeSpendNonParents3m", 
                "inf_HasHospitalization_past3m", 
                "inf_cummulative_NbHospi", 
                "inf_Sex", "SUM", 
                "chao1", "Cluster") 

model_list12 <- list()
AIC_list12 <- list()

### model selection#####
# DIY drop 1
# loop without the technical variable
for (i in 2:length(mod_vars12)) {
  curr_vars <- paste(mod_vars12[-i], collapse = "+")
  tmp_model <- adonis2(as.formula(paste("dist_b5", "~", curr_vars)),
                       meta_b5, 
                       permutations = 99
  )
  tmp_AIC <- AICc.PERMANOVA2(tmp_model, meta_b5)$AICc
  
  # save models to list
  model_list12[[i]] <- tmp_model
  AIC_list12[[i]] <- tmp_AIC
  new_AIC <- Reduce(min, AIC_list12[2:length(AIC_list12)])
  # if new model has AIC of 2 of more smaller, update model and variable list (and use new AIC in next round)
}


AIC_list12
indx <- which(AIC_list12[2:length(AIC_list12)] == min(unlist(AIC_list12[2:length(AIC_list12)]))) 
indx

model_list12[2:length(model_list12)][[indx]]

# use more strict AIC for selection


if (curr_AIC - new_AIC >= 1.5) {
  curr_model <- model_list12[2:length(model_list12)][[indx]]
  curr_AIC <- new_AIC
}

# next round
mod_vars13 <- c("MGI_Plate", "inf_DeliveryMode", "Solids_ever",
                # "Maternal_ab", #"Paternal_ab", #"IAP_substance",
                # "birth_season",
                # "Solids_quantity",
                # "inf_Gestational_term",
                "env_FurOrFeathers",
                "inf_AgeHome", "m_PreviousDeliveries",
                # "inf_TravelAbroad3m",
                "breastfeeding",
                "ExclusiveBF", #"genotype",
                # "m_Education", #"familly_NbCaregivers3m_Cat",
                "Solids", "Formula_1stmonth", "Formula",
                # "ATB_prev6m_corr",
                "inf_SkinContactAge", 
                "Penicillins", "m_AgeDelivery", "m_BMI", 
                # "inf_IsAway3m", 
                "inf_HospitalLocal", 
                "familly_TimeSpendNonParents3m", 
                "inf_HasHospitalization_past3m", 
                "inf_cummulative_NbHospi", 
                "inf_Sex", "SUM", 
                "chao1", "Cluster")

model_list13 <- list()
AIC_list13 <- list()

### model selection#####
# DIY drop 1
# loop without the technical variable
for (i in 2:length(mod_vars13)) {
  curr_vars <- paste(mod_vars13[-i], collapse = "+")
  tmp_model <- adonis2(as.formula(paste("dist_b5", "~", curr_vars)),
                       meta_b5, 
                       permutations = 99
  )
  tmp_AIC <- AICc.PERMANOVA2(tmp_model, meta_b5)$AICc
  
  # save models to list
  model_list13[[i]] <- tmp_model
  AIC_list13[[i]] <- tmp_AIC
  new_AIC <- Reduce(min, AIC_list13[2:length(AIC_list13)])
  # if new model has AIC of 2 of more smaller, update model and variable list (and use new AIC in next round)
}


AIC_list13
indx <- which(AIC_list13[2:length(AIC_list13)] == min(unlist(AIC_list13[2:length(AIC_list13)]))) 
indx

model_list13[2:length(model_list13)][[indx]]

# stricter


if (curr_AIC - new_AIC >= 1.5) {
  curr_model <- model_list13[2:length(model_list13)][[indx]]
  curr_AIC <- new_AIC
}

# next round
mod_vars14 <- c("MGI_Plate", "inf_DeliveryMode", "Solids_ever",
                # "Maternal_ab", #"Paternal_ab", #"IAP_substance",
                # "birth_season",
                # "Solids_quantity",
                # "inf_Gestational_term",
                # "env_FurOrFeathers",
                "inf_AgeHome", "m_PreviousDeliveries",
                # "inf_TravelAbroad3m",
                "breastfeeding",
                "ExclusiveBF", #"genotype",
                # "m_Education", #"familly_NbCaregivers3m_Cat",
                "Solids", "Formula_1stmonth", "Formula",
                # "ATB_prev6m_corr",
                "inf_SkinContactAge", 
                "Penicillins", "m_AgeDelivery", "m_BMI", 
                # "inf_IsAway3m", 
                "inf_HospitalLocal", 
                "familly_TimeSpendNonParents3m", 
                "inf_HasHospitalization_past3m", 
                "inf_cummulative_NbHospi", 
                "inf_Sex", "SUM", 
                "chao1", "Cluster") 

model_list14 <- list()
AIC_list14 <- list()

### model selection#####
# DIY drop 1
# loop without the technical variable
for (i in 2:length(mod_vars14)) {
  curr_vars <- paste(mod_vars14[-i], collapse = "+")
  tmp_model <- adonis2(as.formula(paste("dist_b5", "~", curr_vars)),
                       meta_b5, 
                       permutations = 99
  )
  tmp_AIC <- AICc.PERMANOVA2(tmp_model, meta_b5)$AICc
  
  # save models to list
  model_list14[[i]] <- tmp_model
  AIC_list14[[i]] <- tmp_AIC
  new_AIC <- Reduce(min, AIC_list14[2:length(AIC_list14)])
  # if new model has AIC of 2 of more smaller, update model and variable list (and use new AIC in next round)
}


AIC_list14
indx <- which(AIC_list14[2:length(AIC_list14)] == min(unlist(AIC_list14[2:length(AIC_list14)]))) 
indx

model_list14[2:length(model_list14)][[indx]]

# stricter


if (curr_AIC - new_AIC >= 1.5) {
  curr_model <- model_list14[2:length(model_list14)][[indx]]
  curr_AIC <- new_AIC
}

# next round
mod_vars15 <- c("MGI_Plate", "inf_DeliveryMode", "Solids_ever",
                # "Maternal_ab", #"Paternal_ab", #"IAP_substance",
                # "birth_season",
                # "Solids_quantity",
                # "inf_Gestational_term",
                # "env_FurOrFeathers",
                "inf_AgeHome", "m_PreviousDeliveries",
                # "inf_TravelAbroad3m",
                # "breastfeeding",
                "ExclusiveBF", #"genotype",
                # "m_Education", #"familly_NbCaregivers3m_Cat",
                "Solids", "Formula_1stmonth", "Formula",
                # "ATB_prev6m_corr",
                "inf_SkinContactAge", 
                "Penicillins", "m_AgeDelivery", "m_BMI", 
                # "inf_IsAway3m", 
                "inf_HospitalLocal", 
                "familly_TimeSpendNonParents3m", 
                "inf_HasHospitalization_past3m", 
                "inf_cummulative_NbHospi", 
                "inf_Sex", "SUM", 
                "chao1", "Cluster") 

model_list15 <- list()
AIC_list15 <- list()

### model selection#####
# DIY drop 1
# loop without the technical variable
for (i in 2:length(mod_vars15)) {
  curr_vars <- paste(mod_vars15[-i], collapse = "+")
  tmp_model <- adonis2(as.formula(paste("dist_b5", "~", curr_vars)),
                       meta_b5, 
                       permutations = 99
  )
  tmp_AIC <- AICc.PERMANOVA2(tmp_model, meta_b5)$AICc
  
  # save models to list
  model_list15[[i]] <- tmp_model
  AIC_list15[[i]] <- tmp_AIC
  new_AIC <- Reduce(min, AIC_list15[2:length(AIC_list15)])
  # if new model has AIC of 2 of more smaller, update model and variable list (and use new AIC in next round)
}


AIC_list15
indx <- which(AIC_list15[2:length(AIC_list15)] == min(unlist(AIC_list15[2:length(AIC_list15)]))) 
indx

model_list15[2:length(model_list15)][[indx]]

# stricter


if (curr_AIC - new_AIC >= 1.5) {
  curr_model <- model_list15[2:length(model_list15)][[indx]]
  curr_AIC <- new_AIC
}

# next round
mod_vars16 <-  c("MGI_Plate", "inf_DeliveryMode", "Solids_ever",
                 # "Maternal_ab", #"Paternal_ab", #"IAP_substance",
                 # "birth_season",
                 # "Solids_quantity",
                 # "inf_Gestational_term",
                 # "env_FurOrFeathers",
                 "inf_AgeHome", "m_PreviousDeliveries",
                 # "inf_TravelAbroad3m",
                 # "breastfeeding",
                 "ExclusiveBF", #"genotype",
                 # "m_Education", #"familly_NbCaregivers3m_Cat",
                 "Solids", "Formula_1stmonth", "Formula",
                 # "ATB_prev6m_corr",
                 "inf_SkinContactAge", 
                 "Penicillins", "m_AgeDelivery", "m_BMI", 
                 # "inf_IsAway3m", 
                 "inf_HospitalLocal", 
                 "familly_TimeSpendNonParents3m", 
                 # "inf_HasHospitalization_past3m", 
                 "inf_cummulative_NbHospi", 
                 "inf_Sex", "SUM", 
                 "chao1", "Cluster")  

model_list16 <- list()
AIC_list16 <- list()

### model selection#####
# DIY drop 1
# loop without the technical variable
for (i in 2:length(mod_vars16)) {
  curr_vars <- paste(mod_vars16[-i], collapse = "+")
  tmp_model <- adonis2(as.formula(paste("dist_b5", "~", curr_vars)),
                       meta_b5, 
                       permutations = 99
  )
  tmp_AIC <- AICc.PERMANOVA2(tmp_model, meta_b5)$AICc
  
  # save models to list
  model_list16[[i]] <- tmp_model
  AIC_list16[[i]] <- tmp_AIC
  new_AIC <- Reduce(min, AIC_list16[2:length(AIC_list16)])
  # if new model has AIC of 2 of more smaller, update model and variable list (and use new AIC in next round)
}


AIC_list16
indx <- which(AIC_list16[2:length(AIC_list16)] == min(unlist(AIC_list16[2:length(AIC_list16)]))) 
indx

model_list16[2:length(model_list16)][[indx]]

# done


if (curr_AIC - new_AIC >= 1.5) {
  curr_model <- model_list16[2:length(model_list16)][[indx]]
  curr_AIC <- new_AIC
}

final_ado_b5 <- adonis2(dist_b5 ~ MGI_Plate + inf_DeliveryMode + Solids_ever +
                          # Maternal_ab + 
                          # Paternal_ab + #IAP_substance + 
                          # birth_season + Solids_quantity +
                          # inf_Gestational_term +
                          # env_FurOrFeathers + 
                          inf_AgeHome + m_PreviousDeliveries + 
                          # inf_TravelAbroad3m + #breastfeeding + 
                          ExclusiveBF + #genotype +
                          # m_Education + familly_NbCaregivers3m_Cat + 
                          Solids + Formula_1stmonth + Formula + 
                          # ATB_prev6m_corr + 
                          inf_SkinContactAge + 
                          Penicillins + m_AgeDelivery + m_BMI + 
                          # inf_IsAway3m +
                          inf_HospitalLocal + 
                          familly_TimeSpendNonParents3m + 
                          # inf_HasHospitalization_past3m + 
                          inf_cummulative_NbHospi +
                          inf_Sex + SUM + 
                          chao1 + Cluster,
                        meta_b5,
                        permutations = 999) 


# B7 #####

skim(meta[meta$Sample_type == "B7",])

# B7 ####
meta_b7 <- meta_df[meta_df$Sample_type == "B7",]
meta_b7 <- meta_b7 %>%
  select_if(~ sum(!is.na(.)) > 0.9 * 475)
meta_b7 <- na.omit(meta_b7)

summarytools::dfSummary(meta_b7)

targ_b7 <- targ[grepl("B7", rownames(targ)),]
targ_b7 <- targ_b7[,colSums(targ_b7) > 0]

targ_b7 <- targ_b7[rownames(targ_b7) %in% meta_b7$Sample_ID,]

identical(rownames(targ_b7), meta_df$Sample_ID) #?
all(rownames(targ_b7) == meta_b7$Sample_ID) # OK


dist_b7 <- rtarg_b7 %>%
  vegdist(method = "bray")

ado_b7 <- adonis2(dist_b7 ~ MGI_Plate + 
                    Maternal_ab + Paternal_ab + IAP_substance + 
                    familly_NbCaregivers3m_Cat + inf_DeliveryMode + 
                    inf_Gestational_term + 
                    birth_season + inf_AgeHome + 
                    Formula + Formula_1stmonth + 
                    m_PreviousDeliveries + 
                    m_BMI + 
                    ExclusiveBF_3m + breastfeeding + 
                    genotype +
                    inf_HospitalLocal +
                    m_Education +  
                    familly_TimeSpendNonParents3m + 
                    env_FurOrFeathers + 
                    inf_HasHospitalization_past3m + 
                    env_NbHoursDayCare3m +
                    ATB_prev6m_corr + inf_TravelAbroad3m + 
                    m_AgeDelivery + 
                    inf_SkinContactAge +
                    Penicillins + inf_IsAway3m +
                    Aminoglycosides + Ext_BristolScore +
                    DayCare + SUM + 
                    inf_Sex + 
                    chao1 + Cluster,
                  meta_b7,
                  permutations = 99) 

AICc.PERMANOVA2(ado_b7, meta_b7)$small
# Small, use AICc
curr_AIC <- AICc.PERMANOVA2(ado_b7, meta_b7)$AICc


mod_vars <- c("MGI_Plate", "Maternal_ab", "Paternal_ab",
              "IAP_substance", "familly_NbCaregivers3m_Cat",
              "inf_DeliveryMode", "inf_Gestational_term", "birth_season",
              "inf_AgeHome", "Formula", "Formula_1stmonth", 
              "m_PreviousDeliveries", "m_BMI", "ExclusiveBF_3m",
              "breastfeeding", "genotype",
              "inf_HospitalLocal", "m_Education",
              "familly_TimeSpendNonParents3m", "env_FurOrFeathers",
              "inf_HasHospitalization_past3m", "env_NbHoursDayCare3m",
              "ATB_prev6m_corr", "inf_TravelAbroad3m", "m_AgeDelivery",
              "inf_SkinContactAge", "Penicillins", "inf_IsAway3m",
              "Aminoglycosides", "Ext_BristolScore",
              "DayCare", "SUM", "inf_Sex", "chao1", "Cluster") 
curr_model <- ado_b7


model_list <- list()
AIC_list <- list()

### model selection#####
# DIY drop 1
# loop without the technical variable
for (i in 2:length(mod_vars)) {
  curr_vars <- paste(mod_vars[-i], collapse = "+")
  tmp_model <- adonis2(as.formula(paste("dist_b7", "~", curr_vars)),
                       meta_b7, 
                       permutations = 99
  )
  tmp_AIC <- AICc.PERMANOVA2(tmp_model, meta_b7)$AICc
  
  # save models to list
  model_list[[i]] <- tmp_model
  AIC_list[[i]] <- tmp_AIC
  new_AIC <- Reduce(min, AIC_list[2:length(AIC_list)])
  # if new model has AIC of 2 of more smaller, update model and variable list (and use new AIC in next round)
}


AIC_list
indx <- which(AIC_list[2:length(AIC_list)] == min(unlist(AIC_list[2:length(AIC_list)]))) 
indx

model_list[2:length(model_list)][[indx]]

# selection 

if (curr_AIC - new_AIC >= 2) {
  curr_model <- model_list[2:length(model_list)][[indx]]
  curr_AIC <- new_AIC
}


# next round
mod_vars2 <- c("MGI_Plate", "Maternal_ab", "Paternal_ab",
               "IAP_substance", "familly_NbCaregivers3m_Cat",
               "inf_DeliveryMode", #"inf_Gestational_term",
               "birth_season",
               "inf_AgeHome", "Formula", "Formula_1stmonth", 
               "m_PreviousDeliveries", "m_BMI", "ExclusiveBF_3m",
               "breastfeeding", "genotype",
               "inf_HospitalLocal", "m_Education",
               "familly_TimeSpendNonParents3m", "env_FurOrFeathers",
               "inf_HasHospitalization_past3m", "env_NbHoursDayCare3m",
               "ATB_prev6m_corr", "inf_TravelAbroad3m", "m_AgeDelivery",
               "inf_SkinContactAge", "Penicillins", "inf_IsAway3m",
               "Aminoglycosides", "Ext_BristolScore",
               "DayCare", "SUM", "inf_Sex", "chao1", "Cluster") 

model_list2 <- list()
AIC_list2 <- list()

### model selection#####
# DIY drop 1
# loop without the technical variable
for (i in 2:length(mod_vars2)) {
  curr_vars <- paste(mod_vars2[-i], collapse = "+")
  tmp_model <- adonis2(as.formula(paste("dist_b7", "~", curr_vars)),
                       meta_b7, 
                       permutations = 99
  )
  tmp_AIC <- AICc.PERMANOVA2(tmp_model, meta_b7)$AICc
  
  # save models to list
  model_list2[[i]] <- tmp_model
  AIC_list2[[i]] <- tmp_AIC
  new_AIC <- Reduce(min, AIC_list2[2:length(AIC_list2)])
  # if new model has AIC of 2 of more smaller, update model and variable list (and use new AIC in next round)
}


AIC_list2
indx <- which(AIC_list2[2:length(AIC_list2)] == min(unlist(AIC_list2[2:length(AIC_list2)]))) 
indx

model_list2[2:length(model_list2)][[indx]]

# clear lowest


if (curr_AIC - new_AIC >= 2) {
  curr_model <- model_list2[2:length(model_list2)][[indx]]
  curr_AIC <- new_AIC
}

# next round
mod_vars3 <- c("MGI_Plate", "Maternal_ab", "Paternal_ab",
               # "IAP_substance", Removed for being obsolete
               "familly_NbCaregivers3m_Cat",
               "inf_DeliveryMode", #"inf_Gestational_term", 
               # "birth_season", Removed for being obsolete
               "inf_AgeHome", "Formula", "Formula_1stmonth", 
               "m_PreviousDeliveries", "m_BMI", "ExclusiveBF_3m",
               "breastfeeding", 
               # "genotype", Removed for being obsolete
               "inf_HospitalLocal", 
               # "m_Education", Removed for being obsolete
               "familly_TimeSpendNonParents3m", "env_FurOrFeathers",
               "inf_HasHospitalization_past3m", "env_NbHoursDayCare3m",
               "ATB_prev6m_corr", "inf_TravelAbroad3m", "m_AgeDelivery",
               "inf_SkinContactAge", "Penicillins", "inf_IsAway3m",
               "Aminoglycosides", "Ext_BristolScore",
               "DayCare", "SUM", "inf_Sex", "chao1", "Cluster") 

model_list3 <- list()
AIC_list3 <- list()

### model selection#####
# DIY drop 1
# loop without the technical variable
for (i in 2:length(mod_vars3)) {
  curr_vars <- paste(mod_vars3[-i], collapse = "+")
  tmp_model <- adonis2(as.formula(paste("dist_b7", "~", curr_vars)),
                       meta_b7, 
                       permutations = 99
  )
  tmp_AIC <- AICc.PERMANOVA2(tmp_model, meta_b7)$AICc
  
  # save models to list
  model_list3[[i]] <- tmp_model
  AIC_list3[[i]] <- tmp_AIC
  new_AIC <- Reduce(min, AIC_list3[2:length(AIC_list3)])
  # if new model has AIC of 2 of more smaller, update model and variable list (and use new AIC in next round)
}


AIC_list3
indx <- which(AIC_list3[2:length(AIC_list3)] == min(unlist(AIC_list3[2:length(AIC_list3)]))) 
indx

model_list3[2:length(model_list3)][[indx]]

# select


if (curr_AIC - new_AIC >= 2) {
  curr_model <- model_list3[2:length(model_list3)][[indx]]
  curr_AIC <- new_AIC
}

# next round
mod_vars4 <- c("MGI_Plate", "Maternal_ab", "Paternal_ab",
               # "IAP_substance", Removed for being obsolete
               "familly_NbCaregivers3m_Cat",
               "inf_DeliveryMode", #"inf_Gestational_term", 
               # "birth_season", Removed for being obsolete
               "inf_AgeHome", "Formula", #"Formula_1stmonth", 
               "m_PreviousDeliveries", "m_BMI", "ExclusiveBF_3m",
               "breastfeeding", 
               # "genotype", Removed for being obsolete
               "inf_HospitalLocal", 
               # "m_Education", Removed for being obsolete
               "familly_TimeSpendNonParents3m", "env_FurOrFeathers",
               "inf_HasHospitalization_past3m", "env_NbHoursDayCare3m",
               "ATB_prev6m_corr", "inf_TravelAbroad3m", "m_AgeDelivery",
               "inf_SkinContactAge", "Penicillins", "inf_IsAway3m",
               "Aminoglycosides", "Ext_BristolScore",
               "DayCare", "SUM", "inf_Sex", "chao1", "Cluster") 

model_list4 <- list()
AIC_list4 <- list()

### model selection#####
# DIY drop 1
# loop without the technical variable
for (i in 2:length(mod_vars4)) {
  curr_vars <- paste(mod_vars4[-i], collapse = "+")
  tmp_model <- adonis2(as.formula(paste("dist_b7", "~", curr_vars)),
                       meta_b7, 
                       permutations = 99
  )
  tmp_AIC <- AICc.PERMANOVA2(tmp_model, meta_b7)$AICc
  
  # save models to list
  model_list4[[i]] <- tmp_model
  AIC_list4[[i]] <- tmp_AIC
  new_AIC <- Reduce(min, AIC_list4[2:length(AIC_list4)])
  # if new model has AIC of 2 of more smaller, update model and variable list (and use new AIC in next round)
}


AIC_list4
indx <- which(AIC_list4[2:length(AIC_list4)] == min(unlist(AIC_list4[2:length(AIC_list4)]))) 
indx

model_list4[2:length(model_list4)][[indx]]

# select


if (curr_AIC - new_AIC >= 2) {
  curr_model <- model_list4[2:length(model_list4)][[indx]]
  curr_AIC <- new_AIC
}

# next round
mod_vars5 <- c("MGI_Plate", "Maternal_ab", "Paternal_ab",
               # "IAP_substance", Removed for being obsolete
               "familly_NbCaregivers3m_Cat",
               "inf_DeliveryMode", #"inf_Gestational_term", 
               # "birth_season", Removed for being obsolete
               "inf_AgeHome", "Formula", #"Formula_1stmonth", 
               "m_PreviousDeliveries", "m_BMI", "ExclusiveBF_3m",
               "breastfeeding", 
               # "genotype", Removed for being obsolete
               "inf_HospitalLocal", 
               # "m_Education", Removed for being obsolete
               "familly_TimeSpendNonParents3m", #"env_FurOrFeathers",
               "inf_HasHospitalization_past3m", "env_NbHoursDayCare3m",
               "ATB_prev6m_corr", "inf_TravelAbroad3m", "m_AgeDelivery",
               "inf_SkinContactAge", "Penicillins", "inf_IsAway3m",
               "Aminoglycosides", "Ext_BristolScore",
               "DayCare", "SUM", "inf_Sex", "chao1", "Cluster") 

model_list5 <- list()
AIC_list5 <- list()

### model selection#####
# DIY drop 1
# loop without the technical variable
for (i in 2:length(mod_vars5)) {
  curr_vars <- paste(mod_vars5[-i], collapse = "+")
  tmp_model <- adonis2(as.formula(paste("dist_b7", "~", curr_vars)),
                       meta_b7, 
                       permutations = 99
  )
  tmp_AIC <- AICc.PERMANOVA2(tmp_model, meta_b7)$AICc
  
  # save models to list
  model_list5[[i]] <- tmp_model
  AIC_list5[[i]] <- tmp_AIC
  new_AIC <- Reduce(min, AIC_list5[2:length(AIC_list5)])
  # if new model has AIC of 2 of more smaller, update model and variable list (and use new AIC in next round)
}


AIC_list5
indx <- which(AIC_list5[2:length(AIC_list5)] == min(unlist(AIC_list5[2:length(AIC_list5)]))) 
indx

model_list5[2:length(model_list5)][[indx]]

# stricter selection


if (curr_AIC - new_AIC >= 1.5) {
  curr_model <- model_list5[2:length(model_list5)][[indx]]
  curr_AIC <- new_AIC
}

# next round
mod_vars6 <- c("MGI_Plate", "Maternal_ab", "Paternal_ab",
               # "IAP_substance", Removed for being obsolete
               "familly_NbCaregivers3m_Cat",
               "inf_DeliveryMode", #"inf_Gestational_term", 
               # "birth_season", Removed for being obsolete
               "inf_AgeHome", "Formula", #"Formula_1stmonth", 
               "m_PreviousDeliveries", "m_BMI", "ExclusiveBF_3m",
               "breastfeeding", 
               # "genotype", Removed for being obsolete
               "inf_HospitalLocal", 
               # "m_Education", Removed for being obsolete
               "familly_TimeSpendNonParents3m", #"env_FurOrFeathers",
               "inf_HasHospitalization_past3m", "env_NbHoursDayCare3m",
               "ATB_prev6m_corr", "inf_TravelAbroad3m", "m_AgeDelivery",
               "inf_SkinContactAge", "Penicillins", "inf_IsAway3m",
               "Aminoglycosides", "Ext_BristolScore",
               "DayCare", "SUM", #"inf_Sex", 
               "chao1", "Cluster") 


model_list6 <- list()
AIC_list6 <- list()

### model selection#####
# DIY drop 1
# loop without the technical variable
for (i in 2:length(mod_vars6)) {
  curr_vars <- paste(mod_vars6[-i], collapse = "+")
  tmp_model <- adonis2(as.formula(paste("dist_b7", "~", curr_vars)),
                       meta_b7, 
                       permutations = 99
  )
  tmp_AIC <- AICc.PERMANOVA2(tmp_model, meta_b7)$AICc
  
  # save models to list
  model_list6[[i]] <- tmp_model
  AIC_list6[[i]] <- tmp_AIC
  new_AIC <- Reduce(min, AIC_list6[2:length(AIC_list6)])
  # if new model has AIC of 2 of more smaller, update model and variable list (and use new AIC in next round)
}


AIC_list6
indx <- which(AIC_list6[2:length(AIC_list6)] == min(unlist(AIC_list6[2:length(AIC_list6)]))) 
indx

model_list6[2:length(model_list6)][[indx]]

# stricter selection

if (curr_AIC - new_AIC >= 1.5) {
  curr_model <- model_list6[2:length(model_list6)][[indx]]
  curr_AIC <- new_AIC
}

# next round
mod_vars7 <- c("MGI_Plate", #"Maternal_ab",
               "Paternal_ab",
               # "IAP_substance", Removed for being obsolete
               "familly_NbCaregivers3m_Cat",
               "inf_DeliveryMode", #"inf_Gestational_term", 
               # "birth_season", Removed for being obsolete
               "inf_AgeHome", "Formula", #"Formula_1stmonth", 
               "m_PreviousDeliveries", "m_BMI", "ExclusiveBF_3m",
               "breastfeeding", 
               # "genotype", Removed for being obsolete
               "inf_HospitalLocal", 
               # "m_Education", Removed for being obsolete
               "familly_TimeSpendNonParents3m", #"env_FurOrFeathers",
               "inf_HasHospitalization_past3m", "env_NbHoursDayCare3m",
               "ATB_prev6m_corr", "inf_TravelAbroad3m", "m_AgeDelivery",
               "inf_SkinContactAge", "Penicillins", "inf_IsAway3m",
               "Aminoglycosides", "Ext_BristolScore",
               "DayCare", "SUM", #"inf_Sex", 
               "chao1", "Cluster") 

model_list7 <- list()
AIC_list7 <- list()

### model selection#####
# DIY drop 1
# loop without the technical variable
for (i in 2:length(mod_vars7)) {
  curr_vars <- paste(mod_vars7[-i], collapse = "+")
  tmp_model <- adonis2(as.formula(paste("dist_b7", "~", curr_vars)),
                       meta_b7, 
                       permutations = 99
  )
  tmp_AIC <- AICc.PERMANOVA2(tmp_model, meta_b7)$AICc
  
  # save models to list
  model_list7[[i]] <- tmp_model
  AIC_list7[[i]] <- tmp_AIC
  new_AIC <- Reduce(min, AIC_list7[2:length(AIC_list7)])
  # if new model has AIC of 2 of more smaller, update model and variable list (and use new AIC in next round)
}


AIC_list7
indx <- which(AIC_list7[2:length(AIC_list7)] == min(unlist(AIC_list7[2:length(AIC_list7)]))) 
indx

model_list7[2:length(model_list7)][[indx]]

# selection

if (curr_AIC - new_AIC >= 2) {
  curr_model <- model_list7[2:length(model_list7)][[indx]]
  curr_AIC <- new_AIC
}

# next round
mod_vars8 <- c("MGI_Plate", #"Maternal_ab",
               # "Paternal_ab",
               # "IAP_substance", Removed for being obsolete
               "familly_NbCaregivers3m_Cat",
               "inf_DeliveryMode", #"inf_Gestational_term", 
               # "birth_season", Removed for being obsolete
               "inf_AgeHome", "Formula", #"Formula_1stmonth", 
               "m_PreviousDeliveries", "m_BMI", "ExclusiveBF_3m",
               "breastfeeding", 
               # "genotype", Removed for being obsolete
               "inf_HospitalLocal", 
               # "m_Education", Removed for being obsolete
               "familly_TimeSpendNonParents3m", #"env_FurOrFeathers",
               "inf_HasHospitalization_past3m", "env_NbHoursDayCare3m",
               "ATB_prev6m_corr", "inf_TravelAbroad3m", "m_AgeDelivery",
               "inf_SkinContactAge", "Penicillins", "inf_IsAway3m",
               "Aminoglycosides", "Ext_BristolScore",
               "DayCare", "SUM", #"inf_Sex", 
               "chao1", "Cluster") 

model_list8 <- list()
AIC_list8 <- list()

### model selection#####
# DIY drop 1
# loop without the technical variable
for (i in 2:length(mod_vars8)) {
  curr_vars <- paste(mod_vars8[-i], collapse = "+")
  tmp_model <- adonis2(as.formula(paste("dist_b7", "~", curr_vars)),
                       meta_b7, 
                       permutations = 99
  )
  tmp_AIC <- AICc.PERMANOVA2(tmp_model, meta_b7)$AICc
  
  # save models to list
  model_list8[[i]] <- tmp_model
  AIC_list8[[i]] <- tmp_AIC
  new_AIC <- Reduce(min, AIC_list8[2:length(AIC_list8)])
  # if new model has AIC of 2 of more smaller, update model and variable list (and use new AIC in next round)
}


AIC_list8
indx <- which(AIC_list8[2:length(AIC_list8)] == min(unlist(AIC_list8[2:length(AIC_list8)]))) 
indx

model_list8[2:length(model_list8)][[indx]]

# stricter selection 


if (curr_AIC - new_AIC >= 1.5) {
  curr_model <- model_list8[2:length(model_list8)][[indx]]
  curr_AIC <- new_AIC
}

# next round
mod_vars9 <- c("MGI_Plate", #"Maternal_ab",
               # "Paternal_ab",
               # "IAP_substance", Removed for being obsolete
               "familly_NbCaregivers3m_Cat",
               "inf_DeliveryMode", #"inf_Gestational_term", 
               # "birth_season", Removed for being obsolete
               "inf_AgeHome", "Formula", #"Formula_1stmonth", 
               "m_PreviousDeliveries", "m_BMI", "ExclusiveBF_3m",
               "breastfeeding", 
               # "genotype", Removed for being obsolete
               "inf_HospitalLocal", 
               # "m_Education", Removed for being obsolete
               "familly_TimeSpendNonParents3m", #"env_FurOrFeathers",
               "inf_HasHospitalization_past3m", "env_NbHoursDayCare3m",
               "ATB_prev6m_corr", "inf_TravelAbroad3m", "m_AgeDelivery",
               "inf_SkinContactAge", "Penicillins", "inf_IsAway3m",
               "Aminoglycosides", "Ext_BristolScore",
               # "DayCare",
               "SUM", #"inf_Sex", 
               "chao1", "Cluster")  

model_list9 <- list()
AIC_list9 <- list()

### model selection#####
# DIY drop 1
# loop without the technical variable
for (i in 2:length(mod_vars9)) {
  curr_vars <- paste(mod_vars9[-i], collapse = "+")
  tmp_model <- adonis2(as.formula(paste("dist_b7", "~", curr_vars)),
                       meta_b7, 
                       permutations = 99
  )
  tmp_AIC <- AICc.PERMANOVA2(tmp_model, meta_b7)$AICc
  
  # save models to list
  model_list9[[i]] <- tmp_model
  AIC_list9[[i]] <- tmp_AIC
  new_AIC <- Reduce(min, AIC_list9[2:length(AIC_list9)])
  # if new model has AIC of 2 of more smaller, update model and variable list (and use new AIC in next round)
}


AIC_list9
indx <- which(AIC_list9[2:length(AIC_list9)] == min(unlist(AIC_list9[2:length(AIC_list9)]))) 
indx

model_list9[2:length(model_list9)][[indx]]

# stricter AIC criteria

if (curr_AIC - new_AIC >= 1.5) {
  curr_model <- model_list9[2:length(model_list9)][[indx]]
  curr_AIC <- new_AIC
}


# next round
mod_vars10 <- c("MGI_Plate", #"Maternal_ab",
                # "Paternal_ab",
                # "IAP_substance", Removed for being obsolete
                "familly_NbCaregivers3m_Cat",
                "inf_DeliveryMode", #"inf_Gestational_term", 
                # "birth_season", Removed for being obsolete
                "inf_AgeHome", "Formula", #"Formula_1stmonth", 
                # "m_PreviousDeliveries", 
                "m_BMI", "ExclusiveBF_3m",
                "breastfeeding", 
                # "genotype", Removed for being obsolete
                "inf_HospitalLocal", 
                # "m_Education", Removed for being obsolete
                "familly_TimeSpendNonParents3m", #"env_FurOrFeathers",
                "inf_HasHospitalization_past3m", "env_NbHoursDayCare3m",
                "ATB_prev6m_corr", "inf_TravelAbroad3m", "m_AgeDelivery",
                "inf_SkinContactAge", "Penicillins", "inf_IsAway3m",
                "Aminoglycosides", "Ext_BristolScore",
                # "DayCare",
                "SUM", #"inf_Sex", 
                "chao1", "Cluster") 

model_list10 <- list()
AIC_list10 <- list()

### model selection#####
# DIY drop 1
# loop without the technical variable
for (i in 2:length(mod_vars10)) {
  curr_vars <- paste(mod_vars10[-i], collapse = "+")
  tmp_model <- adonis2(as.formula(paste("dist_b7", "~", curr_vars)),
                       meta_b7, 
                       permutations = 99
  )
  tmp_AIC <- AICc.PERMANOVA2(tmp_model, meta_b7)$AICc
  
  # save models to list
  model_list10[[i]] <- tmp_model
  AIC_list10[[i]] <- tmp_AIC
  new_AIC <- Reduce(min, AIC_list10[2:length(AIC_list10)])
  # if new model has AIC of 2 of more smaller, update model and variable list (and use new AIC in next round)
}


AIC_list10
indx <- which(AIC_list10[2:length(AIC_list10)] == min(unlist(AIC_list10[2:length(AIC_list10)]))) 
indx

model_list10[2:length(model_list10)][[indx]]

# stricter AIC criteria

if (curr_AIC - new_AIC >= 1.5) {
  curr_model <- model_list10[2:length(model_list10)][[indx]]
  curr_AIC <- new_AIC
}

# next round
mod_vars11 <- c("MGI_Plate", #"Maternal_ab",
                # "Paternal_ab",
                # "IAP_substance", Removed for being obsolete
                "familly_NbCaregivers3m_Cat",
                "inf_DeliveryMode", #"inf_Gestational_term", 
                # "birth_season", Removed for being obsolete
                "inf_AgeHome", #"Formula", #"Formula_1stmonth", 
                # "m_PreviousDeliveries", 
                "m_BMI", "ExclusiveBF_3m",
                "breastfeeding", 
                # "genotype", Removed for being obsolete
                "inf_HospitalLocal", 
                # "m_Education", Removed for being obsolete
                "familly_TimeSpendNonParents3m", #"env_FurOrFeathers",
                "inf_HasHospitalization_past3m", "env_NbHoursDayCare3m",
                "ATB_prev6m_corr", "inf_TravelAbroad3m", "m_AgeDelivery",
                "inf_SkinContactAge", "Penicillins", "inf_IsAway3m",
                "Aminoglycosides", "Ext_BristolScore",
                # "DayCare",
                "SUM", #"inf_Sex", 
                "chao1", "Cluster") 

model_list11 <- list()
AIC_list11 <- list()

### model selection#####
# DIY drop 1
# loop without the technical variable
for (i in 2:length(mod_vars11)) {
  curr_vars <- paste(mod_vars11[-i], collapse = "+")
  tmp_model <- adonis2(as.formula(paste("dist_b7", "~", curr_vars)),
                       meta_b7, 
                       permutations = 99
  )
  tmp_AIC <- AICc.PERMANOVA2(tmp_model, meta_b7)$AICc
  
  # save models to list
  model_list11[[i]] <- tmp_model
  AIC_list11[[i]] <- tmp_AIC
  new_AIC <- Reduce(min, AIC_list11[2:length(AIC_list11)])
  # if new model has AIC of 2 of more smaller, update model and variable list (and use new AIC in next round)
}


AIC_list11
indx <- which(AIC_list11[2:length(AIC_list11)] == min(unlist(AIC_list11[2:length(AIC_list11)]))) 
indx

model_list11[2:length(model_list11)][[indx]]

# stricter AIC criteria


if (curr_AIC - new_AIC >= 1.5) {
  curr_model <- model_list11[2:length(model_list11)][[indx]]
  curr_AIC <- new_AIC
}


# next round
mod_vars12 <- c("MGI_Plate", #"Maternal_ab",
                # "Paternal_ab",
                # "IAP_substance", Removed for being obsolete
                "familly_NbCaregivers3m_Cat",
                "inf_DeliveryMode", #"inf_Gestational_term", 
                # "birth_season", Removed for being obsolete
                "inf_AgeHome", #"Formula", #"Formula_1stmonth", 
                # "m_PreviousDeliveries", 
                "m_BMI", "ExclusiveBF_3m",
                "breastfeeding", 
                # "genotype", Removed for being obsolete
                "inf_HospitalLocal", 
                # "m_Education", Removed for being obsolete
                "familly_TimeSpendNonParents3m", #"env_FurOrFeathers",
                "inf_HasHospitalization_past3m", "env_NbHoursDayCare3m",
                "ATB_prev6m_corr",# "inf_TravelAbroad3m", 
                "m_AgeDelivery",
                "inf_SkinContactAge", "Penicillins", "inf_IsAway3m",
                "Aminoglycosides", "Ext_BristolScore",
                # "DayCare",
                "SUM", #"inf_Sex", 
                "chao1", "Cluster") 

model_list12 <- list()
AIC_list12 <- list()

### model selection#####
# DIY drop 1
# loop without the technical variable
for (i in 2:length(mod_vars12)) {
  curr_vars <- paste(mod_vars12[-i], collapse = "+")
  tmp_model <- adonis2(as.formula(paste("dist_b7", "~", curr_vars)),
                       meta_b7, 
                       permutations = 99
  )
  tmp_AIC <- AICc.PERMANOVA2(tmp_model, meta_b7)$AICc
  
  # save models to list
  model_list12[[i]] <- tmp_model
  AIC_list12[[i]] <- tmp_AIC
  new_AIC <- Reduce(min, AIC_list12[2:length(AIC_list12)])
  # if new model has AIC of 2 of more smaller, update model and variable list (and use new AIC in next round)
}


AIC_list12
indx <- which(AIC_list12[2:length(AIC_list12)] == min(unlist(AIC_list12[2:length(AIC_list12)]))) 
indx

model_list12[2:length(model_list12)][[indx]]

# stricter AIC criteria


if (curr_AIC - new_AIC >= 1.5) {
  curr_model <- model_list12[2:length(model_list12)][[indx]]
  curr_AIC <- new_AIC
}

# next round
mod_vars13 <- c("MGI_Plate", #"Maternal_ab",
                # "Paternal_ab",
                # "IAP_substance", Removed for being obsolete
                "familly_NbCaregivers3m_Cat",
                "inf_DeliveryMode", #"inf_Gestational_term", 
                # "birth_season", Removed for being obsolete
                "inf_AgeHome", #"Formula", #"Formula_1stmonth", 
                # "m_PreviousDeliveries", 
                "m_BMI", "ExclusiveBF_3m",
                "breastfeeding", 
                # "genotype", Removed for being obsolete
                "inf_HospitalLocal", 
                # "m_Education", Removed for being obsolete
                "familly_TimeSpendNonParents3m", #"env_FurOrFeathers",
                "inf_HasHospitalization_past3m", "env_NbHoursDayCare3m",
                "ATB_prev6m_corr",# "inf_TravelAbroad3m", 
                "m_AgeDelivery",
                "inf_SkinContactAge", "Penicillins", "inf_IsAway3m",
                # "Aminoglycosides", 
                "Ext_BristolScore",
                # "DayCare",
                "SUM", #"inf_Sex", 
                "chao1", "Cluster") 

model_list13 <- list()
AIC_list13 <- list()

### model selection#####
# DIY drop 1
# loop without the technical variable
for (i in 2:length(mod_vars13)) {
  curr_vars <- paste(mod_vars13[-i], collapse = "+")
  tmp_model <- adonis2(as.formula(paste("dist_b7", "~", curr_vars)),
                       meta_b7, 
                       permutations = 99
  )
  tmp_AIC <- AICc.PERMANOVA2(tmp_model, meta_b7)$AICc
  
  # save models to list
  model_list13[[i]] <- tmp_model
  AIC_list13[[i]] <- tmp_AIC
  new_AIC <- Reduce(min, AIC_list13[2:length(AIC_list13)])
  # if new model has AIC of 2 of more smaller, update model and variable list (and use new AIC in next round)
}


AIC_list13
indx <- which(AIC_list13[2:length(AIC_list13)] == min(unlist(AIC_list13[2:length(AIC_list13)]))) 
indx

model_list13[2:length(model_list13)][[indx]]

# stricter AIC criteria


if (curr_AIC - new_AIC >= 1.5) {
  curr_model <- model_list13[2:length(model_list13)][[indx]]
  curr_AIC <- new_AIC
}


# next round
mod_vars14 <- c("MGI_Plate", #"Maternal_ab",
                # "Paternal_ab",
                # "IAP_substance", Removed for being obsolete
                "familly_NbCaregivers3m_Cat",
                "inf_DeliveryMode", #"inf_Gestational_term", 
                # "birth_season", Removed for being obsolete
                # "inf_AgeHome", #"Formula", #"Formula_1stmonth", 
                # "m_PreviousDeliveries", 
                "m_BMI", "ExclusiveBF_3m",
                "breastfeeding", 
                # "genotype", Removed for being obsolete
                "inf_HospitalLocal", 
                # "m_Education", Removed for being obsolete
                "familly_TimeSpendNonParents3m", #"env_FurOrFeathers",
                "inf_HasHospitalization_past3m", "env_NbHoursDayCare3m",
                "ATB_prev6m_corr",# "inf_TravelAbroad3m", 
                "m_AgeDelivery",
                "inf_SkinContactAge", "Penicillins", "inf_IsAway3m",
                # "Aminoglycosides", 
                "Ext_BristolScore",
                # "DayCare",
                "SUM", #"inf_Sex", 
                "chao1", "Cluster") 

model_list14 <- list()
AIC_list14 <- list()

### model selection#####
# DIY drop 1
# loop without the technical variable
for (i in 2:length(mod_vars14)) {
  curr_vars <- paste(mod_vars14[-i], collapse = "+")
  tmp_model <- adonis2(as.formula(paste("dist_b7", "~", curr_vars)),
                       meta_b7, 
                       permutations = 99
  )
  tmp_AIC <- AICc.PERMANOVA2(tmp_model, meta_b7)$AICc
  
  # save models to list
  model_list14[[i]] <- tmp_model
  AIC_list14[[i]] <- tmp_AIC
  new_AIC <- Reduce(min, AIC_list14[2:length(AIC_list14)])
  # if new model has AIC of 2 of more smaller, update model and variable list (and use new AIC in next round)
}


AIC_list14
indx <- which(AIC_list14[2:length(AIC_list14)] == min(unlist(AIC_list14[2:length(AIC_list14)]))) 
indx

model_list14[2:length(model_list14)][[indx]]

# stricter AIC criteria


if (curr_AIC - new_AIC >= 1.5) {
  curr_model <- model_list14[2:length(model_list14)][[indx]]
  curr_AIC <- new_AIC
}

# next round
mod_vars15 <- c("MGI_Plate", #"Maternal_ab",
                # "Paternal_ab",
                # "IAP_substance", Removed for being obsolete
                "familly_NbCaregivers3m_Cat",
                "inf_DeliveryMode", #"inf_Gestational_term", 
                # "birth_season", Removed for being obsolete
                # "inf_AgeHome", #"Formula", #"Formula_1stmonth", 
                # "m_PreviousDeliveries", 
                "m_BMI", #"ExclusiveBF_3m",
                "breastfeeding", 
                # "genotype", Removed for being obsolete
                "inf_HospitalLocal", 
                # "m_Education", Removed for being obsolete
                "familly_TimeSpendNonParents3m", #"env_FurOrFeathers",
                "inf_HasHospitalization_past3m", "env_NbHoursDayCare3m",
                "ATB_prev6m_corr",# "inf_TravelAbroad3m", 
                "m_AgeDelivery",
                "inf_SkinContactAge", "Penicillins", "inf_IsAway3m",
                # "Aminoglycosides", 
                "Ext_BristolScore",
                # "DayCare",
                "SUM", #"inf_Sex", 
                "chao1", "Cluster") 

model_list15 <- list()
AIC_list15 <- list()

### model selection#####
# DIY drop 1
# loop without the technical variable
for (i in 2:length(mod_vars15)) {
  curr_vars <- paste(mod_vars15[-i], collapse = "+")
  tmp_model <- adonis2(as.formula(paste("dist_b7", "~", curr_vars)),
                       meta_b7, 
                       permutations = 99
  )
  tmp_AIC <- AICc.PERMANOVA2(tmp_model, meta_b7)$AICc
  
  # save models to list
  model_list15[[i]] <- tmp_model
  AIC_list15[[i]] <- tmp_AIC
  new_AIC <- Reduce(min, AIC_list15[2:length(AIC_list15)])
  # if new model has AIC of 2 of more smaller, update model and variable list (and use new AIC in next round)
}


AIC_list15
indx <- which(AIC_list15[2:length(AIC_list15)] == min(unlist(AIC_list15[2:length(AIC_list15)]))) 
indx

model_list15[2:length(model_list15)][[indx]]

# stricter AIC criteria


if (curr_AIC - new_AIC >= 1.5) {
  curr_model <- model_list15[2:length(model_list15)][[indx]]
  curr_AIC <- new_AIC
}

# next round
mod_vars16 <- c("MGI_Plate", #"Maternal_ab",
                # "Paternal_ab",
                # "IAP_substance", Removed for being obsolete
                "familly_NbCaregivers3m_Cat",
                "inf_DeliveryMode", #"inf_Gestational_term", 
                # "birth_season", Removed for being obsolete
                # "inf_AgeHome", #"Formula", #"Formula_1stmonth", 
                # "m_PreviousDeliveries", 
                "m_BMI", #"ExclusiveBF_3m",
                "breastfeeding", 
                # "genotype", Removed for being obsolete
                "inf_HospitalLocal", 
                # "m_Education", Removed for being obsolete
                "familly_TimeSpendNonParents3m", #"env_FurOrFeathers",
                "inf_HasHospitalization_past3m", "env_NbHoursDayCare3m",
                "ATB_prev6m_corr",# "inf_TravelAbroad3m",
                # "m_AgeDelivery",
                "inf_SkinContactAge", "Penicillins", "inf_IsAway3m",
                # "Aminoglycosides", 
                "Ext_BristolScore",
                # "DayCare",
                "SUM", #"inf_Sex", 
                "chao1", "Cluster") 

model_list16 <- list()
AIC_list16 <- list()

### model selection#####
# DIY drop 1
# loop without the technical variable
for (i in 2:length(mod_vars16)) {
  curr_vars <- paste(mod_vars16[-i], collapse = "+")
  tmp_model <- adonis2(as.formula(paste("dist_b7", "~", curr_vars)),
                       meta_b7, 
                       permutations = 99
  )
  tmp_AIC <- AICc.PERMANOVA2(tmp_model, meta_b7)$AICc
  
  # save models to list
  model_list16[[i]] <- tmp_model
  AIC_list16[[i]] <- tmp_AIC
  new_AIC <- Reduce(min, AIC_list16[2:length(AIC_list16)])
  # if new model has AIC of 2 of more smaller, update model and variable list (and use new AIC in next round)
}


AIC_list16
indx <- which(AIC_list16[2:length(AIC_list16)] == min(unlist(AIC_list16[2:length(AIC_list16)]))) 
indx

model_list16[2:length(model_list16)][[indx]]

# done


if (curr_AIC - new_AIC >= 1.5) {
  curr_model <- model_list16[2:length(model_list16)][[indx]]
  curr_AIC <- new_AIC
}

final_ado_b7 <- adonis2(dist_b7 ~ MGI_Plate + 
                          # Maternal_ab + Paternal_ab + IAP_substance + 
                          # familly_NbCaregivers3m_Cat + 
                          inf_DeliveryMode + 
                          # inf_Gestational_term + 
                          # birth_season +
                          # inf_AgeHome + 
                          # Formula +# Formula_1stmonth + 
                          # m_PreviousDeliveries + 
                          m_BMI + 
                          # ExclusiveBF_3m + 
                          breastfeeding + 
                          inf_HospitalLocal +
                          # m_Education +  
                          familly_TimeSpendNonParents3m + 
                          # env_FurOrFeathers + 
                          inf_HasHospitalization_past3m + 
                          env_NbHoursDayCare3m +
                          ATB_prev6m_corr + #inf_TravelAbroad3m + 
                          # m_AgeDelivery + 
                          inf_SkinContactAge +
                          Penicillins + inf_IsAway3m +
                          # Aminoglycosides +
                          Ext_BristolScore +
                          # DayCare +
                          SUM + 
                          # inf_Sex + 
                          chao1 + Cluster,
                        meta_b7,
                        permutations = 999) 


# B9 #####

skim(meta[meta$Sample_type == "B9",])

# B9 ####
meta_b9 <- meta_df[meta_df$Sample_type == "B9",]
meta_b9 <- meta_b9 %>%
  select_if(~ sum(!is.na(.)) > 0.9 * 475)
meta_b9 <- na.omit(meta_b9)

summarytools::dfSummary(meta_b9)

targ_b9 <- targ[grepl("B9", rownames(targ)),]
targ_b9 <- targ_b9[,colSums(targ_b9) > 0]

targ_b9 <- targ_b9[rownames(targ_b9) %in% meta_b9$Sample_ID,]

identical(rownames(targ_b9), meta_df$Sample_ID) #?
all(rownames(targ_b9) == meta_b9$Sample_ID) # OK

dist_b9 <- rtarg_b9 %>%
  vegdist(method = "bray")

ado_b9 <- adonis2(dist_b9 ~ MGI_Plate + 
                    inf_DeliveryMode + 
                    m_PreviousDeliveries + 
                    Formula_1stmonth + 
                    ExclusiveBF_3m + 
                    inf_SkinContactAge + Ext_BristolScore + 
                    inf_AgeHome + inf_HospitalLocal + 
                    m_AgeDelivery + env_FurOrFeathers + 
                    SUM + inf_TravelAbroad3m + 
                    ATB_prev6m_corr + Penicillins + 
                    inf_Sex + familly_NbCaregivers3m_Cat + 
                    inf_HasHospitalization_past3m + 
                    familly_TimeSpendNonParents3m + 
                    m_BMI + inf_IsAway3m + 
                    DayCare + breastfeeding + 
                    env_NbHoursDayCare3m + 
                    Aminoglycosides + 
                    Formula + 
                    chao1 + Cluster,
                  meta_b9,
                  permutations = 99) 

AICc.PERMANOVA2(ado_b9, meta_b9)$small
# Small, use AICc
curr_AIC <- AICc.PERMANOVA2(ado_b9, meta_b9)$AICc


mod_vars <- c("MGI_Plate", 
              "inf_DeliveryMode",
              "m_PreviousDeliveries",
              "Formula_1stmonth", "ExclusiveBF_3m",
              "inf_SkinContactAge", "Ext_BristolScore",
              "inf_AgeHome", "inf_HospitalLocal", 
              "m_AgeDelivery", "env_FurOrFeathers", 
              "SUM", "inf_TravelAbroad3m", 
              "ATB_prev6m_corr", "Penicillins", 
              "inf_Sex", "familly_NbCaregivers3m_Cat", 
              "inf_HasHospitalization_past3m", 
              "familly_TimeSpendNonParents3m", 
              "m_BMI", "inf_IsAway3m", 
              "DayCare", "breastfeeding", 
              "env_NbHoursDayCare3m", "Aminoglycosides", 
              "Formula", "chao1", "Cluster") 
curr_model <- ado_b9

model_list <- list()
AIC_list <- list()

### model selection#####
# DIY drop 1
# loop without the technical variable
for (i in 2:length(mod_vars)) {
  curr_vars <- paste(mod_vars[-i], collapse = "+")
  tmp_model <- adonis2(as.formula(paste("dist_b9", "~", curr_vars)),
                       meta_b9, 
                       permutations = 99
  )
  tmp_AIC <- AICc.PERMANOVA2(tmp_model, meta_b9)$AICc
  
  # save models to list
  model_list[[i]] <- tmp_model
  AIC_list[[i]] <- tmp_AIC
  new_AIC <- Reduce(min, AIC_list[2:length(AIC_list)])
  # if new model has AIC of 2 of more smaller, update model and variable list (and use new AIC in next round)
}


AIC_list
indx <- which(AIC_list[2:length(AIC_list)] == min(unlist(AIC_list[2:length(AIC_list)]))) 
indx

model_list[2:length(model_list)][[indx]]

# Go with selection


if (curr_AIC - new_AIC >= 2) {
  curr_model <- model_list[2:length(model_list)][[indx]]
  curr_AIC <- new_AIC
}

# next round
mod_vars2 <- c("MGI_Plate", 
               "inf_DeliveryMode",
               "m_PreviousDeliveries", 
               "Formula_1stmonth", "ExclusiveBF_3m",
               "inf_SkinContactAge", "Ext_BristolScore",
               "inf_AgeHome", "inf_HospitalLocal", 
               "m_AgeDelivery", "env_FurOrFeathers", 
               "SUM", "inf_TravelAbroad3m", 
               "ATB_prev6m_corr", "Penicillins", 
               "inf_Sex", #"familly_NbCaregivers3m_Cat", 
               "inf_HasHospitalization_past3m", 
               "familly_TimeSpendNonParents3m", 
               "m_BMI", "inf_IsAway3m", 
               "DayCare", "breastfeeding", 
               "env_NbHoursDayCare3m", "Aminoglycosides", 
               "Formula", "chao1", "Cluster") 

model_list2 <- list()
AIC_list2 <- list()

### model selection#####
# DIY drop 1
# loop without the technical variable
for (i in 2:length(mod_vars2)) {
  curr_vars <- paste(mod_vars2[-i], collapse = "+")
  tmp_model <- adonis2(as.formula(paste("dist_b9", "~", curr_vars)),
                       meta_b9, 
                       permutations = 99
  )
  tmp_AIC <- AICc.PERMANOVA2(tmp_model, meta_b9)$AICc
  
  # save models to list
  model_list2[[i]] <- tmp_model
  AIC_list2[[i]] <- tmp_AIC
  new_AIC <- Reduce(min, AIC_list2[2:length(AIC_list2)])
  # if new model has AIC of 2 of more smaller, update model and variable list (and use new AIC in next round)
}


AIC_list2
indx <- which(AIC_list2[2:length(AIC_list2)] == min(unlist(AIC_list2[2:length(AIC_list2)]))) 
indx

model_list2[2:length(model_list2)][[indx]]

# stricter AIC

if (curr_AIC - new_AIC >= 1.5) {
  curr_model <- model_list2[2:length(model_list2)][[indx]]
  curr_AIC <- new_AIC
}

# next round
mod_vars3 <- c("MGI_Plate", 
               "inf_DeliveryMode",
               "m_PreviousDeliveries", 
               "Formula_1stmonth", "ExclusiveBF_3m",
               "inf_SkinContactAge", "Ext_BristolScore",
               "inf_AgeHome", "inf_HospitalLocal", 
               "m_AgeDelivery", "env_FurOrFeathers", 
               "SUM", "inf_TravelAbroad3m", 
               "ATB_prev6m_corr", "Penicillins", 
               "inf_Sex", #"familly_NbCaregivers3m_Cat", 
               "inf_HasHospitalization_past3m", 
               "familly_TimeSpendNonParents3m", 
               "m_BMI", "inf_IsAway3m", 
               "DayCare", "breastfeeding", 
               "env_NbHoursDayCare3m",# "Aminoglycosides", 
               "Formula", "chao1", "Cluster") 

model_list3 <- list()
AIC_list3 <- list()

### model selection#####
# DIY drop 1
# loop without the technical variable
for (i in 2:length(mod_vars3)) {
  curr_vars <- paste(mod_vars3[-i], collapse = "+")
  tmp_model <- adonis2(as.formula(paste("dist_b9", "~", curr_vars)),
                       meta_b9, 
                       permutations = 99
  )
  tmp_AIC <- AICc.PERMANOVA2(tmp_model, meta_b9)$AICc
  
  # save models to list
  model_list3[[i]] <- tmp_model
  AIC_list3[[i]] <- tmp_AIC
  new_AIC <- Reduce(min, AIC_list3[2:length(AIC_list3)])
  # if new model has AIC of 2 of more smaller, update model and variable list (and use new AIC in next round)
}


AIC_list3
indx <- which(AIC_list3[2:length(AIC_list3)] == min(unlist(AIC_list3[2:length(AIC_list3)]))) 
indx

model_list3[2:length(model_list3)][[indx]]

# stricter AIC

if (curr_AIC - new_AIC >= 1.5) {
  curr_model <- model_list3[2:length(model_list3)][[indx]]
  curr_AIC <- new_AIC
}

# next round
mod_vars4 <- c("MGI_Plate", 
               "inf_DeliveryMode",
               "m_PreviousDeliveries", 
               "Formula_1stmonth", "ExclusiveBF_3m",
               "inf_SkinContactAge", "Ext_BristolScore",
               "inf_AgeHome", "inf_HospitalLocal", 
               "m_AgeDelivery", "env_FurOrFeathers", 
               "SUM", "inf_TravelAbroad3m", 
               "ATB_prev6m_corr", #"Penicillins", 
               "inf_Sex", #"familly_NbCaregivers3m_Cat", 
               "inf_HasHospitalization_past3m", 
               "familly_TimeSpendNonParents3m", 
               "m_BMI", "inf_IsAway3m", 
               "DayCare", "breastfeeding", 
               "env_NbHoursDayCare3m",# "Aminoglycosides", 
               "Formula", "chao1", "Cluster") 

model_list4 <- list()
AIC_list4 <- list()

### model selection#####
# DIY drop 1
# loop without the technical variable
for (i in 2:length(mod_vars4)) {
  curr_vars <- paste(mod_vars4[-i], collapse = "+")
  tmp_model <- adonis2(as.formula(paste("dist_b9", "~", curr_vars)),
                       meta_b9, 
                       permutations = 99
  )
  tmp_AIC <- AICc.PERMANOVA2(tmp_model, meta_b9)$AICc
  
  # save models to list
  model_list4[[i]] <- tmp_model
  AIC_list4[[i]] <- tmp_AIC
  new_AIC <- Reduce(min, AIC_list4[2:length(AIC_list4)])
  # if new model has AIC of 2 of more smaller, update model and variable list (and use new AIC in next round)
}


AIC_list4
indx <- which(AIC_list4[2:length(AIC_list4)] == min(unlist(AIC_list4[2:length(AIC_list4)]))) 
indx

model_list4[2:length(model_list4)][[indx]]

# stricter AIC

if (curr_AIC - new_AIC >= 1.5) {
  curr_model <- model_list4[2:length(model_list4)][[indx]]
  curr_AIC <- new_AIC
}

# next round
mod_vars5 <- c("MGI_Plate", 
               "inf_DeliveryMode",
               "m_PreviousDeliveries",
               "Formula_1stmonth", "ExclusiveBF_3m",
               # "inf_SkinContactAge",
               "Ext_BristolScore",
               "inf_AgeHome", "inf_HospitalLocal", 
               "m_AgeDelivery", "env_FurOrFeathers", 
               "SUM", "inf_TravelAbroad3m", 
               "ATB_prev6m_corr", #"Penicillins", 
               "inf_Sex", #"familly_NbCaregivers3m_Cat", 
               "inf_HasHospitalization_past3m", 
               "familly_TimeSpendNonParents3m", 
               "m_BMI", "inf_IsAway3m", 
               "DayCare", "breastfeeding", 
               "env_NbHoursDayCare3m",# "Aminoglycosides", 
               "Formula", "chao1", "Cluster") 

model_list5 <- list()
AIC_list5 <- list()

### model selection#####
# DIY drop 1
# loop without the technical variable
for (i in 2:length(mod_vars5)) {
  curr_vars <- paste(mod_vars5[-i], collapse = "+")
  tmp_model <- adonis2(as.formula(paste("dist_b9", "~", curr_vars)),
                       meta_b9, 
                       permutations = 99
  )
  tmp_AIC <- AICc.PERMANOVA2(tmp_model, meta_b9)$AICc
  
  # save models to list
  model_list5[[i]] <- tmp_model
  AIC_list5[[i]] <- tmp_AIC
  new_AIC <- Reduce(min, AIC_list5[2:length(AIC_list5)])
  # if new model has AIC of 2 of more smaller, update model and variable list (and use new AIC in next round)
}


AIC_list5
indx <- which(AIC_list5[2:length(AIC_list5)] == min(unlist(AIC_list5[2:length(AIC_list5)]))) 
indx

model_list5[2:length(model_list5)][[indx]]

# stricter AIC

if (curr_AIC - new_AIC >= 1.5) {
  curr_model <- model_list5[2:length(model_list5)][[indx]]
  curr_AIC <- new_AIC
}

# next round
mod_vars6 <- c("MGI_Plate", 
               "inf_DeliveryMode",
               "m_PreviousDeliveries",
               "Formula_1stmonth", "ExclusiveBF_3m",
               # "inf_SkinContactAge",
               "Ext_BristolScore",
               "inf_AgeHome", "inf_HospitalLocal", 
               "m_AgeDelivery", "env_FurOrFeathers", 
               "SUM", "inf_TravelAbroad3m", 
               "ATB_prev6m_corr", #"Penicillins", 
               # "inf_Sex", #"familly_NbCaregivers3m_Cat", 
               "inf_HasHospitalization_past3m", 
               "familly_TimeSpendNonParents3m", 
               "m_BMI", "inf_IsAway3m", 
               "DayCare", "breastfeeding", 
               "env_NbHoursDayCare3m",# "Aminoglycosides", 
               "Formula", "chao1", "Cluster") 

model_list6 <- list()
AIC_list6 <- list()

### model selection#####
# DIY drop 1
# loop without the technical variable
for (i in 2:length(mod_vars6)) {
  curr_vars <- paste(mod_vars6[-i], collapse = "+")
  tmp_model <- adonis2(as.formula(paste("dist_b9", "~", curr_vars)),
                       meta_b9, 
                       permutations = 99
  )
  tmp_AIC <- AICc.PERMANOVA2(tmp_model, meta_b9)$AICc
  
  # save models to list
  model_list6[[i]] <- tmp_model
  AIC_list6[[i]] <- tmp_AIC
  new_AIC <- Reduce(min, AIC_list6[2:length(AIC_list6)])
  # if new model has AIC of 2 of more smaller, update model and variable list (and use new AIC in next round)
}


AIC_list6
indx <- which(AIC_list6[2:length(AIC_list6)] == min(unlist(AIC_list6[2:length(AIC_list6)]))) 
indx

model_list6[2:length(model_list6)][[indx]]

# stricter AIC

if (curr_AIC - new_AIC >= 1.5) {
  curr_model <- model_list6[2:length(model_list6)][[indx]]
  curr_AIC <- new_AIC
}

# next round
mod_vars7 <- c("MGI_Plate", 
               "inf_DeliveryMode",
               "m_PreviousDeliveries", 
               "Formula_1stmonth", "ExclusiveBF_3m",
               # "inf_SkinContactAge",
               "Ext_BristolScore",
               "inf_AgeHome", "inf_HospitalLocal", 
               "m_AgeDelivery", "env_FurOrFeathers", 
               "SUM", "inf_TravelAbroad3m", 
               "ATB_prev6m_corr", #"Penicillins", 
               # "inf_Sex", #"familly_NbCaregivers3m_Cat", 
               "inf_HasHospitalization_past3m", 
               "familly_TimeSpendNonParents3m", 
               # "m_BMI", 
               "inf_IsAway3m", 
               "DayCare", "breastfeeding", 
               "env_NbHoursDayCare3m",# "Aminoglycosides", 
               "Formula", "chao1", "Cluster") 

model_list7 <- list()
AIC_list7 <- list()

### model selection#####
# DIY drop 1
# loop without the technical variable
for (i in 2:length(mod_vars7)) {
  curr_vars <- paste(mod_vars7[-i], collapse = "+")
  tmp_model <- adonis2(as.formula(paste("dist_b9", "~", curr_vars)),
                       meta_b9, 
                       permutations = 99
  )
  tmp_AIC <- AICc.PERMANOVA2(tmp_model, meta_b9)$AICc
  
  # save models to list
  model_list7[[i]] <- tmp_model
  AIC_list7[[i]] <- tmp_AIC
  new_AIC <- Reduce(min, AIC_list7[2:length(AIC_list7)])
  # if new model has AIC of 2 of more smaller, update model and variable list (and use new AIC in next round)
}


AIC_list7
indx <- which(AIC_list7[2:length(AIC_list7)] == min(unlist(AIC_list7[2:length(AIC_list7)]))) 
indx

model_list7[2:length(model_list7)][[indx]]

# stricter AIC

if (curr_AIC - new_AIC >= 1.5) {
  curr_model <- model_list7[2:length(model_list7)][[indx]]
  curr_AIC <- new_AIC
}

# next round
mod_vars8 <- c("MGI_Plate", 
               # "Maternal_ab", "Paternal_ab", obsolete
               # "IAP_substance", obsolete
               "inf_DeliveryMode",
               # "birth_season", "inf_Gestational_term", obsolete
               "m_PreviousDeliveries", #"m_Education", obsolete
               "Formula_1stmonth", "ExclusiveBF_3m",
               # "genotype", obsolete
               # "inf_SkinContactAge",
               "Ext_BristolScore",
               "inf_AgeHome", "inf_HospitalLocal", 
               "m_AgeDelivery", "env_FurOrFeathers", 
               "SUM", "inf_TravelAbroad3m", 
               "ATB_prev6m_corr", #"Penicillins", 
               # "inf_Sex", #"familly_NbCaregivers3m_Cat", 
               "inf_HasHospitalization_past3m", 
               "familly_TimeSpendNonParents3m", 
               # "m_BMI", 
               # "inf_IsAway3m", 
               "DayCare", "breastfeeding", 
               "env_NbHoursDayCare3m",# "Aminoglycosides", 
               "Formula", "chao1", "Cluster") 

model_list8 <- list()
AIC_list8 <- list()

### model selection#####
# DIY drop 1
# loop without the technical variable
for (i in 2:length(mod_vars8)) {
  curr_vars <- paste(mod_vars8[-i], collapse = "+")
  tmp_model <- adonis2(as.formula(paste("dist_b9", "~", curr_vars)),
                       meta_b9, 
                       permutations = 99
  )
  tmp_AIC <- AICc.PERMANOVA2(tmp_model, meta_b9)$AICc
  
  # save models to list
  model_list8[[i]] <- tmp_model
  AIC_list8[[i]] <- tmp_AIC
  new_AIC <- Reduce(min, AIC_list8[2:length(AIC_list8)])
  # if new model has AIC of 2 of more smaller, update model and variable list (and use new AIC in next round)
}


AIC_list8
indx <- which(AIC_list8[2:length(AIC_list8)] == min(unlist(AIC_list8[2:length(AIC_list8)]))) 
indx

model_list8[2:length(model_list8)][[indx]]

# stricter AIC

if (curr_AIC - new_AIC >= 1.5) {
  curr_model <- model_list8[2:length(model_list8)][[indx]]
  curr_AIC <- new_AIC
}


final_ado_b9 <- adonis2(dist_b9 ~ MGI_Plate + 
                          # Maternal_ab + Paternal_ab + IAP_substance +
                          inf_DeliveryMode + #birth_season +
                          # inf_Gestational_term + 
                          m_PreviousDeliveries + 
                          # m_Education +  
                          Formula_1stmonth + 
                          ExclusiveBF_3m + 
                          # inf_SkinContactAge +
                          Ext_BristolScore + 
                          inf_AgeHome + inf_HospitalLocal + 
                          m_AgeDelivery + env_FurOrFeathers + 
                          SUM + inf_TravelAbroad3m + 
                          ATB_prev6m_corr +# Penicillins + 
                          # inf_Sex +# familly_NbCaregivers3m_Cat + 
                          inf_HasHospitalization_past3m +
                          familly_TimeSpendNonParents3m + 
                          # m_BMI + inf_IsAway3m + 
                          DayCare + breastfeeding + 
                          env_NbHoursDayCare3m + 
                          # Aminoglycosides + 
                          Formula + 
                          chao1 + Cluster,
                        meta_b9,
                        permutations = 999) 


# write out tables
ado_bact <- final_ado_b4 %>%
  as.data.frame() %>%
  rownames_to_column("Variable") %>%
  mutate(Age = "B4") %>%
  rbind(final_ado_b5 %>%
          as.data.frame() %>%
          rownames_to_column("Variable") %>%
          mutate(Age = "B5")) %>%
  rbind(final_ado_b7 %>%
          as.data.frame() %>%
          rownames_to_column("Variable") %>%
          mutate(Age = "B7")) %>%
  rbind(final_ado_b9 %>%
          as.data.frame() %>%
          rownames_to_column("Variable") %>%
          mutate(Age = "B9")) 


write_csv(ado_bact, "combined_permanova_BC.csv")
