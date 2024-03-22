library(tidyverse)
library(mboost)
library(caret)
library(SummarizedExperiment)

setwd("/path/to/folder")

meta <- read_csv("metadata/metadata.csv")
colnames(meta)

genus <- read_csv("taxonomical_tables/genus_table_relative.csv")
genus <- genus %>%
  column_to_rownames("genus") 
rownames(genus) <- gsub(" ", "_", rownames(genus))

  
family <- read_csv("taxonomical_tables/family_table_relative.csv")
family <- family %>%
  column_to_rownames("family") 

meta_num <- meta %>% 
  as.data.frame() %>%
  dplyr::select(where(is.numeric)) 


## Example #####
meta_num$ARG_Obs_richness <- NULL

meta_num_b4 <- meta_num[meta$Sample_type == "B4",]
meta_num_b4 <- meta_num_b4 %>%
  select_if(~ sum(!is.na(.)) > 0.9 * 475)

# return sample IDs to combine with family data
meta_num_b4 <- cbind(meta[meta$Sample_type == "B4",]$Sample_ID, meta_num_b4) %>%
  dplyr::rename(c("Sample_ID" = 'meta[meta$Sample_type == "B4", ]$Sample_ID')) 

skimr::skim(meta_num_b4) # remove varibles with uneven distribution

family_b4 <- family[, grepl("B4", colnames(family))]
family_b4 <- family_b4[rowSums(family_b4 > 0.00005) > 0.3 * ncol(family_b4),] %>%
  t() %>%
  as.data.frame() %>%
  select(-Unknown) %>%
  rownames_to_column("Sample_ID")

meta_num_b4 <- meta_num_b4 %>%
  left_join(family_b4, by = "Sample_ID") %>%
  select(-c(Sample_ID))

# 70% train and 30% test sets

train_ind <- createDataPartition(meta_num_b4$ARG_load, times = 1, p = 0.7, list = FALSE) 
b4_train <- meta_num_b4[train_ind,]
b4_test <- meta_num_b4[-train_ind,]


fitControl <- trainControl(method = "cv")

# Set the seed for reproducibility
set.seed(123)

# Train GLMBoost model for ARG load
glmBoostModel <-
  caret::train(
    ARG_load ~ .,
    data = b4_train,
    method = "glmboost",
    tuneLength = 5,
    center = TRUE,
    na.action = na.omit
  )
#------------------------------------------------------
# Regular GLM with GLMboost variables

# Merge tables for coefficients
cf <-
  summary(glmBoostModel)$object %>% coef(off2int = TRUE) %>%
  as.list %>%
  as.matrix %>% 
  as.data.frame %>% 
  dplyr::select(-(matches("Intercept")))

selprob <- summary(glmBoostModel)$selprob %>%
  as.data.frame
coefs_table <- merge(cf, selprob, by = 0)

colnames(coefs_table) <-
  c("Variable", "Coefficient", "Selection probablity")

covariates <-
  c(coefs_table$"Variable"[1:nrow(coefs_table)], "ARG_load")

b4_train_glm <- b4_train %>%
  dplyr::select(all_of(covariates)) %>% as.data.frame()
fit_glm_b4 <- glm(ARG_load ~ ., 
                  data = b4_train_glm,
                  family = "gaussian"(link = "log"))
summary_glm <- summary(fit_glm_b4)
summary_glm

#———————————————————
#Cross validation
# Predicted vs observed plot
predictions <- predict(glmBoostModel, newdata = b4_test)
test.target <- subset(b4_test, select = ARG_load)
test.target <- test.target %>%
  t %>% 
  subset(select = rownames(predictions)) %>% 
  t

b4_test <- na.omit(b4_test)
plot(predictions, y = b4_test$ARG_load,
     xlab='Predicted Values',
     ylab='Actual Values',
     main='Predicted vs. Actual Values')
abline(a=0, b=1)

# RMSE
sqrt(mean((test.target - predictions) ^ 2))
# R2
cor(test.target, predictions) ^ 2

glm_b4 <- as.data.frame(summary_glm$coefficients) %>%
  rownames_to_column("variable")
glm_b4$age <- "B4"
glm_b4$RMSE <- sqrt(mean((test.target - predictions) ^ 2))
glm_b4$R2 <- cor(test.target, predictions)[1] ^ 2
