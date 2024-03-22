## Group test for ARGs?
library(tidyverse)
library(DESeq2)

setwd("/path/to/folder/")


meta <- read_csv("metadata/metadata.csv")
meta <- meta %>%
  column_to_rownames("Sample_ID")


# multiply the ARG counts by 10 000 and round to integer

arg <- read_csv("arg_files/normalised_simple_bla.csv")
arg <- column_to_rownames(arg, var = "simplified bla")
arg <- arg[rowSums(arg) > 0,]
meta$blaTEM_load <- unlist(arg[rownames(arg) == "blaTEM" ,])

min(arg[arg > 0]) # 0.001751133 <- Multiply by 10000 to differ form 0 counts
arg <- round(arg * 10000 + 1) # add pseudo count 

arg_class <- read_csv("arg_files/normalised_class.csv")
arg_class <- column_to_rownames(arg_class, var = "Class")
arg_class <- arg_class[rowSums(arg_class) > 0,]

# top classes
identical(colnames(arg_class), rownames(meta))
meta$Tetracycline_load <- unlist(arg_class[rownames(arg_class) == "Tetracycline" ,])
meta$Betalactam_load <- unlist(arg_class[rownames(arg_class) == "Beta-lactam" ,])
meta$Aminoglycoside_load <- unlist(arg_class[rownames(arg_class) == "Aminoglycoside" ,])
meta$erm_gene_load <- unlist(arg_class[rownames(arg_class) == "Macrolide, Lincosamide, Streptogramin B" ,])
meta$mdf_load <- unlist(arg_class[rownames(arg_class) == "Macrolide, Aminoglycoside, Tetracycline, Quinolone, Amphenicol, Rifamycin" ,])
meta$oqx_gene_load <- unlist(arg_class[rownames(arg_class) == "Amphenicol, Quinolone, Quaternary Ammonium Compounds, Folate pathway antagonist" ,])

min(arg_class[arg_class > 0]) # 0.001751133 <- Multiply by 10000 to differ form 0 counts
arg_class <- round(arg_class * 10000 + 1)



identical(colnames(arg), rownames(meta))
identical(colnames(arg_class), rownames(meta))

b4 <- meta$Sample_type == "B4"
b5 <- meta$Sample_type == "B5"
b7 <- meta$Sample_type == "B7"
b9 <- meta$Sample_type == "B9"


## gene names ######

# Delivery ####

arg_b4 <- arg[,b4]
common <- rownames(arg_b4[rowSums(arg_b4 > 1) > 0.02 * ncol(arg_b4),])

arg_b4 <- rbind(arg_b4[rownames(arg_b4) %in% common,], 
                colSums(arg_b4[!rownames(arg_b4) %in% common,]))

meta$Delivery = factor(meta$Delivery, levels = c("Vaginal FALSE", "Vaginal TRUE",
                                                 "C-section TRUE"))
# datasets #####
# B4 #### 
ds_b4 <- DESeqDataSetFromMatrix(countData = arg_b4,
                                colData = meta[b4,],
                                design = ~ MGI_Plate + Delivery)
m_del4 <- DESeq(ds_b4)
res_del4_vagTF <- results(m_del4, contrast = c("Delivery", "Vaginal TRUE", "Vaginal FALSE"))
res_del4_vagFCS <- results(m_del4, contrast = c("Delivery", "C-section TRUE", "Vaginal FALSE"))

# trues
meta$Delivery = factor(meta$Delivery, levels = c("Vaginal TRUE",
                                                 "C-section TRUE", "Vaginal FALSE"))
# datasets #####
# B4 #### 
ds_b4 <- DESeqDataSetFromMatrix(countData = arg_b4,
                                colData = meta[b4,],
                                design = ~ MGI_Plate + Delivery)
m_del4 <- DESeq(ds_b4)
res_del4_vagTCS <- results(m_del4, contrast = c("Delivery", "C-section TRUE", "Vaginal TRUE"))

# convert to data frames
res_del4_vagTF <- as.data.frame(res_del4_vagTF)
res_del4_vagFCS <- as.data.frame(res_del4_vagFCS)
res_del4_vagTCS <- as.data.frame(res_del4_vagTCS)

res_del4_vagTF$contrast <- "Vaginal FALSE vs Vaginal TRUE"
res_del4_vagFCS$contrast <- "Vaginal FALSE vs C-section TRUE"
res_del4_vagTCS$contrast <- "Vaginal TRUE vs C-section TRUE" 

res_del4_vagTF <- res_del4_vagTF %>%
  rownames_to_column("genename")
res_del4_vagFCS <- res_del4_vagFCS %>%
  rownames_to_column("genename")
res_del4_vagTCS <- res_del4_vagTCS %>%
  rownames_to_column("genename")

res_del4 <- rbind(res_del4_vagTF, res_del4_vagFCS, res_del4_vagTCS)
res_del4$timepoint <- "B4"


#### B5 #####
arg_b5 <- arg[,b5]
common <- rownames(arg_b5[rowSums(arg_b5 > 1) > 0.02 * ncol(arg_b5),])

arg_b5 <- rbind(arg_b5[rownames(arg_b5) %in% common,], 
                colSums(arg_b5[!rownames(arg_b5) %in% common,]))

meta$Delivery = factor(meta$Delivery, levels = c("Vaginal FALSE", "Vaginal TRUE",
                                                 "C-section TRUE"))
ds_b5 <- DESeqDataSetFromMatrix(countData = arg_b5,
                                colData = meta[b5,],
                                design = ~ MGI_Plate + Delivery)
m_del5 <- DESeq(ds_b5)
res_del5_vagTF <- results(m_del5, contrast = c("Delivery", "Vaginal TRUE", "Vaginal FALSE"))
res_del5_vagFCS <- results(m_del5, contrast = c("Delivery", "C-section TRUE", "Vaginal FALSE"))

# trues
meta$Delivery = factor(meta$Delivery, levels = c("Vaginal TRUE",
                                                 "C-section TRUE", "Vaginal FALSE"))
ds_b5 <- DESeqDataSetFromMatrix(countData = arg_b5,
                                colData = meta[b5,],
                                design = ~ MGI_Plate + Delivery)
m_del5 <- DESeq(ds_b5)
res_del5_vagTCS <- results(m_del5, contrast = c("Delivery", "C-section TRUE", "Vaginal TRUE"))
# convert to data frames
res_del5_vagTF <- as.data.frame(res_del5_vagTF)
res_del5_vagFCS <- as.data.frame(res_del5_vagFCS)
res_del5_vagTCS <- as.data.frame(res_del5_vagTCS)

res_del5_vagTF$contrast <- "Vaginal FALSE vs Vaginal TRUE"
res_del5_vagFCS$contrast <- "Vaginal FALSE vs C-section TRUE"
res_del5_vagTCS$contrast <- "Vaginal TRUE vs C-section TRUE" 

res_del5_vagTF <- res_del5_vagTF %>%
  rownames_to_column("genename")
res_del5_vagFCS <- res_del5_vagFCS %>%
  rownames_to_column("genename")
res_del5_vagTCS <- res_del5_vagTCS %>%
  rownames_to_column("genename")

res_del5 <- rbind(res_del5_vagTF, res_del5_vagFCS, res_del5_vagTCS)
res_del5$timepoint <- "B5"


#### B7 #####
arg_b7 <- arg[,b7]
common <- rownames(arg_b7[rowSums(arg_b7 > 1) > 0.02 * ncol(arg_b7),])

arg_b7 <- rbind(arg_b7[rownames(arg_b7) %in% common,], 
                colSums(arg_b7[!rownames(arg_b7) %in% common,]))

meta$Delivery = factor(meta$Delivery, levels = c("Vaginal FALSE", "Vaginal TRUE",
                                                 "C-section TRUE"))
ds_b7 <- DESeqDataSetFromMatrix(countData = arg_b7,
                                colData = meta[b7,],
                                design = ~ MGI_Plate + Delivery)
m_del7 <- DESeq(ds_b7)
res_del7_vagTF <- results(m_del7, contrast = c("Delivery", "Vaginal TRUE", "Vaginal FALSE"))
res_del7_vagFCS <- results(m_del7, contrast = c("Delivery", "C-section TRUE", "Vaginal FALSE"))

# trues

meta$Delivery = factor(meta$Delivery, levels = c("Vaginal TRUE",
                                                 "C-section TRUE", "Vaginal FALSE"))
ds_b7 <- DESeqDataSetFromMatrix(countData = arg_b7,
                                colData = meta[b7,],
                                design = ~ MGI_Plate + Delivery)
m_del7 <- DESeq(ds_b7)
res_del7_vagTCS <- results(m_del7, contrast = c("Delivery", "C-section TRUE", "Vaginal TRUE"))

# convert to data frames
res_del7_vagTF <- as.data.frame(res_del7_vagTF)
res_del7_vagFCS <- as.data.frame(res_del7_vagFCS)
res_del7_vagTCS <- as.data.frame(res_del7_vagTCS)

res_del7_vagTF$contrast <- "Vaginal FALSE vs Vaginal TRUE"
res_del7_vagFCS$contrast <- "Vaginal FALSE vs C-section TRUE"
res_del7_vagTCS$contrast <- "Vaginal TRUE vs C-section TRUE" 

res_del7_vagTF <- res_del7_vagTF %>%
  rownames_to_column("genename")
res_del7_vagFCS <- res_del7_vagFCS %>%
  rownames_to_column("genename")
res_del7_vagTCS <- res_del7_vagTCS %>%
  rownames_to_column("genename")

res_del7 <- rbind(res_del7_vagTF, res_del7_vagFCS, res_del7_vagTCS)
res_del7$timepoint <- "B7"


#### B9 #####
arg_b9 <- arg[,b9]
common <- rownames(arg_b9[rowSums(arg_b9 > 1) > 0.02 * ncol(arg_b9),])

arg_b9 <- rbind(arg_b9[rownames(arg_b9) %in% common,], 
                colSums(arg_b9[!rownames(arg_b9) %in% common,]))

meta$Delivery = factor(meta$Delivery, levels = c("Vaginal FALSE", "Vaginal TRUE",
                                                 "C-section TRUE"))
ds_b9 <- DESeqDataSetFromMatrix(countData = arg_b9,
                                colData = meta[b9,],
                                design = ~ MGI_Plate + Delivery)
m_del9 <- DESeq(ds_b9)
res_del9_vagTF <- results(m_del9, contrast = c("Delivery", "Vaginal TRUE", "Vaginal FALSE"))
res_del9_vagFCS <- results(m_del9, contrast = c("Delivery", "C-section TRUE", "Vaginal FALSE"))

# trues

meta$Delivery = factor(meta$Delivery, levels = c("Vaginal TRUE",
                                                 "C-section TRUE", "Vaginal FALSE"))
ds_b9 <- DESeqDataSetFromMatrix(countData = arg_b9,
                                colData = meta[b9,],
                                design = ~ MGI_Plate + Delivery)
m_del9 <- DESeq(ds_b9)
res_del9_vagTCS <- results(m_del9, contrast = c("Delivery", "C-section TRUE", "Vaginal TRUE"))

# convert to data frames
res_del9_vagTF <- as.data.frame(res_del9_vagTF)
res_del9_vagFCS <- as.data.frame(res_del9_vagFCS)
res_del9_vagTCS <- as.data.frame(res_del9_vagTCS)

res_del9_vagTF$contrast <- "Vaginal FALSE vs Vaginal TRUE"
res_del9_vagFCS$contrast <- "Vaginal FALSE vs C-section TRUE"
res_del9_vagTCS$contrast <- "Vaginal TRUE vs C-section TRUE" 

res_del9_vagTF <- res_del9_vagTF %>%
  rownames_to_column("genename")
res_del9_vagFCS <- res_del9_vagFCS %>%
  rownames_to_column("genename")
res_del9_vagTCS <- res_del9_vagTCS %>%
  rownames_to_column("genename")

res_del9 <- rbind(res_del9_vagTF, res_del9_vagFCS, res_del9_vagTCS)
res_del9$timepoint <- "B9"


res_del <- rbind(res_del4, res_del5, res_del7, res_del9)
write_csv(res_del, "DAs_NEGBIN/Delivery_genenames.csv")


# Delivery mode ####
meta$inf_DeliveryMode = factor(meta$inf_DeliveryMode, levels = c("Vaginal", "C-section"))

# datasets #####
# B4 ####
arg_b4 <- arg[,b4]
common <- rownames(arg_b4[rowSums(arg_b4 > 1) > 0.02 * ncol(arg_b4),])

arg_b4 <- rbind(arg_b4[rownames(arg_b4) %in% common,], 
                colSums(arg_b4[!rownames(arg_b4) %in% common,]))

ds_b4 <- DESeqDataSetFromMatrix(countData = arg_b4,
                                colData = meta[b4,],
                                design = ~ MGI_Plate + inf_DeliveryMode)
m_delM4 <- DESeq(ds_b4)
res_delM4 <- results(m_delM4, contrast = c("inf_DeliveryMode", "C-section", "Vaginal"))

# convert to data frames
res_delM4 <- as.data.frame(res_delM4)

res_delM4$contrast <- "Vaginal vs C-section"

res_delM4$timepoint <- "B4"
res_delM4 <- res_delM4 %>%
  rownames_to_column("genename")

#### B5 #####
arg_b5 <- arg[,b5]
common <- rownames(arg_b5[rowSums(arg_b5 > 1) > 0.02 * ncol(arg_b5),])

arg_b5 <- rbind(arg_b5[rownames(arg_b5) %in% common,], 
                colSums(arg_b5[!rownames(arg_b5) %in% common,]))

ds_b5 <- DESeqDataSetFromMatrix(countData = arg_b5,
                                colData = meta[b5,],
                                design = ~ MGI_Plate + inf_DeliveryMode)
m_delM5 <- DESeq(ds_b5)
res_delM5 <- results(m_delM5, contrast = c("inf_DeliveryMode", "C-section", "Vaginal"))

# convert to data frames
res_delM5 <- as.data.frame(res_delM5)

res_delM5$contrast <- "Vaginal vs C-section"

res_delM5$timepoint <- "B5"
res_delM5 <- res_delM5 %>%
  rownames_to_column("genename")

#### B7 #####
arg_b7 <- arg[,b7]
common <- rownames(arg_b7[rowSums(arg_b7 > 1) > 0.02 * ncol(arg_b7),])

arg_b7 <- rbind(arg_b7[rownames(arg_b7) %in% common,], 
                colSums(arg_b7[!rownames(arg_b7) %in% common,]))

ds_b7 <- DESeqDataSetFromMatrix(countData = arg_b7,
                                colData = meta[b7,],
                                design = ~ MGI_Plate + inf_DeliveryMode)
m_delM7 <- DESeq(ds_b7)
res_delM7 <- results(m_delM7, contrast = c("inf_DeliveryMode", "C-section", "Vaginal"))
# convert to data frames
res_delM7 <- as.data.frame(res_delM7)
res_delM7$contrast <- "Vaginal vs C-section"

res_delM7$timepoint <- "B7"
res_delM7 <- res_delM7 %>%
  rownames_to_column("genename")

#### B9 #####
arg_b9 <- arg[,b9]
common <- rownames(arg_b9[rowSums(arg_b9 > 1) > 0.02 * ncol(arg_b9),])

arg_b9 <- rbind(arg_b9[rownames(arg_b9) %in% common,], 
                colSums(arg_b9[!rownames(arg_b9) %in% common,]))

ds_b9 <- DESeqDataSetFromMatrix(countData = arg_b9,
                                colData = meta[b9,],
                                design = ~ MGI_Plate + inf_DeliveryMode)
m_delM9 <- DESeq(ds_b9)
res_delM9 <- results(m_delM9, contrast = c("inf_DeliveryMode", "C-section", "Vaginal"))

# convert to data frames
res_delM9 <- as.data.frame(res_delM9)

res_delM9$contrast <- "Vaginal vs C-section"

res_delM9$timepoint <- "B9"
res_delM9 <- res_delM9 %>%
  rownames_to_column("genename")

res_delM <- rbind(res_delM4, res_delM5, res_delM7, res_delM9)
write_csv(res_delM, "DAs_NEGBIN/DeliveryMode_genenames.csv")


# IAP antibiotic only vd ######

b4VD <- meta$Sample_type == "B4" & meta$inf_DeliveryMode == "Vaginal"
b5VD <- meta$Sample_type == "B5" & meta$inf_DeliveryMode == "Vaginal"
b7VD <- meta$Sample_type == "B7" & meta$inf_DeliveryMode == "Vaginal"
b9VD <- meta$Sample_type == "B9" & meta$inf_DeliveryMode == "Vaginal"

arg_b4 <- arg[,b4VD]
common <- rownames(arg_b4[rowSums(arg_b4 > 1) > 0.02 * ncol(arg_b4),])

arg_b4 <- rbind(arg_b4[rownames(arg_b4) %in% common,], 
                colSums(arg_b4[!rownames(arg_b4) %in% common,]))

meta$reg_AntibioticSubstance <- ifelse(meta$corr_RecievingIAP == FALSE, "None",
                                       meta$reg_AntibioticSubstance)
meta$reg_AntibioticSubstance = factor(meta$reg_AntibioticSubstance, levels = c("None", "Cephalosporin", 
                                                                               "Penicillin", "Other"))
# datasets #####
# B4 #### 
ds_b4 <- DESeqDataSetFromMatrix(countData = arg_b4[, colnames(arg_b4) %in% rownames(meta[b4VD,])],
                                colData = meta[b4VD,],
                                design = ~ MGI_Plate + reg_AntibioticSubstance)
m_iaab4 <- DESeq(ds_b4)
res_iaab4_NoCep <- results(m_iaab4, contrast = c("reg_AntibioticSubstance", "Cephalosporin", "None"))
res_iaab4_NoPen <- results(m_iaab4, contrast = c("reg_AntibioticSubstance", "Penicillin", "None"))
res_iaab4_NoOth <- results(m_iaab4, contrast = c("reg_AntibioticSubstance", "Other", "None"))

# Cephalosporin as ref
meta$reg_AntibioticSubstance = factor(meta$reg_AntibioticSubstance, levels = c("Cephalosporin", 
                                                                               "Penicillin", "Other",
                                                                               "None"))

ds_b4 <- DESeqDataSetFromMatrix(countData = arg_b4[, colnames(arg_b4) %in% rownames(meta[b4VD,])],
                                colData = meta[b4VD & !is.na(meta$reg_AntibioticSubstance),],
                                design = ~ MGI_Plate + reg_AntibioticSubstance)
m_iaab4 <- DESeq(ds_b4)
res_iaab4_CepPen <- results(m_iaab4, contrast = c("reg_AntibioticSubstance", "Penicillin", "Cephalosporin"))
res_iaab4_CepOth <- results(m_iaab4, contrast = c("reg_AntibioticSubstance", "Other", "Cephalosporin"))


# Penicillin as ref
meta$reg_AntibioticSubstance = factor(meta$reg_AntibioticSubstance, levels = c("Penicillin", "Other",
                                                                               "Cephalosporin", "None"))

ds_b4 <- DESeqDataSetFromMatrix(countData = arg_b4[, colnames(arg_b4) %in% rownames(meta[b4VD,])],
                                colData = meta[b4VD,],
                                design = ~ MGI_Plate + reg_AntibioticSubstance)
m_iaab4 <- DESeq(ds_b4)
res_iaab4_PenOth <- results(m_iaab4, contrast = c("reg_AntibioticSubstance", "Other", "Penicillin"))

# convert to data frames
res_iaab4_NoCep <- as.data.frame(res_iaab4_NoCep)
res_iaab4_NoPen <- as.data.frame(res_iaab4_NoPen)
res_iaab4_NoOth <- as.data.frame(res_iaab4_NoOth)
res_iaab4_CepPen <- as.data.frame(res_iaab4_CepPen)
res_iaab4_CepOth <- as.data.frame(res_iaab4_CepOth)
res_iaab4_PenOth <- as.data.frame(res_iaab4_PenOth)

res_iaab4_NoCep$contrast <- "None vs Cephalosporin"
res_iaab4_NoPen$contrast <- "None vs Penicillin"
res_iaab4_NoOth$contrast <- "None vs Other"
res_iaab4_CepPen$contrast <- "Cephalosporin vs Penicillin"
res_iaab4_CepOth$contrast <- "Cephalosporin FALSE vs Other"
res_iaab4_PenOth$contrast <- "Penicillin vs Other"


res_iaab4_NoCep <- res_iaab4_NoCep %>%
  rownames_to_column("genename")
res_iaab4_NoPen <- res_iaab4_NoPen %>%
  rownames_to_column("genename")
res_iaab4_NoOth <- res_iaab4_NoOth %>%
  rownames_to_column("genename")
res_iaab4_CepPen <- res_iaab4_CepPen %>%
  rownames_to_column("genename")
res_iaab4_CepOth <- res_iaab4_CepOth %>%
  rownames_to_column("genename")
res_iaab4_PenOth <- res_iaab4_PenOth %>%
  rownames_to_column("genename")

res_iaab4 <- rbind(res_iaab4_NoCep, res_iaab4_NoPen, res_iaab4_NoOth,
                   res_iaab4_CepPen, res_iaab4_CepOth, res_iaab4_PenOth)
res_iaab4$timepoint <- "B4"


#### B5 #####
arg_b5 <- arg[,b5VD]
common <- rownames(arg_b5[rowSums(arg_b5 > 1) > 0.02 * ncol(arg_b5),])

arg_b5 <- rbind(arg_b5[rownames(arg_b5) %in% common,], 
                colSums(arg_b5[!rownames(arg_b5) %in% common,]))

ds_b5 <- DESeqDataSetFromMatrix(countData = arg_b5[, colnames(arg_b5) %in% rownames(meta[b5VD,])],
                                colData = meta[b5VD,],
                                design = ~ MGI_Plate + reg_AntibioticSubstance)
m_iaab5 <- DESeq(ds_b5)
res_iaab5_NoCep <- results(m_iaab5, contrast = c("reg_AntibioticSubstance", "Cephalosporin", "None"))
res_iaab5_NoPen <- results(m_iaab5, contrast = c("reg_AntibioticSubstance", "Penicillin", "None"))
res_iaab5_NoOth <- results(m_iaab5, contrast = c("reg_AntibioticSubstance", "Other", "None"))

# Cephalosporin as ref
meta$reg_AntibioticSubstance = factor(meta$reg_AntibioticSubstance, levels = c("Cephalosporin", 
                                                                               "Penicillin", "Other",
                                                                               "None"))

ds_b5 <- DESeqDataSetFromMatrix(countData = arg_b5[, colnames(arg_b5) %in% rownames(meta[b5VD,])],
                                colData = meta[b5VD ,],
                                design = ~ MGI_Plate + reg_AntibioticSubstance)
m_iaab5 <- DESeq(ds_b5)
res_iaab5_CepPen <- results(m_iaab5, contrast = c("reg_AntibioticSubstance", "Penicillin", "Cephalosporin"))
res_iaab5_CepOth <- results(m_iaab5, contrast = c("reg_AntibioticSubstance", "Other", "Cephalosporin"))


# Penicillin as ref
meta$reg_AntibioticSubstance = factor(meta$reg_AntibioticSubstance, levels = c("Penicillin", "Other",
                                                                               "Cephalosporin", "None"))

ds_b5 <- DESeqDataSetFromMatrix(countData = arg_b5[, colnames(arg_b5) %in% rownames(meta[b5VD,])],
                                colData = meta[b5VD,],
                                design = ~ MGI_Plate + reg_AntibioticSubstance)
m_iaab5 <- DESeq(ds_b5)
res_iaab5_PenOth <- results(m_iaab5, contrast = c("reg_AntibioticSubstance", "Other", "Penicillin"))

# convert to data frames
res_iaab5_NoCep <- as.data.frame(res_iaab5_NoCep)
res_iaab5_NoPen <- as.data.frame(res_iaab5_NoPen)
res_iaab5_NoOth <- as.data.frame(res_iaab5_NoOth)
res_iaab5_CepPen <- as.data.frame(res_iaab5_CepPen)
res_iaab5_CepOth <- as.data.frame(res_iaab5_CepOth)
res_iaab5_PenOth <- as.data.frame(res_iaab5_PenOth)

res_iaab5_NoCep$contrast <- "None vs Cephalosporin"
res_iaab5_NoPen$contrast <- "None vs Penicillin"
res_iaab5_NoOth$contrast <- "None vs Other"
res_iaab5_CepPen$contrast <- "Cephalosporin vs Penicillin"
res_iaab5_CepOth$contrast <- "Cephalosporin FALSE vs Other"
res_iaab5_PenOth$contrast <- "Penicillin vs Other"


res_iaab5_NoCep <- res_iaab5_NoCep %>%
  rownames_to_column("genename")
res_iaab5_NoPen <- res_iaab5_NoPen %>%
  rownames_to_column("genename")
res_iaab5_NoOth <- res_iaab5_NoOth %>%
  rownames_to_column("genename")
res_iaab5_CepPen <- res_iaab5_CepPen %>%
  rownames_to_column("genename")
res_iaab5_CepOth <- res_iaab5_CepOth %>%
  rownames_to_column("genename")
res_iaab5_PenOth <- res_iaab5_PenOth %>%
  rownames_to_column("genename")

res_iaab5 <- rbind(res_iaab5_NoCep, res_iaab5_NoPen, res_iaab5_NoOth,
                   res_iaab5_CepPen, res_iaab5_CepOth, res_iaab5_PenOth)
res_iaab5$timepoint <- "B5"


#### B7 #####
arg_b7 <- arg[,b7VD]
common <- rownames(arg_b7[rowSums(arg_b7 > 1) > 0.02 * ncol(arg_b7),])

arg_b7 <- rbind(arg_b7[rownames(arg_b7) %in% common,], 
                colSums(arg_b7[!rownames(arg_b7) %in% common,]))

ds_b7 <- DESeqDataSetFromMatrix(countData = arg_b7[, colnames(arg_b7) %in% rownames(meta[b7VD,])],
                                colData = meta[b7VD,],
                                design = ~ MGI_Plate + reg_AntibioticSubstance)
m_iaab7 <- DESeq(ds_b7)
res_iaab7_NoCep <- results(m_iaab7, contrast = c("reg_AntibioticSubstance", "Cephalosporin", "None"))
res_iaab7_NoPen <- results(m_iaab7, contrast = c("reg_AntibioticSubstance", "Penicillin", "None"))
res_iaab7_NoOth <- results(m_iaab7, contrast = c("reg_AntibioticSubstance", "Other", "None"))

# Cephalosporin as ref
meta$reg_AntibioticSubstance = factor(meta$reg_AntibioticSubstance, levels = c("Cephalosporin", 
                                                                               "Penicillin", "Other",
                                                                               "None"))

ds_b7 <- DESeqDataSetFromMatrix(countData = arg_b7[, colnames(arg_b7) %in% rownames(meta[b7VD,])],
                                colData = meta[b7VD,],
                                design = ~ MGI_Plate + reg_AntibioticSubstance)
m_iaab7 <- DESeq(ds_b7)
res_iaab7_CepPen <- results(m_iaab7, contrast = c("reg_AntibioticSubstance", "Penicillin", "Cephalosporin"))
res_iaab7_CepOth <- results(m_iaab7, contrast = c("reg_AntibioticSubstance", "Other", "Cephalosporin"))


# Penicillin as ref
meta$reg_AntibioticSubstance = factor(meta$reg_AntibioticSubstance, levels = c("Penicillin", "Other",
                                                                               "Cephalosporin", "None"))

ds_b7 <- DESeqDataSetFromMatrix(countData = arg_b7[, colnames(arg_b7) %in% rownames(meta[b7VD,])],
                                colData = meta[b7VD,],
                                design = ~ MGI_Plate + reg_AntibioticSubstance)
m_iaab7 <- DESeq(ds_b7)
res_iaab7_PenOth <- results(m_iaab7, contrast = c("reg_AntibioticSubstance", "Other", "Penicillin"))

# convert to data frames
res_iaab7_NoCep <- as.data.frame(res_iaab7_NoCep)
res_iaab7_NoPen <- as.data.frame(res_iaab7_NoPen)
res_iaab7_NoOth <- as.data.frame(res_iaab7_NoOth)
res_iaab7_CepPen <- as.data.frame(res_iaab7_CepPen)
res_iaab7_CepOth <- as.data.frame(res_iaab7_CepOth)
res_iaab7_PenOth <- as.data.frame(res_iaab7_PenOth)

res_iaab7_NoCep$contrast <- "None vs Cephalosporin"
res_iaab7_NoPen$contrast <- "None vs Penicillin"
res_iaab7_NoOth$contrast <- "None vs Other"
res_iaab7_CepPen$contrast <- "Cephalosporin vs Penicillin"
res_iaab7_CepOth$contrast <- "Cephalosporin FALSE vs Other"
res_iaab7_PenOth$contrast <- "Penicillin vs Other"


res_iaab7_NoCep <- res_iaab7_NoCep %>%
  rownames_to_column("genename")
res_iaab7_NoPen <- res_iaab7_NoPen %>%
  rownames_to_column("genename")
res_iaab7_NoOth <- res_iaab7_NoOth %>%
  rownames_to_column("genename")
res_iaab7_CepPen <- res_iaab7_CepPen %>%
  rownames_to_column("genename")
res_iaab7_CepOth <- res_iaab7_CepOth %>%
  rownames_to_column("genename")
res_iaab7_PenOth <- res_iaab7_PenOth %>%
  rownames_to_column("genename")

res_iaab7 <- rbind(res_iaab7_NoCep, res_iaab7_NoPen, res_iaab7_NoOth,
                   res_iaab7_CepPen, res_iaab7_CepOth, res_iaab7_PenOth)
res_iaab7$timepoint <- "B7"

#### B9 #####
arg_b9 <- arg[,b9VD]
common <- rownames(arg_b9[rowSums(arg_b9 > 1) > 0.02 * ncol(arg_b9),])

arg_b9 <- rbind(arg_b9[rownames(arg_b9) %in% common,], 
                colSums(arg_b9[!rownames(arg_b9) %in% common,]))

ds_b9 <- DESeqDataSetFromMatrix(countData = arg_b9[, colnames(arg_b9) %in% rownames(meta[b9VD,])],
                                colData = meta[b9VD,],
                                design = ~ MGI_Plate + reg_AntibioticSubstance)
m_iaab9 <- DESeq(ds_b9)
res_iaab9_NoCep <- results(m_iaab9, contrast = c("reg_AntibioticSubstance", "Cephalosporin", "None"))
res_iaab9_NoPen <- results(m_iaab9, contrast = c("reg_AntibioticSubstance", "Penicillin", "None"))
res_iaab9_NoOth <- results(m_iaab9, contrast = c("reg_AntibioticSubstance", "Other", "None"))

# Cephalosporin as ref
meta$reg_AntibioticSubstance = factor(meta$reg_AntibioticSubstance, levels = c("Cephalosporin", 
                                                                               "Penicillin", "Other",
                                                                               "None"))

ds_b9 <- DESeqDataSetFromMatrix(countData = arg_b9[, colnames(arg_b9) %in% rownames(meta[b9VD,])],
                                colData = meta[b9VD,],
                                design = ~ MGI_Plate + reg_AntibioticSubstance)
m_iaab9 <- DESeq(ds_b9)
res_iaab9_CepPen <- results(m_iaab9, contrast = c("reg_AntibioticSubstance", "Penicillin", "Cephalosporin"))
res_iaab9_CepOth <- results(m_iaab9, contrast = c("reg_AntibioticSubstance", "Other", "Cephalosporin"))


# Penicillin as ref
meta$reg_AntibioticSubstance = factor(meta$reg_AntibioticSubstance, levels = c("Penicillin", "Other",
                                                                               "Cephalosporin", "None"))

ds_b9 <- DESeqDataSetFromMatrix(countData = arg_b9[, colnames(arg_b9) %in% rownames(meta[b9VD,])],
                                colData = meta[b9VD,],
                                design = ~ MGI_Plate + reg_AntibioticSubstance)
m_iaab9 <- DESeq(ds_b9)
res_iaab9_PenOth <- results(m_iaab9, contrast = c("reg_AntibioticSubstance", "Other", "Penicillin"))

# convert to data frames
res_iaab9_NoCep <- as.data.frame(res_iaab9_NoCep)
res_iaab9_NoPen <- as.data.frame(res_iaab9_NoPen)
res_iaab9_NoOth <- as.data.frame(res_iaab9_NoOth)
res_iaab9_CepPen <- as.data.frame(res_iaab9_CepPen)
res_iaab9_CepOth <- as.data.frame(res_iaab9_CepOth)
res_iaab9_PenOth <- as.data.frame(res_iaab9_PenOth)

res_iaab9_NoCep$contrast <- "None vs Cephalosporin"
res_iaab9_NoPen$contrast <- "None vs Penicillin"
res_iaab9_NoOth$contrast <- "None vs Other"
res_iaab9_CepPen$contrast <- "Cephalosporin vs Penicillin"
res_iaab9_CepOth$contrast <- "Cephalosporin FALSE vs Other"
res_iaab9_PenOth$contrast <- "Penicillin vs Other"


res_iaab9_NoCep <- res_iaab9_NoCep %>%
  rownames_to_column("genename")
res_iaab9_NoPen <- res_iaab9_NoPen %>%
  rownames_to_column("genename")
res_iaab9_NoOth <- res_iaab9_NoOth %>%
  rownames_to_column("genename")
res_iaab9_CepPen <- res_iaab9_CepPen %>%
  rownames_to_column("genename")
res_iaab9_CepOth <- res_iaab9_CepOth %>%
  rownames_to_column("genename")
res_iaab9_PenOth <- res_iaab9_PenOth %>%
  rownames_to_column("genename")

res_iaab9 <- rbind(res_iaab9_NoCep, res_iaab9_NoPen, res_iaab9_NoOth,
                   res_iaab9_CepPen, res_iaab9_CepOth, res_iaab9_PenOth)
res_iaab9$timepoint <- "B9"

res_iaab <- rbind(res_iaab4, res_iaab5, res_iaab7, res_iaab9)
write_csv(res_iaab, "DAs_NEGBIN/IAP_ab_VD_genenames.csv")

## Age #####
# genenames ####
inf <- meta$Sample_type %in% c("B4", "B5", "B7", "B9")

arg_inf <- arg[,inf]

common <- rownames(arg_inf[rowSums(arg > 1) > 0.02 * ncol(arg_inf),])

arg_nonrare <- rbind(arg_inf[rownames(arg_inf) %in% common,], 
                colSums(arg_inf[!rownames(arg_inf) %in% common,]))

meta$Age <- ifelse(meta$Sample_type == "M", meta$m_AgeDelivery,
                   ifelse(meta$Sample_type == "F", meta$p_AgeDelivery,
                          ifelse(meta$Sample_type == "B4", 3,
                                 ifelse(meta$Sample_type == "B5", 6,
                                        ifelse(meta$Sample_type == "B7", 12, 24)))))


ds_inf <- DESeqDataSetFromMatrix(countData = arg_nonrare,
                                colData = meta[inf,],
                                design = ~ MGI_Plate + Family_ID + Age)

ds_inf <- estimateSizeFactors(ds_inf, type = "poscounts")

m_age <- DESeq(ds_inf, fitType = "parametric")
res_age <- results(m_age)

# convert to data frames
res_age <- as.data.frame(res_age) %>%
  rownames_to_column("genename")

write_csv(res_age, "DAs_NEGBIN/Age_genenames.csv")

# Class ####
arg_class_inf <- arg_class[,inf]


ds_inf <- DESeqDataSetFromMatrix(countData = arg_class_inf,
                                 colData = meta[inf,],
                                 design = ~ MGI_Plate + Family_ID + Age)

ds_inf <- estimateSizeFactors(ds_inf, type = "poscounts")

m_age <- DESeq(ds_inf, fitType = "parametric")
res_age <- results(m_age)

# convert to data frames
res_age <- as.data.frame(res_age) %>%
  rownames_to_column("class")

write_csv(res_age, "DAs_NEGBIN/Age_class.csv")


