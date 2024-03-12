# Script for Computational Analysis of Lung_data and Kidney Cancer gene expressions and their CNV relationshi project.
# Overview:
# This script analyzes gene expression data for Lung_data and kidney cancer, performs hypothesis testing,
# creates volcano plots, and conducts regression analysis on selected genes.

# Function to remove gene IDs
source("G:/Amany.Project/Func1.R")

# Function to calculate Shapiro p-values
source("G:/Amany.Project/Func2.R")

# Setting the path for the project data
project_data_path = "G:/NileUniversity/Statistical analysis and Visualization/Project/Project_Data"

# Constructing file paths for Lung_data and kidney cancer data
Lung_data_path = file.path(project_data_path, "/lusc-rsem-fpkm-tcga_paired.txt", sep = "")
Lung_tumor_data_path = file.path(project_data_path, "/lusc-rsem-fpkm-tcga-t_paired.txt", sep = "")
Kidney_data_path = file.path(project_data_path, "/Kidney_data-rsem-fpkm-tcga_paired.txt", sep = "")
Kidney_tumor_data_path = file.path(project_data_path, "/Kidney_data-rsem-fpkm-tcga-t_paired.txt", sep = "")

# Loading gene expression data into R data frames for Lung_data and kidney
Lung_data = read.table(Lung_data_path, header = TRUE, sep = "\t", row.names = 1)
Lung_tumor_data = read.table(Lung_tumor_data_path, header = TRUE, sep = "\t", row.names = 1)
Kidney_data = read.table(Kidney_data_path, header = TRUE, sep = "\t", row.names = 1)
Kidney_tumor_data = read.table(Kidney_tumor_data_path, header = TRUE, sep = "\t", row.names = 1)

# Applying the function to remove gene IDs from the data frames
Lung_data = delet.gene.Ids(Lung_data)
Lung_tumor_data = delet.gene.Ids(Lung_tumor_data)
Kidney_data = delet.gene.Ids(Kidney_data)
Kidney_tumor_data = delet.gene.Ids(Kidney_tumor_data)

# Filtration to delete genes with more than or equal to 50% zeroes
rows_to_delete = c()

# Going through kidney cancer data to identify and remove rows with excessive zeroes
rows_to_delete = which(rowSums(kirc_tumor_data == 0) >= ncol(kirc_tumor_data)/2)
kirc_tumor_data = kirc_tumor_data[-rows_to_delete,]
kirc_data = kirc_data[-rows_to_delete,]

# Similar filtration process for Lung_data data
rows_to_delete = which(rowSums(Lung_tumor_data == 0) >= ncol(Lung_tumor_data)/2)
Lung_tumor_data = Lung_tumor_data[-rows_to_delete,]
Lung_data = Lung_data[-rows_to_delete,]

# Removing identified rows with more than or equal to 50% zeroes from kidney data
Kidney_tumor_data = Kidney_tumor_data[-rows_to_delete,]
Kidney_data = Kidney_data[-rows_to_delete,]

Lung_tumor_data = Lung_tumor_data[-rows_to_delete2,]
Lung_data = Lung_data[-rows_to_delete2,]

############ Kidney, Hypothesis Testing in the first case: ################
############ Samples are paired ###############
### Computing differences between Kidney_tumor_data and Kidney_data, checking normality using shapiro test for the difference ###

diff_Kidney_H = sapply(Kidney_tumor_data - Kidney_data , as.numeric)

shapiro.pvalues.diff = shapiro.pvalues(diff_Kidney_H)

# Checking normality and selecting appropriate statistical test
if(length(which(shapiro.pvalues.diff < 0.05)) >= 1) {
  print("welcoxon test")
} else {
  print("T test")
}

# Performing Wilcoxon test and adjusting p-values
Wilcox_P = c()
Wilcox_T = c()
n1 = nrow(Kidney_tumor_data)
for (i in 1:n1){
  Wilcox_Kidney = wilcox.test(x = as.numeric(Kidney_tumor_data[i,]) , y = as.numeric(Kidney_data[i,]), alternative = 'two.sided', paired = TRUE)
  Wilcox_P[i] = Wilcox_Kidney$p.value
  Wilcox_T[i] = Wilcox_Kidney$statistic
}

# Adjusting p-values, calculating Log2FoldChange, and creating a data frame
padj_Kidney_data = p.adjust(Wilcox_P, method = 'fdr')
Log2FoldChange_Kidney_data = log2(rowMeans(Kidney_tumor_data)) - log2(rowMeans(Kidney_data))
Kidney.df <- data.frame(Wilcox_P, padj_Kidney_data, Wilcox_T, Log2FoldChange_Kidney_data, row.names = row.names(Kidney_tumor_data))

### Degs Using Hypothesis
significant_DEGs_Kidney = which(Kidney.df$padj_Kidney_data < 0.05)

### Degs using logFold
significant_DEGs_Kidney_Fold = which(abs(Kidney.df$Log2FoldChange_Kidney_data) > log2(1.5))

#### Volcano plot

par(mfrow=c(1,1))

with(Kidney.df, plot(Log2FoldChange_Kidney_data, -log10(padj_Kidney_data), pch=20, 
                 main = "DEGs for kidney cancer with healthy", col="gray",
                   xlim = c(-3,3), ylim = c(0,15)))

with(subset(Kidney.df, padj_Kidney_data < 0.05 & Log2FoldChange_Kidney_data > log2(1.5)), 
     points(Log2FoldChange_Kidney_data, -log10(padj_Kidney_data), pch=20, 
                        col="red", xlim = c(-3,3), ylim = c(0,15)))

with(subset(Kidney.df, padj_Kidney_data < 0.05 & Log2FoldChange_Kidney_data < -log2(1.5)), 
     points(Log2FoldChange_Kidney_data, -log10(padj_Kidney_data), pch=20, 
                         col="blue", xlim = c(-3,3), ylim = c(0,15)))

legend("topright", legend = c('Increased', 'Decreased'), 
       col = c('red', 'blue'), pch = c(20, 20))
lines(c(0.6,0.6), c(0,15), lwd = 2, lty = 2)
lines(c(-0.6,-0.6), c(0,15), lwd = 2, lty = 2)
lines(c(-3.5,3.5), c(1.3,1.3), lwd = 2, lty = 2)

########################### GSEA Kidney Paired Files ##########################
###############################################################################

Final_Degs_Kidney_C = Kidney_tumor_data[significant_DEGs_Kidney,]
Final_Degs_Kidney_H = Kidney_data[significant_DEGs_Kidney,]

project_data_path_out_kidney_C = file.path(project_data_path, "Kidney_Degs_Paired_Cancer.csv")
write.table(Final_Degs_Kidney_C, project_data_path_out_kidney_C, sep = ",", row.names = TRUE, col.names = NA, quote = FALSE)

project_data_path_out_kidney_H = file.path(project_data_path, "Kidney_Degs_Paired_Healthy.csv")
write.table(Final_Degs_Kidney_H, project_data_path_out_kidney_H, sep = ",", row.names = TRUE, col.names = NA, quote = FALSE)


########## To Pick up the five most differentially expressed genes ############
###############################################################################
print("the five most differentially expressed genes [Kidney:paired]")
Ordered_Wilcox_T_Kidney = order(-Kidney.df$Wilcox_T)
Top_5_Kidney = Kidney.df[Ordered_Wilcox_T_Kidney[1:5],]
print(row.names(Top_5_Kidney))

############ Lung_data, Hypothesis Testing in the First case: ##################
############ Samples are paired ###############
### Computing differences between Lung_tumor_data and Lung_data, checking normality using shapiro test for the difference ###

diff_Lung_data_H = sapply(Lung_tumor_data - Lung_data , as.numeric)

shapiro.pvalues.diff_2 = shapiro.pvalues(diff_Lung_data_H)

# Checking normality and selecting appropriate statistical test
if(length(which(shapiro.pvalues.diff_2 < 0.05)) > 1) {
  print("welcoxon test")
} else {
  print("T test")
}

# Performing Wilcoxon test and adjusting p-values
Wilcox_P_2 = c()
Wilcox_T_2 = c()
n2 = nrow(Lung_tumor_data)
for (i in 1:n2){
  Wilcox_Lung_data = wilcox.test(x = as.numeric(Lung_tumor_data[i,]) , y = as.numeric(Lung_data[i,]), alternative = 'two.sided', paired = TRUE)
  Wilcox_P_2[i] = Wilcox_Lung_data$p.value
  Wilcox_T_2[i] = Wilcox_Lung_data$statistic
}

# Adjusting p-values, calculating Log2FoldChange, and creating a data frame
padj_Lung_data = p.adjust(Wilcox_P_2, method = 'fdr')
Log2FoldChange_Lung_data = log2(rowMeans(Lung_tumor_data)) - log2(rowMeans(Lung_data))
Lung_data.df <- data.frame(Wilcox_P_2, padj_Lung_data, Wilcox_T_2, Log2FoldChange_Lung_data, row.names = row.names(Lung_tumor_data))

### Degs Using Hypothesis
DEGs_Paired_Lung_data = which(Lung_data.df$padj_Lung_data < 0.05)

### Degs using logFold
DEGs_Paired_Lung_data_Fold = which(abs(Lung_data.df$Log2FoldChange_Lung_data) > log2(1.5))

#### Volcano plot

par(mfrow=c(1,1))

with(Lung_data.df, plot(Log2FoldChange_Lung_data, -log10(padj_Lung_data), pch=20, 
                     main = "Diff. Exp. for Lung_data cancer with healthy", 
            col="gray", xlim = c(-3,3), ylim = c(0,15)))

with(subset(Lung_data.df, padj_Lung_data < 0.05 & Log2FoldChange_Lung_data > log2(1.5)), 
     points(Log2FoldChange_Lung_data, -log10(padj_Lung_data), pch=20, 
            col="red", xlim = c(-3,3), ylim = c(0,15)))

with(subset(Lung_data.df, padj_Lung_data < 0.05 & Log2FoldChange_Lung_data < -log2(1.5)), 
     points(Log2FoldChange_Lung_data, -log10(padj_Lung_data), pch=20, 
            col="blue", xlim = c(-3,3), ylim = c(0,15)))

legend("topright", legend = c('Increased', 'Decreased'), 
       col = c('red', 'blue'), pch = c(20, 20))
lines(c(0.6,0.6), c(0,15), lwd = 2, lty = 2)
lines(c(-0.6,-0.6), c(0,15), lwd = 2, lty = 2)
lines(c(-3.5,3.5), c(1.3,1.3), lwd = 2, lty = 2)

############################ GSEA Lung_data Paired Files ##############################
###############################################################################

Final_Degs_Lung_data_C = Lung_tumor_data[DEGs_Paired_Lung_data, ]
Final_Degs_Lung_data_H = Lung_data[DEGs_Paired_Lung_data, ]

project_data_path.out_L_C = paste(project_data_path, "/Lung_data_Degs_Paired_Cancer.csv", sep = "")
write.table(Final_Degs_Lung_data_C, project_data_path.out_L_C, sep = ",", row.names = T, col.names = NA, quote = F)

project_data_path.out_L_H = paste(project_data_path, "/Lung_data_Degs_Paired_Healthy.csv", sep = "")
write.table(Final_Degs_Lung_data_H, project_data_path.out_L_H, sep = ",", row.names = T, col.names = NA, quote = F)

########## To Pick up the five most differentially expressed genes ############
###############################################################################
print("the five most differentially expressed genes [Lung_data:paired]")
Ordered_Wilcox_T_Lung_data = order(-Lung_data.df$Wilcox_T_2)
Top_5_Lung_data = Lung_data.df[Ordered_Wilcox_T_Lung_data[1:5], ]
Name = row.names(Top_5_Lung_data)
print(Name)

############ Kidney, Hypothesis Testing in the second case: ###############
############ Samples are Independent ###############
## Checking normality using shapiro test for each group of Kidney_tumor_data and kir (one by one) ##

shapiro.pvalues.Kidney_data_Inde=shapiro.pvalues(Kidney_data)
shapiro.pvalues.Kidney_tumor_data_Inde=shapiro.pvalues(Kidney_tumor_data)

if (length(which(shapiro.pvalues.Kidney_data_Inde<0.05)) >= 1){
  print("welcoxon test")
}else{ 
  if (length(which(shapiro.pvalues.Kidney_tumor_data_Inde<0.05)) >= 1){
    print("welcoxon test")
  }else{ 
    print("T test")
  }
}

####### Applying t-test/wilcoxon-test according to the p-values of shapiro
#############################################################################

Wilcox_Inde_P=c()
Wilcox_Inde_T=c()
for (i in 1:nrow(Kidney_tumor_data)){
  Wilcox_Inde_Kidney=wilcox.test(x = as.numeric(Kidney_tumor_data[i,]) , y = as.numeric(Kidney_data[i,]), alternative = 'two.sided',paired = FALSE)
  Wilcox_Inde_P[i]=Wilcox_Inde_Kidney$p.value
  Wilcox_Inde_T[i]=Wilcox_Inde_Kidney$statistic
}

padj_Inde_Kidney_data=p.adjust(Wilcox_Inde_P, method = 'fdr')

Kidney_Inde.df <- data.frame(Wilcox_Inde_P,padj_Inde_Kidney_data,Wilcox_Inde_T,row.names = row.names(Kidney_tumor_data))

DEGs_Inde_Kidney=which(Kidney_Inde.df$padj_Inde_Kidney_data < 0.05)

############ Lung_data, the Hypothesis Testing in the second case: ##################
############                 Samples are Independent                ############
## So, for all genes we check the normality using shapiro test for each group ## 
## of GE_lusc_cancer and GE_lusc_healthy (one by one)                ###########
################################################################################
## Calling a function calculates the shapiro p-values 

shapiro.pvalues.Lung_data_Inde=shapiro.pvalues(Lung_data)
shapiro.pvalues.Lung_tumor_data_Inde=shapiro.pvalues(Lung_tumor_data)

if (length(which(shapiro.pvalues.Lung_data_Inde<0.05)) > 1){
  print("welcoxon test")
}else{ 
  if (length(which(shapiro.pvalues.Lung_tumor_data_Inde<0.05)) > 1){
    print("welcoxon test")
  }else{ 
    print("T test")
  }
}

####### Applying t-test/wilcoxon-test according to the p-values of shapiro
###############################################################################

Wilcox_Lung_data_Inde_P=c()
Wilcox_Lung_data_Inde_T=c()
for (i in 1:nrow(Lung_tumor_data)){
  Wilcox_Inde_Lung_data=wilcox.test(x = as.numeric(Lung_tumor_data[i,]) , y = as.numeric(Lung_data[i,]), alternative = 'two.sided',paired = FALSE)
  Wilcox_Lung_data_Inde_P[i]=Wilcox_Inde_Lung_data$p.value
  Wilcox_Lung_data_Inde_T[i]=Wilcox_Inde_Lung_data$statistic
}

padj_Inde_Lung_data=p.adjust(Wilcox_Lung_data_Inde_P, method = 'fdr')
Lung_data_Inde.df <- data.frame(Wilcox_Lung_data_Inde_P,padj_Inde_Lung_data,Wilcox_Lung_data_Inde_T,row.names = row.names(Lung_tumor_data))
DEGs_Inde_Lung_data=which(Lung_data_Inde.df$padj_Inde_Lung_data < 0.05)

################Venn_Diagram###################################
#############################################################################
install.packages("ggVennDiagram")
library(ggVennDiagram)
library(ggplot2)

significant_DEGs_Kidney_list=rownames(Kidney.df[padj_Kidney_data < 0.05, ])
DEGs_Inde_Kidney_list=rownames(Kidney_Inde.df[padj_Inde_Kidney_data < 0.05,])

ggVennDiagram(list(significant_DEGs_Kidney_list, DEGs_Inde_Kidney_list), 
              category.names = c("Kidney paired","Kidney independent"), color = "red")

DEGs_Paired_Lung_data_list=rownames(Lung_data.df[padj_Lung_data < 0.05, ])
DEGs_Inde_Lung_data_list=rownames(Lung_data_Inde.df[padj_Inde_Lung_data < 0.05,])

ggVennDiagram(list(DEGs_Paired_Lung_data_list, DEGs_Inde_Lung_data_list), 
              category.names = c("Lung_data paired","Lung_data independent"), color = "red")


###############################################################################
################################ Regression : Kidney ##########################
###############################################################################
Kidney_data_cnv_path=paste(project_data_path,"/Kidney_data_CNV_core.txt",sep = "")
lusc_cnv_path=paste(project_data_path,"/lusc_CNV_core.txt",sep = "")
Kidney_data_cnv=read.table(Kidney_data_cnv_path,header = TRUE, sep = "\t")
lusc_cnv=read.table(lusc_cnv_path,header = TRUE, sep = "\t")

######## CNV filtration based on 50 % zero threshold and the elimination of any 
# NA columns

cnv_delete3=c()
R_K=nrow(Kidney_data_cnv)/2
for (i in 1:ncol(Kidney_data_cnv)){
  if (length(which(Kidney_data_cnv[,i]==0)) >= R_K ) {
    cnv_delete3=append(cnv_delete3,i)
  }
}
Kidney_data_cnv=Kidney_data_cnv[,-cnv_delete3]
K_Column_NA=names(which(colSums(is.na(Kidney_data_cnv))>0))
Delete_K=which(colnames(Kidney_data_cnv)==K_Column_NA)
Kidney_data_cnv=Kidney_data_cnv[,-Delete_K]

cnv_delete4=c()
R_L=nrow(lusc_cnv)/2
for (i in 1:ncol(lusc_cnv)){
  if (length(which(lusc_cnv[,i]==0)) >= R_L ) {
    cnv_delete4=append(cnv_delete4,i)
  }
}
lusc_cnv=lusc_cnv[,-cnv_delete4]
L_Column_NA=names(which(colSums(is.na(lusc_cnv))>0))
Delete_L=which(colnames(lusc_cnv)==L_Column_NA)
lusc_cnv=lusc_cnv[,-Delete_L]
################################################################################
# Intersection between the most expressed genes and the cancer data so have a data frame
# that  has only the most expressed genes. 

five_most_expressed_genes_Kidney <- Kidney_tumor_data[intersect(row.names(Top_5_Kidney),rownames(Kidney_tumor_data)),]

# Changing the sample names in CNV files to match the sample names in the gene expression 
# as the CNV files had "." instead of "-" in the sample name .
for (i in 1:nrow(Kidney_data_cnv)){
  Kidney_data_cnv[i,1]=gsub("-",".",Kidney_data_cnv[i,1])
}
rownames(Kidney_data_cnv)=Kidney_data_cnv[,1]
Kidney_data_cnv=Kidney_data_cnv[,-1]

# Intersection between the most expressed genes and the CNVs so have the same 
# samples in order for the regression to work

cnv_of_five_most_expressed_Kidney = Kidney_data_cnv[intersect(colnames(five_most_expressed_genes_Kidney),rownames(Kidney_data_cnv)),]
cnv_of_five_most_expressed_Kidney <- as.matrix(cbind(cnv_of_five_most_expressed_Kidney))
five_most_expressed_genes_Kidney  <- five_most_expressed_genes_Kidney[,intersect(colnames(five_most_expressed_genes_Kidney),rownames(cnv_of_five_most_expressed_Kidney)),]
five_most_expressed_genes_Kidney  <- as.matrix(cbind(five_most_expressed_genes_Kidney[,]))

# Creating a Linear model function 

linear_model <-function(x,y){
  model = lm(x ~ y)
  return(model)
}


# For loop on all genes and exporting it into separate files.

Gene_Names_K=rownames(five_most_expressed_genes_Kidney)
for (i in 1:5) {
  sink(paste(project_data_path,"/","Kidney_",Gene_Names_K[i],"_regression.txt",sep = ""))
  Genes_K_Regression=linear_model(five_most_expressed_genes_Kidney[i,],cnv_of_five_most_expressed_Kidney)
  Genes_Coeff_K= data.frame(summary(Genes_K_Regression)$coefficients)
  print("The significant CNV are :")
  print(Genes_Coeff_K[Genes_Coeff_K[,4]<=0.05,])
  print("The Full results of the regression:")
  print(summary(Genes_K_Regression))
  sink()
  closeAllConnections()
}
print("Kidney Regression Done :D")



###############################################################################
################################ Regression : Lung_data ############################
###############################################################################
###############################################################################

# same Thing in Lung_data just like the kidny 
five_most_expressed_genes_Lung_data <- Lung_tumor_data[intersect(Name,rownames(Lung_tumor_data)),]

for (i in 1:nrow(lusc_cnv)){
  lusc_cnv[i,1]=gsub("-",".",lusc_cnv[i,1])
}
rownames(lusc_cnv)=lusc_cnv[,1]
lusc_cnv=lusc_cnv[,-1]

cnv_of_five_most_expressed_Lung_data=lusc_cnv[intersect(colnames(five_most_expressed_genes_Lung_data),rownames(lusc_cnv)),]
cnv_of_five_most_expressed_Lung_data <- as.matrix(cbind(cnv_of_five_most_expressed_Lung_data))
five_most_expressed_genes_Lung_data  <- five_most_expressed_genes_Lung_data[,intersect(colnames(five_most_expressed_genes_Lung_data),rownames(cnv_of_five_most_expressed_Lung_data)),]
five_most_expressed_genes_Lung_data  <- as.matrix(cbind(five_most_expressed_genes_Lung_data[,]))

# Since predictors (CNVs) in the Lung_data case are bigger then the data points we have 
# to first penalize the CNV selection then run the regression.
# so we used package glmnet and for loop on all the genes

library(glmnet)
CNV_Accepted <- list()
Variables_Number <- dim(cnv_of_five_most_expressed_Lung_data)[2]
for ( i in 1:nrow(five_most_expressed_genes_Lung_data)){
  fit_cv <- cv.glmnet(cnv_of_five_most_expressed_Lung_data, five_most_expressed_genes_Lung_data[i,], family="gaussian", alpha=1, standardize=FALSE, nfolds=5)
  lambda <-fit_cv$lambda.min
  model <- glmnet(cnv_of_five_most_expressed_Lung_data, five_most_expressed_genes_Lung_data[i,],alpha=1, lambda=lambda, standardize=FALSE)
  coef_fit <- coef(model, s=lambda)[2:(Variables_Number+1)]
  CNV_Accepted[[i]] <- which(abs(coef_fit) > 0)
}

Gene_Names_L=rownames(five_most_expressed_genes_Lung_data)

# since some genes doesn't have any CNV we only ran the regression on specific 
# genes.

sink(paste(project_data_path,"/","Lung_data_",Gene_Names_L[1],"_regression.txt",sep = ""))
Genes_L_Regression=linear_model(five_most_expressed_genes_Lung_data[1,],cnv_of_five_most_expressed_Lung_data[,CNV_Accepted[[1]]])
Genes_Coeff_L= data.frame(summary(Genes_L_Regression)$coefficients)
print("The CNV are :")
print(Genes_Coeff_L[Genes_Coeff_L[,4]<=0.05,])
print(summary(Genes_L_Regression))
sink()
closeAllConnections()

sink(paste(project_data_path,"/","Lung_data_",Gene_Names_L[3],"_regression.txt",sep = ""))
Genes_L_Regression=linear_model(five_most_expressed_genes_Lung_data[3,],cnv_of_five_most_expressed_Lung_data[,CNV_Accepted[[3]]])
Genes_Coeff_L= data.frame(summary(Genes_L_Regression)$coefficients)
print("The CNV are :")
print(Genes_Coeff_L[Genes_Coeff_L[,4]<=0.05,])
print(summary(Genes_L_Regression))
sink()
closeAllConnections()

sink(paste(project_data_path,"/","Lung_data_",Gene_Names_L[4],"_regression.txt",sep = ""))
Genes_L_Regression=linear_model(five_most_expressed_genes_Lung_data[4,],cnv_of_five_most_expressed_Lung_data[,CNV_Accepted[[4]]])
Genes_Coeff_L= data.frame(summary(Genes_L_Regression)$coefficients)
print("The CNV are :")
print(Genes_Coeff_L[Genes_Coeff_L[,4]<=0.05,])
print(summary(Genes_L_Regression))
sink()
closeAllConnections()

print("Lung_data Regression Done :D")
###############################################################################
###############################################################################
print("Script is Done :D")






