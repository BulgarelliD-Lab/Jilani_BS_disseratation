#############################################################
#
# Ref to the ARTICLE
# 
#  Code to split dataset for microhabitat
#  Revision 11/2020 mjilani@dundee.ac.uk  
#  
# In this code we will do the following things:
# 1) Compute statistical analysis to identify individual bacteria differentiating between microhabitats
# 2) Check proportion of bacteria in the roots that are derived from the bulk soil 
#
#############################################################

#############################################################
# Clean-up the memory and start a new session
#############################################################

rm(list=ls())
dev.off()

#############################################################
# Libraries required
#############################################################

#required packages 
library("phyloseq")
library ("DESeq2")

#retrieve R and package versions and compare to the uploaded file in gitHub for the reproducibility of the code
sessionInfo()

#UoD_JHI_Davide-Aman
setwd("C:/Users/Aman J/University of Dundee/Davide Bulgarelli (Staff) - Aman_Hons/R_Analysis")

#############################################################
#import the RDS object
#############################################################

#Import the .RDS file 
dat_info_42k <- readRDS("JH18_rare_nowater_42K_1020.rds")

#inspect the files 
class(dat_info_42k)
dat_info_42k

#Davide's solution for field effect corrected for microhabitat
dat_info_42k_bulk <- subset_samples(dat_info_42k, Microhabitat=="Bulk")
dat_info_42k_rhizo <- subset_samples(dat_info_42k, Microhabitat=="Rhizosphere")
dat_info_42k_root<- subset_samples(dat_info_42k, Microhabitat=="Root")

#############################################################
#Construct the DESeq object for Bulk microhabitat
#############################################################

#extract count data (usually no issues in exporting this piece of information) 
JH18_42k_counts_integer <- otu_table(dat_info_42k_bulk)
countData = as.data.frame(JH18_42k_counts_integer)
class(countData)
colnames(JH18_42k_counts_integer)

#the design file containing sample information
#If extracted directly from phyloseq objects (e.g., using sample_data) it will return still a phyloseq-class object not a data frame and DESeq won't work
#workaround: subset the sample data, first convert it to matrix and then to data frame
colData = as.data.frame(as.matrix(sample_data(dat_info_42k_bulk)[colnames(JH18_42k_counts_integer), ]))
rownames(colData)
class(colData)

#############################################################
#DESeq calculation: field effect for Bulk microhabitat
#############################################################

#construct a DESeq dataset combining count data and sample information 
JH18_42k_cds <- DESeqDataSetFromMatrix(countData =countData, colData=colData , design= ~Field)

#execute the differential count analysis with the function DESeq 
JH18_42k_cds_test <- DESeq(JH18_42k_cds, fitType="local", betaPrior=FALSE) 

#set the contrast
JH18_42k_contrast <- results(JH18_42k_cds_test, contrast = c("Field",  "Dragoni", "Minzoni")) 

#visualise the data
plotMA(JH18_42k_contrast, main = "Field effect (Bulk)", sub= "Pos=Dragoni, Neg=Minzoni", 
      colNonSig = "gray60", colSig = "blue", colLine = "grey40")

#let's inspect the result file
colnames(JH18_42k_contrast)

# extract  Genera whose adjusted p.value in a given comparison is below 0.05. 
Field_FDR_005 <- JH18_42k_contrast[(rownames(JH18_42k_contrast)[which(JH18_42k_contrast$padj <0.05)]), ]
rownames(Field_FDR_005)

#Identify Genera enriched in dragoni bulk (first term of the comparison, positive fold change)
Field_Dragoni_enriched <-  JH18_42k_contrast[(rownames(JH18_42k_contrast)[which(JH18_42k_contrast$log2FoldChange > 0)]), ]
#Identify Genera enriched in minzoni bulk (second term of the comparison, negative fold change)
Field_Minzoni_enriched <-  JH18_42k_contrast[(rownames(JH18_42k_contrast)[which(JH18_42k_contrast$log2FoldChange < 0)]), ]

#Intersection: this is the list of genera signficanlty enriched in the bulk dragoni field
Field_Dragoni_enriched_FDR005_bulk <- intersect(rownames(Field_FDR_005), rownames(Field_Dragoni_enriched))
length(Field_Dragoni_enriched_FDR005_bulk)
#Intersection: this is the list of genera significantly enriched in the bulk minzoni field
Field_Minzoni_enriched_FDR005_bulk <- intersect(rownames(Field_FDR_005), rownames(Field_Minzoni_enriched))
length(Field_Minzoni_enriched_FDR005_bulk)

#############################################################
#Construct the DESeq object for Rhizosphere microhabitat
#############################################################

#extract count data (usually no issues in exporting this piece of information) 
JH18_42k_counts_integer <- otu_table(dat_info_42k_rhizo)
countData = as.data.frame(JH18_42k_counts_integer)
class(countData)
colnames(JH18_42k_counts_integer)

#the design file containing sample information
#If extracted directly from phyloseq objects (e.g., using sample_data) it will return still a phyloseq-class object not a data frame and DESeq won't work
#workaround: subset the sample data, first convert it to matrix and then to data frame
colData = as.data.frame(as.matrix(sample_data(dat_info_42k_rhizo)[colnames(JH18_42k_counts_integer), ]))
rownames(colData)
class(colData)

#############################################################
#DESeq calculation: field effect for rhizosphere
#############################################################

#construct a DESeq dataset combining count data and sample information
JH18_42k_cds <- DESeqDataSetFromMatrix(countData =countData, colData=colData , design= ~Field)

#execute the differential count analysis with the function DESeq 
JH18_42k_cds_test <- DESeq(JH18_42k_cds, fitType="local", betaPrior=FALSE) 

#set the contrast
JH18_42k_contrast <- results(JH18_42k_cds_test, contrast = c("Field",  "Dragoni", "Minzoni")) 

#visualise the data
plotMA(JH18_42k_contrast, main = "Field effect (Rhizosphere)",  sub= "Pos=Dragoni, Neg=Minzoni",
      colNonSig = "gray60", colSig = "blue", colLine = "grey40")

#let's inspect the result file  
colnames(JH18_42k_contrast)

# extract  Genera whose adjusted p.value in a given comparison is below 0.05. 
Field_FDR_005 <- JH18_42k_contrast[(rownames(JH18_42k_contrast)[which(JH18_42k_contrast$padj <0.05)]), ]
rownames(Field_FDR_005)

#Identify Genera enriched in dragoni (first term of the comparison, pos fold change)
Field_Dragoni_enriched <-  JH18_42k_contrast[(rownames(JH18_42k_contrast)[which(JH18_42k_contrast$log2FoldChange > 0)]), ]
#Identify Genera enriched in minzoni bulk (second term of the comparison, negative fold change)
Field_Minzoni_enriched <-  JH18_42k_contrast[(rownames(JH18_42k_contrast)[which(JH18_42k_contrast$log2FoldChange < 0)]), ]

#Intersection: this is the list of genera significantly enriched in the dragoni rhizosphere
Field_Dragoni_enriched_FDR005_rhizo <- intersect(rownames(Field_FDR_005), rownames(Field_Dragoni_enriched))
length(Field_Dragoni_enriched_FDR005_rhizo)
#Intersection: this is the list of genera significantly enriched in the rhizosphere minzoni field
Field_Minzoni_enriched_FDR005_rhizo <- intersect(rownames(Field_FDR_005), rownames(Field_Minzoni_enriched))
length(Field_Minzoni_enriched_FDR005_rhizo)

#############################################################
#Construct the DESeq object for Root microhabitat
#############################################################

#extract count data (usually no issues in exporting this piece of information) 
JH18_42k_counts_integer <- otu_table(dat_info_42k_root)
countData = as.data.frame(JH18_42k_counts_integer)
class(countData)
colnames(JH18_42k_counts_integer)

#the design file containing sample information
#If extracted directly from phyloseq objects (e.g., using sample_data) it will return still a phyloseq-class object not a data frame and DESeq won't work
#workaround: subset the sample data, first convert it to matrix and then to data frame
colData = as.data.frame(as.matrix(sample_data(dat_info_42k_root)[colnames(JH18_42k_counts_integer), ]))
rownames(colData)
class(colData)

#############################################################
#DESeq calculation: field effect for root
#############################################################

#construct a DESeq dataset combining count data and sample information t
JH18_42k_cds <- DESeqDataSetFromMatrix(countData =countData, colData=colData , design= ~Field)

#execute the differential count analysis with the function DESeq 
JH18_42k_cds_test <- DESeq(JH18_42k_cds, fitType="local", betaPrior=FALSE) 

#set the contrast
JH18_42k_contrast <- results(JH18_42k_cds_test, contrast = c("Field",  "Dragoni", "Minzoni")) 

#visualise the data
plotMA(JH18_42k_contrast, main = "Field effect (Root)", sub = "Pos=Dragoni, Neg=Minzoni", 
       colNonSig = "gray60", colSig = "blue", colLine = "grey40")

#let's inspect the result file
colnames(JH18_42k_contrast)

# extract  Genera whose adjusted p.value in a given comparison is below 0.05. 
Field_FDR_005 <- JH18_42k_contrast[(rownames(JH18_42k_contrast)[which(JH18_42k_contrast$padj <0.05)]), ]
rownames(Field_FDR_005)

#Identify Genera enriched in the root (first term of the comparison, pos fold change)
Field_Dragoni_enriched <-  JH18_42k_contrast[(rownames(JH18_42k_contrast)[which(JH18_42k_contrast$log2FoldChange > 0)]), ]
#Identify Genera enriched in the root (second term of the comparison, neg fold change)
Field_Minzoni_enriched <-  JH18_42k_contrast[(rownames(JH18_42k_contrast)[which(JH18_42k_contrast$log2FoldChange < 0)]), ]

#Intersection: this is the list of genera significantly enriched in the dragoni root
Field_Dragoni_enriched_FDR005_root <- intersect(rownames(Field_FDR_005), rownames(Field_Dragoni_enriched))
length(Field_Dragoni_enriched_FDR005_root)
#Intersection: this is the list of genera significantly enriched in the minzoni root
Field_Minzoni_enriched_FDR005_root <- intersect(rownames(Field_FDR_005), rownames(Field_Minzoni_enriched))
length(Field_Minzoni_enriched_FDR005_root)

#############################################################
#Checking the proportion of bacteria that are derived from the soil 
#############################################################
Dragoni_proportion <- intersect(Field_Dragoni_enriched_FDR005_root, Field_Dragoni_enriched_FDR005_bulk)
Dragoni_proportion

Minzoni_proportion <- intersect(Field_Minzoni_enriched_FDR005_root, Field_Minzoni_enriched_FDR005_bulk)
Minzoni_proportion

#Can't use the soil to predict the outcome of the plant, differences in spinach dragoni and minzoni are not
#just from differences in the soil, plant effect could be identified under agricultural conditions. Look
#at Minute 40 9/12/2020