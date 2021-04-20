#############################################################
#
# Ref to the ARTICLE
# 
#  Code to split dataset for microhabitat corrected for field  
#  Revision 11/2020 mjilani@dundee.ac.uk  
#  
# In this code we will do the following things:
# 1) Compute statistical analysis to identify individual bacteria differentiating between microhabitats (corrected for field effect)
# 2) We will also create venn diagrams for the two separate fields
# 3) We will create and save lists for the next code (creating the UpsetR plot)
# 4) We will create and save phyloseq objects for the next code (visualizing enriched bacterial communities in phyloseq)
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
#import the RDS object and subset for field
#############################################################

#Import the .RDS file 
dat_info_42k <- readRDS("JH18_rare_nowater_42K_1020.rds")

#inspect the files 
class(dat_info_42k)
dat_info_42k

dat_info_42k_dragoni <- subset_samples(dat_info_42k, Field=="Dragoni")
dat_info_42k_dragoni

dat_info_42k_minzoni <- subset_samples(dat_info_42k, Field=="Minzoni")
dat_info_42k_minzoni

#############################################################
#Construct the DESeq object for Dragoni
#############################################################

#extract count data (usually no issues in exporting this piece of information) 
JH18_42k_counts_integer <- otu_table(dat_info_42k_dragoni)
countData = as.data.frame(JH18_42k_counts_integer)
class(countData)
colnames(JH18_42k_counts_integer)

#the design file containing sample information
#If extracted directly from phyloseq objects (e.g., using sample_data) it will return still a phyloseq-class object not a data frame and DESeq won't work
#workaround: subset the sample data, first convert it to matrix and then to data frame
colData = as.data.frame(as.matrix(sample_data(dat_info_42k_dragoni)[colnames(JH18_42k_counts_integer), ]))
rownames(colData)
class(colData)

#############################################################
#DESeq calculation: Dragoni Root vs Bulk 
#############################################################

#construct a DESeq dataset combining count data and sample information 
JH18_42k_cds_dragoni <- DESeqDataSetFromMatrix(countData =countData, colData=colData , design= ~Microhabitat)

#execute the differential count analysis with the function DESeq 
JH18_42k_cds_test_dragoni <- DESeq(JH18_42k_cds_dragoni, fitType="local", betaPrior=FALSE) 

#set the contrast
JH18_42k_contrast_dragoni <- results(JH18_42k_cds_test_dragoni, contrast = c("Microhabitat",  "Root", "Bulk")) 

#visualise the data
plotMA(JH18_42k_contrast_dragoni)
#OK what is this?

#Let's add more info
plotMA(JH18_42k_contrast_dragoni, main = "Root vs Bulk (Dragoni)", sub= "Pos=Root, Neg=Bulk", 
       colNonSig = "gray60", colSig = "blue", colLine = "grey40")

#let's inspect the result file
colnames(JH18_42k_contrast_dragoni)

# extract  Genera whose adjusted p.value in a given comparison is below 0.05. 
Dragoni_FDR_005 <- JH18_42k_contrast_dragoni[(rownames(JH18_42k_contrast_dragoni)[which(JH18_42k_contrast_dragoni$padj <0.05)]), ]
rownames(Dragoni_FDR_005)

#Identify Genera enriched in dragoni root (first term of the comparison, positive fold change)
Dragoni_enriched <-  JH18_42k_contrast_dragoni[(rownames(JH18_42k_contrast_dragoni)
                                                [which(JH18_42k_contrast_dragoni$log2FoldChange > 0)]), ]

#Intersection: this is the list of genera significantly enriched in the root field
Dragoni_enriched_FDR005_root <- intersect(rownames(Dragoni_FDR_005), rownames(Dragoni_enriched))
length(Dragoni_enriched_FDR005_root)

#############################################################
#DESeq calculation: Dragoni Rhizo vs Bulk
#############################################################

#construct a DESeq dataset combining count data and sample information
JH18_42k_cds_dragoni <- DESeqDataSetFromMatrix(countData =countData, colData=colData , design= ~Microhabitat)

#execute the differential count analysis with the function DESeq 
JH18_42k_cds_test_dragoni <- DESeq(JH18_42k_cds_dragoni, fitType="local", betaPrior=FALSE) 

#set the contrast
JH18_42k_contrast_dragoni <- results(JH18_42k_cds_test_dragoni, contrast = c("Microhabitat",  "Rhizosphere", "Bulk")) 

#visualise the data
plotMA(JH18_42k_contrast_dragoni)
#OK what is this?

#Let's add more info
plotMA(JH18_42k_contrast_dragoni, main = "Rhizo vs Bulk (Dragoni)", sub = "Pos=Rhizo, Neg=Bulk",
       colNonSig = "gray60", colSig = "blue", colLine = "grey40")

#let's inspect the result file
colnames(JH18_42k_contrast_dragoni)

#extract  Genera whose adjusted p.value in a given comparison is below 0.05. 
Dragoni_FDR_005_rhizo <- JH18_42k_contrast_dragoni[(rownames(JH18_42k_contrast_dragoni)[which(JH18_42k_contrast_dragoni$padj <0.05)]), ]
rownames(Dragoni_FDR_005_rhizo)

#Identify Genera enriched in dragoni bulk (first term of the comparison, positive fold change)
Dragoni_enriched <-  JH18_42k_contrast_dragoni[(rownames(JH18_42k_contrast_dragoni)
                                                [which(JH18_42k_contrast_dragoni$log2FoldChange > 0)]), ]

#Intersection: this is the list of genera significantly enriched in the bulk field
Dragoni_enriched_FDR005_rhizo <- intersect(rownames(Dragoni_FDR_005_rhizo), rownames(Dragoni_enriched))
length(Dragoni_enriched_FDR005_rhizo)

#############################################################
#Construct the DESeq object for Minzoni
#############################################################

#extract count data (usually no issues in exporting this piece of information) 
JH18_42k_counts_integer <- otu_table(dat_info_42k_minzoni)
countData = as.data.frame(JH18_42k_counts_integer)
class(countData)
colnames(JH18_42k_counts_integer)

#the design file containing sample information
#If extracted directly from phyloseq objects (e.g., using sample_data) it will return still a phyloseq-class object not a data frame and DESeq won't work
#workaround: subset the sample data, first convert it to matrix and then to data frame
colData = as.data.frame(as.matrix(sample_data(dat_info_42k_minzoni)[colnames(JH18_42k_counts_integer), ]))
rownames(colData)
class(colData)

#############################################################
#DESeq calculation: Minzoni Root vs Bulk
#############################################################

#construct a DESeq dataset combining count data and sample information 
JH18_42k_cds_minzoni <- DESeqDataSetFromMatrix(countData =countData, colData=colData , design= ~Microhabitat)

#execute the differential count analysis with the function DESeq 
JH18_42k_cds_test_minzoni <- DESeq(JH18_42k_cds_minzoni, fitType="local", betaPrior=FALSE) 

#set the contrast
JH18_42k_contrast_minzoni <- results(JH18_42k_cds_test_minzoni, contrast = c("Microhabitat",  "Root", "Bulk")) 

#visualise the data
plotMA(JH18_42k_contrast_minzoni)
#OK what is this?

#Let's add more info
plotMA(JH18_42k_contrast_minzoni, main = "Root vs Bulk (Minzoni)", sub= "Pos=Root, Neg=Bulk", 
       colNonSig = "gray60", colSig = "blue", colLine = "grey40")

#let's inspect the result file
colnames(JH18_42k_contrast_minzoni)

#extract  Genera whose adjusted p.value in a given comparison is below 0.05. 
Minzoni_FDR_005 <- JH18_42k_contrast_minzoni[(rownames(JH18_42k_contrast_minzoni)[which(JH18_42k_contrast_minzoni$padj <0.05)]), ]
rownames(Minzoni_FDR_005)

#Identify Genera enriched in dragoni root (first term of the comparison, positive fold change)
Minzoni_enriched <-  JH18_42k_contrast_minzoni[(rownames(JH18_42k_contrast_minzoni)
                                                [which(JH18_42k_contrast_minzoni$log2FoldChange > 0)]), ]

#Intersection: this is the list of genera significantly enriched in the root field
Minzoni_enriched_FDR005_root <- intersect(rownames(Minzoni_FDR_005), rownames(Minzoni_enriched))
length(Minzoni_enriched_FDR005_root)

#############################################################
#DESeq calculation: Minzoni Rhizo vs Bulk
#############################################################

#construct a DESeq dataset combining count data and sample information
JH18_42k_cds_minzoni <- DESeqDataSetFromMatrix(countData =countData, colData=colData , design= ~Microhabitat)

#execute the differential count analysis with the function DESeq 
JH18_42k_cds_test_minzoni <- DESeq(JH18_42k_cds_minzoni, fitType="local", betaPrior=FALSE) 

#set the contrast
JH18_42k_contrast_minzoni <- results(JH18_42k_cds_test_minzoni, contrast = c("Microhabitat",  "Rhizosphere", "Bulk")) 

#visualise the data
plotMA(JH18_42k_contrast_minzoni)
#OK what is this?

#Let's add more info
plotMA(JH18_42k_contrast_minzoni, main = "Rhizo vs Bulk (Minzoni)", sub = "Pos=Rhizo, Neg=Bulk",
       colNonSig = "gray60", colSig = "blue", colLine = "grey40")

#let's inspect the result file
colnames(JH18_42k_contrast_minzoni)

#extract  Genera whose adjusted p.value in a given comparison is below 0.05. 
Minzoni_FDR_005_rhizo <- JH18_42k_contrast_minzoni[(rownames(JH18_42k_contrast_minzoni)[which(JH18_42k_contrast_minzoni$padj <0.05)]), ]
rownames(Minzoni_FDR_005_rhizo)

#Identify Genera enriched in dragoni bulk (first term of the comparison, positive fold change)
Minzoni_enriched <-  JH18_42k_contrast_minzoni[(rownames(JH18_42k_contrast_minzoni)
                                                [which(JH18_42k_contrast_minzoni$log2FoldChange > 0)]), ]

#Intersection: this is the list of genera signficanlty enriched in the bulk field
Minzoni_enriched_FDR005_rhizo <- intersect(rownames(Minzoni_FDR_005_rhizo), rownames(Minzoni_enriched))
length(Minzoni_enriched_FDR005_rhizo)

#############################################################
#Creating the intersections and saving them for an UpsetR plot (Code 6)
#############################################################
#Total amount (use this for sanity checks)
length(Dragoni_enriched_FDR005_rhizo)
length(Dragoni_enriched_FDR005_root)
length(Minzoni_enriched_FDR005_rhizo)
length(Minzoni_enriched_FDR005_root)

########Creating the data frames 
#Dragoni first
Dragoni_enriched_FDR005_root_df <- Dragoni_FDR_005[Dragoni_enriched_FDR005_root, ]
Dragoni_enriched_FDR005_rhizo_df <- Dragoni_FDR_005_rhizo[Dragoni_enriched_FDR005_rhizo, ]

#Now Minzoni
Minzoni_enriched_FDR005_root_df <- Minzoni_FDR_005[Minzoni_enriched_FDR005_root, ]
Minzoni_enriched_FDR005_rhizo_df <- Minzoni_FDR_005_rhizo[Minzoni_enriched_FDR005_rhizo, ]

#######Modifying the data frames
#We are adding the basemean columns using this code below
Dragoni_enriched_FDR005_root_df_1 <- as.data.frame(Dragoni_enriched_FDR005_root_df[ ,1])
Dragoni_enriched_FDR005_rhizo_df_1 <- as.data.frame(Dragoni_enriched_FDR005_rhizo_df[ ,1])

#Now Minzoni
Minzoni_enriched_FDR005_root_df_1 <- as.data.frame(Minzoni_enriched_FDR005_root_df [ ,1])
Minzoni_enriched_FDR005_rhizo_df_1 <- as.data.frame(Minzoni_enriched_FDR005_rhizo_df[ ,1])

######Renaming the rows and columns
#We are renaming the rownames to the ASV names they correspond to 
rownames(Dragoni_enriched_FDR005_root_df_1) <- Dragoni_enriched_FDR005_root
rownames(Dragoni_enriched_FDR005_rhizo_df_1) <- Dragoni_enriched_FDR005_rhizo
#Now Minzoni
rownames(Minzoni_enriched_FDR005_root_df_1) <- Minzoni_enriched_FDR005_root
rownames(Minzoni_enriched_FDR005_rhizo_df_1) <- Minzoni_enriched_FDR005_rhizo

#Now renaming the columns
colnames(Dragoni_enriched_FDR005_root_df_1) <- c("Dragoni_root")
colnames(Dragoni_enriched_FDR005_rhizo_df_1) <- c("Dragoni_rhizo")
#Now Minzoni
colnames(Minzoni_enriched_FDR005_root_df_1) <- c("Minzoni_root")
colnames(Minzoni_enriched_FDR005_rhizo_df_1) <- c("Minzoni_rhizo")

#####Turning the column from abundances (basemean) to boolean
Dragoni_enriched_FDR005_root_df_1[Dragoni_enriched_FDR005_root_df_1 > 1] <- 1
Dragoni_enriched_FDR005_rhizo_df_1[Dragoni_enriched_FDR005_rhizo_df_1 > 1] <- 1

#Now Minzoni
Minzoni_enriched_FDR005_root_df_1[Minzoni_enriched_FDR005_root_df_1 > 1] <- 1
Minzoni_enriched_FDR005_rhizo_df_1[Minzoni_enriched_FDR005_rhizo_df_1 > 1] <- 1

####Creating the list of total unique enriched ASVs
#Union of unique roots first
Enriched_root_ASV <- unique(union(Dragoni_enriched_FDR005_root, Minzoni_enriched_FDR005_root))
#Union of unique rhizo next
Enriched_rhizo_ASV <- unique(union(Dragoni_enriched_FDR005_rhizo, Minzoni_enriched_FDR005_rhizo))
#Final union of unique all 
Enriched_spinach_ASV <- unique(union(Enriched_root_ASV, Enriched_rhizo_ASV))
length(Enriched_spinach_ASV)

###Unifying the dataset with the list of spinach asvs
Dragoni_root_merging <- as.data.frame(Dragoni_enriched_FDR005_root_df_1[Enriched_spinach_ASV, ])
Dragoni_rhizo_merging <- as.data.frame(Dragoni_enriched_FDR005_rhizo_df_1[Enriched_spinach_ASV, ])
#Minzoni next
Minzoni_root_merging <- as.data.frame(Minzoni_enriched_FDR005_root_df_1[Enriched_spinach_ASV, ])
Minzoni_rhizo_merging <- as.data.frame(Minzoni_enriched_FDR005_rhizo_df_1[Enriched_spinach_ASV, ])

#Now we need to rename the row names to the spinach_asv list (order is conserved)
rownames(Dragoni_root_merging) <- as.vector(Enriched_spinach_ASV)
rownames(Dragoni_rhizo_merging) <- as.vector(Enriched_spinach_ASV)
#Now Minzoni
rownames(Minzoni_root_merging) <- as.vector(Enriched_spinach_ASV)
rownames(Minzoni_rhizo_merging) <- as.vector(Enriched_spinach_ASV)

#Now to rename column names to something more appropriate
colnames(Dragoni_root_merging) <- c("Dragoni_root")
colnames(Dragoni_rhizo_merging) <- c("Dragoni_rhizo")
#Now Minzoni
colnames(Minzoni_root_merging) <- c("Minzoni_root")
colnames(Minzoni_rhizo_merging) <- c("Minzoni_rhizo")

##Merging the data frames

dat_spinach_ASVs <- cbind(Dragoni_root_merging, Dragoni_rhizo_merging)
dat_spinach_ASVs <- cbind(dat_spinach_ASVs, Minzoni_root_merging)
dat_spinach_ASVs <- cbind(dat_spinach_ASVs, Minzoni_rhizo_merging)

#Great, now we just need to rename the N/A to 0
dat_spinach_ASVs[is.na(dat_spinach_ASVs)] <- 0


########Now for the intersections
#Bacteria enriched in both Dragoni and Minzoni Root
#Intersect the two fields
root_intersect_enriched <- intersect(Dragoni_enriched_FDR005_root, Minzoni_enriched_FDR005_root)
length(root_intersect_enriched)

#Uniquely enriched bacteria
#You can just subtract from the intersect from the total to find uniquely enriched bacteria,
#but for the sake of simplicity and reproducibility i have included the commands in this code
#Bacteria enriched only in Dragoni root 

root_enriched_dragoni <- setdiff(Dragoni_enriched_FDR005_root, Minzoni_enriched_FDR005_root)
length(root_enriched_dragoni)

#Bacteria enriched only in Minzoni root
root_enriched_minzoni <- setdiff(Minzoni_enriched_FDR005_root, Dragoni_enriched_FDR005_root)
length(root_enriched_minzoni)

###Now moving on to the rhizosphere 
#Bacteria enriched in both Dragoni and Minzoni Rhizo

#Intersect the two fields
rhizo_intersect_enriched <- intersect(Dragoni_enriched_FDR005_rhizo, Minzoni_enriched_FDR005_rhizo)
length(rhizo_intersect_enriched)

#Bacteria enriched only in Dragoni rhizo

rhizo_enriched_dragoni <- setdiff(Dragoni_enriched_FDR005_rhizo, Minzoni_enriched_FDR005_rhizo)
length(rhizo_enriched_dragoni)

#Bacteria enriched only in Minzoni rhizo
rhizo_enriched_minzoni <- setdiff(Minzoni_enriched_FDR005_rhizo, Dragoni_enriched_FDR005_rhizo)
length(rhizo_enriched_minzoni)

#Intersection of the 2 intersections (Core spinach microbiota)
core_enriched <- intersect(root_intersect_enriched, rhizo_intersect_enriched)
length(core_enriched)

#Microbiota enriched only in Dragoni and Minzoni
dragoni_only <- setdiff(unique(union(Dragoni_enriched_FDR005_root, Dragoni_enriched_FDR005_rhizo)), 
                        (unique(union(Minzoni_enriched_FDR005_root, Minzoni_enriched_FDR005_rhizo))))

minzoni_only <- setdiff(unique(union(Minzoni_enriched_FDR005_root, Minzoni_enriched_FDR005_rhizo)), 
                        (unique(union(Dragoni_enriched_FDR005_root, Dragoni_enriched_FDR005_rhizo))))


#We have now what we need to create venn diagrams corrected for the field effect
#Replace this with UpsetR

###Save files for UpsetR analysis in subsequent code

####################REMEMBER TO REPLACE THE HASHTAGS FOR THESE LINES BEFORE EXITING THE CODE####################
#saveRDS(dat_spinach_ASVs, file ="JH18_dat_spinach_ASVs.rds")

#saveRDS(Dragoni_enriched_FDR005_root, file ="JH18_dragoni_enriched_root_all.rds")
#saveRDS(Dragoni_enriched_FDR005_rhizo, file ="JH18_dragoni_enriched_rhizo_all.rds")
#saveRDS(Minzoni_enriched_FDR005_root, file ="JH18_minzoni_enriched_root_all.rds")
#saveRDS(Minzoni_enriched_FDR005_rhizo, file ="JH18_minzoni_enriched_rhizo_all.rds")

#saveRDS(root_intersect_enriched, file ="JH18_enriched_root_intersects.rds")
#saveRDS(root_enriched_dragoni, file ="JH18_only_dragoni_enriched_root.rds" )
#saveRDS(root_enriched_minzoni, file ="JH18_only_minzoni_enriched_root.rds")

#saveRDS(rhizo_intersect_enriched, file ="JH18_enriched_rhizo_intersect.rds")
#saveRDS(rhizo_enriched_dragoni, file ="JH18_only_dragoni_enriched_rhizo.rds")
#saveRDS(rhizo_enriched_minzoni, file ="JH18_only_minzoni_enriched_rhizo.rds")

#saveRDS(core_enriched, file ="Spinach_core_microbiota.rds")

#saveRDS(dragoni_only, file ="JH18_dragoni_only.rds")
#saveRDS(minzoni_only, file ="JH18_minzoni_only.rds")
####################REMEMBER TO REPLACE THE HASHTAGS FOR THESE LINES BEFORE EXITING THE CODE####################

#############################################################
#Generating and saving phyloseq objects for analysis in next code
#############################################################
#Creating the objects
JH18_root_enriched_intersect_phyloseq <- prune_taxa(root_intersect_enriched, dat_info_42k)
JH18_root_enriched_dragoni_phyloseq <- prune_taxa(root_enriched_dragoni, dat_info_42k)
JH18_root_enriched_minzoni_phyloseq <- prune_taxa(root_enriched_minzoni, dat_info_42k)

JH18_rhizo_enriched_intersect_phyloseq <- prune_taxa(rhizo_intersect_enriched, dat_info_42k)
JH18_rhizo_enriched_dragoni_phyloseq <- prune_taxa(rhizo_enriched_dragoni, dat_info_42k)
JH18_rhizo_enriched_minzoni_phyloseq <- prune_taxa(rhizo_enriched_minzoni, dat_info_42k)

JH18_core_enriched_phyloseq <- prune_taxa(core_enriched, dat_info_42k)

JH18_dragoni_only <- prune_taxa(dragoni_only, dat_info_42k)
JH18_minzoni_only <- prune_taxa(minzoni_only, dat_info_42k)

####Now save them

####################REMEMBER TO REPLACE THE HASHTAGS FOR THESE LINES BEFORE EXITING THE CODE####################
#saveRDS(JH18_root_enriched_intersect_phyloseq, file ="JH18_root_enriched_intersect_phyloseq.rds")
#saveRDS(JH18_root_enriched_dragoni_phyloseq, file ="JH18_root_enriched_dragoni_phyloseq.rds")
#saveRDS(JH18_root_enriched_minzoni_phyloseq, file ="JH18_root_enriched_minzoni_phyloseq.rds")

#saveRDS(JH18_root_enriched_intersect_phyloseq, file ="JH18_rhizo_enriched_intersect_phyloseq.rds")
#saveRDS(JH18_rhizo_enriched_dragoni_phyloseq, file ="JH18_rhizo_enriched_dragoni_phyloseq.rds")
#saveRDS(JH18_rhizo_enriched_minzoni_phyloseq, file ="JH18_rhizo_enriched_minzoni_phyloseq.rds")

#saveRDS(JH18_core_enriched_phyloseq, file ="JH18_core_enriched_phyloseq.rds")

#saveRDS(JH18_dragoni_only, file ="JH18_dragoni_only_phyloseq.rds")
#saveRDS(JH18_minzoni_only, file ="JH18_minzoni_only_phyloseq.rds")
####################REMEMBER TO REPLACE THE HASHTAGS FOR THESE LINES BEFORE EXITING THE CODE####################

#In next code we will use these pre-generated objects to visualize the bacteria being differentially enriched

