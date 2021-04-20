#############################################################
#
# Ref to the Flowchart
# 
#  Code to import and pre-process the ASV dataset
#  Revision 11/2020 mjilani@dundee.ac.uk  
#  
# In this code we will do the following things:
# 1) Import the phyloseq object (ASV Dataset) https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0061217 
# 2) Remove potential contaminants such as, chloroplast and mitochondria derived sequences
# 3) Further filter samples that were poorly sequenced and/or low count observations (interfering with the analysis)
# 
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

#retrieve R and package versions and compare to the uploaded file in gtHub for the reproducibility of the code
sessionInfo()

#UoD_JHI_Davide-Aman
setwd("C:/Users/Aman J/University of Dundee/Davide Bulgarelli (Staff) - Aman_Hons/R_Analysis")

#############################################################
#import the RDS object
#############################################################

#Import the .RDS file 
dat_info <- readRDS("JH18_silva138_dada2.rds")

#inspect the file 
class(dat_info)
dat_info

##################################################################
#Pre-processing: remove sequencing reads assigned to Chloroplast and Mitochondria
#################################################################

#remove chloroplasts
dat_info_no_chlor <-subset_taxa(dat_info, (Order!="Chloroplast") | is.na(Order))
dat_info_no_chlor

#remove mitochondrias
dat_info_no_plants <-subset_taxa(dat_info_no_chlor, (Family!="Mitochondria") | is.na(Family))
dat_info_no_plants

##################################################################
#Prune putative contaminant ASVs
##################################################################

#Import the list of contaminant ASVs from JH06 library
Contaminant_ASVs <- read.delim("JH06_contaminant_ASVs_ids.txt", header = FALSE)

#identify the proportion of putative contaminants in the merged object
dat_info_contaminants <- intersect(taxa_names(dat_info_no_plants),Contaminant_ASVs)
dat_info_contaminants

#We have no contaminants, so we are good to proceed to the next step
#########################################################################
# Remove ASVs assigned to NA at phylum level
#########################################################################

dat_info_no_plants_1 <- subset_taxa(dat_info_no_plants, Phylum!= "NA")
dat_info_no_plants_1

#########################################################################
# Abundance threshold 1: remove samples with low counts
#########################################################################

sort(sample_sums(dat_info_no_plants_1))
#No low counts, no need to remove any samples

#########################################################################
#Inspect the experimental design, to work out the minimum amount of 
#samples for which a bacterium could be present
#########################################################################

sample_data(dat_info_no_plants_1)
levels(sample_data(dat_info_no_plants_1)$Description)

#########################################################################
# Abundance threshold: 20 reads, set as 16% the minimum number of samples
# This proportion is given as we have 10 reps for 6 treatments:
# if we have one bacterium in just one treatment this ~16% of the samples
########################################################################

dat_info_no_plants_2 = filter_taxa(dat_info_no_plants_1, function(x) sum(x > 20) > (0.16 *length(x)), TRUE)
dat_info_no_plants_2 
sort(sample_sums(dat_info_no_plants_2))

#filtering is very severe for water samples, because we only have 2 each  
#water samples instead of ten each, we do not expect to able to identify
#water specific 

#ratio filtered reads/total reads
ratio <- sum(sample_sums(dat_info_no_plants_2))/sum(sample_sums(dat_info_no_plants_1))*100
ratio

#########################################################################
# Aggregate samples at genus level (Note: NArm set to false as Phylum NA ASVs pruned above)
#########################################################################

dat_info_genus <- tax_glom(dat_info_no_plants_2, taxrank= "Genus", NArm=FALSE, bad_empty=c(NA, "", " ", "\t"))

#compare the two objects
#ASVs
dat_info_no_plants_2
sort(sample_sums(dat_info_no_plants_2))

#Genera
dat_info_genus 
sort(sample_sums(dat_info_genus))
hist(sample_sums(dat_info_genus), main="Distribution of Reads")

#We will create two phyloseq objects, one with water to see its effect on the
#overall bacterial composition, and another one with plant only samples for
#a more detailed investigation

####################################################################################
# Rarefy at an even sequencing depth (18,000) and "freeze" these objects for downstream analyses
####################################################################################

#ASVs : ignore the warnings, the object will be saved right after

####################REMEMBER TO REPLACE THE HASHTAGS FOR THESE LINES BEFORE EXITING THE CODE####################
#JH18_rare_18K <- rarefy_even_depth(dat_info_genus, 18000)
#saveRDS(JH18_rare_18K, file ="JH18_rare_water_18K_1020.rds")
#Sample file with water created
####################REMEMBER TO REPLACE THE HASHTAGS FOR THESE LINES BEFORE EXITING THE CODE####################

####################################################################################
# Remove water samples prior to creating a w/o water sample file
####################################################################################

dat_info_nowater <- prune_samples(sample_sums(dat_info_genus) > 42000, dat_info_genus)
dat_info_nowater
sort(sample_sums(dat_info_nowater))

####################################################################################
# Rarefy at an even sequencing depth (42,000) and "freeze" these objects for downstream analyses
####################################################################################

#ASVs : ignore the warnings, the object will be saved right after

####################REMEMBER TO REPLACE THE HASHTAGS FOR THESE LINES BEFORE EXITING THE CODE####################
#JH18_rare_42K <- rarefy_even_depth(dat_info_nowater, 42000)
#saveRDS(JH18_rare_42K, file ="JH18_rare_nowater_42K_1020.rds")
####################REMEMBER TO REPLACE THE HASHTAGS FOR THESE LINES BEFORE EXITING THE CODE####################

#End
