#############################################################
#
# Ref to the ARTICLE
# 
#  Code to generate the UpsetR plot and the phyloseq bar plots
#  Revision 11/2020 mjilani@dundee.ac.uk  
#  
# In this code we will do the following things:
# 1) Generate an UpsetR plot to visualize intersections between field microhabitats and their enriched bacteria
# 2) We will use phyloseq to visualize the different enriched bacterial communities
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
library("UpSetR")
library("ggplot2")

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
#Import the files for the UpsetR plot and the phyloseq analysis
#############################################################
#UpsetR plot files
dat_spinach_ASVs <- readRDS("JH18_dat_spinach_ASVs.rds")
Dragoni_enriched_FDR005_root <- readRDS("JH18_dragoni_enriched_root_all.rds")
Dragoni_enriched_FDR005_rhizo <- readRDS("JH18_dragoni_enriched_rhizo_all.rds")
Minzoni_enriched_FDR005_root <- readRDS("JH18_minzoni_enriched_root_all.rds")
Minzoni_enriched_FDR005_rhizo <- readRDS("JH18_minzoni_enriched_rhizo_all.rds")

root_intersect_enriched <- readRDS("JH18_enriched_root_intersects.rds")
root_enriched_dragoni <- readRDS("JH18_only_dragoni_enriched_root.rds" )
root_enriched_minzoni <- readRDS("JH18_only_minzoni_enriched_root.rds")

rhizo_intersect_enriched <- readRDS("JH18_enriched_rhizo_intersect.rds")
rhizo_enriched_dragoni <- readRDS("JH18_only_dragoni_enriched_rhizo.rds")
rhizo_enriched_minzoni <- readRDS("JH18_only_minzoni_enriched_rhizo.rds")

core_enriched <- readRDS("Spinach_core_microbiota.rds")

dragoni_only <- readRDS("JH18_dragoni_only.rds")
minzoni_only <- readRDS("JH18_minzoni_only.rds")
                    
#Phyloseq analysis files
JH18_root_enriched_intersect <- readRDS("JH18_root_enriched_intersect_phyloseq.rds")
JH18_root_enriched_dragoni <- readRDS("JH18_root_enriched_dragoni_phyloseq.rds")
JH18_root_enriched_minzoni <- readRDS("JH18_root_enriched_minzoni_phyloseq.rds")

JH18_root_enriched_intersect <- readRDS("JH18_rhizo_enriched_intersect_phyloseq.rds")
JH18_rhizo_enriched_dragoni <- readRDS("JH18_rhizo_enriched_dragoni_phyloseq.rds")
JH18_rhizo_enriched_minzoni <- readRDS("JH18_rhizo_enriched_minzoni_phyloseq.rds")

JH18_core_enriched <- readRDS("JH18_core_enriched_phyloseq.rds")

JH18_dragoni_only <- readRDS("JH18_dragoni_only_phyloseq.rds")
JH18_minzoni_only <- readRDS("JH18_minzoni_only_phyloseq.rds")

#############################################################
#Generating the UpsetR plot (https://academic.oup.com/bioinformatics/article/33/18/2938/3884387)
#############################################################

#We've done all the pre-processing in the previous code, so this is really just one 
#line of code to get a graphical output
list_all_intersection <- list("Dragoni_root", "Minzoni_root", "Dragoni_rhizo", "Minzoni_rhizo")
list_dragoni_shared <- list("Dragoni_root", "Dragoni_rhizo")
list_minzoni_shared <- list("Minzoni_root", "Minzoni_rhizo")
list_dragoni_rhizo <- list("Dragoni_rhizo")
list_minzoni_rhizo <- list("Minzoni_rhizo")
list_dragoni_root <- list("Dragoni_root")
list_minzoni_root <- list("Minzoni_root")

lists_merged <- list(list_all_intersection, list_dragoni_shared, list_dragoni_rhizo, list_dragoni_root,
                     list_minzoni_shared, list_minzoni_rhizo, list_minzoni_root)

upset(dat_spinach_ASVs, sets = c("Dragoni_root", "Minzoni_root", "Dragoni_rhizo", "Minzoni_rhizo"), 
      sets.bar.color = "#56B4E9", order.by = "freq", intersections = lists_merged, 
      matrix.color = "#000000", main.bar.color = "#000000", mainbar.y.label = "Number of Shared Bacterial Families",
      text.scale = 1.3)  


#############################################################
#Phyloseq analysis
#############################################################

#Starting with Dragoni and Minzoni root
plot_dragoni_root <- plot_bar(JH18_root_enriched_dragoni, "Microhabitat", fill="Family", facet_grid= "Field", title= "Dragoni Root Enriched")
plot_minzoni_root <- plot_bar(JH18_root_enriched_minzoni, "Microhabitat", fill="Family", facet_grid= "Field", title= "Minzoni Root Enriched")

#Now for the Dragoni and Minzoni rhizo
plot_dragoni_rhizo <- plot_bar(JH18_rhizo_enriched_dragoni, "Microhabitat", fill="Family", facet_grid= "Field", title= "Dragoni Rhizo Enriched")
plot_minzoni_rhizo <- plot_bar(JH18_rhizo_enriched_minzoni, "Microhabitat", fill="Family", facet_grid= "Field", title= "Minzoni Rhizo Enriched")

#Now for the main 3 things we are interested in (the first intersection (in UpsetR), the second intersection(Dragoni only),
#and the third intersection (Minzoni only))

#Core enriched bacteria (shared between all microhabitats and fields)
plot_core_enriched <- plot_bar(JH18_core_enriched, "Microhabitat", fill="Family", facet_grid= "Field", title= "Core Enriched")

#Dragoni only enriched
plot_dragoni_only <- plot_bar(JH18_dragoni_only, "Microhabitat", fill="Family", facet_grid= "Field", title= "Dragoni Only Enriched")

#Minzoni only enriched
plot_minzoni_only <- plot_bar(JH18_minzoni_only, "Microhabitat", fill="Family", facet_grid= "Field", title= "Minzoni Only Enriched")

#############################################################
#Finalizing the bar plots
#############################################################
#Create the 30 colors, color-blind friendly palette, courtesy of https://coolors.co/ for generating 
#and https://davidmathlogic.com/colorblind/ for double-checking
cbPalette <- c("#F633F8", "#A9BED7", "#E4CCA2", "#8A705F", "#634649", "#675B59", 
               "#2F4F5C", "#12227A", "#5B4B5F", "#2EC4B6", "#E71D36", "#FF9F1C",
               "#840032", "#E59500", "#E08DAC", "#6A7FDB", "#1ED7F0", "#AF3B6E",
               "#21FA90", "#478978", "#B6EEA6", "#F2FF49", "#355834", "#6E633D",
               "#EEF5DB", "#D5C5C8", "#AA3A84", "#953AE0", "#47E066")
#Please note that this palette is made with Protoanomaly (the type of colorblindness I have) in mind. It is not
#suitable for the other types of colorblindess.

#In the following order, this code is going to:
#Apply the colorblind friendly palette, apply the bw theme from ggplot2, change the facet grid text color and size, as well
#as the facet box color
#This will finalize our graphs 
plot_dragoni_root_final <- plot_dragoni_root + scale_fill_manual(values=cbPalette) + 
  theme_bw() + theme(strip.background =element_rect(fill="black")) + theme(strip.text = element_text(colour = 'white', size = 13))
plot_minzoni_root_final <- plot_minzoni_root + scale_fill_manual(values=cbPalette) + 
  theme_bw() + theme(strip.background =element_rect(fill="black")) + theme(strip.text = element_text(colour = 'white', size = 13))

plot_dragoni_rhizo_final <- plot_dragoni_rhizo + scale_fill_manual(values=cbPalette) + 
  theme_bw() + theme(strip.background =element_rect(fill="black")) + theme(strip.text = element_text(colour = 'white', size = 13))
plot_minzoni_rhizo_final <- plot_minzoni_rhizo + scale_fill_manual(values=cbPalette) + 
  theme_bw() + theme(strip.background =element_rect(fill="black")) + theme(strip.text = element_text(colour = 'white', size = 13))

plot_core_enriched_final <- plot_core_enriched + scale_fill_manual(values=cbPalette) + 
  theme_bw() + theme(strip.background =element_rect(fill="black")) + theme(strip.text = element_text(colour = 'white', size = 13))
plot_dragoni_only_final <- plot_dragoni_only + scale_fill_manual(values=cbPalette) + 
  theme_bw() + theme(strip.background =element_rect(fill="black")) + theme(strip.text = element_text(colour = 'white', size = 13))
plot_minzoni_only_final <- plot_minzoni_only + scale_fill_manual(values=cbPalette) + 
  theme_bw() + theme(strip.background =element_rect(fill="black")) + theme(strip.text = element_text(colour = 'white', size = 13))

#Now run the plots so that we can save them for presentating later
plot_dragoni_root_final
plot_minzoni_root_final

plot_dragoni_rhizo_final
plot_minzoni_rhizo_final

plot_core_enriched_final
plot_dragoni_only_final
plot_minzoni_only_final

#############################################################
#Network analysis
#############################################################

#Plot extremes first
plot_net(dat_info_42k, maxdist= 0.3, color= "Microhabitat", shape= "Field")
plot_net(dat_info_42k, maxdist= 0.7, color= "Microhabitat", shape= "Field")
plot_net(dat_info_42k, maxdist= 0.5, color= "Microhabitat", shape= "Field")

#Seems like the middle plot will be the one we want
p_42k <- plot_net(dat_info_42k, maxdist= 0.5, color= "Microhabitat", shape= "Field")
p_42k + ggtitle("42k Network Plot")
p_42k

#Further Justification
#First remove bulk from samples
JH18_core_enriched_wo_bulk <- subset_samples(JH18_core_enriched, Microhabitat!="Bulk")

plot_net(JH18_core_enriched_wo_bulk, maxdist= 0.3, color= "Microhabitat", shape= "Field")
plot_net(JH18_dragoni_only, maxdist= 0.3, color= "Microhabitat", shape= "Field")
plot_net(JH18_minzoni_only, maxdist= 0.2, color= "Microhabitat", shape= "Field")
