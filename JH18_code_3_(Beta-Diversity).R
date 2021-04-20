#############################################################
#
# Ref to the ARTICLE
# 
#  Code to perform beta diversity calculations
#  Revision 11/2020 mjilani@dundee.ac.uk  
#  
# In this code we will do the following things:
# 1) Visualize beta diversity indexes http://www.metagenomics.wiki/pdf/definition/alpha-beta-diversity 
# 2) Compute statistical analysis to test the hypothesis that fields have a difference in beta diversity
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
library ("ggplot2")
library ("vegan")

#retrieve R and package versions and compare to the uploaded file in gtHub for the reproducibility of the code
sessionInfo()

#UoD_JHI_Davide-Aman
setwd("C:/Users/Aman J/University of Dundee/Davide Bulgarelli (Staff) - Aman_Hons/R_Analysis")

#############################################################
#import the RDS objects
#############################################################

#Import the .RDS file 
dat_info_18k <- readRDS("JH18_rare_water_18K_1020.rds")
dat_info_42k <- readRDS("JH18_rare_nowater_42K_1020.rds")

#inspect the files 
class(dat_info_18k)
dat_info_18k

class(dat_info_42k)
dat_info_42k

#################################################################################
#Data visualisation for 18k (https://joey711.github.io/phyloseq/distance.html)
#################################################################################

#key to any beta diversity calculation is the construction of the so called distance matrix
#in other words we will "transform" differences in microbial composition across samples as "distances"
#there are many options in term of computation, we will be using the Bray-Curtis one which is one of the most common/used
#which is sensitive to ASVs presence/absence and abundances https://en.wikipedia.org/wiki/Bray%E2%80%93Curtis_dissimilarity 
BC_18k <- phyloseq::distance(dat_info_18k, "bray")
#let's inspect the file
BC_18k


#We will be using Canonical Analysis of Principal coordinates (CAP)
dat_info_CAP_BC_18k <- ordinate(dat_info_18k, "CAP","bray", ~ Field)

#Follow the same style as the alpha diversity (Microhabitat=Colored, Field=Shapes, etc.)
p_18k =  plot_ordination(dat_info_18k, dat_info_CAP_BC_18k, shape ="Field", color = "Microhabitat")
p_18k = p_18k + geom_point(size = 5, alpha = 0.75)
p_18k + ggtitle("CAP Field 18k, Bray distance")

#################################################################################
#Data visualization for 42k
#################################################################################

#Create the bray-distance matrix for the 42k dataset
BC_42k <- phyloseq::distance(dat_info_42k, "bray")
#let's inspect the file
BC_42k


#Create the file used for statistical analysis later
dat_info_CAP_BC_42k <- ordinate(dat_info_42k, "CAP","bray", ~ Field)

p_42k = plot_ordination(dat_info_42k, dat_info_CAP_BC_42k, shape ="Field", color = "Microhabitat")
p_42k = p_42k + geom_point(size = 5, alpha = 0.75)
p_42k <- p_42k + ggtitle("CAP Field 42k, Bray Distance Dissimilarity")

#Conform the 42k to final figure guidelines
p_42k <- p_42k + theme_bw() 
p_42k
#################################################################################
#Statistical analysis 
#################################################################################
#for 18k
#anova on the axis (this is to assess the robustness of the ordination)
anova(dat_info_CAP_BC_18k, permutations=5000)
#If this anova result is insignificant, we are not allowed to present this ordination
#because the CAP forces the graph to look nicer, the anova checks if any random presentation
#will be better than the one we have generated. 

#permutational analysis of variance (also called adonis)
adonis(BC_18k ~ Microhabitat * Field , data= as.data.frame(as.matrix(sample_data(dat_info_18k))), permutations = 5000)
#R2 = the amount of variation explained by that variable (convert to percentages)
#Residuals is everything not listed as an independent variable, it may be the variation between PCR machines, 
#and/or people carrying it out, if residual goes above 50% it becomes a problem

#for 42k
#anova on the axis (this is to assess the robustness of the ordination)
anova(dat_info_CAP_BC_42k, permutations=5000)

#permutational analysis of variance (also called adonis)
adonis(BC_42k ~ Microhabitat * Field , data= as.data.frame(as.matrix(sample_data(dat_info_42k))), permutations = 5000)

##################################################################
#We can now conclude that our CAP is solid, the microhabitat is the major effect accounting for 80% of the variance.
#leaving water "inside or outside" the analysis will not have any effect.
##################################################################