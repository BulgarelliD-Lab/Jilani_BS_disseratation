#############################################################
#
# Ref to the ARTICLE
# 
#  Code to perform alpha diversity calculations
#  Revision 11/2020 mjilani@dundee.ac.uk  
#  
# In this code we will do the following things:
# 1) Visualize alpha diversity indexes http://www.metagenomics.wiki/pdf/definition/alpha-beta-diversity 
# 2) Compute statistical analysis to test the hypothesis that field influences alpha diversity
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
library("ggplot2")
library("PMCMR")
library("PMCMRplus")

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
#Data visualization https://joey711.github.io/phyloseq/plot_richness-examples.html 
#################################################################################


#Create plot_richness with water (18k)
p=plot_richness(dat_info_18k, x="Description", color="Microhabitat", shape="Field", measures=c("Observed", "Shannon"))
plot <- p + geom_point(size=5, alpha=0.7)
plot + ggtitle("Alpha Diversity (18k)")

#Create plot_richness without water (42k)
p=plot_richness(dat_info_42k, x="Description", color="Microhabitat", shape="Field", measures=c("Observed", "Shannon"))
plot_42k <- p + geom_point(size=5, alpha=0.7) 
plot_42k <- plot_42k + ggtitle("Alpha Diversity (42k)")

#Change the facet text and box then save
plot_42k <- plot_42k + theme_bw() + theme(strip.background =element_rect(fill="black")) + 
  theme(strip.text = element_text(colour = 'white', size = 13)) +
  theme(axis.text.x = element_text(angle = 90))
plot_42k
##################################################################
#Statistical analysis prep for 18k
##################################################################

#first we need to create a new dataset which will be compatible with statistical ana  lysis
dat_alpha_rare_18k <-  estimate_richness(dat_info_18k, measures = c("Observed", "Shannon"))

#generate a new dataframe for statistical analysis
design_alpha_18k <-  as.data.frame(as.matrix(sample_data(dat_info_18k)[rownames(dat_alpha_rare_18k), ]))
dat_alpha_rare_info_18k <- cbind(design_alpha_18k, dat_alpha_rare_18k)

#check the new dataset: it contains both description of the samples and alpha 
dat_alpha_rare_info_18k

#plot histograms to see if the data is normally distributed 
hist_observed_18k <- hist(dat_alpha_rare_info_18k$Observed)
hist_shannon_18k <- hist(dat_alpha_rare_info_18k$Shannon)

#test the normality of the observed data 
shapiro.test(dat_alpha_rare_info_18k$Observed)
shapiro.test(dat_alpha_rare_info_18k$Shannon)

#name the histograms and add their respective shapiro-wilk results 
hist_observed_18k <- hist(dat_alpha_rare_info_18k$Observed, main='Distribution Observed (18k)', 
                          sub='Shapiro P-value= 7.157e-07')
hist_shannon_18k <- hist(dat_alpha_rare_info_18k$Shannon, main='Distribution Shannon (18k)',
                          sub='Shapiro P-value= 0.001574')

#Neither data is normally distributed
#Need to run non-parametric tests on both

##################################################################
#Statistical analysis prep for 42k
##################################################################

#first we need to create a new dataset which will be compatible with statistical analysis
dat_alpha_rare_42k <-  estimate_richness(dat_info_42k, measures = c("Observed", "Shannon"))

#generate a new dataframe for statistical analysis
design_alpha_42k <-  as.data.frame(as.matrix(sample_data(dat_info_42k)[rownames(dat_alpha_rare_42k), ]))
dat_alpha_rare_info_42k <- cbind(design_alpha_42k, dat_alpha_rare_42k)

#check the new dataset: it contains both description of the samples and alpha 
dat_alpha_rare_info_42k

#plot histograms to see if the data is normally distributed 
hist_observed <- hist(dat_alpha_rare_info_42k$Observed)
hist_shannon <- hist(dat_alpha_rare_info_42k$Shannon)

#test the normality of the observed data 
shapiro.test(dat_alpha_rare_info_42k$Observed)
shapiro.test(dat_alpha_rare_info_42k$Shannon)

#name the histograms and add their respective shapiro-wilk results 
hist_observed_42k <- hist(dat_alpha_rare_info_42k$Observed, main='Distribution Observed (42k)', 
                          sub='Shapiro P-value= 3.84e-07')
hist_shannon_42k <- hist(dat_alpha_rare_info_42k$Shannon, main='Distribution Shannon (42k)',
                         sub='Shapiro P-value= 0.001487')

#Neither data is normally distributed
#Need to run non-parametric tests on both

##################################################################
#Statistical analysis
##################################################################
#We run kruskal tests and not wilcox tests, because kruskal test is
#an extension of the wilcox test. While wilcox tests only work if   
#the variable has 2 factors, kruskal tests can work in the presence 
#of 2 or more variables

#Start with 18k
#first check whether the microhabitat has an effect
kruskal.test(Observed ~ Microhabitat, data = dat_alpha_rare_info_18k)
#signifcant effect of the microhabitat: we can now consider field
kruskal.test(Observed ~ Field, data = dat_alpha_rare_info_18k)
#no significant effect between fields

#Post-hoc Dunn's test, to look at this variance pairwise
#https://www.rdocumentation.org/packages/PMCMR/versions/4.3/topics/posthoc.kruskal.dunn.test

#Turn the microhabitat into a factor
microhabitat_character_18k <- as.factor(dat_alpha_rare_info_18k$Microhabitat)
#Now the posthoc test
posthoc.kruskal.dunn.test (x=dat_alpha_rare_info_18k$Observed, g=microhabitat_character_18k, p.adjust.method="BH")


#repeat for shannon
kruskal.test(Shannon ~ Microhabitat, data= dat_alpha_rare_info_18k)
kruskal.test(Shannon ~ Field, data = dat_alpha_rare_info_18k)
#same results for shannon as well

#Repeat the post-hoc test for shannon
posthoc.kruskal.dunn.test (x=dat_alpha_rare_info_18k$Shannon, g=microhabitat_character_18k, p.adjust.method="BH")

#Moving on to 42k
#first check whether the microhabitat has an effect on observed
kruskal.test(Observed ~ Microhabitat, data = dat_alpha_rare_info_42k)
#signifcant effect of the microhabitat: we can now consider field
kruskal.test(Observed ~ Field, data = dat_alpha_rare_info_42k)
#no significant effect between fields

#Post-hoc Dunn's test, to look at this variance pairwise
#Turn the microhabitat into a factor
microhabitat_character_42k <- as.factor(dat_alpha_rare_info_42k$Microhabitat)
posthoc.kruskal.dunn.test (x=dat_alpha_rare_info_42k$Observed, g=microhabitat_character_42k, p.adjust.method="BH")

#repeat for shannon
kruskal.test(Shannon ~ Microhabitat, data= dat_alpha_rare_info_42k)
kruskal.test(Shannon ~ Field, data = dat_alpha_rare_info_42k)

#same results for shannon as well
posthoc.kruskal.dunn.test (x=dat_alpha_rare_info_42k$Shannon, g=microhabitat_character_42k, p.adjust.method="BH")

##################################################################
#We can now conclude that, in the tested conditions, the fields do  
#not have a significant difference in  alpha diversity. Microhabitat
#has a significant effect on the alpha diversity as expected
##################################################################