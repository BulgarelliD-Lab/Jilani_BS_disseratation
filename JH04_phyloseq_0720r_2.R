#############################################################
#
# Ref to the ARTICLE
# 
# Rodrigo Alegria Terrazas (r.z.alegriaterrazas@dundee.ac.uk) and Davide Bulgarelli (d.bulgarelli@dundee.ac.uk)
# 09/07/20
# script to reproduce calculations and figures presented in the manuscript
# 
# 
#############################################################

#############################################################
# Libraries and functions required
#############################################################

#These initial commands are required to clean-up the memory and start a new session
rm(list=ls())
dev.off()

#setwd("C:/Users/rzalegriaterrazas/University of Dundee/Davide Bulgarelli (Staff) - Davide_lab_manuscripts/Rodrigo_B1K_2019/R_data")
setwd("/cluster/db/ralegriaterrazas/R_data_JH04")


#library location
.libPaths()


#load the required packages 
library("phyloseq")
library ("foreach")
library("ggplot2")
library("vegan")
library ("ape")
library("PMCMR")
library("DESeq2")
library("UpSetR")


#R session info
sessionInfo()

#############################################################
#############################################################

#import the count matrix and the desing file

#OTU table this file has been generated using QIIME 1.9.0. In the OTU ids, OTU abundance information has been removed
#file JH04_otu_table_SILVA132_97_24_nc2_woc_nch_nmi.txt is Supplementary worksheet ws2 in Supplementary Dataset 1
dat_info <- read.delim("JH04_otu_table_SILVA132_97_24_nc2_woc_nch_nmi.txt" ,skip=1, sep = "\t", header=T, row.names= 1, blank.lines.skip = FALSE)
colnames(dat_info)

#inspect the file using the command dim 
dim(dat_info)

#extract the total number of reads clustered at OTU 97% identiy (the number underneath the Hv identifier represents the total number of reads clustered for that sample) 
OTU_97_reads <- sort(colSums(dat_info[, 1:76]))
OTU_97_reads
#total reads  6646864
OTU_97_reads_sum <- sum(colSums(dat_info[, 1:76]))
OTU_97_reads_sum

#design file
#file JH04_otu_table_SILVA132_97_24_nc2_woc_nch_nmi.txt is Supplementary worksheet ws1 in Supplementary Dataset 1
design <- read.delim("Map_JH04_24.txt", sep = "\t", header=TRUE, row.names=1)
design

#rename row names of dat_info as noPlants 
noPlants <- rownames(dat_info)
#and check the results
length(noPlants)# 11212

#save this piece of information for qiime purposes. We need to create a a new OTU table in QIIME to generate the taxa tables
#write(noPlants, "JH04_noPlant_OTUs_id_24.txt")

#create a new count matrix 
dat_count <- dat_info[, rownames(design)]
dat_count_noplants <- dat_info[noPlants, rownames(design)]
dim(dat_count_noplants)

#and a new taxa table
dat_tax_noPlants <- as.data.frame(dat_info[rownames(dat_count_noplants), 77])
rownames(dat_tax_noPlants) <- rownames(dat_count_noplants)
#now we need to save the above file and in excel generate a tax table where column represents a taxonomic rank
#write.table(dat_tax_noPlants, file="JH04_dat_tax_noPlants_24.txt", sep="\t")

##########
#Figure 1
##########

#remove soil samples
rhizo_samples <- rownames(design)[which(design$Description!= "Soil")]
design_rhizo <- design[rhizo_samples, ]

#re-order the factors
design_rhizo$Geo <- ordered(design_rhizo$Geo, levels=c("D1","D2","C1","C2", "N","M"))

#stem dw: Figure 1B
p <- ggplot(design_rhizo, aes(x=Geo, y=Stemdryweight, fill=Geo)) + geom_boxplot()
p + geom_jitter( size=5,shape=21, position=position_jitter(0.2))+ scale_fill_manual(values = c("yellow","orange", "deepskyblue1","blue", "green", "magenta"))
#saved as .eps and integrate in the composite figure in Illustratior

##root-shoot ratio: Figure 1C
p <- ggplot(design_rhizo, aes(x=Geo, y=Rootshoot, fill=Geo)) + geom_boxplot()
p + geom_jitter( size=5,shape=21, position=position_jitter(0.2))+ scale_fill_manual(values = c("yellow","orange", "deepskyblue1","blue", "green", "magenta"))
#saved as .eps and integrate in the composite figure in Illustratior

#statistical analysis: data distribution
shapiro.test(design_rhizo$Stemdryweight)
#dryweight not normally distributed (below o.o5)
shapiro.test(design_rhizo$Rootshoot)
#Rootshoot not normally distributed (below o.o5)

################################

#Non parametric anaova: stem dw
kruskal.test(Stemdryweight ~ Geo, data =design_rhizo )
posthoc.kruskal.dunn.test (x=design_rhizo$Stemdryweight, g=design_rhizo$Geo, p.adjust.method="BH")
#Non parametric anaova: roo:shoot
kruskal.test(Rootshoot ~ Geo, data =design_rhizo )
posthoc.kruskal.dunn.test (x=design_rhizo$Rootshoot, g=design_rhizo$Geo, p.adjust.method="BH")

#####################
#Figure S2
#####################

#DNA_concentration test
#Re-order the factors
design_DNA <-design
design_DNA$Geo <- ordered(design_DNA$Geo, levels=c("Bulk","D1","D2","C1","C2", "N","M"))

#DNA_Con: Figure S2
p <- ggplot(design_DNA, aes(x=Geo, y=DNA_Con, fill=Geo)) + geom_boxplot()
p + geom_jitter( size=5,shape=21, position=position_jitter(0.2))+ scale_fill_manual(values = c("brown","yellow","orange", "deepskyblue1","blue", "green", "magenta"))
#saved as .eps and edit figure in Illustratior

shapiro.test(design_DNA$DNA_Con)
#DNA_Con not normally distributed (below o.o5)

#Non parametric anova: DNA_Con
kruskal.test(DNA_Con ~ Geo, data =design_DNA )
posthoc.kruskal.dunn.test (x=design$DNA_Con, g=design$Geo, p.adjust.method="BH")


################
#Figure 2
#################

#This table was generated in QIIME w/o plant-derived OTUs and biological replicates are already averaged according to the levels of the factor "Geo" in the mapping file
#file Geo_otu_table_L2.txt is Supplementary worksheet ws3 in Supplementary Dataset 1
dat_info_taxa_Phylum <- read.delim("Geo_otu_table_L2.txt", skip=1, sep = "\t", header=T, row.names= 1, blank.lines.skip = FALSE)
colnames(dat_info_taxa_Phylum)
dim(dat_info_taxa_Phylum)

# transform the data in % for visualisation
dat_count_taxa_Phylum <- dat_info_taxa_Phylum
dat_norm_Phylum <- ((dat_count_taxa_Phylum)/colSums(dat_count_taxa_Phylum,na=T)) * 100 
dat_norm_Phylum[1:5, ]

# determine the average % of reads for each phylum
Phylum_mean_sorted <- dat_norm_Phylum[(order(-rowSums(dat_norm_Phylum))), ] 

##Realtive abundance for above 1%
Phylum_mean_topRank <- Phylum_mean_sorted[rownames(Phylum_mean_sorted)[which(rowMeans(Phylum_mean_sorted) > 1)], ]
dim(Phylum_mean_topRank)
colSums(Phylum_mean_topRank)

#overall
mean(colSums(Phylum_mean_topRank))

Phylum_mean_topRank<- as.matrix(Phylum_mean_topRank) 
#first we need to arrange the samples in a coherent way (i.e. according to the experiment)
colnames(Phylum_mean_topRank)


Phylum_mean_topRank_samples <- c("Bulk","D1","D2","C1","C2", "N","M")

#now we will use the order of samples we have generated above to sort the count matrix
Phylum_mean_topRank_ordered <- Phylum_mean_topRank[ ,Phylum_mean_topRank_samples] 

#Let's check the first 5 rows of the two files just to chek that everything is OK
Phylum_mean_topRank_ordered[1:5, ]
Phylum_mean_topRank[1:5, ]

#Now we can plot them
# Stacked Bar Plot with Colors and Legend
barplot(Phylum_mean_topRank_ordered, main="Phylum Distribution",
        xlab="Samples", ylab = "Reads %", ylim = c(0,100), col=c("darkblue","red", "green", "cyan", "orange", "brown", "black", "magenta", "grey", "yellow", "darkgreen"), beside=FALSE,   legend = rownames(Phylum_mean_topRank_ordered))

#OK due to size limits the legend covers part of the graph, save the graph as .eps file and in illustrator you can adjust this (and other)  graphical issues
#barplot_no legend: you can save this as image and then reconstruct the leged using the previous figure. Note: the order of samples can be inferred from the command on line 166
barplot(Phylum_mean_topRank_ordered, main="Phylum Distribution JH04",
        xlab="Samples", ylab = "Reads %", ylim = c(0,100), col=c("darkblue","red", "green", "cyan", "orange", "brown", "black", "magenta", "grey", "yellow", "darkgreen"), beside=FALSE)

#Use the files generated with the commands above to create a figure that combines both barchart and legends in an output readable for the users.


###########################################
#Figure 3; Additional File 1 Figure S3 & S4
###########################################

#Generate a Phyloseq object
#a) The OTU Table counts
JH04_OTU <- otu_table(dat_count_noplants, taxa_are_rows=TRUE)
dim(JH04_OTU)
#b) The taxonomy information
#Note this is a new file generated in excell from the output of the command of lines 93-97
#it is a tab-delimited file with 8 columns, the header names are: OTU id; "Kingdom", "Phylum",  "Class",   "Order",   "Family",  "Genus",   "Species",
#and the empty cells are filled with the term 'Unassigned'
JH04_taxa_ordered <- read.delim ("JH04_dat_tax_noPlants_ordered_24.txt", sep = "\t", row.names=1, header=T, blank.lines.skip = FALSE)
JH04_taxa <- tax_table(as.matrix(JH04_taxa_ordered))
dim(JH04_taxa)

#c) The mapping file 
JH04_map <- sample_data(design)

#d) The phylogenetic tree: the OTU table has been generated using a closed reference approach agains the SILVA database use the corresponding phylogenetic tree 
JH04_tree <- read_tree("97_otus.tre")

#check whether the tree is rooted
is.rooted(JH04_tree)


# Is there any unique class in the database?
unique_OTUs <- unique(JH04_taxa_ordered[,2])
unique_OTUs

#We could use a Bacteroidetes as an outgroup 
outgroup <- JH04_taxa_ordered[grepl("D_1__Bacteroidetes", JH04_taxa_ordered$Phylum), ]

#How many of theseD_1__Bacteroidetes do we have? 
dim(outgroup)

#Identify the top abundant OTU of this group
sort(rowSums(dat_count_noplants[rownames(outgroup), ]))

#bit confusing but with 969 reads, OTU DQ450754.1.1381 is relatively abundant and part of Bacteroidetes, let's pick this OTU as an outgroup
newRoot = c("DQ450754.1.1381")
JH04_tree <- root(JH04_tree,newRoot,resolve.root = TRUE)

#is it rooted now?
is.rooted(JH04_tree)

#merge the files and create the phyloseq object
JH04_data_phyloseq <- merge_phyloseq(JH04_OTU, JH04_taxa, JH04_map,  JH04_tree)

#inspect the generated data
JH04_data_phyloseq
sum(colSums(otu_table(JH04_data_phyloseq)))
dim(dat_count_noplants)
sum(colSums(dat_count_noplants))


##########
#Figure Additional File 1 Figure S3
##########

#generate a rerefied dataset
#JH04_data_phyloseq_rare <- rarefy_even_depth(JH04_data_phyloseq, rngseed=TRUE)

#extract and save the OTU table for reproducibiity of the code
#JH04_data_phyloseq_rare_table <- as.data.frame(otu_table(JH04_data_phyloseq_rare))
##inspect the generated file
#class(JH04_data_phyloseq_rare_table)
#dim(JH04_data_phyloseq_rare_table)

#save the file for the reproducibility of the code
#write.table(JH04_data_phyloseq_rare_table, file="JH04_data_phyloseq_rare_table_counts.txt", sep="\t")

#import the rarefied OTU counts (note file name counts2.txt this file has been generated in excel and includes the #OTU ID as header of the first column) worksheet ws6 in Supplementary Dataset 1
dat_count_rare <- read.delim("JH04_data_phyloseq_rare_table_counts_2.txt", sep = "\t", header=T, row.names= 1, blank.lines.skip = FALSE)

#inspect the generated file, 8744, 76, 849680
dim(dat_count_rare)
sum(colSums(dat_count_rare))

#taxonomy information
JH04_taxa_rare <- tax_table(as.matrix(JH04_taxa_ordered[rownames(dat_count_rare), ]))
dim(JH04_taxa_rare)


#generate a new phyloseq object wich will contain only the rarefied counts and the design file (only these two pieces of information are required for alphadiversity calculation)
JH04_OTU_rare <- otu_table(dat_count_rare, taxa_are_rows=TRUE)
JH04_data_rare_phyloseq <- merge_phyloseq(JH04_OTU_rare, JH04_taxa_rare, JH04_map, JH04_tree)

#Inspect the generated file
JH04_data_rare_phyloseq 

#number of reads per sample: original dataset differeSR01 reads per sample
sample_sums(JH04_data_phyloseq)

#number of reads per sample: rarefied dataset, they are all the same. 11180 reads/sample
sample_sums(JH04_data_rare_phyloseq)

#Index calculation
JH04_alpha_rare <-  estimate_richness(JH04_data_rare_phyloseq, measures = c("Observed", "Chao1", "Shannon"))

#generate a new dataframe for data visualisation and statistical analysis
design2 <- design[rownames(JH04_alpha_rare), ]
JH04_alpha_rare_info <- cbind(design2, JH04_alpha_rare)

#check the new dataset: it contains both description of the samples and alpha 
JH04_alpha_rare_info 

#generate a box plot for data visualisation
#re-order the factors
#JH04_alpha_rare_info$Description <- ordered(JH04_alpha_rare_info$Description, levels=c("Soil", "Desert", "Coast", "North","Modern", "Traditional"))
JH04_alpha_rare_info$Geo <- ordered(JH04_alpha_rare_info$Geo, levels=c("Bulk","D1","D2","C1","C2", "N","M"))

#chao 1                                                                                            
p <- ggplot(JH04_alpha_rare_info, aes(x=Geo, y=Chao1, fill=Geo)) + geom_boxplot()
p + geom_jitter( size=5,shape=21, position=position_jitter(0.2))+ scale_fill_manual(values = c("brown","yellow","orange", "deepskyblue1","blue", "green", "magenta"))


#observed
p <- ggplot(JH04_alpha_rare_info, aes(x=Geo, y=Observed, fill=Geo)) + geom_boxplot()
p + geom_jitter( size=5,shape=21, position=position_jitter(0.2))+ scale_fill_manual(values = c("brown","yellow","orange", "deepskyblue1","blue", "green", "magenta"))


#shannon
p <- ggplot(JH04_alpha_rare_info, aes(x=Geo, y=Shannon, fill=Geo)) + geom_boxplot()
p + geom_jitter( size=5,shape=21, position=position_jitter(0.2))+ scale_fill_manual(values = c("brown","yellow","orange", "deepskyblue1","blue", "green", "magenta"))

#test the normality of the observed data 
shapiro.test(JH04_alpha_rare_info$Observed)#not normal
shapiro.test(JH04_alpha_rare_info$Chao1)#not normal
shapiro.test(JH04_alpha_rare_info$Shannon)#normal


#Non parametric test 
#Observed- significant
wilcox.test(Observed ~ Microhabitat, data = JH04_alpha_rare_info)
#Shannon-significant
wilcox.test(Shannon ~ Microhabitat, data = JH04_alpha_rare_info)
#Chao1- not significant
wilcox.test(Chao1 ~ Microhabitat, data = JH04_alpha_rare_info)

#############
#Removing soil samples: stats within rhizosphere samples
#############

JH04_alpha_rare_info_rhizo<- subset(JH04_alpha_rare_info, Microhabitat=="Rhizosphere")
dim(JH04_alpha_rare_info_rhizo)

#Chao1 (significant )
kruskal.test(Chao1 ~ Geo, data =JH04_alpha_rare_info_rhizo )
#(significant D1 vs D2, C1, C2, N, M)
posthoc.kruskal.dunn.test (x=JH04_alpha_rare_info_rhizo$Chao1, g=JH04_alpha_rare_info_rhizo$Geo, p.adjust.method="BH")

#Shannon (slightly significant)
kruskal.test(Shannon ~ Geo, data = JH04_alpha_rare_info_rhizo)
posthoc.kruskal.dunn.test (x=JH04_alpha_rare_info_rhizo$Shannon, g=JH04_alpha_rare_info_rhizo$Geo, p.adjust.method="BH")

#Observed (not significant)
kruskal.test(Observed ~ Description, data = JH04_alpha_rare_info_rhizo)
posthoc.kruskal.dunn.test (x=JH04_alpha_rare_info_rhizo$Observed, g=JH04_alpha_rare_info_rhizo$Geo, p.adjust.method="BH")

#######################################################################################################################################
#Figure 3: betadiversity calculation (proportional transformation)
###########################################################################################################################################

#abundance filtering
#Remove OTUs not seen more than 5 times in at least 20% of the samples. This protects against an OTU with small mean & trivially large C.V.
JH04_data_phyloseq_3 = filter_taxa(JH04_data_phyloseq, function(x) sum(x > 5) > (0.2*length(x)), TRUE)
JH04_data_phyloseq_3
sum(colSums(otu_table(JH04_data_phyloseq))) #6646864, 11212 taxa and 76 samples
sum(colSums(otu_table(JH04_data_phyloseq_3))) #6253979, 2379 taxa and 76 samples

#Transform the count in relative abundance
JH04_data_phyloseq_prop <- transform_sample_counts(JH04_data_phyloseq_3,  function(x) 1e+06 * x/sum(x))

#identify the levels
sample_data(JH04_data_phyloseq_prop)$Description
#Order the levels according to a desired sequence (e.g., first the soil then the ecotypes)
sample_data(JH04_data_phyloseq_prop)$Description <- ordered(sample_data(JH04_data_phyloseq_prop)$Description, levels=c("Soil", "Desert", "Coast", "North","Modern"))

#CAP WU distance (Figure 3)
JH04_data_phyloseq_prop_CAP_WU <- ordinate(JH04_data_phyloseq_prop, "CAP","wunifrac", ~ Description)
plot_ordination(JH04_data_phyloseq, JH04_data_phyloseq_prop_CAP_WU, color = "Description")

#assign shapes to Soil and color to Ecotype
p=plot_ordination(JH04_data_phyloseq_prop, JH04_data_phyloseq_prop_CAP_WU , shape ="Location", color = "Description")
p = p + geom_point(size = 5, alpha = 0.75)
#this is the same color coding of the the manuscript J . EVOL. BIOL. 26 ( 2013) 163
p = p + scale_colour_manual(values = c("brown","orange","blue","green", "magenta","black","red"))
p + ggtitle("CAP 16S data, WU distance")

#anova on the axis
anova(JH04_data_phyloseq_prop_CAP_WU, permutations=5000)

#CAP bray distance (Additional file 1: Figure S3)
JH04_data_phyloseq_prop_CAP <- ordinate(JH04_data_phyloseq_prop, "CAP", "bray", ~ Description)
plot_ordination(JH04_data_phyloseq, JH04_data_phyloseq_prop_CAP, color = "Description", shape ="Location")

#assign shapes to Soil and color to Ecotype
p=plot_ordination(JH04_data_phyloseq_prop, JH04_data_phyloseq_prop_CAP , shape ="Location", color = "Description")
p = p + geom_point(size = 5, alpha = 0.75)
p = p + scale_colour_manual(values = c("brown","orange","blue","green", "magenta","orange", "black","red"))
p + ggtitle("CAP 16S data, Bray distance")


#Permanova calculation: Table 2

#microhabitat effect
#BC distance
BC <- phyloseq::distance(JH04_data_phyloseq_prop, "bray")
adonis(BC ~ Microhabitat , data= design, permutations = 5000)
#WU distnace
WU <- phyloseq::distance(JH04_data_phyloseq_prop, "unifrac", weighted= TRUE)
adonis(WU ~ Microhabitat, data= design, permutations = 5000)

#ecotype effect (rhizosphere only)
JH04_data_phyloseq_prop_rhizo <- subset_samples(JH04_data_phyloseq_prop, Microhabitat == "Rhizosphere")
design_rhizosphere <- design[colnames(otu_table(JH04_data_phyloseq_prop_rhizo)), ]

#BC distance
BC <- phyloseq::distance(JH04_data_phyloseq_prop_rhizo, "bray")
adonis(BC ~ Geo, data= design_rhizosphere, permutations = 5000)
#WU distnace
WU <- phyloseq::distance(JH04_data_phyloseq_prop_rhizo, "unifrac", weighted= TRUE)
adonis(WU ~ Geo, data= design_rhizosphere, permutations = 5000)

###################################################################################################
#Figures 4 and 5
#Hypothesis testing using DESeq2 package 
###################################################################################################
Map_JH04_DESeq <-design
#otu TABLE FILTERED for low abundance
#extract count data and 
JH04_OTU_counts_integer <- otu_table(JH04_data_phyloseq_3)
countData = as.data.frame(JH04_OTU_counts_integer)
class(countData)
colnames(JH04_OTU_counts_integer)

#the design file containing sample information
colData = Map_JH04_DESeq[colnames(JH04_OTU_counts_integer), ]
class(colData)

#construct a DESeq dataset combining count data and sample information
#A DESeqDataSet object must have an associated design formula  The formula should be a tilde (???) followed by the variables of interest. In this case the column "Description" in the desing file depicts the variable of interest
JH04_cds <- DESeqDataSetFromMatrix(countData =countData, colData=colData, design= ~ Geo)

#execute the differential count analysis with the function DESeq 
JH04_cds_test <- DESeq(JH04_cds, fitType="local", betaPrior=FALSE) 

#define the OTUs significantly enriched in the rhizosphere samples
Soil_Rhizosphere_D1 <- results(JH04_cds_test , contrast = c("Geo",  "Bulk", "D1")) 
Soil_Rhizosphere_D2 <- results(JH04_cds_test , contrast = c("Geo",  "Bulk", "D2")) 
Soil_Rhizosphere_C1 <- results(JH04_cds_test , contrast = c("Geo",  "Bulk", "C1")) 
Soil_Rhizosphere_C2 <- results(JH04_cds_test , contrast = c("Geo",  "Bulk", "C2")) 
Soil_Rhizosphere_N <- results(JH04_cds_test , contrast = c("Geo",  "Bulk", "N"))
Soil_Rhizosphere_M <- results(JH04_cds_test , contrast = c("Geo",  "Bulk", "M"))

#inspect a result file
Soil_Rhizosphere_D1  
mcols(Soil_Rhizosphere_D1 , use.names=TRUE)

# extract  OTUs whose adjusted p.value in a given comparison is below 0.05. 
Soil_Rhizosphere_D1_FDR_005 <- Soil_Rhizosphere_D1[(rownames(Soil_Rhizosphere_D1)[which(Soil_Rhizosphere_D1$padj <0.05)]), ]
Soil_Rhizosphere_D2_FDR_005 <- Soil_Rhizosphere_D2[(rownames(Soil_Rhizosphere_D2)[which(Soil_Rhizosphere_D2$padj <0.05)]), ]
Soil_Rhizosphere_C1_FDR_005 <- Soil_Rhizosphere_C1[(rownames(Soil_Rhizosphere_C1)[which(Soil_Rhizosphere_C1$padj <0.05)]), ]
Soil_Rhizosphere_C2_FDR_005 <- Soil_Rhizosphere_C2[(rownames(Soil_Rhizosphere_C2)[which(Soil_Rhizosphere_C2$padj <0.05)]), ]
Soil_Rhizosphere_N_FDR_005 <- Soil_Rhizosphere_N[(rownames(Soil_Rhizosphere_N)[which(Soil_Rhizosphere_N$padj <0.05)]), ]
Soil_Rhizosphere_M_FDR_005 <- Soil_Rhizosphere_M[(rownames(Soil_Rhizosphere_M)[which(Soil_Rhizosphere_M$padj <0.05)]), ]

#Identify OTUs enriched in the rhizosphere (second term of the comparison, negative fold change)
Soil_Rhizosphere_D1_enriched <-  Soil_Rhizosphere_D1[(rownames(Soil_Rhizosphere_D1)[which(Soil_Rhizosphere_D1$log2FoldChange < 0)]), ]
Soil_Rhizosphere_D2_enriched <-  Soil_Rhizosphere_D2[(rownames(Soil_Rhizosphere_D2)[which(Soil_Rhizosphere_D2$log2FoldChange < 0)]), ]
Soil_Rhizosphere_C1_enriched <-  Soil_Rhizosphere_C1[(rownames(Soil_Rhizosphere_C1)[which(Soil_Rhizosphere_C1$log2FoldChange < 0)]), ]
Soil_Rhizosphere_C2_enriched <-  Soil_Rhizosphere_C2[(rownames(Soil_Rhizosphere_C2)[which(Soil_Rhizosphere_C2$log2FoldChange < 0)]), ]
Soil_Rhizosphere_N_enriched <-  Soil_Rhizosphere_N[(rownames(Soil_Rhizosphere_N)[which(Soil_Rhizosphere_N$log2FoldChange < 0)]), ]
Soil_Rhizosphere_M_enriched <-  Soil_Rhizosphere_M[(rownames(Soil_Rhizosphere_M)[which(Soil_Rhizosphere_M$log2FoldChange < 0)]), ]

#Commands on lines 440/454  provides lists of OTUs fulfilling the imposed criteria. To identify OTUs significantly enriched in the compartment of interest an intersection of these lists is required
Soil_Rhizosphere_D1_enriched_FDR005 <- intersect(rownames(Soil_Rhizosphere_D1_FDR_005), rownames(Soil_Rhizosphere_D1_enriched))
Soil_Rhizosphere_D2_enriched_FDR005 <- intersect(rownames(Soil_Rhizosphere_D2_FDR_005), rownames(Soil_Rhizosphere_D2_enriched))
Soil_Rhizosphere_C1_enriched_FDR005 <- intersect(rownames(Soil_Rhizosphere_C1_FDR_005), rownames(Soil_Rhizosphere_C1_enriched))
Soil_Rhizosphere_C2_enriched_FDR005 <- intersect(rownames(Soil_Rhizosphere_C2_FDR_005), rownames(Soil_Rhizosphere_C2_enriched))
Soil_Rhizosphere_N_enriched_FDR005 <- intersect(rownames(Soil_Rhizosphere_N_FDR_005), rownames(Soil_Rhizosphere_N_enriched))
Soil_Rhizosphere_M_enriched_FDR005 <- intersect(rownames(Soil_Rhizosphere_M_FDR_005), rownames(Soil_Rhizosphere_M_enriched))

#Define the number of OTUs significantly enriched in and differentiating between each plant-associated compartment and unplanted soil
length(Soil_Rhizosphere_D1_enriched_FDR005)# 725
length(Soil_Rhizosphere_D2_enriched_FDR005)# 759
length(Soil_Rhizosphere_C1_enriched_FDR005)# 777
length(Soil_Rhizosphere_C2_enriched_FDR005)# 761
length(Soil_Rhizosphere_N_enriched_FDR005)# 791
length(Soil_Rhizosphere_M_enriched_FDR005)# 787

##################
#GEO
##################
#D1 vs.M (4-53) 16 actinobacteria

D1_M <- results(JH04_cds_test , contrast = c("Geo",  "D1", "M"))

# extract  OTUs whose adjusted p.value in a given comparison is below 0.05. 
D1_M_FDR_005 <- D1_M[(rownames(D1_M)[which(D1_M$padj <0.05)]), ]

#Identify OTUs enriched in Modern (first term of the comparison)
D1_enriched_vs_M <-  D1_M[(rownames(D1_M)[which(D1_M$log2FoldChange > 0)]), ]

#Identify OTUs enriched in Desert (second term of the comparison)
M_enriched_vs_D1 <-  D1_M[(rownames(D1_M)[which(D1_M$log2FoldChange < 0)]), ]

#intersect the datasets
D1_enriched_vs_M_FDR_005 <- intersect(rownames(D1_M_FDR_005), rownames(D1_enriched_vs_M))
M_enriched_vs_D1_FDR_005 <- intersect(rownames(D1_M_FDR_005), rownames(M_enriched_vs_D1))

#further filtering for OTUs significantly enriched vs. soil in the respective comparision
D1_enriched_vs_M_FDR_005_2 <- intersect(D1_enriched_vs_M_FDR_005, Soil_Rhizosphere_D1_enriched_FDR005)
M_enriched_vs_D1_FDR_005_2 <- intersect(M_enriched_vs_D1_FDR_005, Soil_Rhizosphere_M_enriched_FDR005)

#Supplementary worksheet ws 7
M_D1_enriched_taxa<- JH04_taxa_ordered[M_enriched_vs_D1_FDR_005_2, ]
M_D1_enriched_info <- cbind(as.data.frame(D1_M[M_enriched_vs_D1_FDR_005_2, ]), M_D1_enriched_taxa)
#write.table(M_D1_enriched_info, file="ws5.txt")

########################################]
#D2 vs. M (12-104)
D2_M <- results(JH04_cds_test , contrast = c("Geo",  "D2", "M"))

# extract  OTUs whose adjusted p.value in a given comparison is below 0.05. 
D2_M_FDR_005 <- D2_M[(rownames(D2_M)[which(D2_M$padj <0.05)]), ]

#Identify OTUs enriched in Modern (first term of the comparison)
D2_enriched_vs_M <-  D2_M[(rownames(D2_M)[which(D2_M$log2FoldChange > 0)]), ]

#Identify OTUs enriched in Desert (second term of the comparison)
M_enriched_vs_D2 <-  D2_M[(rownames(D2_M)[which(D2_M$log2FoldChange < 0)]), ]

#intersect the datasets
D2_enriched_vs_M_FDR_005 <- intersect(rownames(D2_M_FDR_005), rownames(D2_enriched_vs_M))
M_enriched_vs_D2_FDR_005 <- intersect(rownames(D2_M_FDR_005), rownames(M_enriched_vs_D2))

#further filtering for OTUs significantly enriched vs. soil in the respective comparision
D2_enriched_vs_M_FDR_005_2 <- intersect(D2_enriched_vs_M_FDR_005, Soil_Rhizosphere_D2_enriched_FDR005)
M_enriched_vs_D2_FDR_005_2 <- intersect(M_enriched_vs_D2_FDR_005, Soil_Rhizosphere_M_enriched_FDR005)

#Supplementary worksheet ws 8
M_D2_enriched_taxa<- JH04_taxa_ordered[M_enriched_vs_D2_FDR_005_2, ]
M_D2_enriched_info <- cbind(as.data.frame(D2_M[M_enriched_vs_D2_FDR_005_2, ]), M_D2_enriched_taxa)
#write.table(M_D2_enriched_info, file="ws6.txt")

#################################################################
#C1 vs M (1-8)
C1_M <- results(JH04_cds_test , contrast = c("Geo",  "C1", "M"))

# extract  OTUs whose adjusted p.value in a given comparison is below 0.05. 
C1_M_FDR_005 <- C1_M[(rownames(C1_M)[which(C1_M$padj <0.05)]), ]

#Identify OTUs enriched in Modern (first term of the comparison)
C1_enriched_vs_M <-  C1_M[(rownames(C1_M)[which(C1_M$log2FoldChange > 0)]), ]

#Identify OTUs enriched in Desert (second term of the comparison)
M_enriched_vs_C1 <-  C1_M[(rownames(C1_M)[which(C1_M$log2FoldChange < 0)]), ]

#intersect the datasets
C1_enriched_vs_M_FDR_005 <- intersect(rownames(C1_M_FDR_005), rownames(C1_enriched_vs_M))
M_enriched_vs_C1_FDR_005 <- intersect(rownames(C1_M_FDR_005), rownames(M_enriched_vs_C1))

#further filtering for OTUs significantly enriched vs. soil in the respective comparision
C1_enriched_vs_M_FDR_005_2 <- intersect(C1_enriched_vs_M_FDR_005, Soil_Rhizosphere_C1_enriched_FDR005)
M_enriched_vs_C1_FDR_005_2 <- intersect(M_enriched_vs_C1_FDR_005, Soil_Rhizosphere_M_enriched_FDR005)

#Supplementary worksheet ws 9
M_C1_enriched_taxa<- JH04_taxa_ordered[M_enriched_vs_C1_FDR_005_2, ]
M_C1_enriched_info <- cbind(as.data.frame(C1_M[M_enriched_vs_C1_FDR_005_2, ]), M_C1_enriched_taxa)
#write.table(M_C1_enriched_info, file="ws7.txt")

####################################################################################################
#C2 vs M (1-7)
C2_M <- results(JH04_cds_test , contrast = c("Geo",  "C2", "M"))

# extract  OTUs whose adjusted p.value in a given comparison is below 0.05. 
C2_M_FDR_005 <- C2_M[(rownames(C2_M)[which(C2_M$padj <0.05)]), ]

#Identify OTUs enriched in Modern (first term of the comparison)
C2_enriched_vs_M <-  C2_M[(rownames(C2_M)[which(C2_M$log2FoldChange > 0)]), ]

#Identify OTUs enriched in Desert (second term of the comparison)
M_enriched_vs_C2 <-  C2_M[(rownames(C2_M)[which(C2_M$log2FoldChange < 0)]), ]

#intersect the datasets
C2_enriched_vs_M_FDR_005 <- intersect(rownames(C2_M_FDR_005), rownames(C2_enriched_vs_M))
M_enriched_vs_C2_FDR_005 <- intersect(rownames(C2_M_FDR_005), rownames(M_enriched_vs_C2))

#further filtering for OTUs significantly enriched vs. soil in the respective comparision
C2_enriched_vs_M_FDR_005_2 <- intersect(C2_enriched_vs_M_FDR_005, Soil_Rhizosphere_C2_enriched_FDR005)
M_enriched_vs_C2_FDR_005_2 <- intersect(M_enriched_vs_C2_FDR_005, Soil_Rhizosphere_M_enriched_FDR005)

#Supplementary worksheet ws 10
M_C2_enriched_taxa<- JH04_taxa_ordered[M_enriched_vs_C2_FDR_005_2, ]
M_C2_enriched_info <- cbind(as.data.frame(C2_M[M_enriched_vs_C2_FDR_005_2, ]), M_C2_enriched_taxa)
#write.table(M_C2_enriched_info, file="ws8.txt")

#select enriched Actinobacteria in Modern vs Ecotype comparisons for Fig_Sup_Act
###########################################################################################################
#M vs N (10-6)
M_N <- results(JH04_cds_test , contrast = c("Geo",  "M", "N"))

# extract  OTUs whose adjusted p.value in a given comparison is below 0.05. 
M_N_FDR_005 <- M_N[(rownames(M_N)[which(M_N$padj <0.05)]), ]

#Identify OTUs enriched in Modern (first term of the comparison)
M_enriched_vs_N <-  M_N[(rownames(M_N)[which(M_N$log2FoldChange > 0)]), ]

#Identify OTUs enriched in Desert (second term of the comparison)
N_enriched_vs_M <-  M_N[(rownames(M_N)[which(M_N$log2FoldChange < 0)]), ]

#intersect the datasets
M_enriched_vs_N_FDR_005 <- intersect(rownames(M_N_FDR_005), rownames(M_enriched_vs_N))
N_enriched_vs_M_FDR_005 <- intersect(rownames(M_N_FDR_005), rownames(N_enriched_vs_M))

#further filtering for OTUs significantly enriched vs. soil in the respective comparision
M_enriched_vs_N_FDR_005_2 <- intersect(M_enriched_vs_N_FDR_005, Soil_Rhizosphere_M_enriched_FDR005)
N_enriched_vs_M_FDR_005_2 <- intersect(N_enriched_vs_M_FDR_005, Soil_Rhizosphere_N_enriched_FDR005)

#Supplementary worksheet ws 11
M_N_enriched_taxa<- JH04_taxa_ordered[M_enriched_vs_N_FDR_005_2, ]
M_N_enriched_info <- cbind(as.data.frame(M_N[M_enriched_vs_N_FDR_005_2, ]), M_N_enriched_taxa)
#write.table(M_N_enriched_info, file="ws9.txt")


######################################################################################################################
#SUPPLementary Fig ACT

#Select enriched Actinobacteria in Modern vs Ecotype comparisons 
M_D1_enriched_info_A <- M_D1_enriched_info[grepl("D_1__Actinobacteria", M_D1_enriched_info$Phylum), ]
M_D2_enriched_info_A <- M_D2_enriched_info[grepl("D_1__Actinobacteria", M_D2_enriched_info$Phylum), ]
M_C1_enriched_info_A <- M_C1_enriched_info[grepl("D_1__Actinobacteria", M_C1_enriched_info$Phylum), ]
M_C2_enriched_info_A <- M_C2_enriched_info[grepl("D_1__Actinobacteria", M_C2_enriched_info$Phylum), ]
M_N_enriched_info_A <- M_N_enriched_info[grepl("D_1__Actinobacteria", M_N_enriched_info$Phylum), ]

#Select OTU names
D1_enriched_info_A<-rownames(M_D1_enriched_info_A)
D2_enriched_info_A<-rownames(M_D2_enriched_info_A)
C1_enriched_info_A<-rownames(M_C1_enriched_info_A)
C2_enriched_info_A<-rownames(M_C2_enriched_info_A)
N_enriched_info_A<-rownames(M_N_enriched_info_A)

#Join OTU names
enriched_info_A <- c(N_enriched_info_A, D1_enriched_info_A, D2_enriched_info_A, C2_enriched_info_A, C1_enriched_info_A)  

#select enriched Actinobacteria from low abundance trimmed OTU table
JH04_enriched_info_A<-otu_table(JH04_data_phyloseq_prop)[enriched_info_A, ]
dim(JH04_enriched_info_A)

#total number of reads of enriched Actinobacteria
OTU_97_reads_Act <- colSums(JH04_enriched_info_A)

#now sort per sample 
sort(OTU_97_reads_Act)

#create a dataset to visualise the proportion of reads per sample 
OTU_97_reads_Act_d <- as.data.frame(OTU_97_reads_Act)

#rename the columns in the generated datasets
colnames(OTU_97_reads_Act_d) <- c("Act_reads")

#combine these datasets with the design file
design_Act <- cbind(design, OTU_97_reads_Act_d)

#order levels
design_Act$Geo <- ordered(design_Act$Geo, levels=c("Bulk","D1","D2","C1","C2", "N","M"))

#Plot
p <- ggplot(design_Act, aes(x=Geo, y=Act_reads, fill=Geo)) + geom_boxplot()
p + geom_jitter( size=5,shape=21, position=position_jitter(0.2))+ scale_fill_manual(values = c("brown","yellow","orange", "deepskyblue1","blue", "green", "magenta"))

##################
#D1 vs. D2 (7-3)
D1_D2 <- results(JH04_cds_test , contrast = c("Geo",  "D1", "D2"))

# extract  OTUs whose adjusted p.value in a given comparison is below 0.05. 
D1_D2_FDR_005 <- D1_D2[(rownames(D1_D2)[which(D1_D2$padj <0.05)]), ]

#Identify OTUs enriched in Modern (first term of the comparison)
D1_enriched_vs_D2 <-  D1_D2[(rownames(D1_D2)[which(D1_D2$log2FoldChange > 0)]), ]

#Identify OTUs enriched in Desert (second term of the comparison)
D2_enriched_vs_D1 <-  D1_D2[(rownames(D1_D2)[which(D1_D2$log2FoldChange < 0)]), ]

#intersect the datasets
D1_enriched_vs_D2_FDR_005 <- intersect(rownames(D1_D2_FDR_005), rownames(D1_enriched_vs_D2))
D2_enriched_vs_D1_FDR_005 <- intersect(rownames(D1_D2_FDR_005), rownames(D2_enriched_vs_D1))

#further filtering for OTUs significantly enriched vs. soil in the respective comparision
D1_enriched_vs_D2_FDR_005_2 <- intersect(D1_enriched_vs_D2_FDR_005, Soil_Rhizosphere_D1_enriched_FDR005)
D2_enriched_vs_D1_FDR_005_2 <- intersect(D2_enriched_vs_D1_FDR_005, Soil_Rhizosphere_D2_enriched_FDR005)

#Supplementary worksheet ws 12
D1_D2_enriched_taxa<- JH04_taxa_ordered[D1_enriched_vs_D2_FDR_005_2, ]
D1_D2_enriched_info <- cbind(as.data.frame(D1_D2[D1_enriched_vs_D2_FDR_005_2, ]), D1_D2_enriched_taxa)
#write.table(D1_D2_enriched_info, file="ws10.txt")

#Supplementary worksheet ws 13
D2_D1_enriched_taxa<- JH04_taxa_ordered[D2_enriched_vs_D1_FDR_005_2, ]
D2_D1_enriched_info <- cbind(as.data.frame(D1_D2[D2_enriched_vs_D1_FDR_005_2, ]), D2_D1_enriched_taxa)
write.table(D2_D1_enriched_info, file="ws11.txt")

########################
#C1 vs. C2 (0-0)

C1_C2 <- results(JH04_cds_test , contrast = c("Geo",  "C1", "C2"))

# extract  OTUs whose adjusted p.value in a given comparison is below 0.05. 
C1_C2_FDR_005 <- C1_C2[(rownames(C1_C2)[which(C1_C2$padj <0.05)]), ]

#Identify OTUs enriched in Modern (first term of the comparison)
C1_enriched_vs_C2 <-  C1_C2[(rownames(C1_C2)[which(C1_C2$log2FoldChange > 0)]), ]

#Identify OTUs enriched in Desert (second term of the comparison)
C2_enriched_vs_C1 <-  C1_C2[(rownames(C1_C2)[which(C1_C2$log2FoldChange < 0)]), ]

#intersect the datasets
C1_enriched_vs_C2_FDR_005 <- intersect(rownames(C1_C2_FDR_005), rownames(C1_enriched_vs_C2))
C2_enriched_vs_C1_FDR_005 <- intersect(rownames(C1_C2_FDR_005), rownames(C2_enriched_vs_C1))

#further filtering for OTUs significantly enriched vs. soil in the respective comparision
C1_enriched_vs_C2_FDR_005_2 <- intersect(C1_enriched_vs_C2_FDR_005, Soil_Rhizosphere_C1_enriched_FDR005)
C2_enriched_vs_C1_FDR_005_2 <- intersect(C2_enriched_vs_C1_FDR_005, Soil_Rhizosphere_C2_enriched_FDR005)

#taxa info
C1_C2_enriched_taxa<- JH04_taxa_ordered[C1_enriched_vs_C2_FDR_005_2, ]
C1_C2_enriched_info <- cbind(as.data.frame(C1_C2[C1_enriched_vs_C2_FDR_005_2, ]), C1_C2_enriched_taxa)

#taxa info
C2_C1_enriched_taxa<- JH04_taxa_ordered[C2_enriched_vs_C1_FDR_005_2, ]
C2_C1_enriched_info <- cbind(as.data.frame(C1_C2[C2_enriched_vs_C1_FDR_005_2, ]), C2_C1_enriched_taxa)

###################################################################################
#Data visualisation using UpSetR (figure 4)
###################################################################################

#Prepare the data for UpSetR visualisation

#Merge OTUs enriched in both comparisons
D1_enriched_vs_M_FDR_005_3<-c(D1_enriched_vs_M_FDR_005_2,M_enriched_vs_D1_FDR_005_2)
D2_enriched_vs_M_FDR_005_3<-c(D2_enriched_vs_M_FDR_005_2,M_enriched_vs_D2_FDR_005_2)
C1_enriched_vs_M_FDR_005_3<-c(C1_enriched_vs_M_FDR_005_2,M_enriched_vs_C1_FDR_005_2)
C2_enriched_vs_M_FDR_005_3<-c(C2_enriched_vs_M_FDR_005_2,M_enriched_vs_C2_FDR_005_2)
N_enriched_vs_M_FDR_005_3<-c(N_enriched_vs_M_FDR_005_2,M_enriched_vs_N_FDR_005_2)

#retrieve info from otu name 
D1_enriched_vs_M_FDR_005_4<-D1_M_FDR_005[D1_enriched_vs_M_FDR_005_3,]
D2_enriched_vs_M_FDR_005_4<-D2_M_FDR_005[D2_enriched_vs_M_FDR_005_3,]
C1_enriched_vs_M_FDR_005_4<-C1_M_FDR_005[C1_enriched_vs_M_FDR_005_3,]
C2_enriched_vs_M_FDR_005_4<-C2_M_FDR_005[C2_enriched_vs_M_FDR_005_3,]
N_enriched_vs_M_FDR_005_4<-M_N_FDR_005[N_enriched_vs_M_FDR_005_3,]


#Extract baseMean (1st column)
D1_enriched_vs_M_FDR_005_4_df <- as.data.frame(D1_enriched_vs_M_FDR_005_4[ ,1])
D2_enriched_vs_M_FDR_005_4_df <- as.data.frame(D2_enriched_vs_M_FDR_005_4[ ,1])
C1_enriched_vs_M_FDR_005_4_df <- as.data.frame(C1_enriched_vs_M_FDR_005_4[ ,1])
C2_enriched_vs_M_FDR_005_4_df <- as.data.frame(C2_enriched_vs_M_FDR_005_4[ ,1])
N_enriched_vs_M_FDR_005_4_df <- as.data.frame(N_enriched_vs_M_FDR_005_4[ ,1])

#rename rows
rownames(D1_enriched_vs_M_FDR_005_4_df)<-D1_enriched_vs_M_FDR_005_3  
rownames(D2_enriched_vs_M_FDR_005_4_df)<-D2_enriched_vs_M_FDR_005_3 
rownames(C1_enriched_vs_M_FDR_005_4_df)<-C1_enriched_vs_M_FDR_005_3 
rownames(C2_enriched_vs_M_FDR_005_4_df)<-C2_enriched_vs_M_FDR_005_3 
rownames(N_enriched_vs_M_FDR_005_4_df)<-N_enriched_vs_M_FDR_005_3 

#rename columns
colnames(D1_enriched_vs_M_FDR_005_4_df) <- c("counts_D1_M")
colnames(D2_enriched_vs_M_FDR_005_4_df) <- c("counts_D2_M")
colnames(C1_enriched_vs_M_FDR_005_4_df) <- c("counts_C1_M")
colnames(C2_enriched_vs_M_FDR_005_4_df) <- c("counts_C2_M")
colnames(N_enriched_vs_M_FDR_005_4_df) <- c("counts_N_M")

D1_enriched_vs_M_FDR_005_4_df[D1_enriched_vs_M_FDR_005_4_df > 1] <- 1
D2_enriched_vs_M_FDR_005_4_df[D2_enriched_vs_M_FDR_005_4_df > 1] <- 1
C1_enriched_vs_M_FDR_005_4_df[C1_enriched_vs_M_FDR_005_4_df > 1] <- 1
C2_enriched_vs_M_FDR_005_4_df[C2_enriched_vs_M_FDR_005_4_df > 1] <- 1
N_enriched_vs_M_FDR_005_4_df[N_enriched_vs_M_FDR_005_4_df > 1] <- 1

dim(N_enriched_vs_M_FDR_005_4_df)

#define a list of unique OTUs
OTU_list <- unique(c(c((c(rownames(D1_enriched_vs_M_FDR_005_4_df), rownames(D2_enriched_vs_M_FDR_005_4_df))),
                       (c(rownames(C1_enriched_vs_M_FDR_005_4_df), rownames(C2_enriched_vs_M_FDR_005_4_df)))), 
                     (c(rownames(N_enriched_vs_M_FDR_005_4_df)))))
length(OTU_list)
#Pellet
D1_eriched_merging <- as.data.frame(D1_enriched_vs_M_FDR_005_4_df[OTU_list, ])
D2_eriched_merging <- as.data.frame(D2_enriched_vs_M_FDR_005_4_df[OTU_list, ])
C1_eriched_merging <- as.data.frame(C1_enriched_vs_M_FDR_005_4_df[OTU_list, ])
C2_eriched_merging <- as.data.frame(C2_enriched_vs_M_FDR_005_4_df[OTU_list, ])
N_eriched_merging <- as.data.frame(N_enriched_vs_M_FDR_005_4_df[OTU_list, ])


colnames(D1_eriched_merging) <- c("OTUs_D1")
colnames(D2_eriched_merging) <- c("OTUs_D2")
colnames(C1_eriched_merging) <- c("OTUs_C1")
colnames(C2_eriched_merging) <- c("OTUs_C2")
colnames(N_eriched_merging) <- c("OTUs_N")

row.names(D1_eriched_merging) <- as.vector(OTU_list)
row.names(D2_eriched_merging) <- as.vector(OTU_list)
row.names(C1_eriched_merging) <- as.vector(OTU_list)
row.names(C2_eriched_merging) <- as.vector(OTU_list)
row.names(N_eriched_merging) <- as.vector(OTU_list)


#Merge the dataset
enriched_OTUs <- cbind(D1_eriched_merging, D2_eriched_merging)
enriched_OTUs <- cbind(enriched_OTUs, C1_eriched_merging)
enriched_OTUs <- cbind(enriched_OTUs, C2_eriched_merging)
enriched_OTUs<- cbind(enriched_OTUs, N_eriched_merging)

#set NA to 0
enriched_OTUs[is.na(enriched_OTUs)] <- 0
#visualisation
upset(enriched_OTUs, sets = c("OTUs_D1", "OTUs_D2", "OTUs_C1", "OTUs_C2", "OTUs_N"), sets.bar.color = "#56B4E9",
      order.by = "freq", empty.intersections = "on")

#########################################
#Figure 5
#########################################
############################################
#Preparation for the phylogenetic tree iTOL
#create a phylogenetic tree for visualisation of rhizosphere and root enriched
JH04_data_phyloseq_4 <- prune_taxa(unique(union(M_enriched_vs_D1_FDR_005_2, M_enriched_vs_D2_FDR_005_2)), JH04_data_phyloseq_3)
JH04_data_phyloseq_4
#phylogenetic tree
tree_Modern_enriched = phy_tree(JH04_data_phyloseq_4) 
#ape::write.tree(tree_Modern_enriched, "JH04_Modern_enriched_Desert.tree")
#create the annotation files
#Taxonomic annotation
tree_Modern_enriched_taxonomic_annotation <- JH04_taxa_ordered[unique(union(M_enriched_vs_D1_FDR_005_2, M_enriched_vs_D2_FDR_005_2)), ]
#write.table(tree_Modern_enriched_taxonomic_annotation, file="tree_Modern_enriched_taxonomic_annotation.txt", sep = "\t")
#Comparison annotation
#Modern vs Desert1
Modern_eriched_Desert1 <- as.data.frame(D1_M_FDR_005[M_enriched_vs_D1_FDR_005_2, 1]) 
rownames(Modern_eriched_Desert1) <- M_enriched_vs_D1_FDR_005_2
colnames(Modern_eriched_Desert1) <- c("counts_MD1")
Modern_eriched_Desert1[Modern_eriched_Desert1 > 1] <- 1
Modern_eriched_Desert1
dim(Modern_eriched_Desert1)
#Modern vs Desert2
Modern_eriched_Desert2 <- as.data.frame(D2_M_FDR_005[M_enriched_vs_D2_FDR_005_2, 1]) 
rownames(Modern_eriched_Desert2) <- M_enriched_vs_D2_FDR_005_2
colnames(Modern_eriched_Desert2) <- c("counts_MD2")
Modern_eriched_Desert2[Modern_eriched_Desert2 > 1] <- 1
Modern_eriched_Desert2
dim(Modern_eriched_Desert2)

#combine the datasets: note they have unequal values
OTU_list_Modern_enriched <- unique(union(M_enriched_vs_D1_FDR_005_2, M_enriched_vs_D2_FDR_005_2))
length(OTU_list_Modern_enriched)
#Modern-Desert1
MD1_eriched_merging <- as.data.frame(Modern_eriched_Desert1[OTU_list_Modern_enriched , ])
colnames(MD1_eriched_merging) <- c("OTUs_MD1")
row.names(MD1_eriched_merging) <- as.vector(OTU_list_Modern_enriched)
#Modern-Desert2
MD2_eriched_merging <- as.data.frame(Modern_eriched_Desert2[OTU_list_Modern_enriched , ])
colnames(MD2_eriched_merging) <- c("OTUs_MD2")
row.names(MD2_eriched_merging) <- as.vector(OTU_list_Modern_enriched)

#Merge the dataset
Modern_enriched_OTUs <- cbind(MD1_eriched_merging, MD2_eriched_merging)
#set NA to 0
Modern_enriched_OTUs[is.na(Modern_enriched_OTUs)] <- 0
#write.table(Modern_enriched_OTUs, file="tree_Modern_comparison_annotation_1612.txt", sep = "\t")

#end of the code