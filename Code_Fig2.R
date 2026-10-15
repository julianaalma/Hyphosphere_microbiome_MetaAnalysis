###Figs 2 and S1: Hyphosphere Bacterial communities composition 
###Diversity analysis### 
###libraries needed####
library(phyloseq)
library(ggplot2)
library(tidyverse)
library(readxl)
library(dplyr)    
library(tibble)
library(vegan)
library(ade4)
library(devtools)
library(RSQLite)    
library(RColorBrewer)
library(rstatix)
library(ggpubr)
library(car)
library(psych)
library(corrplot)
library(MASS)
library(agricolae)
library(stats)
library(RColorBrewer)
library(corrplot)
library(psych)
library(cowplot)
library(stringr)
library(adegraphics)
library(DescTools)
#library(metagMisc)
library(ggpubr)
#install.packages("RColorBrewer")
library(RColorBrewer)
install.packages(ggboxplot)
library(ggboxplot)

#set wd
setwd("C:/Users/Lucie/Desktop/Thèse/Revue/Metaanalyse/Analyse_R")

#Importation of excel and mothur files 
samples_df<- read_excel("Metadata_meta_analyse.xlsx")
sharedfile = "combo.list.trim.contigs.good.cuta_16SFR_RV_RV.good.unique.denovo.vsearch.pick.dgc.0.01.abund.0.01.pick.shared"
otu_mat <- import_mothur(mothur_shared_file = sharedfile)
tax_mat <- read_excel("Table_Taxonomy.xlsx")

#Remove the bracket from taxon's names
tax_mat <- as.data.frame(lapply(tax_mat, function(x) gsub("\\(\\d+\\)", "", x)))
head(tax_mat)


samples_df <- samples_df %>% 
  tibble::column_to_rownames("Library_name") 
tax_mat <- tax_mat %>% 
  tibble::column_to_rownames("OTU") 
otu_mat <- as.matrix(otu_mat)
tax_mat <- as.matrix(tax_mat)

colnames(tax_mat)
OTU = otu_table(otu_mat, taxa_are_rows = TRUE)
TAX = tax_table(tax_mat)
samples = sample_data(samples_df)

phy_obj_plant<- phyloseq(OTU, TAX, samples)
phy_obj_plant<- prune_species(speciesSums(phy_obj_plant) > 0, phy_obj_plant)


#transform data in log(RA+1)
##relative abundance for each OTU abundOTU/abundtotOTU
phy_obj_plant_ra<-transform_sample_counts(phy_obj_plant, function(x) x/sum(x) )
phy_obj_plant_ra<-otu_table(phy_obj_plant_ra)+1
phy_obj_plant_log<-log10(otu_table(phy_obj_plant_ra))
###reassemble otu_table with taxonomy and metadata
phy_obj_plant_log<- phyloseq(otu_table(phy_obj_plant_log), TAX, samples)

#############################################################################################################################

#### Beta diversity plotting #####
##### Stackplot, relative abundance taxa at "Phylum" level ########
y_phylum <- tax_glom(phy_obj_plant, taxrank = "Phylum") # agglomerate taxa
y_ra_phylum<-transform_sample_counts(y_phylum, function(x) x/sum(x) )
y_ra_phylum <- phyloseq(otu_table(y_ra_phylum), tax_table(y_phylum), samples)
y4_phylum <- psmelt(y_ra_phylum) # create dataframe from phyloseq object
y4_phylum$Phylum <- as.character(y4_phylum$Phylum) #convert to character
y4_phylum$Abundance<-as.numeric(y4_phylum$Abundance)

#Calculates the total abundance of each phylum, on all samples
phylum_order <- y4_phylum %>%
  group_by(Phylum) %>%
  summarise(total_abundance = sum(Abundance)) %>%
  arrange(desc(total_abundance)) %>%
  pull(Phylum)

#Applies this order to the factor
y4_phylum$Phylum <- factor(y4_phylum$Phylum, levels = phylum_order)

#Abbreviation and order of the factors Fungi and Soil_detail
y4_phylum$abbrevF<- y4_phylum$Fungi
y4_phylum$abbrevF<-dplyr::recode_factor(y4_phylum$abbrevF,Fusarium_oxysprorum="F.o",Penicillium_canescens="P.c",
                                 Termitomyces_albuminosus="T.a",Chlorophyllum_molybdites="C.m", Agaricus_flocculosipes= "A.f",
                                 Amanita_pantherina="A.p",Tylopilus_felleus="T.f",Suillus_placidus="S.p",Rhizophagus_intraradices="R.i")
y4_phylum$abbrevF <- as.character(y4_phylum$abbrevF)
y4_phylum<-rstatix::reorder_levels(y4_phylum,abbrevF,order=c("F.o","P.c","T.a","C.m","A.f","A.p","T.f","S.p","R.i"))
table(y4_phylum$abbrevF, useNA = "ifany")

y4_phylum$abbrevS<- y4_phylum$Soil_detail
y4_phylum$abbrevS<-dplyr::recode_factor(y4_phylum$abbrevS,Artificial="Artificial",Field_M1P1="FM1P1",
                                 Field_N1K1="FN1K1",	
                                 Field_P1N2K2="FP1N2K2", Forest_soil= "Forest",
                                 Field_high_P="F.High.P",Field_arid_environment="F.arid",Field_semiarid_environment="F.semiarid",Field_wet_environment="F.wet")
y4_phylum$abbrevS <- as.character(y4_phylum$abbrevS)
y4_phylum<-rstatix::reorder_levels(y4_phylum,abbrevS,order=c("Artificial","FM1P1","FN1K1","FP1N2K2","Forest","F.High.P","F.arid","F.semiarid","F.wet"))
table(y4_phylum$abbrevS, useNA = "ifany")

#Color of the different phyla in the histogram
library(randomcoloR)
phylum_unique <- unique(y4_phylum$Phylum)
cols <- distinctColorPalette(length(phylum_unique))
names(cols) <- phylum_unique
cols["Pseudomonadota"] <- "green"
cols["Actinomycetota"] <- "red"
cols["Acidobacteriota"] <- "orange"
cols["Chloroflexota"] <- "blue"
cols["Myxococcota"] <- "yellow"


library(ggh4x)

strip_background <- strip_nested(text_x = elem_list_text(face = "bold"),
                                 background_x = elem_list_rect(fill = c("red","olivedrab3", "orange","orange","orange","orange","orange","orange", "gray58",
                                                                        "gold", "slateblue4", "lightskyblue1",  "plum", "seagreen4", "seagreen4","seagreen4",
                                                                        "seagreen4", "seagreen4", "seagreen4","seagreen1",  "red2", "palevioletred1", "rosybrown1")))

#Histogram at the level of the phylum
p_phylum <- ggplot(y4_phylum, aes(x = Sample, y = Abundance, fill = Phylum)) +
  geom_bar(stat = "identity", position = "fill", width = 0.85) +
  scale_fill_manual(values = cols) +
  labs(x = NULL, y = "Relative abundance", fill = "Bacterial phylum") +
  facet_nested(~ abbrevF + abbrevS, scales = "free_x", space = "free_x", strip = strip_background) +
  guides(fill = guide_legend(nrow = 4, byrow = TRUE)) +
  theme(legend.position = "bottom", legend.direction = "horizontal", legend.title = element_blank(), legend.text = element_text(size = 22), legend.key.size = unit(0.3, "cm"), axis.text.x = element_blank(),element_text(
    angle = 90,
    vjust = 0.5,
    hjust = 1), axis.text.y = element_text(size = 9))
p_phylum


ggsave(plot = p_phylum, filename="C:/Users/Lucie/Desktop/Thèse/Revue/Metaanalyse/Analyse_R/RA_Hist_Phylum.pdf", width=32, height=10)
ggsave(plot = p_phylum, filename="C:/Users/Lucie/Desktop/Thèse/Revue/Metaanalyse/Analyse_R/RA_Hist_Phylum.png", width=32, height=10)

####Order####
y_order <- tax_glom(phy_obj_plant, taxrank = "Order") # agglomerate taxa
y_ra_order<-transform_sample_counts(y_order, function(x) x/sum(x) )
y_ra_order<- phyloseq(otu_table(y_ra_order), tax_table(y_order), samples)
y4_order <- psmelt(y_ra_order) # create dataframe from phyloseq object
y4_order$Order <- as.character(y4_order$Order) #convert to character
y4_order$Abundance<-as.numeric(y4_order$Abundance)

#Calculates the total abundance of each orders, on all samples
order_order <- y4_order %>%
  group_by(Order) %>%
  summarise(total_abundance = sum(Abundance)) %>%
  arrange(desc(total_abundance)) %>%
  pull(Order)

#Applies this order to the factor
y4_order$Order <- factor(y4_order$Order, levels = order_order)

#Keep the 30 most abundant orders in the legend  
top30_order <- y4_order %>%
  group_by(Order) %>%
  summarise(total_abundance = sum(Abundance)) %>%
  arrange(desc(total_abundance)) %>%
  slice_head(n = 30) %>%
  pull(Order)
head(top30_order)
#Add break = top30_order in scall_fill_manual for the histogram

#Abbreviation and order of the factors Fungi and Soil_detail
y4_order$abbrevF_order<- y4_order$Fungi
y4_order$abbrevF_order<-dplyr::recode_factor(y4_order$abbrevF,Fusarium_oxysprorum="F.o",Penicillium_canescens="P.c",
                                             Termitomyces_albuminosus="T.a",Chlorophyllum_molybdites="C.m", Agaricus_flocculosipes= "A.f",
                                             Amanita_pantherina="A.p",Tylopilus_felleus="T.f",Suillus_placidus="S.p",Rhizophagus_intraradices="R.i")
y4_order$abbrevF_order <- as.character(y4_order$abbrevF_order)
y4_order<-rstatix::reorder_levels(y4_order,abbrevF_order,order=c("F.o","P.c","T.a","C.m","A.f","A.p","T.f","S.p","R.i"))
table(y4_order$abbrevF_order, useNA = "ifany")

y4_order$abbrevS_order<- y4_order$Soil_detail
y4_order$abbrevS_order<-dplyr::recode_factor(y4_order$abbrevS_order,Artificial="Artificial",Field_M1P1="FM1P1",
                                             Field_N1K1="FN1K1",	
                                             Field_P1N2K2="FP1N2K2", Forest_soil= "Forest",
                                             Field_high_P="F.High.P",Field_arid_environment="F.arid",Field_semiarid_environment="F.semiarid",Field_wet_environment="F.wet")
y4_order$abbrevS_order <- as.character(y4_order$abbrevS_order)
y4_order<-rstatix::reorder_levels(y4_order,abbrevS_order,order=c("Artificial","FM1P1","FN1K1","FP1N2K2","Forest","F.High.P","F.arid","F.semiarid","F.wet"))
table(y4_order$abbrevS_order, useNA = "ifany")


#Color of the different orders in the histogram
library(randomcoloR)
order_unique <- unique(y4_order$Order)
cols <- distinctColorPalette(length(order_unique))
names(cols) <- order_unique


library(ggh4x)

strip_background <- strip_nested(text_x = elem_list_text(face = "bold"),
                                 background_x = elem_list_rect(fill = c("red","olivedrab3", "orange","orange","orange","orange","orange","orange", "gray58",
                                                                        "gold", "slateblue4", "lightskyblue1",  "plum", "seagreen4", "seagreen4","seagreen4",
                                                                        "seagreen4", "seagreen4", "seagreen4","seagreen1",  "red2", "palevioletred1", "rosybrown1")))


#Histogram at the level of the order
p_order <- ggplot(y4_order, aes(x = Sample, y = Abundance, fill = Order)) +
  geom_bar(stat = "identity", position = "fill", width = 0.85) +
  scale_fill_manual(values = cols, breaks = top30_order) +
  labs(x = NULL, y = "Relative abundance", fill = "Bacterial orders") +
  facet_nested(~ abbrevF_order + abbrevS_order, scales = "free_x", space = "free_x", strip = strip_background) +
  guides(fill = guide_legend(nrow = 4, byrow = TRUE)) +
  theme(legend.position = "bottom", legend.direction = "horizontal", legend.title = element_blank(), legend.text = element_text(size = 22), legend.key.size = unit(0.3, "cm"), axis.text.x = element_blank(),element_text(
    angle = 90,
    vjust = 0.5,
    hjust = 1), axis.text.y = element_text(size = 9))
p_order

ggsave(plot = p_order, filename="C:/Users/Lucie/Desktop/Thèse/Revue/Metaanalyse/Analyse_R/RA_Hist_Order.pdf", width=32, height=10)
ggsave(plot = p_order, filename="C:/Users/Lucie/Desktop/Thèse/Revue/Metaanalyse/Analyse_R/RA_Hist_Order.png", width=32, height=10)


####Family####
y_family <- tax_glom(phy_obj_plant, taxrank = "Family") # agglomerate taxa
y_ra_family<-transform_sample_counts(y_family, function(x) x/sum(x) )
y_ra_family<- phyloseq(otu_table(y_ra_family), tax_table(y_family), samples)
y4_family <- psmelt(y_ra_family) # create dataframe from phyloseq object
y4_family$Family <- as.character(y4_family$Family) #convert to character
y4_family$Abundance<-as.numeric(y4_family$Abundance)

#Calculates the total abundance of each family, on all samples
family_order <- y4_family %>%
  group_by(Family) %>%
  summarise(total_abundance = sum(Abundance)) %>%
  arrange(desc(total_abundance)) %>%
  pull(Family)

#Applies this order to the factor
y4_family$Family <- factor(y4_family$Family, levels = family_order)

#Keep the 30 most abundant families in the legend 
top30_families <- y4_family %>%
  group_by(Family) %>%
  summarise(total_abundance = sum(Abundance)) %>%
  arrange(desc(total_abundance)) %>%
  slice_head(n = 30) %>%
  pull(Family)
head(top30_families)
#Add break = top30_families in scall_fill_manual for the histogram


#Abbreviation and order of the factors Fungi and Soil_detail
y4_family$abbrevF_family<- y4_family$Fungi
y4_family$abbrevF_family<-dplyr::recode_factor(y4_family$abbrevF_family,Fusarium_oxysprorum="F.o",Penicillium_canescens="P.c",
                                               Termitomyces_albuminosus="T.a",Chlorophyllum_molybdites="C.m", Agaricus_flocculosipes= "A.f",
                                               Amanita_pantherina="A.p",Tylopilus_felleus="T.f",Suillus_placidus="S.p",Rhizophagus_intraradices="R.i")
y4_family$abbrevF_family <- as.character(y4_family$abbrevF_family)
y4_family<-rstatix::reorder_levels(y4_family,abbrevF_family,order=c("F.o","P.c","T.a","C.m","A.f","A.p","T.f","S.p","R.i"))
table(y4_family$abbrevF_family, useNA = "ifany")

y4_family$abbrevS_family<- y4_family$Soil_detail
y4_family$abbrevS_family<-dplyr::recode_factor(y4_family$abbrevS_family,Artificial="Artificial",Field_M1P1="FM1P1",
                                               Field_N1K1="FN1K1",	
                                               Field_P1N2K2="FP1N2K2", Forest_soil= "Forest",
                                               Field_high_P="F.High.P",Field_arid_environment="F.arid",Field_semiarid_environment="F.semiarid",Field_wet_environment="F.wet")
y4_family$abbrevS_family <- as.character(y4_family$abbrevS_family)
y4_family<-rstatix::reorder_levels(y4_family,abbrevS_family,order=c("Artificial","FM1P1","FN1K1","FP1N2K2","Forest","F.High.P","F.arid","F.semiarid","F.wet"))
table(y4_family$abbrevS_family, useNA = "ifany")

#Color of the different families in the histogram
library(randomcoloR)
family_unique <- unique(y4_family$Family)
cols <- distinctColorPalette(length(family_unique))
names(cols) <- family_unique
cols["Burkholderiaceae"] <- "red"
cols["Oxalobacteraceae"] <- "royalblue1"
cols["Comamonadaceae"] <- "orange"
cols["Pseudomonadaceae"] <- "green"
cols["Propionibacteriaceae"] <- "yellow"
cols["Microccocaceae"] <- "deeppink1"


library(ggh4x)

strip_background <- strip_nested(text_x = elem_list_text(face = "bold"),
                                 background_x = elem_list_rect(fill = c("red","olivedrab3", "orange","orange","orange","orange","orange","orange", "gray58",
                                                                        "gold", "slateblue4", "lightskyblue1",  "plum", "seagreen4", "seagreen4","seagreen4",
                                                                        "seagreen4", "seagreen4", "seagreen4","seagreen1",  "red2", "palevioletred1", "rosybrown1")))


##Histogram at the level of the family with 30 most abundant families in legend 
p_family <- ggplot(y4_family, aes(x = Sample, y = Abundance, fill = Family)) +
  geom_bar(stat = "identity", position = "fill", width = 0.85) +
  scale_fill_manual(values = cols, breaks = top30_families) +
  labs(x = NULL, y = "Relative abundance", fill = NULL) +
  facet_nested(~ abbrevF_family + abbrevS_family, scales = "free_x", space = "free_x", strip = strip_background) +
  guides(fill = guide_legend(nrow = 4, byrow = TRUE)) +
  theme(legend.position = "bottom", legend.direction = "horizontal", legend.title = element_blank(), legend.text = element_text(size = 22), legend.key.size = unit(0.3, "cm"), legend.key.spacing.x = unit(0.4, "cm"),axis.text.x = element_blank(),element_text(
    angle = 90,
    vjust = 0.5,
    hjust = 1), axis.text.y = element_text(size = 9))
p_family

ggsave(plot = p_family, filename="C:/Users/Lucie/Desktop/Thèse/Revue/Metaanalyse/Analyse_R/RA_Hist_Family.png", width=32, height=10)
ggsave(plot = p_family, filename="C:/Users/Lucie/Desktop/Thèse/Revue/Metaanalyse/Analyse_R/RA_Hist_Family.pdf", width=32, height=10)

##Histogram at the level of the family with the entire legend 
p_family_legend <- ggplot(y4_family, aes(x = Sample, y = Abundance, fill = Family)) +
  geom_bar(stat = "identity", position = "fill", width = 0.85) +
  scale_fill_manual(values = cols) +
  labs(x = NULL, y = "Relative abundance", fill = NULL) +
  facet_nested(~ abbrevF_family + abbrevS_family, scales = "free_x", space = "free_x", strip = strip_background) +
  guides(fill = guide_legend(ncol = 5, byrow = TRUE)) +
  theme(legend.position = "bottom", legend.direction = "horizontal", legend.title = element_blank(), legend.text = element_text(size = 20), legend.key.size = unit(0.3, "cm"),axis.text.x = element_blank(),element_text(
    angle = 90,
    vjust = 0.5,
    hjust = 1), axis.text.y = element_text(size = 9))
p_family_legend

ggsave(plot = p_family_legend, filename="C:/Users/Lucie/Desktop/Thèse/Revue/Metaanalyse/Analyse_R/RA_Hist_Family_Legend.png", width=35, height=35)
ggsave(plot = p_family_legend, filename="C:/Users/Lucie/Desktop/Thèse/Revue/Metaanalyse/Analyse_R/RA_Hist_Family_Legend.pdf", width=35, height=35)

####Genus####
y_genus <- tax_glom(phy_obj_plant, taxrank = "Genus") # agglomerate taxa
y_ra_genus<-transform_sample_counts(y_genus, function(x) x/sum(x) )
y_ra_genus<- phyloseq(otu_table(y_ra_genus), tax_table(y_genus), samples)
y4_genus <- psmelt(y_ra_genus) # create dataframe from phyloseq object
y4_genus$Genus <- as.character(y4_genus$Genus) #convert to character
y4_genus$Abundance<-as.numeric(y4_genus$Abundance)

#Calculates the total abundance of each genera, on all samples
genus_order <- y4_genus %>%
  group_by(Genus) %>%
  summarise(total_abundance = sum(Abundance)) %>%
  arrange(desc(total_abundance)) %>%
  pull(Genus)

#Applies this order to the factor
y4_genus$Genus <- factor(y4_genus$Genus, levels = genus_order)

#Keep the 30 most abundant genera in the legend 
top30_genera <- y4_genus %>%
  group_by(Genus) %>%
  summarise(total_abundance = sum(Abundance)) %>%
  arrange(desc(total_abundance)) %>%
  slice_head(n = 30) %>%
  pull(Genus)
head(top30_genera)
#Add break = top30_genera in scall_fill_manual for the histogram

#Abbreviation and order of the factors Fungi and Soil_detail
y4_genus$abbrevF_genus<- y4_genus$Fungi
y4_genus$abbrevF_genus<-dplyr::recode_factor(y4_genus$abbrevF_genus,Fusarium_oxysprorum="F.o",Penicillium_canescens="P.c",
                                             Termitomyces_albuminosus="T.a",Chlorophyllum_molybdites="C.m", Agaricus_flocculosipes= "A.f",
                                             Amanita_pantherina="A.p",Tylopilus_felleus="T.f",Suillus_placidus="S.p",Rhizophagus_intraradices="R.i")
y4_genus$abbrevF_genus <- as.character(y4_genus$abbrevF_genus)
y4_genus<-rstatix::reorder_levels(y4_genus,abbrevF_genus,order=c("F.o","P.c","T.a","C.m","A.f","A.p","T.f","S.p","R.i"))
table(y4_genus$abbrevF_genus, useNA = "ifany")

y4_genus$abbrevS_genus<- y4_genus$Soil_detail
y4_genus$abbrevS_genus<-dplyr::recode_factor(y4_genus$abbrevS_genus,Artificial="Artificial",Field_M1P1="FM1P1",
                                             Field_N1K1="FN1K1",	
                                             Field_P1N2K2="FP1N2K2", Forest_soil= "Forest",
                                             Field_high_P="F.High.P",Field_arid_environment="F.arid",Field_semiarid_environment="F.semiarid",Field_wet_environment="F.wet")
y4_genus$abbrevS_genus <- as.character(y4_genus$abbrevS_genus)
y4_genus<-rstatix::reorder_levels(y4_genus,abbrevS_genus,order=c("Artificial","FM1P1","FN1K1","FP1N2K2","Forest","F.High.P","F.arid","F.semiarid","F.wet"))
table(y4_genus$abbrevS_genus, useNA = "ifany")

#Color of the different genera in the histogram
library(randomcoloR)
genus_unique <- unique(y4_genus$Genus)
cols <- distinctColorPalette(length(genus_unique))
names(cols) <- genus_unique

library(ggh4x)

strip_background <- strip_nested(text_x = elem_list_text(face = "bold"),
                                 background_x = elem_list_rect(fill = c("red","olivedrab3", "orange","orange","orange","orange","orange","orange", "gray58",
                                                                        "gold", "slateblue4", "lightskyblue1",  "plum", "seagreen4", "seagreen4","seagreen4",
                                                                        "seagreen4", "seagreen4", "seagreen4","seagreen1",  "red2", "palevioletred1", "rosybrown1")))


#Histogram at the level of the genus
p_genus <- ggplot(y4_genus, aes(x = Sample, y = Abundance, fill = Genus)) +
  geom_bar(stat = "identity", position = "fill", width = 0.85) +
  scale_fill_manual(values = cols, breaks = top30_genera) +
  labs(x = NULL, y = "Relative abundance", fill = "Bacterial genus") +
  facet_nested(~ abbrevF_genus + abbrevS_genus, scales = "free_x", space = "free_x", strip = strip_background) +
  guides(fill = guide_legend(nrow = 4, byrow = TRUE)) +
  theme(legend.position = "bottom", legend.direction = "horizontal", legend.title = element_blank(), legend.text = element_text(size = 22), legend.key.size = unit(0.3, "cm"), axis.text.x = element_blank(),element_text(
    angle = 90,
    vjust = 0.5,
    hjust = 1), axis.text.y = element_text(size = 9))
p_genus

ggsave(plot = p_genus, filename="C:/Users/Lucie/Desktop/Thèse/Revue/Metaanalyse/Analyse_R/RA_Hist_Genus.pdf", width=32, height=10)
ggsave(plot = p_genus, filename="C:/Users/Lucie/Desktop/Thèse/Revue/Metaanalyse/Analyse_R/RA_Hist_Genus.png", width=32, height=10)

###Species###
y_species <- tax_glom(phy_obj_plant, taxrank = "Species") # agglomerate taxa
y_ra_species<-transform_sample_counts(y_species, function(x) x/sum(x) )
y_ra_species<- phyloseq(otu_table(y_ra_species), tax_table(y_species), samples)
y4_species <- psmelt(y_ra_species) # create dataframe from phyloseq object
y4_species$Species <- as.character(y4_species$Species) #convert to character
y4_species$Abundance<-as.numeric(y4_species$Abundance)

#Calculates the total abundance of each species, on all samples
species_order <- y4_species %>%
  group_by(Species) %>%
  summarise(total_abundance = sum(Abundance)) %>%
  arrange(desc(total_abundance)) %>%
  pull(Species)

#Applies this order to the factor
y4_species$Species <- factor(y4_species$Species, levels = species_order)

#Keep the 30 most abundant species in the legend 
top30_species <- y4_species %>%
  group_by(Species) %>%
  summarise(total_abundance = sum(Abundance)) %>%
  arrange(desc(total_abundance)) %>%
  slice_head(n = 30) %>%
  pull(Species)
head(top30_species)
#Add break = top30_species in scall_fill_manual for the histogram

#Abbreviation and order of the factors Fungi and Soil_detail
y4_species$abbrevF_species<- y4_species$Fungi
y4_species$abbrevF_species<-dplyr::recode_factor(y4_species$abbrevF_species,Fusarium_oxysprorum="F.o",Penicillium_canescens="P.c",
                                                 Termitomyces_albuminosus="T.a",Chlorophyllum_molybdites="C.m", Agaricus_flocculosipes= "A.f",
                                                 Amanita_pantherina="A.p",Tylopilus_felleus="T.f",Suillus_placidus="S.p",Rhizophagus_intraradices="R.i")
y4_species$abbrevF_species <- as.character(y4_species$abbrevF_species)
y4_species<-rstatix::reorder_levels(y4_species,abbrevF_species,order=c("F.o","P.c","T.a","C.m","A.f","A.p","T.f","S.p","R.i"))
table(y4_species$abbrevF_species, useNA = "ifany")

y4_species$abbrevS_species<- y4_species$Soil_detail
y4_species$abbrevS_species<-dplyr::recode_factor(y4_species$abbrevS_species,Artificial="Artificial",Field_M1P1="FM1P1",
                                                 Field_N1K1="FN1K1",	
                                                 Field_P1N2K2="FP1N2K2", Forest_soil= "Forest",
                                                 Field_high_P="F.High.P",Field_arid_environment="F.arid",Field_semiarid_environment="F.semiarid",Field_wet_environment="F.wet")
y4_species$abbrevS_species <- as.character(y4_species$abbrevS_species)
y4_species<-rstatix::reorder_levels(y4_species,abbrevS_species,order=c("Artificial","FM1P1","FN1K1","FP1N2K2","Forest","F.High.P","F.arid","F.semiarid","F.wet"))
table(y4_species$abbrevS_species, useNA = "ifany")

#Color of the different species in the histogram
library(randomcoloR)
species_unique <- unique(y4_species$Species)
cols <- distinctColorPalette(length(species_unique))
names(cols) <- species_unique

library(ggh4x)

strip_background <- strip_nested(text_x = elem_list_text(face = "bold"),
                                 background_x = elem_list_rect(fill = c("red","olivedrab3", "orange","orange","orange","orange","orange","orange", "gray58",
                                                                        "gold", "slateblue4", "lightskyblue1",  "plum", "seagreen4", "seagreen4","seagreen4",
                                                                        "seagreen4", "seagreen4", "seagreen4","seagreen1",  "red2", "palevioletred1", "rosybrown1")))


#Histogram at the level of the species
p_species <- ggplot(y4_species, aes(x = Sample, y = Abundance, fill = Species)) +
  geom_bar(stat = "identity", position = "fill", width = 0.85) +
  scale_fill_manual(values = cols, breaks = top30_species) +
  labs(x = NULL, y = "Relative abundance", fill = "Bacterial species") +
  facet_nested(~ abbrevF_species + abbrevS_species, scales = "free_x", space = "free_x", strip = strip_background) +
  guides(fill = guide_legend(nrow = 4, byrow = TRUE)) +
  theme(legend.position = "bottom", legend.direction = "horizontal", legend.title = element_blank(), legend.text = element_text(size = 22), legend.key.size = unit(0.3, "cm"), axis.text.x = element_blank(),element_text(
    angle = 90,
    vjust = 0.5,
    hjust = 1), axis.text.y = element_text(size = 9))
p_species

ggsave(plot = p_species, filename="C:/Users/Lucie/Desktop/Thèse/Revue/Metaanalyse/Analyse_R/RA_Hist_Species.pdf", width=32, height=10)
ggsave(plot = p_species, filename="C:/Users/Lucie/Desktop/Thèse/Revue/Metaanalyse/Analyse_R/RA_Hist_Species.png", width=32, height=10)

#############################################################################################################################

####Phylum binary table####
#Create a factor corresponding to the Phylum
gentab_phylum = factor(tax_table(phy_obj_plant)[, "Phylum"])
#Tabulate the counts for each phyla in each sample
gentab_phylum = apply(otu_table(phy_obj_plant), MARGIN = 2, function(x) {
  tapply(x, INDEX = gentab_phylum, FUN = sum, na.rm = TRUE, simplify = TRUE)
})
gentab_phylum <- as.data.frame(gentab_phylum)
otu_binary_phylum <- ifelse(gentab_phylum > 0, 1, 0)
otu_binary_phylum <- as.data.frame(otu_binary_phylum)
write.csv2(otu_binary_phylum,"Binary_table_phylum2.csv")

####Order binary table####
#Create a factor corresponding to the Phylum
gentab_order = factor(tax_table(phy_obj_plant)[, "Order"])
#Tabulate the counts for each phyla in each sample
gentab_order = apply(otu_table(phy_obj_plant), MARGIN = 2, function(x) {
  tapply(x, INDEX = gentab_order, FUN = sum, na.rm = TRUE, simplify = TRUE)
})
gentab_order <- as.data.frame(gentab_order)
otu_binary_order <- ifelse(gentab_order > 0, 1, 0)
otu_binary_order <- as.data.frame(otu_binary_order)
write.csv2(otu_binary_order,"Binary_table_order.csv")

####Family binary table####
# Create a factor corresponding to the Family
gentab_family = factor(tax_table(phy_obj_plant)[, "Family"])
# Tabulate the counts for each families in each sample
gentab_family = apply(otu_table(phy_obj_plant), MARGIN = 2, function(x) {
  tapply(x, INDEX = gentab_family, FUN = sum, na.rm = TRUE, simplify = TRUE)
})
gentab_family <- as.data.frame(gentab_family)
otu_binary_family <- ifelse(gentab_family > 0, 1, 0)
otu_binary_family<- as.data.frame(otu_binary_family)
write.csv2(otu_binary_family,"Binary_table_family.csv")

####Genus binary table####
gentab_genus <- factor(tax_table(phy_obj_plant)[, "Genus"])
gentab_genus <- apply(otu_table(phy_obj_plant), 2, function(x) {
  res <- tapply(x,INDEX = gentab_genus,FUN = sum,na.rm = TRUE)
  # keep all the levels
  res <- res[levels(gentab_genus)]
  res[is.na(res)] <- 0
  return(res)})
gentab_genus <- as.data.frame(gentab_genus)

otu_binary_genus <- ifelse(gentab_genus > 0, 1, 0)
otu_binary_genus <- as.data.frame(otu_binary_genus)
write.csv2(otu_binary_genus,"Binary_table_genus.csv")

####Species binary table####
gentab_species <- factor(tax_table(phy_obj_plant)[, "Species"])
gentab_species <- apply(otu_table(phy_obj_plant), 2, function(x) {
  res <- tapply(x,INDEX = genfac_species,FUN = sum,na.rm = TRUE)
  # keep all the levels
  res <- res[levels(gentab_species)]
  res[is.na(res)] <- 0
  return(res)})
gentab_species <- as.data.frame(gentab_species)

otu_binary_species <- ifelse(gentab_species > 0, 1, 0)
otu_binary_species <- as.data.frame(otu_binary_species)
write.csv2(otu_binary_species,"Binary_table_species.csv")


###OTU binary table###
otu_binary <- ifelse(otu_table(phy_obj_plant) > 0, 1, 0)
otu_binary <- as.data.frame(otu_binary)
write.csv2(otu_binary, "Binary_table_OTU.csv")
