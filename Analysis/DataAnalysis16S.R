######################################
#Data Analysis Buried Microbes Project
#####################################
#import phyloseq
library(phyloseq)
packageVersion("phyloseq")
library(ggplot2)
packageVersion("ggplot2")
library(vegan)
library(dplyr)
library(reshape2)
library(tidyverse)
library(ggpubr)
library(ggrepel)
library(viridis)
library(qiime2R)
library(reshape2)
library(plotly)


getwd()
setwd("/Fecal_compare/exported-feature-table")


###Make physeq object
FSG<-qza_to_phyloseq(
  feature= "Featuretable.qza",
  metadata = "map.txt",
  taxonomy = "Taxonomy.qza"
)

#rooted tree is not merging correctly; maybe because this is a fragment insertion tree; just build unifrac matrices outside of R.
tree_b<-read_qza("sepptree.qza")$data
tree<-phy_tree(tree2)
#merge taxonomy table with exisitng phyloseq object
FSGt <- merge_phyloseq(otu_table(FSG), sample_data(FSG),phy_tree(tree2))
FSGt<- merge_phyloseq(FSGt, tax_table(FSG))


#View physeq object
FSG
#change NAs to Unlcassified to avoid taxa from dropping in the future
tax_table(FSG)[is.na(tax_table(FSG))]<-"Unclassified"

summary(sample_data(FSG))
#change column names in taxonomy table from rank1, rank2, etc
colnames(tax_table(FSG))<-c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")

#need sample summaries
FSGGrave<-subset_samples(FSG, Project=="Grave")

#########################
##look at taxa prevalence and filter rare taxa; save for later
#########################
#Compute prevalence of each feature, store as data.frame
prevdf = apply(X = otu_table(FSG),
               MARGIN = ifelse(taxa_are_rows(FSG), yes = 1, no = 2),
               FUN = function(x){sum(x > 0)})
# Add taxonomy and total read counts to this data.frame
prevdf = data.frame(Prevalence = prevdf,
                    TotalAbundance = taxa_sums(FSG),
                    tax_table(FSG))
# Subset to the remaining phyla
library(ggplot2)
prevdf1 = subset(prevdf, Phylum %in% get_taxa_unique(FSG, "Phylum"))
ggplot(prevdf1, aes(TotalAbundance, Prevalence / nsamples(FSG),color=Phylum)) +
  # Include a guess for parameter
  geom_hline(yintercept = 0.01, alpha = 0.5, linetype = 2) +  geom_point(size = 2, alpha = 0.7) +
  scale_x_log10() +  xlab("Total Abundance") + ylab("Prevalence [Frac. Samples]") +
  facet_wrap(~Phylum) + theme(legend.position="none")

#get rid of phyla not seen across 0.5% of taxa; 2 samples
prevalenceThreshold = 0.005 * nsamples(FSG) 
prevalenceThreshold
# Execute prevalence filter, using `prune_taxa()` function
keepTaxa = rownames(prevdf1)[(prevdf1$Prevalence >= prevalenceThreshold)]
FSGf = prune_taxa(keepTaxa, FSG)

##############################
#prevalence for grave data set; to be used for Simper and corncob
##############################
FSGBone<-subset_samples(FSG, Type == "Bone")
prevdf = apply(X = otu_table(FSGBone),
               MARGIN = ifelse(taxa_are_rows(FSGBone), yes = 1, no = 2),
               FUN = function(x){sum(x > 0)})
# Add taxonomy and total read counts to this data.frame
prevdf = data.frame(Prevalence = prevdf,
                    TotalAbundance = taxa_sums(FSGBone),
                    tax_table(FSGBone))
prevdf1 = subset(prevdf, Phylum %in% get_taxa_unique(FSGBone, "Phylum"))
ggplot(prevdf1, aes(TotalAbundance, Prevalence / nsamples(FSGBone),color=Phylum)) +
  # Include a guess for parameter
  geom_hline(yintercept = 0.01, alpha = 0.5, linetype = 2) +  geom_point(size = 2, alpha = 0.7) +
  scale_x_log10() +  xlab("Total Abundance") + ylab("Prevalence [Frac. Samples]") +
  facet_wrap(~Phylum) + theme(legend.position="none")

#get rid of phyla not seen across 1% of taxa, 4 samples
prevalenceThreshold = 0.01 * nsamples(FSGBone) #~4
prevalenceThreshold
# Execute prevalence filter, using `prune_taxa()` function
keepTaxa = rownames(prevdf1)[(prevdf1$Prevalence >= prevalenceThreshold)]
FSGBone2 = prune_taxa(keepTaxa,FSGBone)

#########################
#look at read distribution
#########################
sample_sum_df <- data.frame(Sample = sample_names(FSGf), sum = sample_sums(FSGf))

# Histogram of sample read counts
ggplot(sample_sum_df, aes(x = sum)) + 
  geom_histogram(color = "black", fill = "indianred", binwidth = 2500) +
  ggtitle("Distribution of sample sequencing depth") + 
  xlab("Read counts") +
  theme(axis.title.y = element_blank())

sample_sum_df1<- sample_sum_df %>%
  arrange(sum)

#Turn your 'treatment' column into a character vector
sample_sum_df1$Sample <- as.character(sample_sum_df1$Sample)
#Then turn it back into an ordered factor
sample_sum_df1$Sample<- factor(sample_sum_df1$Sample, levels=unique(sample_sum_df1$Sample))

ggplot(sample_sum_df1, aes(x = Sample, y=sum)) + 
  geom_bar(stat= "Identity") +
  ggtitle("Distribution of sample sequencing depth 16S rRNA") + 
  xlab("library") +
  ylab("read counts")+
  theme(axis.text.x = element_text(angle = 90, vjust=0.5, hjust = 1))

#prune samples less than 7,000 reads and 10000
FSGf2p<-prune_samples(sample_sums(FSGf)>=7000, FSGf)
pruned<-prune_samples(sample_sums(FSGf)<=10000, FSGf)
FSGf2p10<-prune_samples(sample_sums(FSGf)>=10000, FSGf)

# mean, max and min of sample read counts
smin <- min(sample_sums(FSGf2p10))
smean <- mean(sample_sums(FSGf2p))
smax <- max(sample_sums(FSGf2p))

print(smin) 
print(smean)
print(smax)

###########################
#look at relative abundance plot
##########################
#filter out AG_project
sample_data(FSGf2p)$Project<-as.factor(sample_data(FSGf2p)$Project)
sample_data(FSGf2p)$Type<-as.factor(sample_data(FSGf2p)$Type)
sample_data(FSGf2p)$Body.Site<-as.factor(sample_data(FSGf2p)$Body.Site)

MGR<- subset_samples(FSGf2p, Project != "AmGut" & Type != "PosControl"
                       & Body.Site != "Tooth")

#combine at phylum level and transform abundance to relative abundance
MGR_phylum <- MGR %>%
  tax_glom(taxrank = "Phylum") %>%                     # agglomerate at Phylum level
  transform_sample_counts(function(x) {x/sum(x)} ) %>% # Transform to rel. abundance
  psmelt() %>%                                         # Melt to long format
  
  arrange(Type,Body.Site, Bone)

MGR_phylumG<-MGR_phylum %>% filter(Project=="Grave")

#Average relative abundance of phyla by bone type. I have three replicates per treatment.need to first sum abundance by class to be divided by number of samples per group by merging data sets
MGR_phylum_1 <- MGR_phylum %>%
  dplyr::arrange(OTU, Project, Body.Site, Bone) %>%
  dplyr::group_by(Project, Body.Site, Bone, Phylum) %>%
  dplyr::summarise(avg.rel.abund = mean(Abundance))

#average by bone; get sample info
map2<-data.frame(sample_data(MGR))
map_Samp<-map2 %>% group_by(Project, Bone) %>% dplyr::summarise(n=n())
merged2<-merge(MGR_phylum_1, map_Samp)


#look at dominating taxa
taxabone<-merged2 %>% filter(Project == "Surface") %>%
  group_by(Phylum) %>%
  summarize(mean= mean(avg.rel.abund),max = max(avg.rel.abund), min = min(avg.rel.abund))

#looking at grave bones
taxaGrave<-merged2 %>%  
  filter(Project =="Grave" & Body.Site != "Grave"  
           & Body.Site != "ControlG" 
           & Body.Site != "ControlS" &
           Body.Site != "OG0.5 m") %>%
  group_by(Phylum) %>%
  summarize(mean= mean(avg.rel.abund),max = max(avg.rel.abund), min = min(avg.rel.abund))

#looking at soils
taxaSoil<-merged2 %>% filter(Project == "Grave_Soil") %>% 
group_by(Phylum) %>%
  summarize(mean= mean(avg.rel.abund),max = max(avg.rel.abund), min = min(avg.rel.abund))

#for graphing relabel anything with less than 0.01 relative abundance rare 
#labeling rare taxa
merged2$Phylum<-as.character(merged2$Phylum)
merged2$Phylum[merged2$avg.rel.abund <0.01] <- as.character("RareTaxa")
merged2$Phylum<-as.factor(merged2$Phylum)


library(randomcoloR)
n <- 21
palette <- distinctColorPalette(n)
#use these same colors for the same tatxa
myColors1 <- palette
names(myColors1)<-as.character(levels(merged2$Phylum))


#reorder Body site factor
merged2$Body.Site<-factor(merged2$Body.Site, levels = c("ControlG","ControlS", "OG0.5 m", "Grave", 
                           "Skull", "Upper Torso", "Arm", "Hand", "Lower Torso",
                           "Leg", "Foot"))
#change project levels so that grave soils and bones are from the same project. The dataframe was originally set up this way and I changed it
levels(merged2$Project)<-c("Grave", "Grave", "Surface")

#reset working directory for new data
setwd("./fixedNAanalysis/")
# Plot #need to fix grid facet headers so that they fit; originally only one project label grave
ggplot(merged2, aes(x = Bone, y = avg.rel.abund, fill = Phylum)) + 
  facet_grid(Project~Body.Site, scales="free", space="free") +
  geom_bar(stat = "identity") +
  scale_fill_manual(values= myColors1) +
  theme_classic()+
  # Remove x axis title
  theme(axis.title.x = element_blank()) + 
  guides(fill = guide_legend(reverse = FALSE, keywidth = 1, keyheight = 1)) +
  ylab("Average Relative Abundance") +
  #ggtitle("Phylum Level Composition of Bone-Associated Bacterial Communities Averaged Across Individuals") +
  theme(axis.text.x = element_text(angle = 90, vjust=0.5, hjust = 1, size=8), 
        strip.text.x = element_text(face="bold", angle=90, vjust=0.5, size=8)) +
  theme(legend.position = "bottom") +
  theme(legend.text=element_text(size=8)) +
  theme(axis.text.y=element_text(size=8),
        axis.title.y= element_text(size=8), legend.title=element_text(size=8))

ggsave("RelAbund.png", width=7.5, height=7.5, units="in", dpi=500)

####################
#relative abbundance by individual
####################
levels(MGR_phylum$Project)<-c("Grave", "Grave", "Surface")
JGrave<-filter(MGR_phylum, Project=="Grave" & Type != "Soil")


MGR_ind <- JGrave %>%
  dplyr::arrange(OTU, Individual, Body.Site, Bone) %>%
  dplyr::group_by(Individual, Body.Site, Bone, Phylum) %>%
  dplyr::summarise(avg.rel.abund = mean(Abundance)) 

#%>%dplyr::filter(avg.rel.abund > 0.01)
MGR_ind$Phylum<-as.character(MGR_ind$Phylum)
MGR_ind$Phylum[MGR_ind$avg.rel.abund <0.01] <- as.character("RareTaxa")
MGR_ind$Phylum<-as.factor(MGR_ind$Phylum)


library(randomcoloR)
n <- 3
palette3 <- distinctColorPalette(n)

saveRDS(myColors2, file = "myColors2.rds")

myColors2<-c(myColors1,palette3)
names(myColors2)[22]<-"Unclassified"
names(myColors2)[23]<-"Spirochaetes"
names(myColors2)[24]<-"Hydrogenedentes"
#reorder Body site factor
MGR_ind$Body.Site<-factor(MGR_ind$Body.Site, levels = c( 
                                                        "Skull", "Upper Torso", "Arm", "Hand", "Lower Torso",
                                                        "Leg", "Foot"))

# Plot #need to fix grid facet headers so that they fit
ggplot(MGR_ind, aes(x = Bone, y = avg.rel.abund, fill = Phylum)) + 
  facet_grid(Individual~Body.Site, scales="free", space="free") +
  geom_bar(stat = "identity") +
  scale_fill_manual(values= myColors2) +
  theme_classic()+
  # Remove x axis title
  theme(axis.title.x = element_blank()) + 
  guides(fill = guide_legend(reverse = TRUE, keywidth = 1, keyheight = 1)) +
  ylab("Average Relative Abundance") +
  #ggtitle("Phylum Level Composition of Bone-Associated Bacterial Communities Averaged Across Individuals") +
  theme(axis.text.x = element_text(angle = 90, vjust=0.5, hjust = 1, size=8), 
        strip.text.x = element_text(face="bold", angle=90, vjust=0.5, size=8)) +
  theme(legend.position = "bottom") +
  theme(legend.text=element_text(size=8)) +
  theme(axis.text.y=element_text(size=8),
        axis.title.y= element_text(size=8), legend.title=element_text(size=8))

ggsave("RelAbundbyInd.png", width=7.5, height=7.5, units="in", dpi=300)

#####################
#look at specific groups
#######################look at chloroflexi in grave and chlymidiae
CC<- subset_samples(MGR, Project=="Surface")
CC<-subset_taxa(CC, Phylum == "Chloroflexi")
CC<-prune_taxa(taxa_sums(CC) > 50, CC)
CC1<-tax_glom(CC, taxrank = "Genus")
# Phylum =="Chlamydiae" | 

plot_bar(CC1, x="Body.Site", y="Abundance", fill="Family")+
  scale_fill_manual(values=paletteC)

n <- 39
paletteC <- distinctColorPalette(n)

#####################
#rarefy data; justification Weiss et al.2017
######################
#scale to even depth at 10000 reads
normFGS<-rarefy_even_depth(FSGf, sample.size=10000, rngseed=3) #30 samples removed 123 OTUs removed

#save this object for metabolomics data
saveRDS(normFGS, file = "normFGS.rds")


#check to make sure this looks right
smin3 <- min(sample_sums(normFGS))
smean3 <- mean(sample_sums(normFGS))
smax3 <- max(sample_sums(normFGS))

print(smin3)
print(smean3)
print(smax3)

#look at 7000
normFGS7<-rarefy_even_depth(FSGf, sample.size=7000, rngseed=3) #229 OTUs removed

##############################
#View rarefaction curve
##############################
library(devtools)
#install_github("gauravsk/ranacapa")

library(ranacapa)
#remove mock controls
normfilt<- subset_samples(normFGS, Individual != "Seq")
normfilt1<- subset_samples(normFGS7, Individual != "Seq")

#remove teeth
normfilt <- normfilt %>%
  subset_samples(Body.Site != "Tooth")

teeth<-normFGS %>%
  subset_samples(Body.Site == "Tooth")

p <- ggrare(normfilt, step = 1000, color = "Project", se = FALSE)+
  theme(axis.text.x= element_text(angle=90))
p <- p + facet_wrap(~Individual) + theme_classic()+
  theme(axis.text.x= element_text( size=8))+
  theme(axis.text.y= element_text(size=8), axis.title.x= element_text(size=8),
        axis.title.y= element_text(size=8), legend.title= element_blank(),
        strip.text.x= element_text(size=8),
        legend.position="bottom", legend.text = element_text(size=8))

ggsave("bactrarefy.png", height=5, width=5.5, units="in", dpi=300)

############################
#Betadiversity
#NMDS and Bray-curtis
##combine at genus level to view all samples together
#############################
normfiltg <- normfilt %>%
  tax_glom(taxrank = "Genus") 


#ordinate nmds bray-curtis
nmds1<- phyloseq::ordinate(
  physeq = normfiltg, 
  method = "NMDS", 
  distance = "bray", 
  trymax=1000,
  k=3,
  set.seed(5) 
)
#solution reached stress=0.102; 453 samples, 1790 taxa; teeth removed

#rename levels so that AG_project has a label
sample_data(normfiltg)$Individual<-as.factor(sample_data(normfiltg)$Individual)
levels(sample_data(normfiltg)$Individual)<-c("AmerGut", "A", "B", "C", "SA", 
                                             "SB", "SC", "Soil")
sample_data(normfiltg)$Body.Site<-as.factor(sample_data(normfiltg)$Body.Site)

#colorblind package
library(viridis)
#text repel package
library(ggrepel)

#make character variables
sample_data(normfiltg)$Bone<-as.character(sample_data(normfiltg)$Bone)
sample_data(normfiltg)$Body.Site<-as.character(sample_data(normfiltg)$Body.Site)


#save color palette
cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", 
                "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

#make unifrac distance matrices with qiime2
unif<- read_qza("../qiime_exports/qiime2files/unifrac/unweighted_unifrac_BuriedBone.qza")$data
unifNMDS <- metaMDS(unif, k=3, trymax=200)
#stress 0.075


#change code for A and B; make sure to change axes (B includes axes= 2:3)
B<-plot_ordination(
  physeq = normfiltg,
  ordination = nmds1,
  axes= 3:2,
  color = "Individual",
  shape="Type"
  
) + 
  theme_classic()+
  geom_point(aes(color = Individual,shape=Type), alpha = 0.7, size = 0.5) +
  #geom_point(colour = "grey90", size = 1.5) +
  #scale_shape_manual(values=c(15,16,17,3,18,8,6,5))+ 
  #facet_grid(~Individual)+
  #stat_ellipse() +
  theme(legend.position = "bottom")+
  #geom_text_repel(aes(label = ifelse(Individual == "AmerGut", Body.Site, "")), show.legend=FALSE) +
  #geom_label_repel(aes(label = ifelse(Type == "Soil", Bone, "")),
                   #size=2, show.legend=FALSE)+
  #guides(guide_legend(override.aes = list(shape = 22)))+
  scale_color_manual(values=cbbPalette)+
  theme(axis.title.x = element_text(size=8), axis.text.x= element_text(size=8)) + 
  guides(shape = guide_legend(keywidth = 0.5, keyheight = 0.5, nrow=3), color = guide_legend(keywidth = 0.5, keyheight = 0.5, nrow=3)) +
  # ggtitle("Class Level Composition of Bone-Associated Bacterial Communities Averaged Across Individuals") +
  theme(axis.text.y = element_text(size=8), axis.title.y = element_text(size=8)) +
  theme(legend.text=element_text(size=8), legend.title=element_blank()) 


#combine bar plot with nmds
library(ggpubr)
#png("figrf.png", width=6, height=5.5, units= "in", res=300)
bray1<-ggarrange(A,B, nrow=2, common.legend=TRUE, legend="bottom")

braycomb<-ggarrange(bray1, bar3, nrow =1,ncol=2, labels=c("A","B"), heights = c(2, 1))

ggsave("figbraycomb.png", units="in", width = 6.5, height = 5,  dpi=300, device="png")


########################
#Relative abundance by individual to accompany beta diversity
########################
Ind_phylum <- normfiltg %>%
  tax_glom(taxrank = "Phylum") %>%                     # agglomerate at Phylum level
  transform_sample_counts(function(x) {x/sum(x)} ) %>% # Transform to rel. abundance
  psmelt() %>%  # Melt to long format
  filter(Type != "human-associated" & Type != "human-oral") %>%
  group_by(Individual)

merged2 <- Ind_phylum %>%
  dplyr::arrange(OTU, Project,Individual) %>%
  dplyr::group_by(Project, Individual, Phylum) %>%
  dplyr::summarise(avg.rel.abund = mean(Abundance)) 
#%>% dplyr::filter(avg.rel.abund > 0.01)

merged2$Phylum<-as.character(merged2$Phylum)
merged2$Phylum[merged2$avg.rel.abund <0.01] <- as.character("RareTaxa")
merged2$Phylum<-as.factor(merged2$Phylum)

n <- 20
palette <- distinctColorPalette(n)

# Plot #need to fix grid facet headers so that they fit
bar3<-ggplot(merged2, aes(x = Individual, y = avg.rel.abund, fill = Phylum)) + 
  #facet_grid(Project~Body.Site, scales="free", space="free") +
  geom_bar(stat = "identity") +
  scale_fill_manual(values= myColors2) +
  theme_classic()+
  # Remove x axis title
  theme(axis.title.x = element_blank()) + 
  guides(fill = guide_legend(reverse = TRUE, keywidth = 1, keyheight = 1)) +
  ylab("Avg. Rel. Abundance") +
  #ggtitle("Phylum Level Composition of Bone-Associated Bacterial Communities Averaged Across Individuals") +
  theme(axis.text.x = element_text(angle = 90, vjust=0.5, hjust = 1, size=8), 
        strip.text.x = element_text(face="bold", angle=90, vjust=0.5, size=8)) +
  theme(legend.position = "right") +
  theme(legend.text=element_text(size=8)) +
  theme(axis.text.y=element_text(size=8),
        axis.title.y= element_text(size=8), legend.title=element_text(size=8))+
  guides(fill = guide_legend(keywidth = 0.5, keyheight = 0.5, nrow=5,byrow=TRUE))

ggsave("ProjectRelAbund.png", width=5, height=3.5, units="in", dpi=300)

####################################
#Beta diversity
#weighted unifrac
####################################
#table of excluded samples
excluded<-dplyr::anti_join(data.frame(sample_data(FSG)), data.frame(sample_data(normfilt)))
excluded<-excluded[,c(1:11)]
write.csv(excluded,"exclude.csv")


#calculate weighted and unweighted distance matrices in qiime2 to use rooted tree exported from qiime2 or sepp tree.
#there have been reported inconsistencies with the UniFrac() function from phyloseq
#Because of sepp tree can use ASV level
write.table(otu_table(normfiltg),"rarefiedOTUtable10k.txt",sep="\t",row.names=TRUE,col.names=TRUE)
write.table(otu_table(normfilt),"rarefiedOTUtable10k_notG.txt",sep="\t",row.names=TRUE,col.names=TRUE)

set.seed(4)
#import distance matrix
wunif<- read_qza("../qiime_exports/qiime2files/unifrac/weighted_unifrac_BuriedBone.qza")$data
wunifNMDS <- metaMDS(wunif, k=3, trymax=200, set.seed(4))
#stress 0.104
wunif2<- read_qza("../qiime_exports/qiime2files/unifrac/weighted_unifrac_all_asv.qza")$data
wunifNMDS2 <- metaMDS(wunif2, k=3, trymax=200, set.seed(10))
#stress  0.105

R<-plot_ordination(
  physeq = normfilt,
  ordination = wunifNMDS2,
  axes=3:2,
  color = "Individual",
  shape = "Type"
  
) + 
  theme_classic()+
  geom_point(aes(color = Individual, shape=Type), alpha = 0.7, size = 0.5) +
  #geom_point(colour = "grey90", size = 1.5) +
  #scale_shape_manual(values=c(15,16,17,3,18,8,6,5))+ 
  #facet_grid(~Individual)+
  #stat_ellipse() +
  theme(legend.position = "bottom")+
  #geom_text(aes(label= Residence))+
  #geom_text(aes(label = ifelse(Body.Site == "Foot", "F", "")), show.legend=FALSE) +
  #geom_label_repel(aes(label = ifelse(Type == "Soil", Bone, "")),
  #size=2, show.legend=FALSE)+
  guides(guide_legend(override.aes = list(shape = 22)))+
  scale_color_manual(values=cbbPalette)+
  theme(axis.title.x = element_text(size=8), axis.text.x= element_text(size=8)) + 
  guides(shape = guide_legend(keywidth = 0.5, keyheight = 0.5, nrow=3), color = guide_legend(keywidth = 0.5, keyheight = 0.5, nrow=3)) +
  theme(axis.text.y = element_text(size=8), axis.title.y = element_text(size=8)) +
  theme(legend.text=element_text(size=8), legend.title=element_blank())

wuni1<-ggarrange(Q,R, nrow=2, common.legend=TRUE, legend="bottom")

wunicomb<-ggarrange(wuni1, bar3, nrow =1,ncol=2, labels=c("A","B"), heights = c(2, 1))

ggsave("figwunicomb.png", units="in", width = 6.5, height = 5,  dpi=300, device="png")


#combine bray and unweighted unifrac
BWU<-ggarrange(A,Q,B,R,ncol=2,nrow=2, common.legend=TRUE, legend="bottom", labels=c("A", "B"))

ABC<-ggarrange(BWU, bar3, nrow =2, labels=c("","C"), heights=c(1.75, 1.25))

ggsave("braywuni10000relabund.png", units="in", width = 6.5, height = 8,  dpi=500, device="png")

###############################
#ONLY GRAVE DATA
#beta diversity by individual; back to ASV level rather than genus level
##############################
sample_data(normfilt)$Project<-as.factor(sample_data(normfilt)$Project)
sample_data(normfiltg)$Project<-as.factor(sample_data(normfiltg)$Project)
levels(sample_data(normfilt)$Project)<-c("AmGut","Grave", "Grave", "Surface")
levels(sample_data(normfiltg)$Project)<-c("AmGut","Grave", "Grave", "Surface")

#get subset data sett
Grave<-subset_samples(normfilt, Project=="Grave")
Grave<-prune_taxa(taxa_sums(Grave) > 0, Grave) 

#get table for qiime2
write.table(otu_table(Grave),"rarefiedOTUtable10kGrave.txt",sep="\t",row.names=TRUE,col.names=TRUE)


#distance matrix from qiime2
gwunif_ind<-read_qza("./../qiime_exports/qiime2files/unifrac/weighted_unifrac_BuriedBone_grave.qza")$data

#ordinate using nmds
set.seed(27)
gwunifNMDS_ind <- metaMDS(gwunif_ind, k=2, trymax=400, set.seed(82)) #stress = 0.142
pcoa_unif<-ape::pcoa(gwunif_ind)


sample_data(Grave)$Bone<-as.character(sample_data(Grave)$Bone)
sample_data(Grave)$SampleID<-row.names(sample_data(Grave))
sample_data(Grave)$SampleID<-as.character(sample_data(Grave)$SampleID)
sample_data(Grave)$Body.Site<-as.character(sample_data(Grave)$Body.Site)

set.seed(4)
nmdsgrave<- ordinate(
  physeq = Grave, 
  method = "NMDS", 
  distance = "bray",
  trymax = 200,
  k=2,
  set.seed(4)
)


#Run 20 stress 0.128
plot_ordination(
  physeq = Grave,
  ordination = gwunifNMDS_ind,
  axes=1:2,
  color = "Individual"
  
) + 
  theme_classic()+
  geom_point(aes(color = Individual, shape = ifelse(Body.Site == "Foot" , Body.Site, ""))
             , alpha = 0.7, size = 2.0) +
  #geom_point(aes()) +
  scale_shape_manual(values=c(1,8))+ 
  #facet_grid(~Individual)+
  #stat_ellipse() +
  theme(legend.position = "bottom")+
  #geom_text(aes(label = ifelse(Type == "Soil", Depth, "")) +
  geom_text_repel(aes(label = ifelse(Type == "Soil", Bone, "")),
                  size=2, show.legend=FALSE, color="black")+
  #guides(guide_legend(override.aes = list(shape = 22)))
  scale_color_manual(values=c("#E69F00","#56B4E9","#009E73","#CC79A7"))+
  theme(axis.title = element_text(size=8), axis.text= element_text(size=8)) + 
  guides(color = guide_legend(keywidth = 0.5, keyheight = 0.5, nrow=2),
         shape = guide_legend(keywidth = 0.5, keyheight = 0.5)) +
  theme(legend.text=element_text(size=8), legend.title=element_text(size=8))

ggsave("figNMDSGravebray.png", width=5.5, height=3.5, units="in", dpi=300)

#############
#look at sequencing runs for batch effects
#############
G1<- subset_samples(Grave, Run=="G1")
G2<- subset_samples(Grave, Run=="G2")
G3<- subset_samples(Grave, Run=="G3")
G1ord<- ordinate(
  physeq = G1, 
  method = "NMDS", 
  distance = "bray",
  trymax = 300,
  k=2,
  set.seed(2)
)
#stress  0.16
#"Date.Sampled"
CD_G1<-plot_ordination(
  physeq = G1,
  ordination = G1ord,
  axes=1:2,
  color = "Date.Sampled"
  
) + 
  theme_classic()+
  geom_point(aes(color = Date.Sampled, shape = ifelse(Body.Site == "Foot" , Body.Site, ""))
             , alpha = 0.7, size = 2.0) +
  #geom_point(aes()) +
  scale_shape_manual(values=c(1,8))+ 
  facet_grid(~Run)+
  #stat_ellipse() +
  theme(legend.position = "bottom")+
  #geom_text(aes(label = ifelse(Type == "Soil", Depth, "")) +
  geom_text_repel(aes(label = ifelse(Type == "Soil", Bone, "")),
                  size=2, show.legend=FALSE, color="black")+
  #guides(guide_legend(override.aes = list(shape = 22)))
  #scale_color_manual(values=c("#E69F00","#56B4E9","#009E73","#CC79A7"))+
  theme(axis.title = element_text(size=8), axis.text= element_text(size=8)) + 
  guides(color = guide_legend(keywidth = 0.5, keyheight = 0.5, nrow=6,title=NULL),
         shape = guide_legend(keywidth = 0.5, keyheight = 0.5)) +
  theme(legend.text=element_text(size=8), legend.title=element_text(size=8))

Ind_G1<-plot_ordination(
  physeq = G1,
  ordination = G1ord,
  axes=1:2,
  color = "Individual"
  
) + 
  theme_classic()+
  geom_point(aes(color = Individual, shape = ifelse(Body.Site == "Foot" , Body.Site, ""))
             , alpha = 0.7, size = 2.0) +
  #geom_point(aes()) +
  scale_shape_manual(values=c(1,8))+ 
  facet_grid(~Run)+
  #stat_ellipse() +
  theme(legend.position = "bottom")+
  #geom_text(aes(label = ifelse(Type == "Soil", Depth, "")) +
  geom_text_repel(aes(label = ifelse(Type == "Soil", Bone, "")),
                  size=2, show.legend=FALSE, color="black")+
  #guides(guide_legend(override.aes = list(shape = 22)))
  scale_color_manual(values=c("#E69F00","#56B4E9","#009E73","#CC79A7"))+
  theme(axis.title = element_text(size=8), axis.text= element_text(size=8)) + 
  guides(color = guide_legend(keywidth = 0.5, keyheight = 0.5, nrow=2),
         shape = guide_legend(keywidth = 0.5, keyheight = 0.5)) +
  theme(legend.text=element_text(size=8), legend.title=element_text(size=8))

G2ord<- ordinate(
  physeq = G2, 
  method = "NMDS", 
  distance = "bray",
  trymax = 300,
  k=2,
  set.seed(2)
)
#stress  0.17
#"Date.Sampled"
CD_G2<-plot_ordination(
  physeq = G2,
  ordination = G2ord,
  axes=1:2,
  color = "Date.Sampled"
  
) + 
  theme_classic()+
  geom_point(aes(color = Date.Sampled, shape = ifelse(Body.Site == "Foot" , Body.Site, ""))
             , alpha = 0.7, size = 2.0) +
  #geom_point(aes()) +
  scale_shape_manual(values=c(1,8))+ 
  facet_grid(~Run)+
  #stat_ellipse() +
  theme(legend.position = "bottom")+
  #geom_text(aes(label = ifelse(Type == "Soil", Depth, "")) +
  geom_text_repel(aes(label = ifelse(Type == "Soil", Bone, "")),
                  size=2, show.legend=FALSE, color="black")+
  #guides(guide_legend(override.aes = list(shape = 22)))
  #scale_color_manual(values=c("#E69F00","#56B4E9","#009E73","#CC79A7"))+
  theme(axis.title = element_text(size=8), axis.text= element_text(size=8)) + 
  guides(color = guide_legend(keywidth = 0.5, keyheight = 0.5, nrow=6, title=NULL),
         shape = guide_legend(keywidth = 0.5, keyheight = 0.5)) +
  theme(legend.text=element_text(size=8), legend.title=element_text(size=8))

Ind_G2<-plot_ordination(
  physeq = G2,
  ordination = G2ord,
  axes=1:2,
  color = "Individual"
  
) + 
  theme_classic()+
  geom_point(aes(color = Individual, shape = ifelse(Body.Site == "Foot" , Body.Site, ""))
             , alpha = 0.7, size = 2.0) +
  #geom_point(aes()) +
  scale_shape_manual(values=c(1,8))+ 
  facet_grid(~Run)+
  #stat_ellipse() +
  theme(legend.position = "bottom")+
  #geom_text(aes(label = ifelse(Type == "Soil", Depth, "")) +
  geom_text_repel(aes(label = ifelse(Type == "Soil", Bone, "")),
                  size=2, show.legend=FALSE, color="black")+
  #guides(guide_legend(override.aes = list(shape = 22)))
  scale_color_manual(values=c("#56B4E9","#009E73","#CC79A7"))+
  theme(axis.title = element_text(size=8), axis.text= element_text(size=8)) + 
  guides(color = guide_legend(keywidth = 0.5, keyheight = 0.5, nrow=2),
         shape = guide_legend(keywidth = 0.5, keyheight = 0.5)) +
  theme(legend.text=element_text(size=8), legend.title=element_text(size=8))

G3ord<- ordinate(
  physeq = G3, 
  method = "NMDS", 
  distance = "bray",
  trymax = 300,
  k=2,
  set.seed(2)
)
#stress  0.107
#"Date.Sampled"
CD_G3<-plot_ordination(
  physeq = G3,
  ordination = G3ord,
  axes=1:2,
  color = "Date.Sampled"
  
) + 
  theme_classic()+
  geom_point(aes(color = Date.Sampled, shape = ifelse(Body.Site == "Foot" , Body.Site, ""))
             , alpha = 0.7, size = 2.0) +
  #geom_point(aes()) +
  scale_shape_manual(values=c(1,8))+ 
  facet_grid(~Run)+
  #stat_ellipse() +
  theme(legend.position = "bottom")+
  #geom_text(aes(label = ifelse(Type == "Soil", Depth, "")) +
  geom_text_repel(aes(label = ifelse(Type == "Soil", Bone, "")),
                  size=2, show.legend=FALSE, color="black")+
  #guides(guide_legend(override.aes = list(shape = 22)))
  #scale_color_manual(values=c("#E69F00","#56B4E9","#009E73","#CC79A7"))+
  theme(axis.title = element_text(size=8), axis.text= element_text(size=8)) + 
  guides(color = guide_legend(keywidth = 0.5, keyheight = 0.5, nrow=6,title=NULL),
         shape = guide_legend(keywidth = 0.5, keyheight = 0.5), nrow=2) +
  theme(legend.text=element_text(size=8), legend.title=element_text(size=8))

Ind_G3<-plot_ordination(
  physeq = G3,
  ordination = G3ord,
  axes=1:2,
  color = "Individual"
  
) + 
  theme_classic()+
  geom_point(aes(color = Individual, shape = ifelse(Body.Site == "Foot" , Body.Site, ""))
             , alpha = 0.7, size = 2.0) +
  #geom_point(aes()) +
  scale_shape_manual(values=c(1,8))+ 
  facet_grid(~Run)+
  #stat_ellipse() +
  theme(legend.position = "bottom")+
  #geom_text(aes(label = ifelse(Type == "Soil", Depth, "")) +
  geom_text_repel(aes(label = ifelse(Type == "Soil", Bone, "")),
                  size=2, show.legend=FALSE, color="black")+
  #guides(guide_legend(override.aes = list(shape = 22)))
  scale_color_manual(values=c("#E69F00","#56B4E9","#009E73","#CC79A7"))+
  theme(axis.title = element_text(size=8), axis.text= element_text(size=8)) + 
  guides(color = guide_legend(keywidth = 0.5, keyheight = 0.5, nrow=2),
         shape = guide_legend(keywidth = 0.5, keyheight = 0.5)) +
  theme(legend.text=element_text(size=8), legend.title=element_text(size=8))


Batch_ind<-ggarrange(Ind_G1,Ind_G2,Ind_G3,ncol=3, legend="bottom", labels=c("A", "B","C"))

Batch_run<-ggarrange(CD_G1, CD_G2,CD_G3, ncol =3)

Batch_comb<-ggarrange(Batch_ind,Batch_run,nrow=2)

ggsave("batcheffects.png", units="in", width = 10, height =7.5,  dpi=500, device="png")

###########################
#Within bone variation and replicate variation
##########################
Grave2<-subset_samples(normFGS, Project=="Grave" & Type == "Bone" |Type== "PosControl")
sample_data(Grave2)$SampleID<-rownames(sample_data(Grave2))
#grab otu table for wunifrac calculation
write.table(otu_table(Grave2),"rarefiedOTUtable10kGraveWpos.txt",sep="\t",row.names=TRUE,col.names=TRUE)

sample_data(Grave2)$Bone<-as.character(sample_data(Grave2)$Bone)
sample_data(Grave2)$Location<-as.character(sample_data(Grave2)$Location)
sample_data(Grave2)$Individual<-as.factor(sample_data(Grave2)$Individual)
sample_data(Grave2)$Type<-as.character(sample_data(Grave2)$Type)


Grave3<-Grave2 %>% subset_samples(Bone == "Femur" | Bone == "Humerus" | Bone == "Tibia" |Type == "PosControl")

#make bar plots
gwunif_ind<-read_qza("./qiime_exports/qiime2files/unifrac/weighted_unifrac_BuriedBone_grave_w_pos.qza")$data

#melt distance matrix
m_w<-melt(as.matrix(gwunif_ind))

# remove self-comparisons
m_w = m_w %>%
  filter(as.character(Var1) != as.character(Var2)) %>%
  mutate_if(is.factor,as.character)

dft <- data.frame(t(apply(m_w, 1, sort)))
dftu<-unique(dft)
colnames(dftu)<-c("value","Var1","Var")

# get sample data (S4 error OK and expected)
sd = data.frame(sample_data(Grave2)) %>%
  select(SampleID,Bone,Location,Body.Site,Individual) %>%
  mutate_if(is.factor,as.character)

dft <- data.frame(t(apply(m_w, 1, sort)))
dftu<-unique(dft)
colnames(dftu)<-c("value","Var1","Var2")

# combined distances with sample data
colnames(sd) = c("Var1", "Bone1","Location1","Body.Site1","Individual1")
w_sd = left_join(dftu, sd, by = "Var1")

colnames(sd) = c("Var2", "Bone2","Location2","Body.Site2","Individual2")
w_sd = left_join(w_sd, sd, by = "Var2")

w_sd$value<-as.numeric(w_sd$value)

#look at femur, humerus, tibia, posControl
#Grave3<-Grave2 %>% subset_samples(Bone == "Femur" | Bone == "Humerus" | Bone == "Tibia" |Type == "PosControl")
#A,B, and C by bone
femurA<-w_sd %>% filter(Individual1=="A" & Individual2=="A" & Bone1=="Femur" & Bone2=="Femur")
femurA$Comp<-"Femur"
femurA$Individual<-"A"
femurB<-w_sd %>% filter(Individual1=="B" & Individual2=="B" & Bone1=="Femur" & Bone2=="Femur")
femurB$Comp<-"Femur"
femurB$Individual<-"B"
femurC<-w_sd %>% filter(Individual1=="C" & Individual2=="C" & Bone1=="Femur" & Bone2=="Femur")
femurC$Comp<-"Femur"
femurC$Individual<-"C"
f_c<-rbind(femurA,femurB,femurC)

humerusA<-w_sd %>% filter(Individual1=="A" & Individual2=="A" & Bone1=="Humerus" & Bone2=="Humerus")
humerusA$Comp<-"Humerus"
humerusA$Individual<-"A"
humerusB<-w_sd %>% filter(Individual1=="B" & Individual2=="B" & Bone1=="Humerus" & Bone2=="Humerus")
humerusB$Comp<-"Humerus"
humerusB$Individual<-"B"
humerusC<-w_sd %>% filter(Individual1=="C" & Individual2=="C" & Bone1=="Humerus" & Bone2=="Humerus")
humerusC$Comp<-"Humerus"
humerusC$Individual<-"C"
h_c<-rbind(humerusA,humerusB,humerusC)
  
tibiaA<-w_sd %>% filter(Individual1=="A" & Individual2=="A" & Bone1 == "Tibia" & Bone2=="Tibia")
tibiaA$Comp<-"Tibia"
tibiaA$Individual<-"A"
tibiaB<-w_sd %>% filter(Individual1=="B" & Individual2=="B" & Bone1 == "Tibia" & Bone2=="Tibia")
tibiaB$Comp<-"Tibia"
tibiaB$Individual<-"B"
tibiaC<-w_sd %>% filter(Individual1=="C" & Individual2=="C" & Bone1 == "Tibia" & Bone2=="Tibia")
tibiaC$Comp<-"Tibia"
tibiaC$Individual<-"C"
t_c<-rbind(tibiaA,tibiaB,tibiaC)

PosC<-w_sd %>% filter(Bone1=="Control" & Bone2=="Control")
PosC$Comp<-"ZymoControl"
PosC$Individual<-"ZymoControl"
summary(PosC)

 mult<-rbind(f_c,h_c,t_c,PosC)
mult$Site1 <- paste(mult$Bone1,mult$Location1)
mult$Site2 <- paste(mult$Bone2,mult$Location2)

mult$Individual<-as.factor(mult$Individual)
levels(mult$Individual)<-c("A","B","C","Zymo")
hft<-ggplot(mult, aes(x = Comp, y =value)) +
  geom_boxplot(aes(color = ifelse(Site1 == Site2, "Within Sample Site", "Between Sample Sites"))) +
  scale_color_manual(values=c("black","red")) +
  ylab("Distance (weighted UniFrac)")+
  facet_grid(~Individual, scales = "free_x") +
  scale_y_continuous(limits=c(0,1.75))+
  theme_classic()+
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5), axis.title.x = element_blank(),
        legend.position="bottom", legend.title=element_blank())+
  guides(color=guide_legend(title="Comparison",nrow=2))
ggplot(mult, aes(x = Site1, y =value)) +
  geom_boxplot(aes(color = ifelse(Site1 == Site2, "Within Sample Site", "Between Sample Sites"))) +
  scale_color_manual(values=c("black","red")) +
  ylab("Distance (weighted UniFrac)")+
  facet_grid(~Individual, scales = "free_x") +
  scale_y_continuous(limits=c(0,1.75))+
  theme_classic()+
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5), axis.title.x = element_blank(),
        legend.position="bottom", legend.title=element_blank())+
  guides(color=guide_legend(title="Comparison",nrow=2))

ggsave("NegativeAnalysisDistBoxplotWuni.png", height=6.5, width=5, device="png", dpi=500)

###look at all bones
w_sd$Between_Within<-ifelse(w_sd$Bone1 == w_sd$Bone2, "Within Bones", "Between Bones")
cA<-w_sd %>% filter(Individual1=="A" & Individual2=="A")
cA$Individual<-"A"
cB<-w_sd %>% filter(Individual1=="B" & Individual2=="B")
cB$Individual<-"B"
cC<-w_sd %>% filter(Individual1=="C" & Individual2=="C")
cC$Individual<-"C"

allb<-rbind(cA,cB,cC)


#look at within bone values
within<-allb %>% filter(Between_Within=="Within Bones") %>% select_all() %>% arrange(value) %>% group_by(Individual) 

allb_fig<-ggplot(allb, aes(x = Between_Within, y =value, color=Between_Within)) +
  geom_boxplot() +
  scale_color_manual(values=c("black","red")) +
  ylab("Distance (weighted UniFrac)")+
  facet_grid(~Individual, scales = "free_x") +
  scale_y_continuous(limits=c(0,1.75))+
  theme_classic()+
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5), axis.title.x = element_blank(),
        legend.position="bottom",legend.title=element_blank())+
  guides(color=guide_legend(nrow=2))
ggplot(allb, aes(x = Bone1, y =value, color=Between_Within)) +
  geom_boxplot() +
  scale_color_manual(values=c("black","red")) +
  ylab("Distance (weighted UniFrac)")+
  facet_grid(Between_Within~Individual, scales = "free_x") +
  scale_y_continuous(limits=c(0,1.75))+
  theme_classic()+
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5), axis.title.x = element_blank(),
        legend.position="bottom",legend.title=element_blank())+
  guides(color=guide_legend(nrow=2))

multsite<-ggarrange(allb_fig,hft,ncol=2,nrow=1)
library(cowplot)
multsite<-plot_grid(allb_fig,hft,labels="AUTO",align="h",rel_widths = c(1.5,2))

ggsave("BoneSitesBoxplotWuni.png", height=5.0, width=6, device="png", dpi=500)

##############################
#look at permanova and statistics
#############################
# Calculate bray curtis distance matrix
projbray <- phyloseq::distance(normfiltg, method = "bray")

# make a data frame from the sample_data
sampledf <- data.frame(sample_data(normfiltg))

# Adonis test; also use wunif --use for Project and Type
set.seed(1)
adonis(wunif2 ~ Project, data = sampledf, permutations=999)
set.seed(1)
adonis(projbray ~ Type , data = sampledf, permutations=999)


##look only at bone samples
BoneSamp<- normfiltg %>% subset_samples(Type == "Bone")
BoneSamp1<- normfilt %>% subset_samples(Type == "Bone")

#bray and unifrac
B_projbray <- phyloseq::distance(BoneSamp, method = "bray")
#make a unifrac matrix
write.table(otu_table(BoneSamp1),"rarefiedOTUtable10kAllbone_asv.txt",sep="\t",row.names=TRUE,col.names=TRUE)
subsurfunif<-read_qza("../qiime_exports/qiime2files/unifrac/weighted_unifrac_allbones_asv.qza")$data
#new sample set
sampledf1 <- data.frame(sample_data(BoneSamp))
#Adonis test with interactions only bone samples; use same code for bray and wunif
set.seed(1)
adonis(B_projbray ~ Individual, data= sampledf1,  permutations=999)
adonis(subsurfunif ~ Project, data= sampledf1,  permutations=999)


#homogeneity of dispersion test
beta <- betadisper(wunif, sampledf$Type)
permutest(beta)

beta <- betadisper(projbray, sampledf$Type)
permutest(beta)

#substitute Project(surface vs subsurface) and Individual
beta <- betadisper(subsurfunif, sampledf1$Project)
permutest(beta)

beta <- betadisper(B_projbray, sampledf1$Individual)
permutest(beta)

###look only at bone samples from the grave
BoneSamp2<- normfilt %>% subset_samples(Project == "Grave" & Type == "Bone") 
write.table(otu_table(BoneSamp2),"rarefiedOTUtable10kAllbonegrave.txt",sep="\t",row.names=TRUE,col.names=TRUE)
gunifB<-read_qza("../qiime_exports/qiime2files/unifrac/weighted_unifrac_allbones_grave.qza")$data
#bray and unifrac
BG_projbray <- phyloseq::distance(BoneSamp2, method = "bray")
#new sample set
sampledf2 <- data.frame(sample_data(BoneSamp2))
#Adonis test with interactions only bone samples; use same code for bray and wunif
set.seed(1)
adonis(BG_projbray~ Individual*Body.Site, data= sampledf2,  permutations=999)
adonis(gunifB~ Individual*Body.Site, data= sampledf2,  permutations=999)

#Substitute individual and body site
beta <- betadisper(gunifB, sampledf2$Body.Site)
permutest(beta)

beta <- betadisper(BG_projbray, sampledf2$Body.Site)
permutest(beta)

############################################
#ordinations by indivdiual to look at body regions
############################################
mapNMDS<-sample_data((normfilt))
A_NMDS<-subset_samples(normfilt, Individual=="A")
B_NMDS<-subset_samples(normfilt, Individual=="B" )
C_NMDS<-subset_samples(normfilt, Individual=="C")
Soil_NMDS<-subset_samples(normfilt, Individual=="Soil")

A_nmds <- ordinate(
  physeq = A_NMDS, 
  method = "NMDS", 
  distance = "bray",
  trymax = 200,
  k=2,
  set.seed(5)
)

#stress 0.161 

#make bone not an integer
sample_data(A_NMDS)$Bone<-as.character(sample_data(A_NMDS)$Bone)
sample_data(A_NMDS)$Location<-as.character(sample_data(A_NMDS)$Location)

Ap<-plot_ordination(
  physeq = A_NMDS,
  ordination = A_nmds,
  axes= 1:2,
  color = "Body.Site"
  
) + 
  theme_classic()+
  geom_point(aes(color = Body.Site), alpha = 0.7, size = 2.0) +
  stat_ellipse(aes(color = Body.Site, group = Body.Site))+
  theme(legend.position = "bottom")+
  ggtitle("A")+
  scale_color_manual(values=cbbPalette)+
  #geom_label_repel(aes(label = ifelse(Body.Site == "Leg"|
  #                                      Body.Site == "Lower Torso", Bone, "")),
  #                 show.legend=FALSE, size=2)+
  guides(color=guide_legend(title="Body Site", nrow=4))+
  theme(legend.text = element_text(colour="black", size=6, 
                                   face="bold"), legend.title=element_text(size=6))+
  theme(legend.key.size = unit(2, "mm"), axis.text.x = element_text(size=8, color="black"), axis.text.y = element_text(size=8, color="black"))


plot_net(A_NMDS, color="Body.Site", distance="bray",  maxdist=0.4)

Bones_brayA <- phyloseq::distance(A_NMDS, method = "bray")


ggsave("figregionA.png", width=5, height=3.5, units="in", dpi=300)


B_nmds<- ordinate(
  physeq = B_NMDS, 
  method = "NMDS", 
  distance = "bray",
  trymax = 200,
  k=2,
  set.seed(6)
)
#stress 0.154

sample_data(B_NMDS)$Bone<-as.character(sample_data(B_NMDS)$Bone)
sample_data(B_NMDS)$Location<-as.character(sample_data(B_NMDS)$Location)

Bp<-plot_ordination(
  physeq = B_NMDS,
  ordination = B_nmds,
  axes= 1:2,
  color = "Body.Site"
  
) + 
  theme_classic()+
  geom_point(aes(color = Body.Site),alpha = 0.7, size = 2.0) +
  theme(legend.position = "bottom")+
  stat_ellipse(aes(color = Body.Site, group = Body.Site))+
  ggtitle("B")+
  scale_color_manual(values=cbbPalette)+
  guides(color=guide_legend(title="Body Site", nrow=4))+
  theme(legend.text = element_text(colour="black", size=6, 
                                   face="bold"), legend.title=element_text(size=6))+
  theme(legend.key.size = unit(2, "mm"), axis.text.x = element_text(size=8, color="black"), axis.text.y = element_text(size=8, color="black"))

ggsave("figregionB.png", width=5, height=3.5, units="in", dpi=300)

plot_net(B_NMDS, color="Body.Site", distance="bray",  maxdist=0.4)


set.seed(6)
c_nmds<- ordinate(
  physeq = C_NMDS, 
  method = "NMDS", 
  distance = "bray",
  trymax = 200,
  k=2
)
#stress 0.175
sample_data(C_NMDS)$Bone<-as.character(sample_data(C_NMDS)$Bone)
sample_data(C_NMDS)$Location<-as.character(sample_data(C_NMDS)$Location)


Cp<-plot_ordination(
  physeq = C_NMDS,
  ordination = c_nmds,
  axes= 1:2,
  color = "Body.Site"
  
) + 
  theme_classic()+
  geom_point(aes(color = Body.Site),alpha = 0.7, size = 2.0) +
  stat_ellipse(aes(color = Body.Site, group = Body.Site))+
  theme(legend.position = "bottom")+
  ggtitle("C")+
  scale_color_manual(values=cbbPalette)+
  #geom_label_repel(aes(label = ifelse(Body.Site == "Leg"|
   #                                     Body.Site == "Lower Torso", Bone, "")),
    #               show.legend=FALSE, size=2)+
  guides(color=guide_legend(title="Body Site",nrow=4))+
  theme(legend.text = element_text(colour="black", size=6, 
                                   face="bold"), legend.title=element_text(size=6))+
  theme(legend.key.size = unit(2, "mm"), axis.text.x = element_text(size=8, color="black"), axis.text.y = element_text(size=8, color="black"))

ggsave("figregionC.png", width=5, height=3.5, units="in", dpi=300)

plot_net(C_NMDS, color="Body.Site", distance="bray",  maxdist=0.5)

Regions<-ggarrange(Ap,Bp,Cp,ncol=1,nrow=3, common.legend=TRUE, legend="bottom")

ggsave("Regions.png", width=3.5, height=8, units="in", dpi=300)

#################################
#make venn diagrams
################################
#used code from https://rpubs.com/dillmcfarlan/R_microbiotaSOP
library(VennDiagram)
#subset and transpose and get rid of 0 OTUs; only bone samples from grave
Anaerobic<-normfiltg %>% subset_samples(Project == "Grave" & Type == "Bone" & Individual != "C") 
Anaerobic1<-prune_taxa(taxa_sums(Anaerobic) > 0, Anaerobic)
Anaerobict<-colnames(t(otu_table(Anaerobic1))) 

allbone<-normfiltg %>% subset_samples(Project == "Grave" & Type == "Bone" | Project=="Surface" & Type=="Bone") 
#subset surface samples and transpose (remove teeth)
Aerobic<-subset_samples(normfiltg, Project=="Surface" & Body.Site != "Tooth")
Aerobic1<-prune_taxa(taxa_sums(Aerobic) > 0, Aerobic)
Aerobict<-colnames(t(otu_table(Aerobic1)))

#subset C from grave
GraveC<-subset_samples(normfiltg, Individual == "C")
GraveC1<-prune_taxa(taxa_sums(GraveC) > 0, GraveC)
GraveCt<-colnames(t(otu_table(GraveC1)))

#separate soil samples; remove samples that clustered with individuals
soil.list<-c("MGR3170031", "MGR3170028", "MGR3170026", "MGR3170007", "MGR3170006")
Soil<-subset_samples(normfiltg, Individual == "Soil" & !SampleID %in%  soil.list)
Soil1<-prune_taxa(taxa_sums(Soil) > 0, Soil)
Soilt<-colnames(t(otu_table(Soil1)))

#separate american gut samples only the stool samples
AmGut<-subset_samples(normfiltg,Project == "AmGut" & Body.Site == "Stool")
AmGut1<-prune_taxa(taxa_sums(AmGut) > 0, AmGut)
AmGutT<-colnames(t(otu_table(AmGut1)))

#AmGut1<-subset_samples(normfiltg,Project == "AmGut" & Body.Site != "Stool")
#AmGut2<-prune_taxa(taxa_sums(AmGut1) > 0, AmGut1)
#AmGutT1<-colnames(t(otu_table(AmGut2)))

#And plot
library(gplots)
venn(list(Anaerobict,GraveCt, Aerobict, Soilt, AmGutT))

###use online tool http://bioinformatics.psb.ugent.be/webtools/Venn/
write.table(Anaerobict, "GraveAB2.csv", sep=",", row.names=FALSE, col.names=FALSE)
write.table(GraveCt, "GraveC2.csv", sep=",", row.names=FALSE, col.names=FALSE)
write.table(AmGutT, "AmGut2.csv", sep=",", row.names=FALSE, col.names=FALSE)
write.table(Soilt, "Soil2.csv", sep=",", row.names=FALSE, col.names=FALSE)
write.table(Aerobict, "Surface2.csv", sep=",", row.names=FALSE, col.names=FALSE)
#write.table(Hum, "OralHair.csv", sep=",", row.names=FALSE, col.names=FALSE)

###############
#look at sequence counts from fasta file (BoneBurialProject_Sequences.fasta) exported from rep seqs file
#used this code to get feature id and lengths from fasta 
#awk '$0 ~ ">" {if (NR > 1) {print c;} c=0;printf substr($0,2,100) "\t"; } $0 !~ ">" {c+=length($0);} END { print c; }'
###############
counts<-read.table("../SequenceCountsBBProject.txt", sep="\t")
colnames(counts)<-c("feature","count")
tcountnames<-taxa_names(normfilt)  
counts2<-filter(counts, feature %in% tcountnames)  
gnames<-taxa_names(AmGut1)
countsg<-filter(counts2, feature %in% gnames)  
surfnames<-taxa_names(Aerobic1)
countsSurf<-filter(counts2, feature %in% surfnames)  
gravnames<-taxa_names(Anaerobic1)
countsGrav<-filter(counts2, feature %in% gravnames)  
##################
#look at tentative gut genera
##################
#look at acinteobacter and stenotrophomonas; recalculate AmGut1 at ASV level above
#BLAST ASVs using fasta file
steno<-AmGut1 %>% subset_taxa(Genus=="Stenotrophomonas")

###look at taxa shared with gut
gut<-read.table("guttaxa.txt", sep="\t",header=FALSE)
gut1<-c(as.character(gut$V1))

gravegut<- normfiltg %>% transform_sample_counts(function(x) {x/sum(x)} )
gutsubset <- subset(otu_table(gravegut), rownames(otu_table(gravegut)) %in% gut1)
gutphy <- merge_phyloseq(gutsubset, tax_table(gravegut), sample_data(gravegut), phy_tree(gravegut))

#prune
gut2<-prune_samples(sample_sums(gutphy)>  0, gutphy)



gravegut<-subset_samples(gut2, Project == "Grave")
gravegut<-prune_samples(sample_sums(gravegut)>  0, gravegut)

gravegut1 <- gravegut  %>%                    
  psmelt() 



#Average relative abundance by bodysite. 
gravegut2 <- gravegut1 %>%
  dplyr::arrange(OTU, Individual, Body.Site) %>%
  dplyr::group_by(Individual, Body.Site, Phylum,Class, Order, Family, Genus, OTU) %>%
  dplyr::summarise(avg.rel.abund = mean(Abundance)) %>%
  dplyr::filter(avg.rel.abund > 0)

gravegut2$Genus<-as.character(as.factor(gravegut2$Genus))
gravegut2$Family<-as.character(as.factor(gravegut2$Family))
gravegut2$Genus[is.na(gravegut2$Genus)] <- gravegut2$Order[is.na(gravegut2$Genus)]
gravegut2$Genus[gravegut2$avg.rel.abund <0.005] <- as.character("RareTaxa")

n <- 10
palette2 <- distinctColorPalette(n)


# Plot #need to fix grid facet headers so that they fit
gutplot<-ggplot(gravegut2, aes(x = Body.Site, y = avg.rel.abund, color=OTU,fill = Genus)) + 
  #facet_grid(Project~Body.Site, scales="free", space="free") +
  facet_grid(~Individual, space = "free", scales="free")+
  geom_bar(stat = "identity", color='black') +
  scale_fill_manual(values= palette2) +
  theme_classic()+
  # Remove x axis title
  guides(fill = guide_legend(reverse = TRUE, keywidth = 1, keyheight = 1)) +
  ylab("Avg. Rel. Abundance") +
  #ggtitle("Phylum Level Composition of Bone-Associated Bacterial Communities Averaged Across Individuals") +
  theme(axis.text.x = element_text(angle = 90, vjust=0.5, hjust = 1, size=8), 
        strip.text.x = element_text(face="bold", angle=90, vjust=0.5, size=8)) +
  theme(legend.position = "bottom") +
  theme(legend.text=element_text(size=8)) +
  theme(axis.text.y=element_text(size=8),
        axis.title= element_text(size=8), legend.title=element_blank())+
  guides(fill = guide_legend(keywidth = 0.5, keyheight = 0.5, ncol=2,byrow=TRUE),
         color=FALSE)
ggsave("gutplot.png", width=5.5, height=5.5, units="in", dpi=300)


gtree<-plot_tree(gravegut, color = "Body.Site", shape="Phylum", label.tips = "Family", 
          size = "abundance", plot.margin = 0.5, ladderize = TRUE)+
  scale_color_manual(values=cbbPalette)


#try a heatmap
plot_heatmap(gravegut, sample.label="Bone", 
             taxa.label="Family", taxa.order="Phylum", low="white", high="purple", 
             na.value="grey") + facet_grid(Individual~Body.Site, space = "free", scales= "free")


n <- 11
palette <- distinctColorPalette(n)

###########################################
#Alpha Diversity
############################################
# Initialize matrices to store richness and evenness estimates
nsamp = nsamples(FSGf2p10) #REASON WE ARE DOING 7000?
trials = 100

richness <- matrix(nrow = nsamp, ncol = trials)
row.names(richness) <- sample_names(FSGf2p10)

evenness <- matrix(nrow = nsamp, ncol = trials)
row.names(evenness) <- sample_names(FSGf2p10)

# It is always important to set a seed when you subsample so your result is replicable 
set.seed(3)

for (i in 1:100) {
  # Subsample
  r <- rarefy_even_depth(FSGf2p10, verbose = FALSE, replace = TRUE)
  
  # Calculate richness
  rich <- as.numeric(as.matrix(estimate_richness(r, measures = "Observed")))
  richness[ ,i] <- rich
  
  # Calculate evenness
  even <- as.numeric(as.matrix(estimate_richness(r, measures = "InvSimpson")))
  evenness[ ,i] <- even
}
saveRDS(richness, file = "richness.rds")

#Let's calculate the mean and standard deviation per sample for observed richness and inverse simpson's index and store those values in a dataframe.

# Create a new dataframe to hold the means and standard deviations of richness estimates
SampleID <- row.names(richness)
mean <- apply(richness, 1, mean)
sd <- apply(richness, 1, sd)
measure <- rep("Richness", nsamp)
rich_stats <- data.frame(SampleID, mean, sd, measure)

# Create a new dataframe to hold the means and standard deviations of diversity estimates
SampleID <- row.names(evenness)
mean <- apply(evenness, 1, mean)
sd <- apply(evenness, 1, sd)
measure <- rep("Inverse Simpson", nsamp)
even_stats <- data.frame(SampleID, mean, sd, measure)

alpha <- rbind(rich_stats, even_stats)

s <- data.frame(sample_data(FSGf2p10))
rownames(s)
s$SampleID<-rownames(s)

library(data.table)
#setnames(s, "SampleID", "SampleID")
alphadiv <- merge(alpha, s, by = "SampleID") 

alphadiv1 <- alphadiv %>% dplyr::filter(Body.Site != "Tooth" & Body.Site != "Mouth"&
                                          Body.Site != "Control" &
                                          Body.Site != "Hair") %>%
  group_by(Body.Site) 

alphadiv1<-droplevels(alphadiv1)

#Turn your 'treatment' column into a character vector
alphadiv1$Bone <- as.character(alphadiv1$Bone)
#Then turn it back into an ordered factor
alphadiv1$Bone<- factor(alphadiv1$Bone, levels=unique(alphadiv1$Bone))

#make measure a factor
alphadiv1$measure<-as.factor(alphadiv1$measure)

#reorder Body site factor
levels(as.factor(alphadiv1$Body.Site))
alphadiv1$Body.Site<-factor(alphadiv1$Body.Site, levels = c("ControlS", "OG0.5 m", "Grave","Stool",
                                                        "Skull", "Upper Torso", "Arm", "Hand", "Lower Torso",
                                                        "Leg", "Foot"))

#name AGP samples
alphadiv1$Individual<-as.factor(alphadiv1$Individual)
levels(alphadiv1$Individual)<- c("AmerGut", "A", "B", "C", "SA", "SB", "SC", "Soil")

ggplot(alphadiv1, aes(x = Bone, y =mean, group= Individual, color = Project, shape = Individual)) +
  geom_point(size = 1) +
  #geom_line(aes(group=Individual), colour="black")+
  facet_grid(measure~Body.Site, scales = "free", space = "free_x") +
  scale_color_manual(values = c("blue", "springgreen3", "purple", "orange")) +
  scale_shape_manual(values = c(15,16,17,3,18,8,6,5,7))+
  theme_classic()+
  theme(axis.title.x = element_blank(), axis.title.y =element_blank(),legend.title = element_text(size=8),
        legend.text =element_text(size=8),
        legend.position="bottom") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5, size=8),
        strip.text.x = element_text(angle =90, size=8),
        axis.text.y = element_text(size=7), strip.text.y= element_text(size=8))
  geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=.2,
                position=position_dodge(0.1), colour="black")+
    guides(shape = guide_legend(keywidth = 0.1,keyheight = 0.1 , nrow=3,alpha=1,byrow=TRUE),
           color = guide_legend(keywidth = 0.1,keyheight = 0.1,byrow=TRUE )) 
ggsave("Alpha.png", width=6.5, height=5, units="in", dpi=300)


#divide back into two datasets with sample info
richness<- filter(alphadiv1, measure == "Richness")
invsimp<- filter(alphadiv1, measure == "Inverse Simpson")
invsimp<-mutate(invsimp, transHuman = log10(Human..ng.gbp.))
richness<-mutate(richness, transHuman = log10(Human..ng.gbp.))


#diversity does not relate to human DNA
invsimp %>% filter(Type == "Bone") %>% 
  ggplot(aes(y=transHuman, x=mean))+
  geom_point()+
  facet_grid(Individual~Body.Site,  space="free", scales = "free")+
  geom_smooth(formula=y~x, method="lm")


#divide diversity by individual
t2<-filter(test, p.adj < 0.05
)

#load gpubr
library(ggpubr)
#compare means, requires ggpubr, continuously changed values for mutliple figures
compare_means(mean ~ Individual,  data = BoneIS, p.adjust.method="fdr")

GraveRich<-filter(richness, Project=="Grave" & Type == "Bone")
BoneRich<-filter(richness, Type == "Bone")
GraveIS<-filter(invsimp, Project=="Grave" & Type == "Bone")
BoneIS<-filter(invsimp, Type == "Bone")

kruskal.test(mean~Individual, data=BoneRich)

IA<-filter(invsimp, Individual=="A")
IB<- filter(invsimp, Individual=="B")
IC<- filter(invsimp, Individual=="C")

###invsimpson by individual
ggplot(BoneIS, aes(x = Individual, y =mean, fill=Individual)) +
  geom_boxplot() +
  ylab("Inverse Simpson Index")+
  xlab("Individual")+
  #geom_point(size = 2) +
  #geom_line(aes(group=Individual), colour="black")+
  facet_wrap(~measure, ncol = 1, scales = "free") +
  scale_color_manual(values = c("#a65628", "red", "#ffae19", "#4daf4a", "#1919ff", "darkorchid3", "magenta", "darkmagenta")) +
  theme_classic()+
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5))+
  #stat_compare_means(comparisons= my_comparisons, symnum.args = list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), symbols = c("****", "***", "**", "*", "ns")))+
  #stat_compare_means(label.x = "B", label.y= 160)+
  theme(legend.position="right")


###richness by individual
ggplot(BoneRich, aes(x = Individual, y =mean, fill=Individual)) +
  geom_boxplot() +
  #geom_point(size = 2) +
  #geom_line(aes(group=Individual), colour="black")+
  ylab("Richness (Observed)")+
  xlab("Individual")+
  facet_wrap(~measure, ncol = 1, scales = "free") +
  scale_color_manual(values = c("#a65628", "red", "#ffae19", "#4daf4a", "#1919ff", "darkorchid3", "magenta", "darkmagenta")) +
  theme_classic()+
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5))+
  #stat_compare_means(comparisons=my_comparisons, symnum.args = list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), symbols = c("****", "***", "**", "*", "ns")))+
  #stat_compare_means(label.x = "B", label.y= 2500)+
  theme(legend.position="right")

###Richness by body region
RichRegion<-ggplot(richness, aes(x = Body.Site, y =mean, fill=Body.Site)) +
  geom_boxplot() +
  ylab("Richness (Observed)")+
  xlab("Body Region")+
  facet_grid(measure~Individual,scales = "free", space = "free") +
  scale_color_manual(values = c("#a65628", "red", "#ffae19", "#4daf4a", "#1919ff", "darkorchid3", "magenta", "darkmagenta")) +
  theme_classic()+
  theme(axis.text.x = element_blank(), legend.title=element_text(size=8),
        legend.text = element_text(size=8), axis.title.x=element_blank())+
  guides(fill = guide_legend(title="Anatomical Region / Site"))+
  #stat_compare_means(label.x = "leg", label.y= 3000)+
  theme(legend.position="none")

ggsave("alphaRegionRich.png", height=3.5, width=5, units="in", dpi=300)


InvRegion<-ggplot(invsimp, aes(x = Body.Site, y =mean, fill=Body.Site)) +
  geom_boxplot() +
  ylab("Inverse Simpson (Diversity)")+
  xlab("Body Region")+
  facet_grid(measure~Individual, scales = "free", space= "free") +
  scale_color_manual(values = c("#a65628", "red", "#ffae19", "#4daf4a", "#1919ff", "darkorchid3", "magenta", "darkmagenta")) +
  theme_classic()+
  theme(axis.text.x = element_blank(), legend.title=element_text(size=8),
        legend.text = element_text(size=8), axis.title.x=element_blank())+
  guides(fill = guide_legend(title="Anatomical Region / Site"))+
  theme(legend.position="none")

ggsave("alphaRegionISimp.png", height=5, width=5, units="in", dpi=300)

ggarrange(RichRegion, InvRegion,ncol=1,nrow=2, common.legend=TRUE, legend="bottom")
ggsave("alphaRegion.png", height=6, width=5, units="in", dpi=300)

#put alpha diversity for body region in a table
alphatab<-alphadiv1 %>% 
  select(measure, Individual, Body.Site, mean) %>%
  group_by(measure, Individual, Body.Site) %>%
  
  summarise( sd= sd(mean),mean1= mean(mean), max= max(mean), min = min(mean))
#individual summaries
alphatabind<-alphadiv1 %>% 
  select(measure, Individual, mean) %>%
  group_by(measure, Individual) %>%
  
  summarise( sd= sd(mean),mean1= mean(mean),max= max(mean), min = min(mean))


write.csv(alphatabind,"alpha.csv")

BoneIS %>% 
  select(measure, Project, mean) %>%
  group_by(measure, Project) %>%
  
  summarise( sd= sd(mean),mean1= mean(mean), max= max(mean), min = min(mean))

#plot together

png("figAlphaRegion.png", height=7, width=, units="in", res=300)
grid.arrange(InvRegion, RichRegion, cols=2)
dev.off()


plot_richness(Bones_scale1, x="Body.region", measures=c("Observed"))+
  theme_classic()+
  facet_grid(~Individual)+
  geom_boxplot(aes(group=Body.region))+
  xlab("Body Region")+
  ylab("Richness(Observed)")+
  theme(legend.title= element_blank())+
  theme(axis.text.x= element_text(angle = 90, vjust=0.5, hjust = 1))

########################
#SIMPER
########################
source("/Volumes/GoogleDrive/My Drive/BoneSurfaceProject/16S_Data/Phyloseq2018/simper_pretty.r")
source("./R_krusk.r")
#simper by individual

#use FSGBone2; dropped taxa not seen across 1% of samples
FSGBoneR<-rarefy_even_depth(FSGBone2, sample.size=10000, rngseed=3) #5 samples removed

#filter data set to include only grave samples
Simper<-FSGBoneR %>% subset_samples(Project == "Grave" & Type == "Bone")
Simper<-prune_taxa(taxa_sums(Simper) > 0, Simper)

otuAll<-otu_table(Simper)
simper(otuAll, data.frame(sample_data(Simper))$Individual, permutations=100)
simper.pretty(otuAll, data.frame(sample_data(Simper)), c('Individual'), perc_cutoff=1, low_cutoff = 'y', low_val=0.01, 'IndividualSimper')

simper.results = data.frame(read.csv("IndividualSimper_clean_simper.csv"))


kruskal.pretty(otuAll, data.frame(sample_data(Simper)), simper.results, c('Individual'), 'Individual')


krusk<-read.table("Individual_krusk_simper.csv", header=TRUE, sep=",")
#conservative alpha of 0.01 used for significance
krusk1<-krusk%>% filter(fdr_krusk_p.val<=0.01)

#get taxonomic info to match to the table
sigtaxa<- c(as.character(krusk1$OTU))


simpdat<-Simper %>% subset_taxa(rownames(tax_table(Simper)) %in% sigtaxa)
a<-data.frame(tax_table(simpdat))
a$OTU<-rownames(a)



krusk2<-merge(krusk1, a, by="OTU")
write.csv(krusk2, "SimpInd.csv")

Simtrans <- Simper  %>%                     # agglomerate at Phylum level
  transform_sample_counts(function(x) {x/sum(x)} ) %>% # Transform to rel. abundance
  psmelt() 


#Average relative abundance by bodysite. 
Indsimper <- Simtrans %>%
  dplyr::arrange(OTU, Individual) %>%
  dplyr::group_by(Individual, Phylum,Class, Order, Family, Genus,Species, OTU) %>%
  dplyr::summarise(avg.rel.abund = mean(Abundance)) %>%
  dplyr::filter(OTU %in% sigtaxa & avg.rel.abund > 0)

Indsimper$Genus<-as.character(as.factor(Indsimper$Genus))
Indsimper$Family<-as.character(as.factor(Indsimper$Family))
Indsimper$Genus[Indsimper$Genus=="Unclassified"] <- Indsimper$Family[Indsimper$Genus=="Unclassified"]


n <- 10
palette <- distinctColorPalette(n)


# Plot #need to fix grid facet headers so that they fit
isim<-ggplot(Indsimper, aes(x = Individual, y = avg.rel.abund, color=OTU,fill = Genus)) + 
  #facet_grid(Project~Body.Site, scales="free", space="free") +
  geom_bar(stat = "identity", color='black') +
  scale_fill_manual(values= palette) +
  theme_classic()+
  # Remove x axis title
  guides(fill = guide_legend(reverse = TRUE, keywidth = 1, keyheight = 1)) +
  ylab("Avg. Rel. Abundance") +
  #ggtitle("Phylum Level Composition of Bone-Associated Bacterial Communities Averaged Across Individuals") +
  theme(axis.text.x = element_text(angle = 90, vjust=0.5, hjust = 1, size=8), 
        strip.text.x = element_text(face="bold", angle=90, vjust=0.5, size=8)) +
  theme(legend.position = "right") +
  theme(legend.text=element_text(size=8)) +
  theme(axis.text.y=element_text(size=8),
        axis.title= element_text(size=8), legend.title=element_blank())+
  guides(fill = guide_legend(keywidth = 0.5, keyheight = 0.5, ncol=2,byrow=TRUE),
         color=FALSE)

ggsave("SimperInd.png", height=3.5, width=5, units="in", dpi=300)


###########################
#look at soils and soil geochemistry
###########################
#look at data rarefied to 7,000 to include more soils samples
Soil<-normFGS7 %>% subset_samples(Type == "Soil")
Soil<-prune_taxa(taxa_sums(Soil) > 0, Soil) ##22 samples and 3,090 samples; prevalence filtered was ~2 samples

set.seed(1)
SoilNMDS<- phyloseq::ordinate(
  physeq = Soil, 
  method = "NMDS", 
  distance = "bray",
  trymax=100,
  k=2
)
sample_data(Soil)$Bone<-as.character(sample_data(Soil)$Bone)
sample_data(Soil)$Body.Site<-as.character(sample_data(Soil)$Body.Site)

S<-plot_ordination(
  physeq = Soil,
  ordination = SoilNMDS,
  axes= 1:2,
  color = "Bone",
  shape = "Body.Site"
  
) + 
  theme_classic()+
  geom_point(aes(color = Bone, shape=Body.Site), alpha = 0.7, size = 0.5) +
  #geom_point(colour = "grey90", size = 1.5) +
  #scale_shape_manual(values=c(15,16,17,3,18,8,6,5))+ 
  #facet_grid(~Individual)+
  #stat_ellipse() +
  theme(legend.position = "bottom")+
  geom_text_repel(aes(label = Depth), show.legend=FALSE) +
  #geom_label_repel(aes(label = ifelse(Type == "Soil", Bone, "")),
  #size=2, show.legend=FALSE)+
  guides(guide_legend(override.aes = list(shape = 22)))+
  scale_color_manual(values=cbbPalette)+
  theme(axis.title.x = element_text(size=8), axis.text.x= element_text(size=8)) + 
  guides(fill = guide_legend(keywidth = 0.5, keyheight = 0.5)) +
  # ggtitle("Class Level Composition of Bone-Associated Bacterial Communities Averaged Across Individuals") +
  theme(axis.text.y = element_text(size=8), axis.title.y = element_text(size=8)) +
  theme(legend.text=element_text(size=8), legend.title=element_blank()) 
cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", 
                "#F0E442", "#0072B2", "#D55E00", "#CC79A7")


#calculate weighted unifrac
#get otu_table for qiime2
write.table(otu_table(Soil),"rarefiedOTUtable7kSoil.txt",sep="\t",row.names=TRUE,col.names=TRUE)
#distance matrix from qiime2
gwunif_soil<-read_qza("./../qiime_exports/qiime2files/unifrac/weighted_unifrac_gravesoils.qza")$data

#ordinate using nmds
set.seed(27)
Soilnmds<- metaMDS(gwunif_soil, k=2, trymax=1000, set.seed(27)) #0.046

Soilu<-plot_ordination(
  physeq = Soil,
  ordination = Soilnmds,
  axes= 1:2,
  color = "Bone",
  shape = "Body.Site"
  
) + 
  theme_classic()+
  geom_point(aes(color = Bone, shape=Body.Site), alpha = 0.7, size = 0.5) +
  #geom_point(colour = "grey90", size = 1.5) +
  #scale_shape_manual(values=c(15,16,17,3,18,8,6,5))+ 
  #facet_grid(~Individual)+
  #stat_ellipse() +
  theme(legend.position = "bottom")+
  geom_text_repel(aes(label= Depth), show.legend=FALSE)+
  #geom_text(aes(label = ifelse(Body.Site == "Foot", "F", "")), show.legend=FALSE) +
  #geom_label_repel(aes(label = ifelse(Type == "Soil", Bone, "")),
  #size=2, show.legend=FALSE)+
  guides(guide_legend(override.aes = list(shape = 22)))+
  scale_color_manual(values=cbbPalette)+
  theme(axis.title.x = element_text(size=8), axis.text.x= element_text(size=8)) + 
  guides(fill = guide_legend(keywidth = 0.5, keyheight = 0.5)) +
  # ggtitle("Class Level Composition of Bone-Associated Bacterial Communities Averaged Across Individuals") +
  theme(axis.text.y = element_text(size=8), axis.title.y = element_text(size=8)) +
  theme(legend.text=element_text(size=8), legend.title=element_blank())





Soil2<-sample_data(Soil)[ ,c("Bacteria..copies.gbp.", "Fungi..copies.gbp.", "X..soil.moisture","pH", "Ammonium..ug.g.dry.weight.soil.", "DON..mg.N.L.gdw.", "Nitrification.potential.mg.NO2.g.dry.weight..day","PDE.nmols.hr.gdw", "G.nmols.hr.gdw","CB.nmols.hr.gdw", "LAP.umols.hr.gdw","Conductivity","Nitrate.inorga.I1.T46nic..ug.g.dry.weight.","DOC..ug.C.gdw.", "Microbial.respiration..CO2.released.gdw.day.")]
Soil2$SampleID<-as.character(row.names(Soil2))
Soil2<-as.tibble(Soil2)


#gene abundances missing so load separate soil meta data and merge
Soil3<-read.csv("../../../MGR_soilchem_meta.csv", header = TRUE)
Soil3$SampleID<-as.character(Soil3$SampleID)
Soil3$SampleID<-(gsub("BN","", Soil3$SampleID))
Soil3$SampleID<-(gsub("MGR0","MGR", Soil3$SampleID))



Soil4<-semi_join(Soil3, Soil2, by="SampleID")

#remove a few columns
Soil5<-Soil4[ ,-c(1,3)]


library(vegan)
library(ggplot2)
library(grid)

scrs <- as.data.frame(scores(Soilnmds, display = "sites"))

#save groups
group1<-c(as.character(sample_data(Soil)$Body.Site))
group2<-c(sample_data(Soil)$Depth)
scrs <- cbind(scrs, Type = group1, Depth=group2)



#fit continuous variables
fit.wUF = envfit(Soilnmds, Soil5, na.rm=TRUE)
fit.wUF

spp.scrs <- as.data.frame(scores(fit.wUF, display = "vectors"))
spp.scrs <- cbind(spp.scrs, GeoChem = rownames(spp.scrs))

#make depth a factor
scrs$Depth<-as.factor(as.character(scrs$Depth))

#rename geochemistry variables
spp.scrs$GeoChem<-as.factor(spp.scrs$GeoChem)
levels(spp.scrs$GeoChem)<-c("Ammonium", "CB","COL", "Conductivity", "DOC",
                            "DON", "LAP", "Respiration", "NAG", "Nitrate", "Nitrification",
                            "PDE", "pH", "Total 16S rRNA", "Total ITS","Total DNA", "Soil Moisture")

p <- ggplot(scrs) +
  geom_point(mapping = aes(x = NMDS1, y = NMDS2, colour = Type, shape=Depth)) +
  coord_fixed() + ## need aspect ratio of 1!
  geom_segment(data = spp.scrs,
               aes(x = 0, xend = NMDS1, y = 0, yend = NMDS2),
               arrow = arrow(length = unit(0.1, "cm")), colour = "black", size=0.1) +
  geom_text_repel(data = spp.scrs, aes(x = NMDS1, y = NMDS2, label = GeoChem),
            size = 2)+
  theme_classic()+
  scale_color_manual(values=cbbPalette)+
  theme(axis.title = element_text(size=8), axis.text= element_text(size=8)) + 
  guides(color = guide_legend(keywidth = 0.5, keyheight = 0.5, ncol=1),
         shape = guide_legend(keywidth = 0.5, keyheight = 0.5, ncol=1)) +
  theme(legend.text=element_text(size=8), legend.title=element_text(size=8), legend.position="bottom")

ggsave("SoilNMDS.png", width=4, height=5, units="in", dpi=300)

##combine with figure 4
SoilG<-ggarrange(ind2, p, ncol =2, labels=c("A","B"), align="hv", widths=c(1.75, 1.25))
ggsave("CombGraveNMDS.png", width=6.5, height=5, units="in", dpi=300)

############################
#Corncob
############################
#corncob
library(corncob)

#differential abundance tests are compared by factors and the first factor level is the baseline.

#filter data set to include only grave samples
CCdf<-FSGBone2 %>% subset_samples(Project == "Grave" & Type == "Bone")


#rearrange individuals so that the baseline is individual C instead of A
sample_data(CCdf)$Individual = factor(sample_data(CCdf)$Individual,levels=c("C","B","A"))
#make body site a factor and rearrange the factor levels so that bones from the feet are the baseline.
sample_data(CCdf)$Body.Site = factor(sample_data(CCdf)$Body.Site)
sample_data(CCdf)$Body.Site = factor(sample_data(CCdf)$Body.Site,levels(sample_data(CCdf)$Body.Site)[c(2,1,3,4)])

#CCdf_g<-CCdf %>% tax_glom(taxrank = "Genus")

sample_data(CCdf)$Oxygen<-sample_data(CCdf)$Individual
sample_data(CCdf)$Oxygen<-as.factor(sample_data(CCdf)$Oxygen)
levels(sample_data(CCdf)$Oxygen)<-c("C","A/B","A/B")
set.seed(1) #from vignette tutorial; set seed and run differential abundance test controlling for the effect of individual on dispersion
da_analysis <- differentialTest(formula = ~ Individual,
                                phi.formula = ~ Individual,
                                formula_null = ~ 1,
                                phi.formula_null = ~ Individual,
                                test = "Wald", boot = FALSE,
                                data = CCdf,
                                fdr_cutoff = 0.001)

da_analysis
da_analysis$significant_taxa

#plot and view only species and strain levels
dfccob<-plot(da_analysis,data_only = TRUE) 

#filter dfccob
library(stringr)
dfccob$OTU <- str_extract(dfccob$taxa, "\\(.*\\)")
dfccob$OTU <- str_remove_all(dfccob$OTU, "[\\(\\)]")


f.dfccob<-filter(dfccob, OTU %in% sigtaxa) %>% separate(taxa,into=c("Kingdom","Phylum","Class","Order","Family","Genus","Species"), sep="_") %>% 
  unite(Taxon,c("Family","Genus","OTU"),sep=";")
ggplot(f.dfccob,aes(x=x,y=Taxon))+
  geom_point(size=0.5) +
  geom_pointrange(aes(xmin=xmin, xmax=xmax))+
  theme_bw()+
  facet_grid(~variable) +
  geom_vline(xintercept = 0, linetype="dashed")+
  xlab("Coefficient")+
  theme(axis.title.y=element_blank(), axis.title.x=element_text(size=6),axis.text=element_text(size=6))
ggsave("CorncobIndv.png", width=6.5, height=5, units="in", dpi=300)



######corncob A
#make body site a factor and rearrange the factor levels so that bones from the feet are the baseline.
sample_data(ASimp)$Body.Site = factor(sample_data(ASimp)$Body.Site)
sample_data(ASimp)$Body.Site = factor(sample_data(ASimp)$Body.Site,levels(sample_data(ASimp)$Body.Site)[c(2,1,3,4,5,6,7)])


set.seed(2) #from vignette tutorial; set seed and run differential abundance test controlling for the effect of individual on dispersion
da_analysisA <- differentialTest(formula = ~ Body.Site,
                                phi.formula = ~ Body.Site,
                                formula_null = ~ 1,
                                phi.formula_null = ~ Body.Site,
                                test = "Wald", boot = FALSE,
                                data = ASimp,
                                fdr_cutoff = 0.001)

da_analysisA
da_analysisA$significant_taxa

#plot and view only species and strain levels
dfccobA<-plot(da_analysisA,data_only = TRUE) 
dfccobA$OTU <- str_extract(dfccobA$taxa, "\\(.*\\)")
dfccobA$OTU <- str_remove_all(dfccobA$OTU, "[\\(\\)]")


f.dfccobA<-filter(dfccobA, OTU %in% sigtaxaA) %>% separate(taxa,into=c("Kingdom","Phylum","Class","Order","Family","Genus","Species"), sep="_") %>% 
  unite(Taxon,c("Family","Genus","OTU"),sep=";")
Acc<-ggplot(f.dfccobA,aes(x=x,y=Taxon))+
  geom_point(size=0.5) +
  geom_pointrange(aes(xmin=xmin, xmax=xmax))+
  theme_bw()+
  facet_grid(~variable) +
  geom_vline(xintercept = 0, linetype="dashed")+
  xlab("Coefficient")+
  theme(axis.title.y=element_blank(), axis.title.x=element_text(size=6),axis.text=element_text(size=6))

###########B
sample_data(BSimp)$Body.Site = factor(sample_data(BSimp)$Body.Site)
sample_data(BSimp)$Body.Site = factor(sample_data(BSimp)$Body.Site,levels(sample_data(BSimp)$Body.Site)[c(2,1,3,4,5,6,7)])


set.seed(2) #from vignette tutorial; set seed and run differential abundance test controlling for the effect of individual on dispersion
da_analysisB <- differentialTest(formula = ~ Body.Site,
                                 phi.formula = ~ Body.Site,
                                 formula_null = ~ 1,
                                 phi.formula_null = ~ Body.Site,
                                 test = "Wald", boot = FALSE,
                                 data = BSimp,
                                 fdr_cutoff = 0.001)

da_analysisB
da_analysisB$significant_taxa

#plot and view only species and strain levels
dfccobB<-plot(da_analysisB,data_only = TRUE) 
dfccobB$OTU <- str_extract(dfccobB$taxa, "\\(.*\\)")
dfccobB$OTU <- str_remove_all(dfccobB$OTU, "[\\(\\)]")


f.dfccobB<-filter(dfccobB, OTU %in% sigtaxaB) %>% separate(taxa,into=c("Kingdom","Phylum","Class","Order","Family","Genus","Species"), sep="_") %>% 
  unite(Taxon,c("Family","Genus","OTU"),sep=";")
Bcc<-ggplot(f.dfccobB,aes(x=x,y=Taxon))+
  geom_point(size=0.5) +
  geom_pointrange(aes(xmin=xmin, xmax=xmax))+
  theme_bw()+
  facet_grid(~variable) +
  geom_vline(xintercept = 0, linetype="dashed")+
  xlab("Coefficient")+
  theme(axis.title.y=element_blank(), axis.title.x=element_text(size=6),axis.text=element_text(size=6))


#######C
sample_data(CSimp)$Body.Site = factor(sample_data(CSimp)$Body.Site)
sample_data(CSimp)$Body.Site = factor(sample_data(CSimp)$Body.Site,levels(sample_data(CSimp)$Body.Site)[c(2,1,3,4,5,6,7)])


set.seed(2) #from vignette tutorial; set seed and run differential abundance test controlling for the effect of individual on dispersion
da_analysisC <- differentialTest(formula = ~ Body.Site,
                                 phi.formula = ~ Body.Site,
                                 formula_null = ~ 1,
                                 phi.formula_null = ~ Body.Site,
                                 test = "Wald", boot = FALSE,
                                 data = CSimp,
                                 fdr_cutoff = 0.001)
da_analysisC
da_analysisC$significant_taxa

#plot and view only species and strain levels
dfccobC<-plot(da_analysisC,data_only = TRUE) 
dfccobC$OTU <- str_extract(dfccobC$taxa, "\\(.*\\)")
dfccobC$OTU <- str_remove_all(dfccobC$OTU, "[\\(\\)]")


f.dfccobC<-filter(dfccobC, OTU %in% sigtaxaC) %>% separate(taxa,into=c("Kingdom","Phylum","Class","Order","Family","Genus","Species"), sep="_") %>% 
  unite(Taxon,c("Family","Genus","OTU"),sep=";")
Ccc<-ggplot(f.dfccobC,aes(x=x,y=Taxon))+
  geom_point(size=0.5) +
  geom_pointrange(aes(xmin=xmin, xmax=xmax))+
  theme_bw()+
  facet_grid(~variable) +
  geom_vline(xintercept = 0, linetype="dashed")+
  xlab("Coefficient")+
  theme(axis.title.y=element_blank(), axis.title.x=element_text(size=6),axis.text=element_text(size=6))


###################
#source tracker figure
###################

props<-read.table("sourcetracker2/stAllgenus/mixing_proportions.txt",header=TRUE, sep='\t',check.names=FALSE)
props<-props[,c(1:5)]
colnames(props)[1]<-"SampleID"
sinfo<-read.table("sourcetracker2/st_metadata_forfigs.txt", header=TRUE,comment.char = "", sep='\t',check.names=FALSE)
colnames(sinfo)[1] <- "SampleID"

props<-arrange(props, SampleID)

#gather
props2<-melt(props, value.name="SourceProportion")
#rename variable
colnames(props2)[2] <- "Source"

#merge datasets
props3<-left_join(props2, sinfo, by="SampleID")

#get mean proportion by Environemnt
sumP<-props3 %>% filter(!is.na(Individual)) %>% filter(Source != "Soil") %>% group_by(Individual, Source)%>% summarise(SourceProportion= mean(SourceProportion))%>%
  mutate(SourcePercent = (SourceProportion*100))

#mean source contribution across all buried bones
props3 %>% filter(!is.na(Individual))  %>% group_by(Source)%>% summarise(SourceProportion= mean(SourceProportion))%>%
  mutate(SourcePercent = (SourceProportion*100))

# Create a basic bar
pie = ggplot(sumP, aes(x="", y=SourcePercent, fill=Source)) + geom_bar(stat="identity", width=1,color="black")

# Convert to pie (polar coordinates) and add labels
pie = pie + coord_polar("y", start=0) + geom_text(aes(label = paste0(round(SourcePercent), "%")), position = position_stack(vjust = 0.5),color="darkgoldenrod2",
                                                        size=3) + facet_grid(Individual~.)
newpalette<-c("#CC79A7","#000000", "white")
# Add color scale (hex colors)
pie = pie + scale_fill_manual(values=newpalette) 

# Remove labels and add title
pie = pie + labs(x = NULL, y = NULL, fill = NULL, title = NULL,color="black")

# Tidy up the theme
pie = pie + theme_classic() + theme(axis.line = element_blank(),
                                    axis.text = element_blank(),
                                    axis.ticks = element_blank(),
                                    plot.title = element_text(hjust = 0.5, color = "#666666"))


ggsave("avgsourcebyindividual.png", width=3, height=5, units="in", dpi=500)


