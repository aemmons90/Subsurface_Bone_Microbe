#######################
#Combining metabolomics and 16S rRNA; Only using known metabolites 
#######################
#load libraries
library(ggplot2)
library(dplyr)
library(phyloseq)
library("genefilter")
library("PMA")
library(viridis)
library(vegan)
library(PMA)
library(factoextra)
library(ade4)
library(ggsci)
library(ggrepel)

#####################################################CCA to PCA
#include 16S, metabolomics
#load metabolite data
metab<-read.csv("../Data/Metabolites/MetabolitesKnown_R.csv", header=T)
metab$Sample.name<- as.character(metab$Sample.name)
metab$Replicate<-as.character(metab$Replicate)

#average by sample
metab1<- metab %>% group_by(Sex,Individual,BodySite,Bone,ID,Sample.name, Content) %>% select_all() %>%summarise_if(is.numeric,"mean", na.rm = TRUE)

#this is to make sure the above worked. for some reason i had to quote "mean" above
metabt<- metab %>% group_by(Sex,Individual,BodySite,Bone,ID,Sample.name, Content) %>% select_all() %>% 
  summarise(mean=mean(Alanine.Sarcosine,na.rm = TRUE))
boneID<- c(as.character(metab1$ID))

#rownames should be sample names and column names are metabolites;formatting
metab1<-as.data.frame(metab1)
metab1$ID<-gsub("bone","", metab1$ID)
metab1$ID
metab1<-metab1 %>% arrange(ID)
metabf<- metab1 %>% arrange(ID)  %>% select(-Sex,-Individual, -BodySite, -Bone, -Sample.name, -mass, -Content)


###make metab table for A/B and C
metabAB<-metab1 %>% arrange(ID)  %>% filter(Individual=="A"|Individual=="B") %>%
  select(-Sex,-Individual, -BodySite, -Bone, -Sample.name, -mass, -Content)
rownames(metabAB)<-metabAB$ID
metabAB<-metabAB[,-1]
metabC<-metab1 %>% arrange(ID)  %>% filter(Individual=="C") %>%
  select(-Sex,-Individual, -BodySite, -Bone, -Sample.name, -mass, -Content)
rownames(metabC)<-metabC$ID
metabC<-metabC[,-1]

#metabm<- metab1 %>% select(ID, Individual, BodySite, Bone) 

#remove ID column
rownames(metabf)<-metabf$ID
metabf<-metabf[,-1]

#get normFGS; normalized 16S rRNA table; rarefied at 10,000 reads
normFGS<-readRDS(file = "../Data/Physeq/normFGS.rdsnormFGS.rds")
sample_data(normFGS)$SampleID <-rownames(sample_data(normFGS))

#metabolites on only 41 samples; rarefied to even depth
BBmetab<- normFGS %>% subset_samples(SampleID %in% boneID)

#combine taxa at the genus level
BBmetab1 <- BBmetab %>%
  tax_glom(taxrank = "Genus") #%>%

####perform some additional filtering
#filter features with a frequency less than 10
BBmetab1 <- prune_taxa(taxa_sums(BBmetab1) > 10, BBmetab1)

#filter table so that taxa have more than 2 reads in at least two samples 
BBmetab1<-filter_taxa(BBmetab1, flist=(kOverA(2, 2)), TRUE) 

#look at metadata
map.data<-data.frame(sample_data(BBmetab1))

#get a table for each individual
taxaAB <- BBmetab1 %>% subset_samples(Individual=="A"|Individual=="B")
taxaC <- BBmetab1 %>% subset_samples(Individual=="C")
taxaAB <- prune_taxa(taxa_sums(taxaAB) > 0, taxaAB)
taxaC <- prune_taxa(taxa_sums(taxaC) > 0, taxaC)


#get otu table of BBmetab and transpose so that sample names are the rownames
metabOTU1<- as.data.frame(t(otu_table(BBmetab1)))
row.names(metabOTU1)<-gsub("bone","", row.names(metabOTU1))
row.names(metabOTU1)
metabOTUAB<- as.data.frame(t(otu_table(taxaAB)))
row.names(metabOTUAB)<-gsub("bone","", row.names(metabOTUAB))
row.names(metabOTUAB)
metabOTUC<- as.data.frame(t(otu_table(taxaC)))
row.names(metabOTUC)<-gsub("bone","", row.names(metabOTUC))
row.names(metabOTUC)

#check dimensions
dim(metabOTU1)
dim(metabf)

#log transform metabolites
metabf <- log(1 + metabf, base=10)
metabAB <- log(1 + metabAB, base=10)
metabC <- log(1 + metabC, base=10)

#plot human DNA for just these samples
raw<-ggplot(map.data,aes(x=Bone, y=Degradation_I_Human)) +
  #ggtitle("B")+
  theme(plot.title = element_text(hjust = 0))+
  geom_point(aes(color=Individual), size=1) +
  theme_bw()+
  scale_color_viridis(discrete=TRUE)+
  facet_grid(~Body.Site, scales="free_x",space="free", labeller = labeller(Body.Site = label_wrap_gen(width = 5)))+
  theme(strip.background = element_blank(), strip.placement = "outside")+
  ylab("Human DNA ng/gbp")+
  theme(axis.title.x = element_blank(), legend.title=element_blank(), legend.position= "bottom", 
        axis.title.y = element_text(size=10))+
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5, size=10))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(axis.text.y=element_text(size=10), 
        strip.text.y=element_text(size=10,face="bold"), 
        strip.text.x=element_text(size=10))

#get mapping data
AB<-map.data%>%filter(Individual=="A"| Individual=="B")
C<-map.data%>%filter(Individual=="C")

######################
#distance matrix for known metabolites
dmat<-otu_table(metabf,taxa_are_rows = FALSE)
row.names(map.data)<-gsub("bone","", row.names(map.data))
row.names(map.data)
dmetab<-distance(dmat,"bray")
#make a physeq
metabPhy<-phyloseq(sample_data(map.data),dmat)
metabNMDS<- ordinate(
  physeq =metabPhy, 
  method = "NMDS", 
  distance = "bray",
  set.seed(4)
)
#Stress = 0.19, k=2

bc_met<-plot_ordination(
  physeq = metabPhy,
  ordination = metabNMDS,
  axes= 1:2,
  color = "Individual"
  
) + 
  theme_classic()+
  geom_point(aes(color = Individual), alpha = 0.7, size = 2.0) +
  #geom_point(colour = "grey90", size = 1.5) +
  #scale_shape_manual(values=c(15,16,17,3,18,8,6,5))+ 
  #facet_grid(~Individual)+
  stat_ellipse() +
  ggtitle("(A) NMDS Metabolites") +
  theme(legend.position = "right",axis.text=element_text(size=8), title=element_text(size=8),axis.title=element_text(size=8))+
  #geom_text(aes(label= Residence))+
  #geom_text(aes(label = ifelse(Project=="Grave" & Body.Site=="Foot", Bone, "")), show.legend=FALSE) +
  #geom_label_repel(aes(label = ifelse(Individual == "A" & Depth > 45, Depth, "")))+
  #size=2, show.legend=FALSE)+
  guides(guide_legend(override.aes = list(shape = 22)))+
  scale_color_manual(values=c("#E69F00","#56B4E9","#009E73"))

#ggsave("braytaxa.png", height=6.5, width=5.5, units="in", dpi=300)

adonis(dmetab ~ Individual*Body.Site, data= map.data,  permutations=999)
asvmat<-distance(otu_table(BBmetab1),"bray")
adonis(asvmat ~ Individual, data=map.data,  permutations=999)


mat<-as.data.frame.matrix(otu_table(BBmetab1))
mat<-matrix(mat)
map<-data.frame(sample_data(BBmetab1))


####################run CCA for known metabolites
#reuse the below code and substitute  in the correct data frame
cca_res <- CCA(metabOTUAB, metabAB, penaltyx = .15, penaltyz = .15)
print(cca_res)


combined <- cbind((metabOTUAB[,cca_res$u != 0]),
                  (metabAB[, cca_res$v != 0]))

#get taxonomy to add to graph
colnames(combined)
names<-c(colnames(combined)[1:15])

my_subset <- subset(otu_table(taxaAB), rownames(otu_table(taxaAB)) %in% names)
new_physeq <- merge_phyloseq(my_subset, tax_table(taxaAB), sample_data(taxaAB))

a<-as.list.data.frame(tax_table(new_physeq))
a<-as.data.frame((a))
a$Family<-as.character(a$Family)
names2<-rownames(a)

###for C
#a$Genus[c(5:6,10:11,15:17,19,21:23,25)]<-a$Family[c(5:6,10:11,15:17,19,21:23,25)]

#do these match
identical(names,names2)

#replace features with taxonomy
colnames(combined)[1:15]<-as.character(a$Genus)


#merge data sets directly to use all data
#make metabolite names nice
colnames(combined)[16:18]<-c("N-Acetylglutamine","N-Acetylglutamate","2-Dehydro-D-gluconate")

pca_res <- dudi.pca(combined, scannf = F, nf = 3)

#can also look solely at metabolites
A<-fviz_pca_biplot(pca_res,  repel=TRUE,
                # Individuals
                geom.ind = "point",
                fill.ind = AB$Individual,
                pointshape = 21, pointsize = 2,labelsize=2,arrowsize=0.10,
                #palette = "jco",
                invisible="quali",
                #addEllipses = TRUE,
                # Variables,
                #alpha.var = "contrib",
                col.var = "contrib", 
                gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07")
)+
  labs(fill = "Individual", color = "Contribution (%)") +# Change legend title 
  geom_point(aes(shape=AB$Body.Site),size=2.5) +
  scale_shape_manual(
    name = "Body Region",
                     #labels = c("","Arm", "Foot", "Leg", "Lower Torso","Skull", "Upper Torso"),
                     values = c(0, 1, 2, 3, 4,5,6))+
  theme_classic()+
    ggtitle("(B) PCA A and B")+
  scale_fill_manual(values=c("#E69F00","#56B4E9","#009E73"))+
  theme(axis.text=element_text(size=8), title=element_text(size=8), legend.text=element_text(size=8))+
  guides(fill = FALSE, shape= guide_legend(keywidth = 0.1, keyheight = 0.1))+
  geom_text_repel(label = ifelse(AB$Degradation_I_Human >=3,round(AB$Degradation_I_Human,digits=1), ""),
                  show.legend=FALSE, size=2)

B<-fviz_pca_biplot(pca_res,  repel=TRUE,
                   # Individuals
                   geom.ind = "point",
                   fill.ind = C$Body.Site,
                   pointshape = 21, pointsize = 2,labelsize=2,arrowsize=0.10,
                   palette = "jco",
                   invisible="quali",
                   #addEllipses = TRUE,
                   # Variables,
                   #alpha.var = "contrib",
                   col.var = "contrib", select.var = list(contrib = 20),
                   gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07")
)+
  labs(fill = "Body Region", color = "Contribution (%)") +# Change legend title 
  #geom_point(aes(shape=C$Body.Site),size=2.5) +
  scale_shape_manual(
    name = "Body Region",
    #labels = c("","Arm", "Foot", "Leg", "Lower Torso","Skull", "Upper Torso"),
    values = c(0, 1, 2, 3, 4,5,6))+
  theme_classic()+
  ggtitle("(C) PCA C")+
  #geom_text_repel(xlim=5)+
  theme(axis.text=element_text(size=8), title=element_text(size=8), legend.text=element_text(size=8))+
  guides(fill = guide_legend(keywidth = 0.1, keyheight = 0.1), shape= guide_legend(keywidth = 0.1, keyheight = 0.1))+
  geom_text_repel(label = ifelse(C$Degradation_I_Human >= 2,round(C$Degradation_I_Human,digits=1), ""),
                  show.legend=FALSE, size=2)

C2<-fviz_pca_biplot(pca_res,  repel=TRUE,
                   # Individuals
                   geom.ind = "point",
                   fill.ind = map.data$Individual,
                   pointshape = 21, pointsize = 2,labelsize=1,arrowsize=0.25,
                   palette = "jco",
                   invisible="quali",
                   #addEllipses = TRUE,
                   # Variables,
                   #alpha.var = "contrib",
                   col.var = "contrib", 
                   gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07")
)+
  labs(fill = "Individual", color = "Contribution (%)") +# Change legend title 
  geom_point(aes(shape=map.data$Body.Site),size=2.5) +
  scale_shape_manual(
    name = "Body Region",
    #labels = c("","Arm", "Foot", "Leg", "Lower Torso","Skull", "Upper Torso"),
    values = c(0, 1, 2, 3, 4,5,6))+
  ggtitle("(A) PCA A, B, and C")+
  theme(axis.text=element_text(size=8), title=element_text(size=8), legend.text=element_text(size=8))+
  guides(fill = guide_legend(keywidth = 0.1, keyheight = 0.1), shape= guide_legend(keywidth = 0.1, keyheight = 0.1))+
  geom_text_repel(label = ifelse(map.data$Degradation_I_Human<= 1.5,round(map.data$Degradation_I_Human,digits=1), ""),
                  show.legend=FALSE, size=2)


library(ggpubr)
figPCAknown<-ggarrange(C2,A, nrow=2, common.legend=TRUE, legend="right")
figPCAknown2<-ggarrange(figPCAknown,B, nrow=2,  heights=c(1,0.5), legend=
                         'right')
ggsave("knownPCAbiplot2.0.png", units="in", width =7, height = 6,  dpi=400, device="png")

library(cowplot)
plot_grid(bc_met,A,B, align = "hv",axis=c("r","b"),ncol=1,
          rel_heights=c(1,1,1),rel_widths=c(0.25,1,1))
ggsave("knownPCAbiplot2.png", units="in", width = 6, height =8.5,  dpi=300, device="png")


# Extract shape legend. Returns a gtable
q<-ggplot(map.data, aes(x=Human..ng.gbp., y=Bacteria..copies.gbp., shape=Body.Site))+
  geom_point()+
  theme_classic()+
  scale_shape_manual(name = "Body Region",
    labels = c("Arm", "Foot", "Head","Leg", "Lower Torso", "Upper Torso"),
    values = c(1, 2, 3, 4,5,6))+
  theme(legend.position="bottom", text=element_text(size=8))+
  guides(
         shape = guide_legend(keywidth = 0.5, keyheight = 0.5, nrow=1))

leg <- get_legend(q)

# Convert to a ggplot and print
C<-as_ggplot(leg)

