# T-cell STING activation analysis 
# Authors: G.K,  N.K

# Clean the environment
rm(list = ls())

# Run this script second

# Perform filtering for good quality cells This script uses outputs from ~/Data/2019.07.18.
# Load libraries
library(cowplot)
library(ggsci)
library(data.table)
library(Matrix)
library(tidyverse)

# Load Data
PBMCs_DGE=readRDS("/Home Office/Data/2020.10.26._Andreas_BRB-exp/PBMCs/Analysis/PBMCs_brb_ex.rds")
cell_annotation=readRDS("/Home Office/Data/2020.10.26._Andreas_BRB-exp/PBMCs/Analysis/PBMCs_brb_anno.rds")

# Filtering Cells
# Filter Plots for mapping stats
# Filter thresholds
map_cutoff <- 0.80
min_reads <- 1000000
max_reads <- 12000000
in_cutoff <- 0.30
ex_cutoff <- 0.35
UMI_max <- 430000
UMI_min <- 220000
Gene_Min <- 14000
unmapped_cutoff <- 0.2
# Filter for doublets
########################################################
plot(density(cell_annotation$UMIs),main="UMIs downsampled")
abline(v=UMI_max,lty="dashed")
abline(v=UMI_min,lty="dashed")
########################################################
# cell_annotation=cell_annotation[order(cell_annotation$Set),]
cell_annotation$Stimulus=as.factor(cell_annotation$Stimulus)
cell_annotation$Well=as.factor(cell_annotation$Well)
head(cell_annotation)
cell_annotation$Info=cell_annotation$Stimulus

a=ggplot(cell_annotation,aes(x=TotalReadsPerCell,y=FractionMapped*100))+
  geom_point(shape=21,alpha=0.75,aes(fill=Info))+
  geom_density2d(color='black',alpha=0.75)+
  scale_x_log10()+
  geom_vline(xintercept = max_reads, linetype="dashed")+
  geom_vline(xintercept = min_reads, linetype="dashed")+
  geom_hline(yintercept = 100*map_cutoff, linetype="dashed")+
  annotation_logticks(sides="b")+
  theme_classic()+
  theme(legend.position = "none")+
  xlab("Sequenced reads")+ ylab("% Mapped")+
  ylim((min(cell_annotation$FractionMapped)*100-5),(max(cell_annotation$FractionMapped)*100+5))+
  scale_fill_rickandmorty()

b=ggplot(cell_annotation,aes(x=Fraction_ex*100,y=Fraction_in*100)) + 
  geom_point(shape=21,aes(fill=Info),alpha=0.75)+ 
  geom_density2d(color='black',alpha=0.75) + 
  geom_vline(xintercept = ex_cutoff*100, linetype="dashed") + 
  geom_hline(yintercept = 100*in_cutoff, linetype="dashed") + 
  theme_classic() + 
  theme(legend.position = "none") + 
  xlab("% Exon")+ ylab("% Intron") + 
  ylim((min(cell_annotation$Fraction_in)*100-5),(max(cell_annotation$Fraction_in)*100+5))+ 
  xlim((min(cell_annotation$Fraction_ex)*100-5),(max(cell_annotation$Fraction_ex)*100+5))+ 
  scale_fill_rickandmorty()

c=ggplot(cell_annotation,aes(x=Fraction_ex*100,y=Fraction_Unmapped*100)) +
  geom_point(shape=21,aes(fill=Info),alpha=0.75)+ 
  geom_density2d(color='black',alpha=0.75) + 
  geom_vline(xintercept = ex_cutoff*100, linetype="dashed") + 
  geom_hline(yintercept = 100*unmapped_cutoff, linetype="dashed") + 
  theme_classic()+ xlab("% Exon")+ ylab("% Unmapped") + 
  ylim((min(cell_annotation$Fraction_Unmapped)*100-5),(max(cell_annotation$Fraction_Unmapped)*100+5))+ 
  xlim((min(cell_annotation$Fraction_ex)*100-5),(max(cell_annotation$Fraction_ex)*100+5))+ 
  scale_fill_rickandmorty()

d=ggplot(cell_annotation,aes(x=TotalReadsPerCell,y=UMIs)) + 
  geom_point(shape=21,aes(fill=Info),alpha=0.75)+ 
  geom_density2d(color='black',alpha=0.75) + 
  scale_x_log10(breaks=c(1000,10000,100000,1000000))  + 
  scale_y_log10(breaks=c(100,1000,10000,100000,500000))+ 
  theme_classic() + 
  theme(legend.position = "none") + 
  xlab("Sequenced reads")+ 
  ylab("UMIs detected")   + 
  scale_fill_rickandmorty() + 
  geom_vline(xintercept = max_reads, linetype="dashed")+
  geom_vline(xintercept = min_reads, linetype="dashed")+
  annotation_logticks(sides="bl")+ 
  geom_hline(yintercept = UMI_min, linetype="dashed")+
  geom_hline(yintercept = UMI_max, linetype="dashed")

e=ggplot(cell_annotation,aes(x=Genes,y=UMIs)) + 
  geom_point(shape=21,aes(fill=Info),alpha=0.75)+ 
  geom_density2d(color='black',alpha=0.75) +
  scale_y_log10()+ 
  theme_classic() +
  xlab("Genes detected")+ 
  ylab("UMIs detected")   +
  theme(legend.position = "none") + 
  annotation_logticks(sides="b") +
  geom_hline(yintercept = UMI_min,linetype="dashed")+
  geom_hline(yintercept = UMI_max,linetype="dashed")+
  geom_vline(xintercept = Gene_Min, linetype="dashed") + 
  scale_fill_rickandmorty()

plot_grid(a,b,d,e,c, ncol = 2)

# QC pass assignment
# Filter thresholds
# Look above 

cell_annotation$QCpass <- cell_annotation$TotalReadsPerCell>min_reads & 
  cell_annotation$TotalReadsPerCell<max_reads & 
  cell_annotation$FractionMapped>map_cutoff & 
  cell_annotation$Fraction_ex>ex_cutoff & 
  cell_annotation$Fraction_in<in_cutoff &
  cell_annotation$Fraction_Unmapped<unmapped_cutoff &
  cell_annotation$Genes>=Gene_Min &
  cell_annotation$UMIs>=UMI_min&
  cell_annotation$UMIs<=UMI_max

cell_annotation$QCpass

# QC Assginment OV
head(cell_annotation)
(patient_QCsumm <- cell_annotation %>% group_by(Stimulus) %>% summarise(good_cells=sum(QCpass)))
patient_QCsumm2<- aggregate(cell_annotation$QCpass, by =list( "Stimulus"=cell_annotation$Stimulus), FUN=sum)

ggplot(patient_QCsumm2,aes(x=reorder(Stimulus, x),y=x,fill=Stimulus)) + 
  geom_bar(stat = "identity") + 
  theme_classic() +  
  scale_fill_rickandmorty() +
  xlab("Time") + ylab("Number of good cells") + 
  coord_flip() + 
  theme(legend.position="none")

# Filtering Samples
head(cell_annotation)
head(PBMCs_DGE)
# Keeping the QC passed cells 
PBMCs_DGE_filtered=PBMCs_DGE[,cell_annotation[which(cell_annotation$QCpass==T),"XC"]]
dim(PBMCs_DGE)
dim(PBMCs_DGE_filtered)
# Keep the annot. only for the cells I keep on working
cell_annotation_filtered=cell_annotation[which(cell_annotation$QCpass==T),]
dim(cell_annotation)
dim(cell_annotation_filtered)

#####################################
###### Perform GENE filtering #######
#####################################
df.filt<-data.frame(Gene=rownames(PBMCs_DGE_filtered), NUM_UMI=rowSums(PBMCs_DGE_filtered), 
                    NUM_samp=rowSums(PBMCs_DGE_filtered>1),Mean_UMI=rowMeans(PBMCs_DGE_filtered))
dim(df.filt)
# Minimum number of UMIs to keep 
umi_gene=5 # at least x UMIS per gene
sample_gene=2 # should be an integer : how many samples should have that many umi present
average_gene=0.1 # average gene expression

A<-ggplot()+
  geom_density(data=df.filt, aes(x=NUM_UMI))+
  geom_vline(xintercept = umi_gene , linetype="dashed", colour="red")+
  scale_x_log10()

B<-ggplot()+
  geom_density(data=df.filt, aes(x=NUM_samp))+
  geom_vline(xintercept =sample_gene , linetype="dashed", colour="red")+
  scale_x_log10()

C<-ggplot()+
  geom_density(data=df.filt, aes(x=Mean_UMI))+
  geom_vline(xintercept = average_gene , linetype="dashed", colour="red")+
  scale_x_log10()

plot_grid(A, B,C, nrow = 1, rel_widths = c(1,1,1))

# umi_gene=2
# sample_gene=1.55
# average_gene=0.01
df.filt=df.filt[df.filt$NUM_UMI >average_gene & df.filt$NUM_samp >=sample_gene & df.filt$Mean_UMI >=average_gene,]
dim(PBMCs_DGE_filtered)[1]
dim(df.filt)[1]

PBMCs_DGE_filtered=as.matrix(PBMCs_DGE_filtered[row.names(df.filt),])
dim(PBMCs_DGE_filtered)
dim(cell_annotation_filtered)

# Save Data
saveRDS(PBMCs_DGE_filtered,"/Home Office//Data/2020.10.26._Andreas_BRB-exp/PBMCs/Analysis/PBMCs_DGE_filtered.rds")
saveRDS(cell_annotation_filtered, "/Home Office/Data/2020.10.26._Andreas_BRB-exp/PBMCs/Analysis/PBMCs_annot_filtered.rds")

