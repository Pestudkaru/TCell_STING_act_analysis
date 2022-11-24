# T-cell STING activation analysis 
# Authors: G.K,  N.K

# Clean the environment
rm(list = ls())

# Run this script first
# Generate stats for further filtering and analysis. 
# Filtering is the next script

library(Matrix)
library(gtools)
library(cowplot)
library(ggsci)
library(data.table)
library(Matrix)
library(tidyverse)

# Read in data
samples=read.csv("/Home Office/Data/2020.10.26._Andreas_BRB-exp/T-cells/Sample_table.csv")
head(samples)
barcodes=read.csv("/Home Office/PhD/Experiments/SCRP-Seq_Barcodes/My_1st_set/1st_set.csv",header = F)
head(barcodes)
samples$Well_BC=barcodes[match(samples$Well,barcodes$V1),"V3"]
samples$Full_BC=paste0(samples$Illu_BC,samples$Well_BC)
samples$XC=samples$Well_BC
head(samples)
rownames(samples)=NULL

# Read in zUMIs output
(zDir<-"/Home Office/Data/2020.10.26._Andreas_BRB-exp/T-cells/zUMIs_output/")
(Tcell_brb=list.files(paste0(zDir,"expression/")))
(Tcell_brb=Tcell_brb[grep("counts",Tcell_brb)])
Tcell_brb=readRDS(paste0(zDir,"expression/",Tcell_brb))
gene_names=read.delim("/Home Office/Data/2020.10.26._Andreas_BRB-exp/T-cells/zUMIs_output/expression/N228.gene_names.txt")
head(gene_names)

###########################################################################################################################
# Ex In-Ex or downSampled counts
# If on average we have >60% of UMI counts in the exon mapped reads we will continue with exonic reads only.
# All reads
summary(colSums(as.matrix(Tcell_brb$umicount$exon$all))/colSums(as.matrix(Tcell_brb$umicount$inex$all)))
# Downsampled reads
summary(colSums(as.matrix(Tcell_brb$umicount$exon$downsampling$downsampled_))/colSums(as.matrix(Tcell_brb$umicount$inex$downsampling$downsampled_)))
# KEEP All In-Ex reads because nothing was downsampled and not enough exonic reads!
Tcell_brb_ex=as.data.frame(as.matrix(Tcell_brb$umicount$exon$all))
saveRDS(Tcell_brb_ex,"/Home Office/Data/2020.10.26._Andreas_BRB-exp/T-cells/Analysis/Tcell_brb_ex.rds")
###################################################################################
# STATS
###################################################################################
(stats=list.files(paste0(zDir,"stats")))
(stats=stats[grep("*.txt",stats)])

(stats=list.files(paste0(zDir,"stats")))
(stats=stats[grep("*.rds",stats)])
Stats_Tcell_brb=readRDS(paste0(zDir,"stats/",stats))
head(Stats_Tcell_brb)
(reads=list.files(zDir))
(reads=reads[grep("*barcodes_binned.txt",reads)])
readspercell=read.csv(paste0(zDir,reads),header = T)
head(readspercell)
dim(readspercell)
stats=Stats_Tcell_brb
stats=stats %>% dplyr::filter(RG %in% samples$XC)
head(stats)

# sort stats by barcode order of reads per cell
stats$TotalReadsPerCell <- readspercell[match(stats$RG, readspercell$XC), "n"]
head(stats)
head(stats[order(stats$RG),],10)

# Create Mapping stats DF - fraction of features
head(stats)
stats$FractionReads <- stats$N/stats$TotalReadsPerCell
map_stats <- stats %>% group_by(RG) %>% 
  dplyr::filter(type=="Unmapped") %>% 
  dplyr::mutate(FractionMapped=1-FractionReads) %>% 
  dplyr::select(RG,TotalReadsPerCell,FractionMapped)
head(stats[order(stats$RG),],10)
head(map_stats[order(map_stats$RG),],10)
unique(stats$type)
map_stats$Fraction_ex <- (stats %>% group_by(RG) %>% dplyr::filter(type=="Exon"))$FractionReads
head(map_stats)
map_stats$Fraction_in <- (stats %>% group_by(RG) %>% dplyr::filter(type=="Intron"))$FractionReads
head(map_stats)
map_stats$Fraction_Unmapped <- (stats %>% group_by(RG) %>% dplyr::filter(type=="Unmapped"))$FractionReads
head(map_stats)
# Making annotation DF
colnames(map_stats)[1]="XC"
cell_annotation <- plyr::join(as.data.frame(map_stats),samples,by="XC", type = "inner")
head(cell_annotation)
dim(cell_annotation)
# Getting UMI info
length(cell_annotation$XC)
length(colnames(Tcell_brb_ex))
# Kicking out samples from annot.table that didnt make it to DGE file
cell_annotation=cell_annotation[cell_annotation$XC%in%colnames(Tcell_brb_ex),]
dim(cell_annotation)
# Colsumming for umis
cell_annotation$UMIs=colSums(Tcell_brb_ex[,cell_annotation$XC])
# Getting Gene info
# Colsumming true or false 
cell_annotation$Genes=colSums(Tcell_brb_ex[,cell_annotation$XC]>0)
head(cell_annotation)
saveRDS(cell_annotation, "/Home Office/Data/2020.10.26._Andreas_BRB-exp/T-cells/Analysis/Tcell_brb_anno.rds")

