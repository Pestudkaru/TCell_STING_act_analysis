# T-cell STING activation analysis 
# Authors: G.K,  N.K

# Clean the environment
rm(list = ls())

# Run this script third
# T-cells low input bulk rna-seq analysis:
options(stringsAsFactors = F)
library(cowplot)
library(data.table)
library(gplots)
library(Hmisc)
library(mclust)
library(pheatmap)
library(sctransform)
library(Seurat)
library(tidyverse)
library(ggrepel)
library(Rsamtools)
library(GenomicAlignments)
library(BiocParallel)
library(DESeq2)
library(apeglm)
library(gtools)
library(IHW)
library(ashr)
library(vsn)
library(pheatmap)
library(RColorBrewer)
library(plotly)
library(extrafont)
library(ComplexHeatmap)
library(biomaRt)
library(GO.db)
library(circlize)

# Make some colors
col_green <- "#2CAB44"
col_grey <- "#808080"
col_blue <- "#2A6AAA"
col_purple <- "#803989"
col_red <- "#9E0003"
my_cc <- c(col_green, col_purple, col_blue)
cc1 <- c(brewer.pal(name = "YlOrRd",n=9))[6]#red
cc2 <- c(brewer.pal(name = "YlOrBr",n=9))[5]#or
cc3 <- c(brewer.pal(name = "YlOrRd",n=9))[3]#yellow
cc4 <- c(brewer.pal(name = "YlGn",n=9))[6]#green
cc5 <- c(brewer.pal(name = "YlGnBu",n=9))[5]#lightblue
cc6 <- c(brewer.pal(name = "Blues",n=9))[6]#blue
cc7 <- c(brewer.pal(name = "Purples",n=9))[6]#viol
cc8 <- c(brewer.pal(name = "RdPu",n=9))[7]#purp
cc=c(brewer.pal(name = "Set1",n=9))


# Filtering and stats
Tcell_DGE=readRDS("/Users/niklaskuhl/Library/Mobile Documents/com~apple~CloudDocs/Documents/Promotion/Revision/RNAseq/2020.10.26._Andreas_BRB-exp/T-cells/Analysis/Tcell_DGE_filtered.rds")
head(Tcell_DGE)
cell_annotation=readRDS("/Users/niklaskuhl/Library/Mobile Documents/com~apple~CloudDocs/Documents/Promotion/Revision/RNAseq/2020.10.26._Andreas_BRB-exp/T-cells/Analysis/Tcell_annot_filtered.rds")
#######
comp_gene=read.delim("/Users/niklaskuhl/Library/Mobile Documents/com~apple~CloudDocs/Documents/Promotion/Revision/RNAseq/Complete_gene_names_ENS.tsv",col.names =c("gene_id","gene_name","description") )
# comp_gene[comp_gene$gene_name=="CSF2",]
# comp_gene[grepl("CSF",comp_gene$gene_name),]
# Throw away ribos and MT
# comp_gene=comp_gene[!grepl("^MT-",comp_gene$gene_name),]
# comp_gene=comp_gene[!grepl("ribosomal",comp_gene$description),]
comp_gene=comp_gene[comp_gene$gene_id%in%rownames(Tcell_DGE),]
comp_gene=comp_gene[!duplicated(comp_gene$gene_name),]
nrow(Tcell_DGE)
Tcell_DGE=Tcell_DGE[rownames(Tcell_DGE)%in%comp_gene$gene_id,]
nrow(Tcell_DGE)
rownames(Tcell_DGE)=comp_gene$gene_name[match(rownames(Tcell_DGE),comp_gene$gene_id)]

#######################################################
# tcell_dots <- CreateSeuratObject(counts = Tcell_DGE, project = "tcell_dots cGAMP", min.cells = 3, min.features = 200)
# tcell_dots@meta.data<- cbind(tcell_dots@meta.data,cell_annotation[match(colnames(Tcell_DGE),cell_annotation$XC),])
# 
# tcell_dots[["percent.mt"]] <- PercentageFeatureSet(tcell_dots, pattern = "^MT-")
# VlnPlot(tcell_dots, features = c("percent.mt","nFeature_RNA", "nCount_RNA"), ncol = 3)
# 
# FeatureScatter(tcell_dots, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
# tcell_dots <- NormalizeData(tcell_dots, normalization.method = "LogNormalize", scale.factor = 10000)
# tcell_dots <- FindVariableFeatures(tcell_dots, selection.method = "vst", nfeatures = 2000)
# # Identify the 10 most highly variable genes
# top10 <- head(VariableFeatures(tcell_dots), 20)
# # plot variable features with and without labels
# plot1 <- VariableFeaturePlot(tcell_dots)
# plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
# plot2
# #  Scaling
# all.genes <- rownames(tcell_dots)
# tcell_dots <- ScaleData(tcell_dots, features = all.genes)
# tcell_dots <- RunPCA(tcell_dots, features = VariableFeatures(object = tcell_dots),npcs = 23)
# # DimPlot(tcell_dots, reduction = "pca",group.by = "Stimulus",pt.size = 2)
# #
# # tcell_dots <- ScaleData(tcell_dots, vars.to.regress = "UMIs")
# # tcell_dots <- RunPCA(tcell_dots, features = VariableFeatures(object = tcell_dots),npcs = 15)
# # DimPlot(tcell_dots, reduction = "pca",group.by = "Stimulus",pt.size = 2)
# DimHeatmap(tcell_dots, dims = 1, cells = 500, balanced = TRUE)
# tcell_dots <- FindNeighbors(tcell_dots, dims = 1:14,k.param = 14)
# tcell_dots <- FindClusters(tcell_dots, resolution = 1)
# Idents(tcell_dots)
# # tcell_dots <- RunUMAP(tcell_dots, dims = 1:14,n.neighbors = 2)
# # tcell_dots <- RunTSNE(tcell_dots, dims = 1:14,n.neighbors = 2,perplexity=5)
# # DimPlot(tcell_dots, reduction = "umap")
# # DimPlot(tcell_dots, reduction = "tsne")
# p2=DimPlot(tcell_dots, reduction = "pca",pt.size = 2)
# p1=DimPlot(tcell_dots, reduction = "pca",group.by = "Stimulus",pt.size = 2)
# plot_grid(p1,p2)
# p1
# pca_data=as.data.frame(as.matrix(tcell_dots@reductions$pca@cell.embeddings))
# head(pca_data)
# meta=tcell_dots@meta.data
# head(meta)
# pca_data$Stim=meta$Stimulus
# pca_data$KO=meta$Type
# pca_data$Sample=meta$Well_BC
# ggplot(pca_data,aes(x=PC_1,y = PC_2,color=Stim,shape=KO))+geom_point(size=2)+theme_minimal()
# 
# 
# # HEAT
# (topgenes=VariableFeatures(tcell_dots)[1:100])
# # Lets first make a heatmap with top 2000 HVG
# tcell_dotss_data=as.data.frame(as.matrix(tcell_dots@assays$RNA@data))
# tcell_dotss_data[1:5,1:5]
# colnames(tcell_dotss_data)=tcell_dots@meta.data$Well
# heat_data=tcell_dotss_data[topgenes,]
# # (need=names(sort(rowSums(heat_data),decreasing = T))[1:5])
# # heat_data=heat_data[!rownames(heat_data)%in%need,]
# heat_data=na.omit(heat_data)
# max(heat_data)
# colSums(heat_data)
# meta=tcell_dots@meta.data
# mat_col=data.frame(Genotype = meta$Type, Stimulus=meta$Stimulus)
# meta$sample=paste(meta$Type,meta$Stimulus,meta$Well,sep="_")
# rownames(mat_col)=meta$sample
# colnames(heat_data)=meta$sample
# 
# pheatmap(heat_data, border_color=NA, show_colnames=T, show_rownames=T, annotation_col=mat_col,
#          drop_levels=T, main = "Heatmap top 100 HVGs")
# 
# # ## Sample distance matrix (Eucledian)
# # library(PoiClaClu)
# # meta=tcell_dots@meta.data
# # meta$sample=paste0(meta$Well,"_",meta$Stimulus)
# normdata=as.data.frame(as.matrix(tcell_dots@assays$RNA@data))
# normdata2=normdata[VariableFeatures(tcell_dots),]
# # poisd=PoissonDistance(t(normdata2))
# # samplePoisDistMatrix=as.matrix(poisd$dd)
# # rownames(samplePoisDistMatrix)=paste(meta$sample)
# # colnames(samplePoisDistMatrix)=NULL
# # c3=colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
# # pheatmap(samplePoisDistMatrix,
# #          clustering_distance_rows = poisd$dd,
# #          clustering_distance_cols = poisd$dd,
# #          col = c3)

###############################
#### DE does not work with so few samples!
# Have to use DEseq2
# tcell_dots@meta.data$XC==colnames(tcell_dots)
# tcell_dots@meta.data$Stimulus
# tcell_dots@meta.data$contour=c(0,1,0,0,0,2,0,0,0,1,1,0,2,0,0,0)
# meta=tcell_dots@meta.data
# Idents(tcell_dots)=tcell_dots@meta.data$contour
# cluster.markers <- FindMarkers(tcell_dots, ident.1 = 1, ident.2 = 2, min.pct = 0.25)

"################# DESEQ2 #######################"
cell_annotation$Stimulus=factor(cell_annotation$Stimulus,levels = unique(cell_annotation$Stimulus)[c(3,1,2,4)])
cell_annotation$Type=factor(cell_annotation$Type,levels = rev(unique(cell_annotation$Type)))

cell_annotation$Name=paste(cell_annotation$Stimulus,cell_annotation$Type,cell_annotation$Well,sep="_")
###
colnames(Tcell_DGE)=cell_annotation$Name[match(colnames(Tcell_DGE),cell_annotation$XC)]
rownames(cell_annotation)=cell_annotation$Name


write.table(Tcell_DGE,file="/Users/niklaskuhl/Library/Mobile Documents/com~apple~CloudDocs/Documents/Promotion/Revision/RNAseq/filtered raw read counts.csv", row.names = TRUE) # keeps the rownames

###
dds <- DESeqDataSetFromMatrix(countData = Tcell_DGE, colData = cell_annotation ,design = ~ Stimulus + Type)
##############################################################
# dds$Treatment <- relevel(dds$Treatment, ref = "Unstimulated")
dds$Stimulus <- relevel(dds$Stimulus, ref = "Unstimulated")
# dds$Stimulus
# Pre-filtering
keep <- rowSums(counts(dds)) >= 10
paste("number of genes kept -",sum(keep))
dds <- dds[keep,]
# Doing the DESEQ
dds <- DESeq(dds)
##################################
norm=counts(dds,normalized=T)
unnorm=counts(dds,normalized=F)
par(mfrow=c(1,2))
barplot(colSums(Tcell_DGE))
barplot(colSums(nii))
########################################################################################
resultsNames(dds)
# resLFC <- lfcShrink(dds, coef="Treatment_IL.16_vs_Unstimulated", type="apeglm")
# resLFC <- lfcShrink(dds, coef="Stimulus_IL.16.500.ng.ml_vs_Unstimulated", type="apeglm")
# resOrdered <- resLFC[order(resLFC$pvalue),]
# summary(resLFC)
# sum(resLFC$padj < 0.1, na.rm=TRUE)
# res05 <- results(dds, alpha=0.05)
# summary(res05)
# # Independent hypothesis weighting
# resIHW <- results(dds, filterFun=ihw)
# summary(resIHW)
# sum(resIHW$padj < 0.01, na.rm=TRUE)
# metadata(resIHW)$ihwResult
# # LOOKING AT RESULTS
# plotMA(resLFC, ylim=c(-8,8))
# resNorm <- lfcShrink(dds, coef=2, type="normal")
# resAsh <- lfcShrink(dds, coef=2, type="ashr")
# par(mfrow=c(1,3), mar=c(4,4,2,1))
# xlim <- c(1,1e5); ylim <- c(-8,8)
# plotMA(resLFC, xlim=xlim, ylim=ylim, main="apeglm")
# plotMA(resNorm, xlim=xlim, ylim=ylim, main="normal")
# plotMA(resAsh, xlim=xlim, ylim=ylim, main="ashr")
# # datz weird
# par(mfrow=c(1,1))
# plotCounts(dds, gene=which.min(resLFC$padj), intgroup="Stimulus")
## Data transformations and visualization
rld <- rlog(dds, blind=FALSE)
head(assay(rld), 3)
meanSdPlot(assay(rld))

deseq=assay(rld)
write.table(deseq,file="/Users/niklaskuhl/Library/Mobile Documents/com~apple~CloudDocs/Documents/Promotion/Revision/RNAseq/normalized log transformed read counts.csv", row.names = TRUE) # keeps the rownames




norm_data=assay(rld)
pheatmap(cor(norm_data,method = "pearson"),color = rev(redgreen(100)))

# Heatmap
#highlighted_genes <- as.data.frame(nii) %>% filter(row.names(nii) %in% c("IFIT1"))
col_order <- c(15,3,6,17,9,11, 20,24,23,10,12,22, 4,2,18,14,7,16, 19,5,21,1,8,13)
select <- order(rowVars(assay(rld,normalized=TRUE)),decreasing=TRUE)[1:100]
df <- as.data.frame(colData(rld)[,c("Stimulus","Type")])
nii=assay(rld)[select,]
rownames(nii)
nii=nii[order(rowMeans(nii))[1:100],]
nii=nii[,col_order]
# write.csv(nii, "/Users/niklaskuhl/Downloads/RNAseq/dataforjochen.csv")



col_fun = colorRamp2(c(-2,0,4), c(col_blue, "white", col_red)) #for color gradient



#ISGs
de_genes_isg <- read.csv("/Users/niklaskuhl/Library/Mobile Documents/com~apple~CloudDocs/Documents/Promotion/Revision/RNAseq/de_genes_isg.csv", header = T) #import isg gene list
de_genes_isg <- as.vector(t(de_genes_isg)) #convert to vector 
# de_genes_isg_highlight <- read.csv("/Users/niklaskuhl/Library/Mobile Documents/com~apple~CloudDocs/Documents/Promotion/Revision/RNAseq/de_genes_isg_highlight.csv", header = T)
# de_genes_isg_highlight <- as.vector(t(de_genes_isg_highlight))

labels <- rownames(nii)[rownames(nii) %in% de_genes_isg] #extract isg genes from dge genes for labels 
genes_isg <- nii[rownames(nii) %in% de_genes_isg,]
labels_highlight <- rownames(genes_isg)[rownames(genes_isg) %in% doubleISGs]

# ifna <- read.csv("/Users/niklaskuhl/Library/Mobile Documents/com~apple~CloudDocs/Documents/Promotion/Revision/RNAseq/ifna.csv", header = F)
# ifna <- as.vector(t(ifna))
# ifna <- toupper(ifna)

###########NFKB

# ensembl = useMart("ensembl",dataset="hsapiens_gene_ensembl") #uses human ensembl annotations
# 
# nfkb <- getBM(
#   attributes= "hgnc_symbol",
#   filters="go",
#   values=GOBPOFFSPRING[["GO:0038061"]], 
#   mart=ensembl)
# nfkb <- as.vector(t(nfkb)) #convert to vector 
# nfkb <- nii[rownames(nii) %in% nfkb,]
# 
# Heatmap(nfkb, #heatmap
#         show_row_names=T,
#         cluster_columns = F,
#         show_heatmap_legend = T,
#         show_column_names = T,
#         height = nrow(nfkb)*unit(7, "pt"),
#         width = ncol(nfkb)*unit(7, "pt"),
#         col = colorRamp2(c(2,6,10), c(col_blue, "white", col_red)),
#         row_dend_reorder = F,
#         row_dend_width = unit(0.5, "cm"),
#         row_names_gp = gpar(fontsize = 7),
#         column_names_gp = gpar(fontsize = 7)
#         # clustering_distance_rows = robust_dist
#         #column_split = rep(c("A","B","c","d"), 6), 
# )


###############


genes_of_interest <- nii[rownames(nii_scale) %in% ifna,]
Heatmap(genes_of_interest, 
        cluster_rows = F,
        cluster_columns = F)



group = kmeans(t(nii), centers = 4)$cluster

Heatmap(nii, 
        show_row_names=T, 
        cluster_columns = F,
        show_heatmap_legend = T,
        show_column_names = T,
        # column_km = 4,
        column_names_gp = gpar(fontsize = 6),
        height = nrow(nii)*unit(5, "pt"),
        width = ncol(nii)*unit(1.2, "mm"),
        col = colorRamp2(c(2,6,10), c(col_blue, "white", col_red)),
        # # row_dend_reorder = F,
        # row_dend_gp = gpar(lwd = .3),
        # row_dend_width = unit(0.5, "cm"),
        # column_dend_gp = gpar(lwd = .5),
        row_names_gp = gpar(fontsize = 5),
       # clustering_distance_rows = robust_dist
        #column_split = rep(c("A","B","c","d"), 6), 
        ) #+
  # rowAnnotation(link = anno_mark(at = which(rownames(nii_scale) %in% de_genes_isg), 
  #                                labels = labels, 
  #                                labels_gp = gpar(fontsize= 5),
  #                                link_width = unit(10,"mm"), 
  #                                link_gp = gpar(lwd = .3)), 
  #               width = unit(.1, "mm") 
  #               + max_text_width(labels))

Heatmap(genes_isg, #heatmap
        show_row_names=T,
        cluster_columns = F,
        show_heatmap_legend = T,
        show_column_names = F,
        height = nrow(genes_isg)*unit(3, "pt"),
        width = ncol(genes_isg)*unit(1.56, "mm"),
        col = colorRamp2(c(2,6,10), c(col_blue, "white", col_red)),
        row_dend_reorder = T,
        row_dend_gp = gpar(lwd = .5),
        row_dend_width = unit(0.5, "cm"),
        row_names_gp = gpar(fontsize = 5),
        heatmap_legend_param = list(
          title = "", labels_gp = gpar(fontsize = 7), legend_height = unit(0.3, "cm"), grid_width = unit(3, "mm"), direction = "horizontal"
        ),
        # bottom_annotation = HeatmapAnnotation(boxplot = anno_link(align_to = c(list(1:3, 4:6, 7:9, 10:12, 13:15, 16:18, 19:21, 22:24)), 
        #                                                           which = "column", 
        #                                                           panel_fun = 
        #                                                             function(index, nm) {
        #                                                               pushViewport(viewport(
        #                                                                 yscale = c(8,2), 
        #                                                                 xscale = c(0.9,1.1), 
        #                                                                 width = unit(15, "pt"), 
        #                                                                 y = unit(12,"pt"), 
        #                                                                 just = "centre"))
        #                                                               grid.boxplot(
        #                                                                 genes_isg[,index], 
        #                                                                 pos = 1, 
        #                                                                 direction = "vertical", 
        #                                                                 outline = FALSE, 
        #                                                                 box_width = 0.06, 
        #                                                                 size= unit(4,"pt"), 
        #                                                                 pch = 20, 
        #                                                                 gp = gpar(lwd = .5, fill = rep(c("white", "white", col_green, "white", "white", "white", col_green, "white"), each = 3)[index])
        #                                                                 )
        #                                                               popViewport()
        #                                                             }, 
        #                                                           gap = unit(0, "pt"), height = unit(40, "pt"), link_gp = gpar(lty = 0)),
        #                                       boxplot2 = anno_link(align_to = c(list(22:24)), which = "column", panel_fun = 
        #                                                              function(index, nm) {
        #                                                                pushViewport(viewport(yscale = c(8,2), xscale = c(0.9,1.1), width = unit(15, "pt"), y = unit(53,"pt"), x = unit(5,"pt"),  just = "centre"))
        #                                                                grid.yaxis(main = TRUE, gp = gpar(lwd = 0.7, fontsize = 7))
        #                                                                popViewport()
        #                                                              }, 
        #                                                            gap = unit(0, "pt"), height = unit(40, "pt"), link_gp = gpar(lty = 0)),
        #                                       gap = unit(0, "mm")),
        #clustering_distance_rows = robust_dist
        #column_split = rep(c("A","B","c","d"), 6), 
 ) +
  rowAnnotation(link = anno_mark(at = which(rownames(genes_isg) %in% doubleISGs),
                                 labels = labels_highlight,
                                 labels_gp = gpar(fontsize= 7),
                                 link_width = unit(7,"mm"),
                                 link_gp = gpar(lwd = .5)),
                width = unit(.1, "mm")
                + max_text_width(labels))

write.csv(genes_isg, "/Users/niklaskuhl/Library/Mobile Documents/com~apple~CloudDocs/Documents/Promotion/Revision/RNAseq/isgs.csv")
#APOPTOSIS
# apoptotic_genes <- read.csv("/Users/niklaskuhl/Documents/Promotion/Revision/RNAseq/apoptotic_genes_go.csv", header = F) #import
# apoptotic_genes <- as.vector(t(apoptotic_genes)) #transform to vector
# 
# labels_apop <- rownames(nii)[rownames(nii) %in% apoptotic_genes] #define labels
# genes_apop <- nii[rownames(nii) %in% apoptotic_genes,]
# 
# Heatmap(nii_scale, #heatmap
#         show_row_names=F,
#         cluster_columns = F,
#         show_heatmap_legend = T,
#         show_column_names = F,
#         height = nrow(nii_scale)*unit(.4, "mm"),
#         width = ncol(nii_scale)*unit(2, "mm"),
#         col = colorRamp2(c(2,6,10), c(col_blue, "white", col_red)),
#         row_dend_reorder = F,
#         row_dend_gp = gpar(lwd = .3),
#         row_dend_width = unit(0.5, "cm")
#         # clustering_distance_rows = robust_dist
#         #column_split = rep(c("A","B","c","d"), 6), 
#           ) +
#   rowAnnotation(link = anno_mark(at = which(rownames(nii_scale) %in% apoptotic_genes), 
#                                  labels = labels_apop, 
#                                  labels_gp = gpar(fontsize= 5),
#                                  link_width = unit(10,"mm"), 
#                                  link_gp = gpar(lwd = .3)), 
#                 width = unit(.1, "mm") 
#                 + max_text_width(labels))
# 
# Heatmap(genes_apop, #heatmap
#         show_row_names=T,
#         cluster_columns = F,
#         show_heatmap_legend = T,
#         show_column_names = F,
#         height = nrow(genes_apop)*unit(1.8, "mm"),
#         width = ncol(genes_apop)*unit(1.2, "mm"),
#         col = colorRamp2(c(2,6,10), c(col_blue, "white", col_red)),
#         row_dend_reorder = F,
#         row_dend_gp = gpar(lwd = .3),
#         row_dend_width = unit(0.5, "cm"),
#         row_names_gp = gpar(fontsize = 5)
#         # clustering_distance_rows = robust_dist
#         #column_split = rep(c("A","B","c","d"), 6), 
# )
# 
# #GLYCOLYSIS
# glyco_genes <- read.csv("/Users/niklaskuhl/Downloads/RNAseq/genes_glycolysis.csv", header = F) #import
# glyco_genes <- as.vector(t(glyco_genes)) #transform to vector
# 
# labels_glyco <- rownames(nii_scale)[rownames(nii_scale) %in% glyco_genes] #define labels
# genes_glyco <- nii_scale[rownames(nii_scale) %in% glyco_genes,]
# 
# Heatmap(nii_scale, #heatmap
#         show_row_names=F,
#         cluster_columns = F,
#         show_heatmap_legend = T,
#         show_column_names = F,
#         height = nrow(nii_scale)*unit(.4, "mm"),
#         width = ncol(nii_scale)*unit(2, "mm"),
#         col = col_fun,
#         row_dend_reorder = F,
#         row_dend_gp = gpar(lwd = .3),
#         row_dend_width = unit(0.5, "cm")
#         # clustering_distance_rows = robust_dist
#         #column_split = rep(c("A","B","c","d"), 6), 
# ) +
#   rowAnnotation(link = anno_mark(at = which(rownames(nii_scale) %in% apoptotic_genes), 
#                                  labels = labels_apop, 
#                                  labels_gp = gpar(fontsize= 5),
#                                  link_width = unit(10,"mm"), 
#                                  link_gp = gpar(lwd = .3)), 
#                 width = unit(.1, "mm") 
#                 + max_text_width(labels))
# 
# Heatmap(genes_glyco, #heatmap
#         show_row_names=T,
#         cluster_columns = F,
#         show_heatmap_legend = T,
#         show_column_names = F,
#         height = nrow(genes_glyco)*unit(1.8, "mm"),
#         width = ncol(genes_glyco)*unit(1.2, "mm"),
#         col = col_fun,
#         row_dend_reorder = F,
#          row_dend_width = unit(0.5, "cm"),
#         row_names_gp = gpar(fontsize = 5)
#         # clustering_distance_rows = robust_dist
#         #column_split = rep(c("A","B","c","d"), 6), 
# )
# 
# #MITO
# mito_genes <- read.csv("/Users/niklaskuhl/Downloads/RNAseq/genes_mito.csv", header = F) #import
# mito_genes <- as.vector(t(mito_genes)) #transform to vector
# 
# labels_mito <- rownames(nii_scale)[rownames(nii_scale) %in% mito_genes] #define labels
# genes_mito <- nii_scale[rownames(nii_scale) %in% mito_genes,]
# 
# Heatmap(nii_scale, #heatmap
#         show_row_names=F,
#         cluster_columns = F,
#         show_heatmap_legend = T,
#         show_column_names = F,
#         height = nrow(nii_scale)*unit(.4, "mm"),
#         width = ncol(nii_scale)*unit(2, "mm"),
#         col = col_fun,
#         row_dend_reorder = F,
#         row_dend_gp = gpar(lwd = .3),
#         row_dend_width = unit(0.5, "cm")
#         # clustering_distance_rows = robust_dist
#         #column_split = rep(c("A","B","c","d"), 6), 
# ) +
#   rowAnnotation(link = anno_mark(at = which(rownames(nii_scale) %in% apoptotic_genes), 
#                                  labels = labels_apop, 
#                                  labels_gp = gpar(fontsize= 5),
#                                  link_width = unit(10,"mm"), 
#                                  link_gp = gpar(lwd = .3)), 
#                 width = unit(.1, "mm") 
#                 + max_text_width(labels))
# 
# Heatmap(genes_mito, #heatmap
#         show_row_names=T,
#         cluster_columns = F,
#         show_heatmap_legend = T,
#         show_column_names = F,
#         height = nrow(genes_mito)*unit(1.8, "mm"),
#         width = ncol(genes_mito)*unit(1.2, "mm"),
#         col = col_fun,
#         row_dend_reorder = F,
#         row_dend_gp = gpar(lwd = .3),
#         row_dend_width = unit(0.5, "cm"),
#         row_names_gp = gpar(fontsize = 5)
#         # clustering_distance_rows = robust_dist
#         #column_split = rep(c("A","B","c","d"), 6), 
# )
# 
# #PROLIFERATION
# pro_genes <- read.csv("/Users/niklaskuhl/Downloads/RNAseq/genes_proliferation.csv", header = F) #import
# pro_genes <- as.vector(t(pro_genes)) #transform to vector
# 
# labels_pro <- rownames(nii_scale)[rownames(nii_scale) %in% pro_genes] #define labels
# genes_pro <- nii_scale[rownames(nii_scale) %in% pro_genes,]
# 
# Heatmap(nii_scale, #heatmap
#         show_row_names=F,
#         cluster_columns = F,
#         show_heatmap_legend = T,
#         show_column_names = F,
#         height = nrow(nii_scale)*unit(.4, "mm"),
#         width = ncol(nii_scale)*unit(2, "mm"),
#         col = col_fun,
#         row_dend_reorder = F,
#         row_dend_gp = gpar(lwd = .3),
#         row_dend_width = unit(0.5, "cm")
#         # clustering_distance_rows = robust_dist
#         #column_split = rep(c("A","B","c","d"), 6), 
# ) +
#   rowAnnotation(link = anno_mark(at = which(rownames(nii_scale) %in% apoptotic_genes), 
#                                  labels = labels_apop, 
#                                  labels_gp = gpar(fontsize= 5),
#                                  link_width = unit(10,"mm"), 
#                                  link_gp = gpar(lwd = .3)), 
#                 width = unit(.1, "mm") 
#                 + max_text_width(labels))
# 
# Heatmap(genes_pro, #heatmap
#         show_row_names=T,
#         cluster_columns = F,
#         show_heatmap_legend = T,
#         show_column_names = F,
#         height = nrow(genes_pro)*unit(1.8, "mm"),
#         width = ncol(genes_pro)*unit(1.2, "mm"),
#         col = col_fun,
#         row_dend_reorder = F,
#         row_dend_gp = gpar(lwd = .3),
#         row_dend_width = unit(0.5, "cm"),
#         row_names_gp = gpar(fontsize = 5)
#         # clustering_distance_rows = robust_dist
#         #column_split = rep(c("A","B","c","d"), 6), 
# )
hist(nii)
pdf(heatmap1)
de_genes <- row.names(nii)
de_genes <- as.factor(de_genes)
de_genes
write.csv(de_genes, "/Users/niklaskuhl/Downloads/RNAseq/de_genes.csv", row.names = F, quote = F)
# PCA !!!
pcaData <- plotPCA(rld, intgroup=c("Stimulus"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
pcaData$KO=cell_annotation$Type[match(pcaData$name,cell_annotation$Name)]
# colnames(pcaData)[6]="ID"
p1=ggplot(pcaData, aes(PC1, PC2, color=Stimulus,shape=KO)) +
  geom_point(size=3) + theme_light()+ggtitle("PCA colored by Stimulus")+
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  coord_fixed()+scale_color_manual(values = cc)
ggplotly(p1)

pcaData$KO2="KO"
pcaData$KO2[pcaData$KO=="WT"]="WT"

set.seed(1)
p1=ggplot(pcaData, aes(PC1, PC2, color=Stimulus,label=KO2)) +
  geom_point(size=6) + theme_classic()+ggtitle("PCA colored by Stimulus")+
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  coord_fixed()+scale_color_manual(values = cc)+geom_text(color="white",size=2, aes(fontface = "bold"))

set.seed(1)
p2=ggplot(pcaData, aes(PC1, PC2, color=Stimulus,label=KO2)) +
  geom_point(size=0.1,alpha=0) + theme_classic()+ggtitle("PCA colored by Stimulus")+
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  coord_fixed()+scale_color_manual(values = cc)+geom_text(size=4,hjust = .5,position = position_jitter(0.5,0.5), 
                                                      aes(color=Stimulus,fontface = "bold"))

set.seed(1)
p3=ggplot(pcaData, aes(PC1, PC2, color=Stimulus,label=KO2)) +
  geom_point(size=1.5,alpha=1) + theme_classic()+ggtitle("PCA colored by Stimulus")+
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  coord_fixed()+scale_color_manual(values = cc)+geom_text(size=4,hjust = 0.5,position = position_jitter(0.5,0.5), 
                                                          aes(color=Stimulus,fontface = "bold"))
plot_grid(p1,p2,p3)


library(mclust)
set.seed(123)
clust=kmeans(pcaData[,1:2],centers = 4)
pcaData$Kmeans=as.character(clust$cluster)
ggplot(pcaData, aes(PC1, PC2, color=Kmeans,shape=KO,label=Stimulus)) +
  geom_point(size=3) + theme_light()+
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  coord_fixed()+scale_color_manual(values = cc)#[c(1,2,4,3)])

pcaData$Cluster="Unstim."
pcaData$Cluster[pcaData$Kmeans==4]="cGAMP"
pcaData$Cluster[pcaData$Kmeans==2]="cGAMP_CD3/CD28"
pcaData$Cluster[pcaData$Kmeans==3]="CD3/CD28"
(levels=unique(pcaData$Cluster)[c(2,3,1,4)])
pcaData$Cluster=factor(pcaData$Cluster,levels = levels)
p2=ggplot(pcaData, aes(PC1, PC2, color=Cluster,shape=KO,label=Stimulus)) +
  geom_point(size=3) + theme_light()+ggtitle("PCA colored Kmeans cluster")+
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  coord_fixed()+scale_color_manual(values = cc)#[c(2,1,3,4)])

ggplotly(p2)
plot_grid(p1,p2)

#####################################################
cell_annotation$Kmeans=pcaData$Kmeans[match(cell_annotation$Name,pcaData$name)]
cell_annotation$Clusters=pcaData$Cluster[match(cell_annotation$Name,pcaData$name)]

cell_annotation$Group=paste(cell_annotation$Stimulus,cell_annotation$Type,sep = "_")
cell_annotation$Group=factor(cell_annotation$Group)
head(cell_annotation$Group)
(1:8)[cell_annotation$Group]
pcaData$group=cell_annotation$Group
gp1=ggplot(pcaData, aes(PC1, PC2, color=group,shape=KO,label=Stimulus)) +
  geom_point(size=3) + theme_light()+ggtitle("PCA colored by unique group")+
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  coord_fixed()+scale_color_manual(values = cc)#[c(2,1,3,4)])

ggplotly(gp1)

########################################################################################
########################################################################################
dds <- DESeqDataSetFromMatrix(countData = Tcell_DGE, colData = cell_annotation ,design = ~ Group)
dds$Group <- relevel(dds$Group, ref = "Unstimulated_WT")
dds <- DESeq(dds)
########################################################################################
########################################################################################
head(cell_annotation)
dds <- DESeqDataSetFromMatrix(countData = Tcell_DGE, colData = cell_annotation ,design = ~ Clusters)
dds$Group <- relevel(dds$Clusters, ref = "Unstim.")
dds <- DESeq(dds)
########################################################################################
########################################################################################
#  DE and volcano! 
res <- results(dds,contrast = c("Group","Unstimulated_STING_KO","Unstimulated_WT"))
res <- results(dds,contrast = c("Group","CD3/CD28_STING_KO","CD3/CD28_WT"))
res <- results(dds,contrast = c("Group","CD3/CD28+cGAMP_STING_KO","CD3/CD28+cGAMP_WT"))
res <- results(dds,contrast = c("Group","cGAMP_STING_KO","cGAMP_WT"))


########################################################################################

res <- results(dds,contrast = c("Clusters","CD3.CD28","Unstim."))
res <- results(dds,contrast = c("Clusters","cGAMP_CD3.CD28","Unstim."))
res <- results(dds,contrast = c("Clusters","cGAMP_CD3.CD28","CD3.CD28"))
res <- results(dds,contrast = c("Clusters","cGAMP","Unstim."))
res <- results(dds,contrast = c("Clusters","cGAMP_CD3.CD28","cGAMP"))
res <- results(dds,contrast = c("Clusters","CD3.CD28","cGAMP"))

res <- results(dds,contrast = c("Clusters","cGAMP","cGAMP_CD3/CD28"))


volcano_data=as.data.frame(res)
dim(volcano_data)
volcano_data=na.omit(volcano_data)
dim(volcano_data)
head(volcano_data)
# add a grouping column; default value is "not significant"
volcano_data["group"] <- "NotSignificant"
# for our plot, we want to highlight FDR < 0.05 (significance level) Fold Change > 1.5
# change the grouping for the entries with significance but not a large enough Fold change
#volcano_data[base::which(volcano_data['padj'] < 0.05 & abs(volcano_data['log2FoldChange']) < 2 ),"group"] <- "Significant"
# change the grouping for the entries a large enough Fold change but not a low enough p value
#volcano_data[which(volcano_data['padj'] > 0.05 & abs(volcano_data['log2FoldChange']) > 2 ),"group"] <- "FoldChange"
# change the grouping for the entries with both significance and large enough fold change
volcano_data[which(volcano_data['padj'] < 0.05 & abs(volcano_data['log2FoldChange']) > 1 ),"group"] <- "Significant&FoldChange"
# to manualy set colors
pal <- c("blue", "black", "purple","red")
head(volcano_data)
volcano_data$Symbol=rownames(volcano_data)
# volcano_data$padj[volcano_data$Symbol=="CSF2"]=1e-300
#plot_ly(volcano_data, x=~log2FoldChange, y=~-log10(padj),
#          color=~group,colors=pal,
#          text=~Symbol,type="scatter",size = ~-log10(padj),
#          mode="markers") %>%
#  layout(title = 'Unstimulated vs. cGAMP + CD3/CD28',
#         yaxis = list(zeroline = T),
#         xaxis = list(zeroline = T)) %>% 
#  add_segments(x=1.5,xend = 1.5,y=0,yend = max(-log10(volcano_data$padj)),color="Euler cut-off",line=list(color="black",width =0.5,dash = 'dash'))

# add_lines(x=1.5,color="Euler cut-off",line=list(color="black",width =0.5,dash = 'dash'))

volcano_data_3.28 <- filter(volcano_data, group == "Significant&FoldChange")

volcano_data[which(volcano_data['padj'] < 0.05 & abs(volcano_data['log2FoldChange']) > 1 ),"group"] <- "Significant&FoldChange"
volcano_data_all[which(volcano_data_all$Symbol %in% volcano_data_3.28$Symbol), "group"] <- "3.28"
volcano_data_all_deleted <- filter(volcano_data_all, group != "3.28")

volcano_data_all$top=NA
volcano_data_all$group[volcano_data_all$group == "Significant&FoldChange" & volcano_data_all$Symbol %in% de_genes_isg] = "ISG"
volcano_data_all$top[volcano_data_all$group == "ISG"]=volcano_data_all$Symbol[volcano_data_all$group == "ISG"]

doubleISGs <- volcano_data_all$top
doubleISGs <- na.omit(doubleISGs)

volcano_data_all <- volcano_data
ggplot(volcano_data_all,aes(x=log2FoldChange,y=-log10(padj),label=top,color=group)) +
  geom_point(aes(), size = 2) +
  geom_text_repel(size = 7, max.overlaps = 20, box.padding = unit(0.2, 'lines'), point.padding = unit(0.2, 'lines')) +
  theme_classic() + 
  theme(axis.text = element_text(size = 20, color = "black"), 
        axis.line = element_line(size = 0.6), 
        axis.ticks = element_line(size = 0.6),
        axis.ticks.length = unit(.4, "cm"),
        axis.title = element_text(size = 20),
        text = element_text(family = "Helvetica")) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = col_grey) +
  geom_rect(aes(xmin = 0, xmax = max(abs(volcano_data_all$log2FoldChange)), ymin = 0, ymax = 40 ), fill = NA, linetype = 2) +
  #xlim(-5,5) + 
  xlab("log2FC") +
  ylab("-log10(padj)") +
  scale_color_manual(values = ggplot2::alpha(c(col_green, col_grey, col_red),1)) +
  ggtitle("DE in double vs all other groups") +
  theme(legend.position="none", plot.title = element_blank()) +
  xlim(-max(abs(volcano_data_all$log2FoldChange)),max(abs(volcano_data_all$log2FoldChange))) 
  # xlim(0,max(abs(volcano_data_all$log2FoldChange))) +
  # ylim(0, 40)


#####################################
plot_grid(unst_cd3,unst_cgamp,
unst_cgamp_cd3,cd3_vs_cgamp,
cd3_vs_cgamp_cd3,cgamp_vs_doubl,ncol = 2)

##########################################################################################################################################

############################################################
# one minus other
############################################################
sel=mixedsort(unique(c("IFNG","IL2","IL4","IL5","IL17RA","CSF2","GZMB","GATA3","TBX21","NLRP3","CASP1","RORC","STAT5A","STAT5B","FOXP3")))
exprss_matr=assay(dds)
exprss_matr=as.data.frame(log2(exprss_matr+1))
# exprss_matr=as.data.frame(assay(rld))
head(exprss_matr)
sel=sel[sel%in%rownames(exprss_matr)]
sel_bar=exprss_matr[sel,]
sel_bar$Gene=rownames(sel_bar)
head(sel_bar)
sel_bar=gather(sel_bar,key = Name,value = Exprs,-Gene)
head(sel_bar)
sel_data=left_join(sel_bar,cell_annotation,by="Name")
head(sel_data)
WT_unstim=sel_data[sel_data$Group=="Unstimulated_WT",]
WT_unstim=WT_unstim[order(WT_unstim$XC,WT_unstim$Gene),]
WT_CD3=sel_data[sel_data$Group=="CD3/CD28_WT",]
WT_CD3=WT_CD3[order(WT_CD3$XC,WT_CD3$Gene),]
WT_norm=WT_CD3
WT_norm$Exprs=WT_CD3$Exprs-WT_unstim$Exprs
WT_norm$Exprs
#####
Sting_unst=sel_data[sel_data$Group=="Unstimulated_STING_KO",]
Sting_unst=Sting_unst[order(Sting_unst$XC,Sting_unst$Gene),]
Sting_CD3=sel_data[sel_data$Group=="CD3/CD28_STING_KO",]
Sting_CD3=Sting_CD3[order(Sting_CD3$XC,Sting_CD3$Gene),]
Sting_Norm=Sting_CD3
Sting_Norm$Exprs=Sting_CD3$Exprs-Sting_unst$Exprs
Sting_Norm$Exprs
###########################
head(Sting_Norm)
head(WT_norm)
Compl_norm=rbind(WT_norm,Sting_Norm)
###########################
Compl_norm$Exprs=2^Compl_norm$Exprs
SEM <- function(x) sd(x)/sqrt(length(x))
head(Compl_norm)
mean_sel_data=Compl_norm %>% group_by(Gene, Type)%>%summarise(Value = mean(Exprs),SEM=SEM(Exprs))
head(mean_sel_data)
mean_sel_data$Gene=factor(mean_sel_data$Gene)
p1=ggplot(mean_sel_data, aes(x=Type, y=Value,fill=Type)) + geom_bar(stat = "identity", width=0.7,alpha=0.8,position = position_dodge(0.8))+facet_grid(~Gene)+
  ggtitle("Rlog data")+xlab("")+ylab('Fold Change (Stim. vs. Unstim.)')+theme_classic()+scale_fill_manual(values=cc) + 
  theme(legend.position = "none",axis.text.x = element_text(angle = 45, hjust = 1,color="black"))+labs(fill="")
p2
p1
plot_grid(p1,p2)
##############################################################################################################################


# VENN DIAGRAM

cgampvsunst=volcano_data
up_in_cgamp=cgampvsunst
up_in_cgamp=up_in_cgamp[up_in_cgamp$log2FoldChange>2,]
# The one's I'm kicking away
up_in_cgamp[!up_in_cgamp$padj<0.05,]
up_in_cgamp=up_in_cgamp[up_in_cgamp$padj<0.05,]
up_in_cgamp$from="up_in_cgamp"

doublevsunst=volcano_data
up_in_double=doublevsunst
up_in_double=up_in_double[up_in_double$log2FoldChange>2,]
up_in_double[!up_in_double$padj<0.05,]
up_in_double=up_in_double[up_in_double$padj<0.05,]
up_in_double$from="up_in_double"

# cgamp_cd3_vs_cd3=volcano_data
# up_in_double=cgamp_cd3_vs_unst[cgamp_cd3_vs_unst$log2FoldChange>1.5,]
# up_in_double[!up_in_double$padj<0.01,]
# up_in_double=up_in_double[up_in_double$padj<0.01,]
# up_in_double

cd3vsunst=volcano_data
up_in_cd3=cd3vsunst
up_in_cd3=up_in_cd3[up_in_cd3$log2FoldChange>2,]
up_in_cd3[!up_in_cd3$padj<0.05,]
up_in_cd3=up_in_cd3[up_in_cd3$padj<0.05,]
up_in_cd3$from="up_in_cd3"

# EULERRRRRRRRRR
eulerplot=data.frame(Genes=unique(c(up_in_cgamp$Symbol,up_in_double$Symbol,up_in_cd3$Symbol)),cGAMP=FALSE,CD3CD28=FALSE,Double=FALSE)
eulerplot[eulerplot$Genes%in%up_in_cgamp$Symbol,]$cGAMP=TRUE
eulerplot[eulerplot$Genes%in%up_in_double$Symbol,]$Double=TRUE
eulerplot[eulerplot$Genes%in%up_in_cd3$Symbol,]$CD3CD28=TRUE
eulerplot
together=rbind(up_in_cd3,up_in_cgamp,up_in_double)
together$FDR=-log10(together$padj)
together=together %>% group_by(Symbol)%>%summarise(Value = mean(FDR))#,SEM=SEM(Exprs))
together=together[order(together$Value,decreasing = T),]
together$Value[1]=300
eulerplot$pval=together$Value[match(eulerplot$Genes,together$Symbol)]
write.csv(eulerplot,"/Users/niklaskuhl/Downloads/RNAseq/eulerrrr.csv")

eulerplot=read.csv("/Users/niklaskuhl/Downloads/RNAseq/eulerrrr.csv")
eulerplot=eulerplot[-1]
library(RColorBrewer)
cc=c(brewer.pal(name = "Set3",n=9))
library(eulerr)
head(eulerplot)
# plot(venn(eulerplot[, 2:4]),col=cc[2],fill=cc,labels = list(font = 4,col="black"),
#      quantities = list(font = 4,col="black"),main=list(col="black", label="Overlapping TFs"))
plot(euler(eulerplot[, 2:4]),col = cc,lwd = 0, fill=my_cc,labels = list(col="black", family = "Lato"),
     quantities = list(col="black", family = "Lato"), alpha = 0.7)

# SAD PLOT! 
options(stringsAsFactors = F)
sad=read.csv("/Users/Ashto/Desktop/Sting-exp/euler to sad plot.csv")
sad=sad[,-1]
head(sad)
sad=sad[order(sad$pval,decreasing = T),]
sad$Genes=factor(sad$Genes,levels = sad$Genes)
head(sad)
sad$cGAMP
head(sleep)
#
cgamp=data.frame(gene=sad$Genes,level="cgamp",taide="puudu")
cgamp$taide[sad$cGAMP]="olemas"
cd3=data.frame(gene=sad$Genes,level="cd3",taide="puudu")
cd3$taide[sad$CD3CD28]="olemas"
double=data.frame(gene=sad$Genes,level="double",taide="puudu")
double$taide[sad$Double]="olemas"
#####
tot=rbind(cgamp,cd3,double)
head(tot)
tot$liin=NA
tot$liin=ifelse(tot$taide=="olemas","fill",NA)
bottom=ggplot(tot,aes(x=gene,y=level,fill=taide))+
  geom_point(size=3,pch = 21,color="white")+ theme_classic() +
  theme(legend.position = "none",axis.text.x = element_text(angle=45,hjust = 1))+ylab("")+
  scale_fill_manual(values = c("#000000","#e0e0e0"))+theme(panel.grid.major.x = element_line(color="gray",linetype = "dashed"))

library(ggplot2)
library(cowplot)

top=ggplot(sad,aes(x=Genes,y=pval))+geom_bar(stat="identity")+theme_classic()+xlab("")+theme(axis.text.x = element_blank(),axis.ticks = element_blank())+
  ylab("-log10(FDR)")

plot_grid(top,bottom,ncol = 1)

mid=data.frame(nimi=names(colSums(sad[,c(2:4)])),number=colSums(sad[,c(2:4)]),row.names = NULL)
side0=ggplot(mid,aes(x=nimi,y=number))+geom_bar(stat="identity")+theme_classic()+xlab("")+theme(axis.text.x = element_blank(),axis.ticks = element_blank())+
  ylab("Overlap size")+coord_flip()

side1=ggplot(mid,aes(x=nimi,y=number))+theme_classic()+xlab("")+theme(axis.text.y=element_blank(),axis.text.x = element_blank(),axis.ticks = element_blank())+
  ylab("")+coord_flip()

plot_grid(side1,top,side0,bottom,ncol = 2)
############# ############# ############# ############# ############# ############# ############# 
############# WORD CLOUD ############# 
############# ############# ############# ############# ############# ############# ############# 
library("wordcloud")
library("wordcloud2")

data=sad
data$pval=round(data$pval)
only_double=data[data$Double,]
only_double=only_double[!only_double$cGAMP,]
only_double=only_double[!only_double$CD3CD28,]
rownames(only_double)=NULL
only_double=only_double[-c(2,3,4)]
only_double

Only_cd3=data[data$CD3CD28,]
Only_cd3=Only_cd3[!Only_cd3$Double,]
Only_cd3=Only_cd3[!Only_cd3$cGAMP,]
rownames(Only_cd3)=NULL
Only_cd3=Only_cd3[-c(2,3,4)]
Only_cd3

double_and_cd3=data[data$Double,]
double_and_cd3=double_and_cd3[double_and_cd3$CD3CD28,]
rownames(double_and_cd3)=NULL
double_and_cd3=double_and_cd3[-c(2,3,4)]
double_and_cd3

double_and_cgamp=data[data$Double,]
double_and_cgamp=double_and_cgamp[double_and_cgamp$cGAMP,]
rownames(double_and_cgamp)=NULL
double_and_cgamp=double_and_cgamp[-c(2,3,4)]
double_and_cgamp

# 1 only double
set.seed(7)
pdf(file = "/Users/Ashto/Desktop/test.pdf",width=12,height=12)
par(mfrow=c(2,2))
# pdf(file = "/Users/Ashto/Desktop/test.pdf",width=6,height=6)
wordcloud(words = only_double$Genes, freq = only_double$pval, min.freq = 1,max.words=200, random.order=FALSE, rot.per=0.35,colors=brewer.pal(8, "Dark2"),scale = c(3,1))
# 2 only cd3
set.seed(7)
wordcloud(words = Only_cd3$Genes, freq = Only_cd3$pval, min.freq = 1,max.words=200, random.order=FALSE, rot.per=0.35,colors=brewer.pal(8, "Dark2"),scale = c(3,1))
# 3 cd3 and double
set.seed(7)
wordcloud(words = double_and_cd3$Genes, freq = double_and_cd3$pval, min.freq = 1,max.words=200, random.order=F, rot.per=0.35,colors=brewer.pal(8, "Dark2"),scale = c(3,1),use.r.layout = T)
# 4 cgamp and double
set.seed(7)
wordcloud(words = double_and_cgamp$Genes, freq = double_and_cgamp$pval, min.freq = 1,max.words=200, random.order=F, rot.per=0.35,colors=brewer.pal(8, "Dark2"),scale = c(3,1),use.r.layout = T)
dev.off()

#########################################################################################################
# set.seed(7)
# wordcloud2(only_double, size = 0.5)
# set.seed(7)
# wordcloud2(only_double, size = 0.3, shape = 'circle')
# wordcloud2(double_and_cd3, size = 0.7, shape = 'star')
# wordcloud2(only_double, size = 0.7, shape = 'star')
# wordcloud2(double_and_cgamp, size = 0.4, shape = 'star')
# wordcloud2(demoFreq, figPath = "/Users/Ashto/Desktop/peaceAndLove.jpg", size = 1.5, color = "skyblue", backgroundColor="black")
# letterCloud(demoFreq, word = "C", color='random-light' , backgroundColor="black")
# head(demoFreq)
# head(double_and_cgamp)
# data1=data[-c(2,3,4)]
# rownames(data1)=NULL
# wordcloud2(data1, size=1.6, color='random-dark',shape = "cardioid")
################################################################################################
################################################################################################
################################################################################################
# UNDER THE RADAR!!!
################################################################################################
################################################################################################
take=c("IL8","IL2","IL22","IL17A","IL9","IL13","IL5","IL4","TNF","IFNG","CSF2","GZMB","IL17RA")
exprss_matr=assay(dds)
exprss_matr=as.data.frame(log2(exprss_matr+1))
# exprss_matr=as.data.frame(assay(rld))
head(exprss_matr)
sel=take[take%in%rownames(exprss_matr)]
sel_bar=exprss_matr[sel,]
sel_bar$Gene=rownames(sel_bar)
head(sel_bar)
sel_bar=gather(sel_bar,key = Name,value = Exprs,-Gene)
head(sel_bar)
sel_data=left_join(sel_bar,cell_annotation,by="Name")
head(sel_data)
# sel_data$Exprs=2^sel_data$Exprs
# SEM <- function(x) sd(x)/sqrt(length(x))
mean_sel_data=sel_data %>% group_by(Gene, Group)%>%summarise(Value = mean(Exprs))#,SEM=SEM(Exprs))
head(mean_sel_data)
# mean_sel_data$Gene=factor(mean_sel_data$Gene)
# ggplot(mean_sel_data, aes(x=Clusters, y=Value,fill=Clusters)) + geom_bar(stat = "identity", width=0.7,alpha=0.8,position = position_dodge(0.8))+
#   facet_grid(~Gene)+ggtitle("Rlog data")+xlab("")+theme_classic()+scale_fill_manual(values=cc) + 
#   theme(legend.position = "none",axis.text.x = element_text(angle = 45, hjust = 1,color="black"))+labs(fill="")
mean_sel_data
wide_data=spread(mean_sel_data, key = Gene, value = Value)
wide_data

fig <- plot_ly(
  type = 'scatterpolar',
  fill = 'toself'
) 
fig <- fig %>%
  add_trace(
    r = as.numeric(wide_data[1,2:11]),
    theta = colnames(wide_data)[2:11],
    name =  as.character(wide_data$Clusters[1]) ) 
fig <- fig %>%
  add_trace(
    r = as.numeric(wide_data[2,2:11]),
    theta = colnames(wide_data)[2:11],
    name =  as.character(wide_data$Clusters[2]) )
fig <- fig %>%
  add_trace(
    r = as.numeric(wide_data[3,2:11]),
    theta = colnames(wide_data)[2:11],
    name =  as.character(wide_data$Clusters[3]) )
fig <- fig %>%
  add_trace(
    r = as.numeric(wide_data[4,2:11]),
    theta = colnames(wide_data)[2:11],
    name =  as.character(wide_data$Clusters[4]) )

fig
#####################
library(fmsb)
wide_data=as.data.frame(wide_data)
rownames(wide_data)=wide_data$Clusters
wide_data3=wide_data[2:ncol(wide_data)]
rownames(wide_data3)=wide_data$Group
wide_data_head=wide_data3
rownames(wide_data_head)[1:2]=c(1,2)
wide_data_head[2,]=0
wide_data_head[1,]=round(max(wide_data_head))
wide_data3=rbind(wide_data_head[1:2,],wide_data3)
wide_data3
order=c("IL2","CSF2","GZMB","IFNG","TNF","IL4","IL5","IL13","IL9","IL17A","IL17RA","IL22")
order=order[order%in%colnames(wide_data3)]
wide_data3=wide_data3[rev(order)]

cc=c(brewer.pal(name = "Set1",n=9))
par(mfrow=c(1,1))
max=round(max(wide_data_head))
radarchart(wide_data3, axistype=1 , 
            #custom polygon
            pcol=cc , plwd=2 , plty=c(1:4),
            #custom the grid
            cglcol="grey", cglty=1, axislabcol="grey", caxislabels=c(0,round(max*0.25,2),round(max*0.5,2),round(max*0.75,2),round(max*1,2)), cglwd=0.8,
            #custom labels
            vlcex=0.8
)
legend(x=1.2, y=1.2, legend = rownames(wide_data3)[-c(1,2)], bty = "n", pch=20 , col=cc , text.col = "black", cex=0.7, pt.cex=2)


# Library
library(fmsb)
# Create data: note in High school for several students
set.seed(99)
data <- as.data.frame(matrix( sample( 0:20 , 15 , replace=F) , ncol=5))
colnames(data) <- c("math" , "english" , "biology" , "music" , "R-coding" )
rownames(data) <- paste("mister" , letters[1:3] , sep="-")
# To use the fmsb package, I have to add 2 lines to the dataframe: the max and min of each variable to show on the plot!
data <- rbind(rep(20,5) , rep(0,5) , data)
data
# plot with default options:
radarchart(data)
