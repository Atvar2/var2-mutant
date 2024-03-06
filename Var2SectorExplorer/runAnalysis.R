suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(patchwork))
suppressPackageStartupMessages(library(clustree))
suppressPackageStartupMessages(library(ggsci))
suppressPackageStartupMessages(library(ggplot2)) 
suppressPackageStartupMessages(library(DoubletFinder))
suppressPackageStartupMessages(library(clusterProfiler))
suppressPackageStartupMessages(library(org.At.tair.db))
suppressPackageStartupMessages(library(monocle))
options(future.globals.maxSize = 20000 * 1024^2)  # //Globals of the maximum allowed size in bytes

setwd("./")

# =============================  1. 10 x data quality control =====================
Sys.time()
WT1 <- Read10X(data.dir = "./01.QC/cellRenger/WT1/outs/filtered_feature_bc_matrix")
Sys.time()

WT1 <- CreateSeuratObject(counts =WT1,  project = "WT1", min.cells = 3, min.features = 200)
WT1[["percent.mito"]] <- PercentageFeatureSet(WT1, pattern = "^mt-")   # Arabidopsis thaliana Gene ID conversion gene symbol
VlnPlot( WT1, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3 )
WT1 <- subset(WT1, subset = nFeature_RNA > 200 & nFeature_RNA < 4000 & nCount_RNA > 300 & nCount_RNA < 15000 & percent.mito < 15)
WT1 <- NormalizeData(WT1, assay = "RNA", normalization.method = "LogNormalize", scale.factor = 10000 )
WT1 <- FindVariableFeatures(WT1, selection.method = "vst", nfeatures = 2000, verbose = FALSE)
# WT1 <- SCTransform(WT1, return.only.var.genes = FALSE,assay = "RNA", verbose = FALSE)
WT1 <- ScaleData( WT1, features=rownames(WT1) )
WT1 <- RunPCA( object=WT1, features=VariableFeatures(WT1), npcs=20 )
WT1 <- FindNeighbors(WT1, dim=1:20)
WT1 <- FindClusters(WT1, resolution = 0.3)
WT1 <- RunTSNE( object = WT1, dims=1:20 )
WT1 <- RunUMAP( object = WT1, dims = 1:20 )
#DimPlot( object = WT1, reduction="umap" )
saveRDS( WT1, file = "./WT1.rds" )

##  remove doubletfinde
pcaN <- 20
train <- WT1
DoubletRate <- 0.075   # Manually adjusted basing on cell number in certain samples.

#pK Identification
sweep.res.list <- paramSweep_v3(train, PCs = 1:pcaN)
sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
bcmvn <- find.pK(sweep.stats)
pK_value <- as.numeric( as.character( bcmvn$pK[bcmvn$BCmetric == max(bcmvn$BCmetric)] ) )

#Homotypic Doublet Proportion Estimate
annotations <- train@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotations)
nExp_poi <- round( DoubletRate * length(train@meta.data$orig.ident) )
nExp_poi.adj <- round( nExp_poi * (1-homotypic.prop) )

#Run DoubletFinder with varying classification stringencies
pN_value <- 0.25
dropLable <-  paste("DF.classifications", pN_value, pK_value, nExp_poi.adj, sep="_")
#train <- doubletFinder(train, PCs = 1:pcaN, pN = pN_value, pK = pK_value, nExp = nExp_poi, reuse.pANN = FALSE)
train <- doubletFinder_v3(train, PCs = 1:pcaN, pN = pN_value, pK = pK_value, nExp = nExp_poi.adj, reuse.pANN = FALSE)

setwd("/szrmyy/Drwang/scRNA/chenjunhui/P00017_liver_INH/02.analysis")
pdf( "WT1.DF.cowplot.pdf", width = 6, height = 5 )
cowplot::plot_grid(ncol = 2, DimPlot(train, group.by = "orig.ident") + NoAxes(),
                   DimPlot(train, group.by = dropLable) + NoAxes())
dev.off()
pdf( "WT1.DF.VlnPlot.pdf", width = 6, height = 5 )
VlnPlot(train, features = "nFeature_RNA", group.by = dropLable, pt.size = 0.1)
dev.off()
train = train[, train@meta.data[, dropLable] == "Singlet"]

pdf( "WT1.tsne.DF.pdf" )
DimPlot( object = train, reduction="tsne" )
dev.off()

pdf( "WT1.umap.DF.pdf" )
DimPlot( object = train, reduction="umap", label = T )
dev.off()

# save Seurat object
saveRDS( train, file = "./WT1_DF.rds" )



Sys.time()
WT2 <- Read10X(data.dir = "./01.QC/cellRenger/WT2/outs/filtered_feature_bc_matrix")
Sys.time()

WT2 <- CreateSeuratObject(counts =WT2,  project = "WT2", min.cells = 3, min.features = 200)
WT2[["percent.mito"]] <- PercentageFeatureSet(WT2, pattern = "^mt-")   # Arabidopsis thaliana Gene ID conversion gene symbol
VlnPlot( WT2, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3 )
WT2 <- subset(WT2, subset = nFeature_RNA > 200 & nFeature_RNA < 4000 & nCount_RNA > 300 & nCount_RNA < 15000 & percent.mito < 15)
WT2 <- NormalizeData(WT2, assay = "RNA", normalization.method = "LogNormalize", scale.factor = 10000 )
WT2 <- FindVariableFeatures(WT2, selection.method = "vst", nfeatures = 2000, verbose = FALSE)
# WT2 <- SCTransform(WT2, return.only.var.genes = FALSE,assay = "RNA", verbose = FALSE)
WT2 <- ScaleData( WT2, features=rownames(WT2) )
WT2 <- RunPCA( object=WT2, features=VariableFeatures(WT2), npcs=20 )
WT2 <- FindNeighbors(WT2, dim=1:20)
WT2 <- FindClusters(WT2, resolution = 0.3)
WT2 <- RunTSNE( object = WT2, dims=1:20 )
WT2 <- RunUMAP( object = WT2, dims = 1:20 )
#DimPlot( object = WT2, reduction="umap" )
saveRDS( WT2, file = "./WT2.rds" )

##  remove doubletfinde
pcaN <- 20
train <- WT2
DoubletRate <- 0.075   # Manually adjusted basing on cell number in certain samples.

#pK Identification
sweep.res.list <- paramSweep_v3(train, PCs = 1:pcaN)
sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
bcmvn <- find.pK(sweep.stats)
pK_value <- as.numeric( as.character( bcmvn$pK[bcmvn$BCmetric == max(bcmvn$BCmetric)] ) )

#Homotypic Doublet Proportion Estimate
annotations <- train@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotations)
nExp_poi <- round( DoubletRate * length(train@meta.data$orig.ident) )
nExp_poi.adj <- round( nExp_poi * (1-homotypic.prop) )

#Run DoubletFinder with varying classification stringencies
pN_value <- 0.25
dropLable <-  paste("DF.classifications", pN_value, pK_value, nExp_poi.adj, sep="_")
#train <- doubletFinder(train, PCs = 1:pcaN, pN = pN_value, pK = pK_value, nExp = nExp_poi, reuse.pANN = FALSE)
train <- doubletFinder_v3(train, PCs = 1:pcaN, pN = pN_value, pK = pK_value, nExp = nExp_poi.adj, reuse.pANN = FALSE)

setwd("/szrmyy/Drwang/scRNA/chenjunhui/P00017_liver_INH/02.analysis")
pdf( "WT2.DF.cowplot.pdf", width = 6, height = 5 )
cowplot::plot_grid(ncol = 2, DimPlot(train, group.by = "orig.ident") + NoAxes(),
                   DimPlot(train, group.by = dropLable) + NoAxes())
dev.off()
pdf( "WT2.DF.VlnPlot.pdf", width = 6, height = 5 )
VlnPlot(train, features = "nFeature_RNA", group.by = dropLable, pt.size = 0.1)
dev.off()
train = train[, train@meta.data[, dropLable] == "Singlet"]

pdf( "WT2.tsne.DF.pdf" )
DimPlot( object = train, reduction="tsne" )
dev.off()

pdf( "WT2.umap.DF.pdf" )
DimPlot( object = train, reduction="umap", label = T )
dev.off()

# save Seurat object
saveRDS( train, file = "./WT2_DF.rds" )



Sys.time()
mutant1 <- Read10X(data.dir = "./01.QC/cellRenger/mutant1/outs/filtered_feature_bc_matrix")
Sys.time()

mutant1 <- CreateSeuratObject(counts =mutant1,  project = "mutant1", min.cells = 3, min.features = 200)
mutant1[["percent.mito"]] <- PercentageFeatureSet(mutant1, pattern = "^mt-")   # Arabidopsis thaliana Gene ID conversion gene symbol
VlnPlot( mutant1, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3 )
mutant1 <- subset(mutant1, subset = nFeature_RNA > 200 & nFeature_RNA < 4000 & nCount_RNA > 300 & nCount_RNA < 15000 & percent.mito < 15)
mutant1 <- NormalizeData(mutant1, assay = "RNA", normalization.method = "LogNormalize", scale.factor = 10000 )
mutant1 <- FindVariableFeatures(mutant1, selection.method = "vst", nfeatures = 2000, verbose = FALSE)
# mutant1 <- SCTransform(mutant1, return.only.var.genes = FALSE,assay = "RNA", verbose = FALSE)
mutant1 <- ScaleData( mutant1, features=rownames(mutant1) )
mutant1 <- RunPCA( object=mutant1, features=VariableFeatures(mutant1), npcs=20 )
mutant1 <- FindNeighbors(mutant1, dim=1:20)
mutant1 <- FindClusters(mutant1, resolution = 0.3)
mutant1 <- RunTSNE( object = mutant1, dims=1:20 )
mutant1 <- RunUMAP( object = mutant1, dims = 1:20 )
#DimPlot( object = mutant1, reduction="umap" )
saveRDS( mutant1, file = "./mutant1.rds" )

##  remove doubletfinde
pcaN <- 20
train <- mutant1
DoubletRate <- 0.075   # Manually adjusted basing on cell number in certain samples.

#pK Identification
sweep.res.list <- paramSweep_v3(train, PCs = 1:pcaN)
sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
bcmvn <- find.pK(sweep.stats)
pK_value <- as.numeric( as.character( bcmvn$pK[bcmvn$BCmetric == max(bcmvn$BCmetric)] ) )

#Homotypic Doublet Proportion Estimate
annotations <- train@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotations)
nExp_poi <- round( DoubletRate * length(train@meta.data$orig.ident) )
nExp_poi.adj <- round( nExp_poi * (1-homotypic.prop) )

#Run DoubletFinder with varying classification stringencies
pN_value <- 0.25
dropLable <-  paste("DF.classifications", pN_value, pK_value, nExp_poi.adj, sep="_")
#train <- doubletFinder(train, PCs = 1:pcaN, pN = pN_value, pK = pK_value, nExp = nExp_poi, reuse.pANN = FALSE)
train <- doubletFinder_v3(train, PCs = 1:pcaN, pN = pN_value, pK = pK_value, nExp = nExp_poi.adj, reuse.pANN = FALSE)

setwd("/szrmyy/Drwang/scRNA/chenjunhui/P00017_liver_INH/02.analysis")
pdf( "mutant1.DF.cowplot.pdf", width = 6, height = 5 )
cowplot::plot_grid(ncol = 2, DimPlot(train, group.by = "orig.ident") + NoAxes(),
                   DimPlot(train, group.by = dropLable) + NoAxes())
dev.off()
pdf( "mutant1.DF.VlnPlot.pdf", width = 6, height = 5 )
VlnPlot(train, features = "nFeature_RNA", group.by = dropLable, pt.size = 0.1)
dev.off()
train = train[, train@meta.data[, dropLable] == "Singlet"]

pdf( "mutant1.tsne.DF.pdf" )
DimPlot( object = train, reduction="tsne" )
dev.off()

pdf( "mutant1.umap.DF.pdf" )
DimPlot( object = train, reduction="umap", label = T )
dev.off()

# save Seurat object
saveRDS( train, file = "./mutant1_DF.rds" )



Sys.time()
mutant2 <- Read10X(data.dir = "./01.QC/cellRenger/mutant2/outs/filtered_feature_bc_matrix")
Sys.time()

mutant2 <- CreateSeuratObject(counts =mutant2,  project = "mutant2", min.cells = 3, min.features = 200)
mutant2[["percent.mito"]] <- PercentageFeatureSet(mutant2, pattern = "^mt-")   # Arabidopsis thaliana Gene ID conversion gene symbol
VlnPlot( mutant2, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3 )
mutant2 <- subset(mutant2, subset = nFeature_RNA > 200 & nFeature_RNA < 4000 & nCount_RNA > 300 & nCount_RNA < 15000 & percent.mito < 15)
mutant2 <- NormalizeData(mutant2, assay = "RNA", normalization.method = "LogNormalize", scale.factor = 10000 )
mutant2 <- FindVariableFeatures(mutant2, selection.method = "vst", nfeatures = 2000, verbose = FALSE)
# mutant2 <- SCTransform(mutant2, return.only.var.genes = FALSE,assay = "RNA", verbose = FALSE)
mutant2 <- ScaleData( mutant2, features=rownames(mutant2) )
mutant2 <- RunPCA( object=mutant2, features=VariableFeatures(mutant2), npcs=20 )
mutant2 <- FindNeighbors(mutant2, dim=1:20)
mutant2 <- FindClusters(mutant2, resolution = 0.3)
mutant2 <- RunTSNE( object = mutant2, dims=1:20 )
mutant2 <- RunUMAP( object = mutant2, dims = 1:20 )
#DimPlot( object = mutant2, reduction="umap" )
saveRDS( mutant2, file = "./mutant2.rds" )

##  remove doubletfinde
pcaN <- 20
train <- mutant2
DoubletRate <- 0.075   # Manually adjusted basing on cell number in certain samples.

#pK Identification
sweep.res.list <- paramSweep_v3(train, PCs = 1:pcaN)
sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
bcmvn <- find.pK(sweep.stats)
pK_value <- as.numeric( as.character( bcmvn$pK[bcmvn$BCmetric == max(bcmvn$BCmetric)] ) )

#Homotypic Doublet Proportion Estimate
annotations <- train@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotations)
nExp_poi <- round( DoubletRate * length(train@meta.data$orig.ident) )
nExp_poi.adj <- round( nExp_poi * (1-homotypic.prop) )

#Run DoubletFinder with varying classification stringencies
pN_value <- 0.25
dropLable <-  paste("DF.classifications", pN_value, pK_value, nExp_poi.adj, sep="_")
#train <- doubletFinder(train, PCs = 1:pcaN, pN = pN_value, pK = pK_value, nExp = nExp_poi, reuse.pANN = FALSE)
train <- doubletFinder_v3(train, PCs = 1:pcaN, pN = pN_value, pK = pK_value, nExp = nExp_poi.adj, reuse.pANN = FALSE)

setwd("/szrmyy/Drwang/scRNA/chenjunhui/P00017_liver_INH/02.analysis")
pdf( "mutant2.DF.cowplot.pdf", width = 6, height = 5 )
cowplot::plot_grid(ncol = 2, DimPlot(train, group.by = "orig.ident") + NoAxes(),
                   DimPlot(train, group.by = dropLable) + NoAxes())
dev.off()
pdf( "mutant2.DF.VlnPlot.pdf", width = 6, height = 5 )
VlnPlot(train, features = "nFeature_RNA", group.by = dropLable, pt.size = 0.1)
dev.off()
train = train[, train@meta.data[, dropLable] == "Singlet"]

pdf( "mutant2.tsne.DF.pdf" )
DimPlot( object = train, reduction="tsne" )
dev.off()

pdf( "mutant2.umap.DF.pdf" )
DimPlot( object = train, reduction="umap", label = T )
dev.off()

# save Seurat object
saveRDS( train, file = "./mutant2_DF.rds" )

# ======================= 2. Integrating datasets of samples and removing  batch 
WT1<-readRDS("./WT1_DF.rds")
WT2<-readRDS("./WT2_DF.rds")

mutant1<-readRDS("./mutant1_DF.rds")
mutant2<-readRDS("./mutant2_DF.rds")
object_list <- object_list = list(WT1, WT2, mutant1, mutant2)
selfeatures <- SelectIntegrationFeatures(object.list = object_list, nfeatures = 2000)
scc.list <- PrepSCTIntegration(object.list = object_list, anchor.features = selfeatures, verbose = FALSE)
scc.anchors <- FindIntegrationAnchors(object.list = scc.list, normalization.method = "SCT",anchor.features = selfeatures, verbose = FALSE)
scc_integrated <- IntegrateData(anchorset = scc.anchors, normalization.method = "SCT",verbose = FALSE)

data<-scc_integrated
dims=1:50
DefaultAssay(data) <-"integrated"
data <- RunPCA(data, verbose = FALSE)
data <- FindNeighbors(data,dim=dims)
data <- FindClusters(data,resolution = 0.5)
data <- RunUMAP (data,reduction="pca", dims = dims)
#data <- RunTSNE(data,dims = dims)
pdf("1-Dimplot.pdf")
DimPlot(data,label = TRUE,reduction = "umap")
DimPlot(data,label = TRUE,reduction = "umap", group.by = "orig.ident")
dev.off()
saveRDS(data,file="ATvar2_mt0.5_SCT_combine.samples.rds")

# =================== 3. Cell annotation and functional analysis ==================
DefaultAssay(data) <-"integrated"
markers <- FindAllMarkers(object = data, only.pos = TRUE, min.pct = 0.25, thresh.use = 0.25)
write.csv(markers,paste("1-seurat_clusters","Allmakers.csv",sep=""))

ids=bitr(markers$gene,'SYMBOL','ENTREZID','org.At.tair.db')
markers=merge(markers,ids,by.x='gene',by.y='SYMBOL')
gcSample=split(markers$gene, markers$cluster)
kegg <- compareCluster(gcSample, fun="enrichKEGG",organism="ath",pvalueCutoff=0.05)
pdf("3-All-cluster-Kegg.pdf",width=12,height=10)
dotplot(kegg)
dev.off()
bp <- compareCluster(gcSample,fun="enrichGO", OrgDb="org.At.tair.db", keyType="TAIR",ont= "BP")
pdf("3-All-cluster-Go_BP.pdf",width=20,height=10)
dotplot(bp)
dev.off()
mf <- compareCluster(gcSample,fun="enrichGO", OrgDb="org.At.tair.db", keyType="TAIR",ont= "MF")
pdf("3-All-cluster-Go_MF.pdf",width=20,height=10)
dotplot(mf)
dev.off()
cc <- compareCluster(gcSample,fun="enrichGO", OrgDb="org.At.tair.db", keyType="TAIR",ont= "CC")
pdf("3-All-cluster-Go_CC.pdf",width=20,height=10)
dotplot(cc)
dev.off()
saveRDS(c(kegg,bp,mf,cc),"All-cluster-KeggGo.RDS")


data@meta.data$Celltype="t"
data@meta.data$Celltype[data$seurat_clusters %in% c(0,1,2,3,4,5,6,13)]="Mesophyll cell"
data@meta.data$Celltype[data$seurat_clusters %in% c(7,12)]="Epidermal cell"
data@meta.data$Celltype[data$seurat_clusters %in% c(10,14)]="Guard cell"
data@meta.data$Celltype[data$seurat_clusters %in% c(8)]="Phloem parenchyma"
data@meta.data$Celltype[data$seurat_clusters %in% c(9)]="Xylem cell"
data@meta.data$Celltype[data$seurat_clusters %in% c(11,15)]="Compansion cell"
data@mata.data$Celltype = factor(data@mata.data$Celltype,levels = c("Mesophyll cell",
"Phloem parenchyma","Epidermal cell","Xylem cell","Guard cell","Compansion cell"))
saveRDS(data,file="ATvar2_mt0.5_SCT_combine.samplesAnno.rds")
# Gene Difference Analysis
data <- NormalizeData(object=data)
markers <- FindAllMarkers(data, assay = "RNA",slot= "data", only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
markers <-subset(markers, p_val_adj<0.05)
markers <-write.csv(file="2-Celltype_markers.csv",markers)

Idents(data)<-"Groups"
DefaultAssay(data)<-"RNA"

Idents(data)<-"Celltype"
celltype = as.character(levels(data))
dat = data.frame(gene="",Difference=0,logFC=0,p=0,adj.P=0,Celltype="",cellNumber=0, State="")
for (i in celltype) {
  subdata=subset(data,idents=i)
  print(unique(subdata$Celltype))
  Idents(subdata)<-"Groups"
  cellNumber <- length(colnames(subdata))
  marker<-FindMarkers(subdata,ident.1="INH",idents.2="Ctrl", assay="RNA", only.pos=F, min.pct=0.1,logfc.threshold=0.25)
  marker<- marker %>% mutate(Difference = pct.1-pct.2)
  submarker = data.frame(gene=row.names(marker),Difference=marker$Difference, logFC=marker$avg_log2FC,p=marker$p_val,adj.P=marker$p_val_adj,Celltype=i,cellNumber=cellNumber,State=ifelse(marker$p_val_adj<0.05,ifelse(marker$avg_log2FC>0.25,"Up",ifelse(marker$avg_log2FC< -0.25,"Down","No")),"No"))
  
  dat=rbind(dat,submarker)
}
dat<-dat[-1,]
write.csv(dat,"3-Celltype_INHvsCtrl.csv",quote=F,row.names=F)

# =================== 4.1 sub type analysis ========================================
data <- readRDS("./ATvar2_mt0.5_SCT_combine.samplesAnno.rds")
data<- subset(data, celltype == "Mesophyll cell")

dims=1:48
DefaultAssay(data) <-"integrated"
data <- RunPCA(data, verbose = FALSE)
data <- FindNeighbors(data,dim=dims)
data <- FindClusters(data,resolution = 0.8)
data <- RunUMAP (data,reduction="pca", dims = dims)
#data <- RunTSNE(data,dims = dims)
pdf("1-Dimplot.pdf")
DimPlot(data,label = TRUE,reduction = "umap")
DimPlot(data,label = TRUE,reduction = "umap", group.by = "orig.ident")
dev.off()
saveRDS(data,file="./scc_integrated_48_0.8_leaf.rds")
p2<-FeaturePlot(data, features = gene,pt.size=0.1,reduction="umap")
data$cluster<-"Mc1"
data$cluster[data$seurat_clusters %in% c(0,2,6,12)]<-"Mc2"
data$cluster[data$seurat_clusters %in% c(4,8)]<-"Mc3"
data$cluster[data$seurat_clusters %in% c(3,11)]<-"Mc4"
data$cluster[data$seurat_clusters %in% c(9)]<-"Mc5"
data$cluster[data$seurat_clusters %in% c(15)]<-"Mc6"

saveRDS(data,file="./scc_integrated_48_0.8_leafAnno.rds")
markers <- FindAllMarkers(data, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
markers <-subset(markers, p_val_adj<0.05)
top10 <- markers  %>%  group_by(cluster)  %>%  top_n(n = 10, wt = avg_log2FC)
write.csv(markers,"1-cluster_markers.csv",quote=F,row.names=F)

bp <- compareCluster(gene~cluster,data=markers,fun="enrichGO", OrgDb="org.At.tair.db",keyType="TAIR", ont= "BP")
bp <-simplify(bp,cutoff=0.7,by="p.adjust",select_fun=min)
pdf(paste("1-At-Annocluster-Go_BP_",top,".pdf",sep=""),width=20,height=20)
dotplot(bp,showCategory=top)
dev.off()
mf <- compareCluster(gene~cluster,data=markers,fun="enrichGO", OrgDb="org.At.tair.db",keyType="TAIR" ,ont= "MF")
mf <-simplify(mf,cutoff=0.7,by="p.adjust",select_fun=min)
pdf(paste("1-At-Annocluster-Go_MF_",top,".pdf",sep=""),width=20,height=20)
dotplot(mf,showCategory=top)
dev.off()
cc <- compareCluster(gene~cluster,data=markers,fun="enrichGO", OrgDb="org.At.tair.db",keyType="TAIR" ,ont= "CC")
cc <-simplify(cc,cutoff=0.7,by="p.adjust",select_fun=min)
pdf(paste("1-At-Annocluster-Go_CC_",top,".pdf",sep=""),width=20,height=20)
dotplot(cc,showCategory=top)
dev.off()
path<-rbind(kegg_formula_res@compareClusterResult,bp@compareClusterResult,cc@compareClusterResult, mf@compareClusterResult)
write.csv(path,"1-At-Annocluster-path.csv",quote=F)


# Deg
Idents(data)<-"cluster"
celltype = as.character(unique(data$cluster))
dat = data.frame(gene="",logFC=0,adj.P=0,Celltype="",State="")
for (i in celltype) {
  subdata=subset(data,idents=i)
  Idents(subdata)<-"orig.ident"
  marker<-FindMarkers(subdata,ident.1=c("mutant1","mutant2"),min.pct=0.25,logfc.threshold=0.25)
  submarker = data.frame(gene=row.names(marker),logFC=marker$avg_log2FC,adj.P=marker$p_val_adj,Celltype=i,State=ifelse(marker$p_val_adj<0.05,ifelse(marker$avg_log2FC>0.25,"Up",ifelse(marker$avg_log2FC< -0.25,"Down","No")),"No"))
  dat=rbind(dat,submarker)
}
dat<-dat[-1,]
write.csv(dat,"3-annoClustermutantVSWT.csv",quote=F,row.names=F)

deg=subset(datdeg,State!="No")
celltype<-unique(data$cluster)
dat = data.frame(celltype="",ID="",Description="",GeneRatio="",BgRatio="",pvalue=0,p.adjust=0,qvalue=0,geneID="",Count=0)
for (i in celltype) {
  ddeg<-subset(deg,Celltype==i)
  x<- ddeg$gene
  ego_MF <- enrichGO(OrgDb="org.Mm.eg.db", gene = eg[,2], pvalueCutoff = 0.05, ont = "MF", readable=TRUE)
  ego_CC <- enrichGO(OrgDb="org.Mm.eg.db", gene = eg[,2], pvalueCutoff = 0.05, ont = "CC", readable=TRUE)
  ego_BP <- enrichGO(OrgDb="org.Mm.eg.db", gene = eg[,2], pvalueCutoff = 0.05, ont = "BP", readable=TRUE)
  ego_BP  <- enrichGO(x, OrgDb = "org.At.tair.db", keyType="TAIR",pvalueCutoff = 0.05 ,ont="BP")
  if (!is.null(ego_BP)) {
    subType<-data.frame(celltype=i,as.data.frame(ego_BP@result))
  }
  dat=rbind(dat,subType)
}
dat<-dat[-1,]
write.csv(dat,"4-Celltype_DEGsmutantVSWTDegEnrichment.csv",quote=F,row.names=F)

# ===================== 4.2 Pseudotime analysis on subtypes of leaf ================
data<- readRDS("./leafAnnocluster.RDS")
markers<-read.csv("1-cluster_markers.csv",header=T)
subdata<-subset(data, cluster %in% c("Mc1","Mc2","Mc3"))
markers<-subset(markers, cluster %in% c("Mc1","Mc2","Mc3"))
gene<-unique(markers$gene)
gene<-unique(c(gene, "AT3G48430","AT3G44680","AT1G79000","AT4G34060",
               "AT5G13170","AT1G66580","AT3G10985","AT4G02380","AT2G29350","AT5G14930","AT1G71190","AT1G20620","AT4G23810","AT1G62300","AT2G40750","AT4G01250","AT5G13080","AT3G56400","AT1G34180","AT1G69490","AT3G29035","AT5G39610","AT2G43000","AT1G32640","AT5G46760","AT4G17880","AT4G16430","AT1G01260","AT4G00870","AT2G46510","AT3G18520","AT5G03740","AT4G13250","AT5G04900","AT4G22920","AT5G13800","AT2G23140","AT2G03670")) # // add aging  genes
subdata<-subset(x=subdata,downsample=10000)
subdata <- data[,rownames(pData(cds))]    #// 取出子集
gene <- gene[gene %in% rownames(subdata)]
expr<-as.matrix(subdata@assays$RNA@counts)[gene,] #// add aging gene
expr<-as.matrix(subdata@assays$RNA@counts)[gene,]
gene<-data.frame(gene_short_name=gene,row.names=gene)
met<-subdata@meta.data

pd <- new("AnnotatedDataFrame", data = met)
fd <- new("AnnotatedDataFrame", data = gene)
cds <- newCellDataSet(expr, phenoData = pd, featureData = fd,expressionFamily=negbinomial.size())
cds <- estimateSizeFactors(cds)
cds <- estimateDispersions(cds)

cds <- reduceDimension(cds, max_components = 2, method = 'DDRTree')
cds <- orderCells(cds)
save(cds,file="4-ClusterMarkersGenes.rData")




