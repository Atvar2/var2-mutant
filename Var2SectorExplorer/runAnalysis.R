suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(patchwork))
suppressPackageStartupMessages(library(clustree))
suppressPackageStartupMessages(library(ggsci))
suppressPackageStartupMessages(library(ggplot2)) 
suppressPackageStartupMessages(library(clusterProfiler))
suppressPackageStartupMessages(library(org.At.tair.db))
suppressPackageStartupMessages(library(monocle))
suppressPackageStartupMessages(library(SCP))
options(future.globals.maxSize = 20000 * 1024^2)  # //Globals of the maximum allowed size in bytes

setwd("./")

# global colors used in the project
Cellcolors <- c(
  "#C5D395", # MC
  "#8DC3AC", # EC
  "#7D8BAB", # GC
  "#965797", #PP
  "#3C8E9D", #XC
  "#5284B3"  #CC
)
Groupcolors <- c("#68A28D", "#C09A60")
Mccolors <- c(
	"#4A813D", "#C9C791", "#B37889",
	"#A4CAE1", "#4984AB", "#EDD6DD"
) 

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

# ======================= 2. Integrating datasets of samples and removing  batch 
WT1<-readRDS("./WT1.rds")
WT2<-readRDS("./WT2.rds")

mutant1<-readRDS("./mutant1.rds")
mutant2<-readRDS("./mutant2.rds")
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
pdf("1-VlnPlotbyorig.ident.pdf", width=6,height=3)
VlnPlot(Atdata, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3,pt.size = 0,cols = c("#cccde1","#9fb5d0", "#e88441", "#cf5219"), group.by = "orig.ident")
dev.off()

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

colors <-c("#E5D4DE","#51AA5D","#F1BC74","#F3B3A2","#D7E8A5",
           "#50AED6","#446C88","#EC5C57","#E59DC5","#1E4129",
           "#AC2E83","#BD966A","#8D529B","#A0A3A9","#E2D5CD",
           "#5D3967","#C6E0BD","#E4C753")                         # // colors for seurat_clusters
DimPlot(data,label=F,group.by="seurat_clusters",reduction="umap",pt.size=0.1,cols = colors)		   
gene=c("AT5G38430","AT2G39470","AT1G70760",# Mesophyll
       "AT2G26250", "AT2G39400",           # 表皮
       "AT3G24140","AT2G46070","AT3G26744", #保卫细胞
       "AT3G48740","AT5G23660", #Phloemparenchyma
       "AT3G25710",  #Xylem cell
       "AT1G22710", "AT1G79430" #Companion cell
)
library(MySeuratWrappers)
DefaultAssay(data)<-"RNA"
VlnPlot(data, features = gene,group.by="seurat_clusters",stacked=T,pt.size=0)+theme(axis.ticks.y = element_blank(),axis.text.y = element_blank())+labs(x="",y="") + 
  scale_color_manual(values=colors)+
  scale_fill_manual(values=colors)

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
markers <- FindAllMarkers(data, assay = "RNA", only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
markers <-subset(markers, p_val_adj<0.05)
markers <-write.csv(file="3-Celltype_markers.csv",markers)

VlnPlot(data, features = gene,group.by="Celltype",stacked=T,pt.size=0)+theme(axis.ticks.y = element_blank(),axis.text.y = element_blank())+labs(x="",y="") + 
  scale_color_manual(values=Cellcolors)+
  scale_fill_manual(values=Cellcolors)  # // classic Markers 

# // marker gene featureplot
umap_theme <- theme(
  plot.title =element_blank(),
  axis.line=element_blank(),
  axis.text.x=element_blank(),
  axis.text.y=element_blank(),
  axis.ticks=element_blank(),
  axis.title.x=element_blank(),
  axis.title.y=element_blank(),
  panel.background=element_blank(),
  panel.border=element_blank(),
  panel.grid.major=element_blank(),
  panel.grid.minor=element_blank()
)
alpha.use=0.8
library(ggplot2)
library(RColorBrewer)
pbmc<-subset(x = data, downsample = 6000)

library(viridis)
col3<-colorRampPalette(c('#F9F8CC','#F18A3D','#7B1823'))(30)
col3<-c("#dcdcdc",col3)
plot_list <-list()
plot_list <- FeaturePlot(
  pbmc,reduction="umap",
  features=gene,
  combine=FALSE, cols=col3,
  min.cutoff = 'q1',
  max.cutoff='q90'
)
# apply theme to each feature plot
for(i in 1:length(plot_list)){
  plot_list[[i]] <- plot_list[[i]] + scale_alpha_continuous(range = alpha.use, guide = F)+umap_theme + NoLegend()
}
pdf('3-basic_canonical_marker_featurePlotUd.pdf', width=12, height=10)
CombinePlots(plot_list)
dev.off()

# // DEGs of cell types
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
  marker<-FindMarkers(subdata,ident.1="mutant",ident.2="WT", assay="RNA", only.pos=F, min.pct=0.1,logfc.threshold=0.25)
  marker<- marker %>% mutate(Difference = pct.1-pct.2)
  submarker = data.frame(gene=row.names(marker),Difference=marker$Difference, logFC=marker$avg_log2FC,p=marker$p_val,adj.P=marker$p_val_adj,Celltype=i,cellNumber=cellNumber,State=ifelse(marker$p_val_adj<0.05,ifelse(marker$avg_log2FC>0.25,"Up",ifelse(marker$avg_log2FC< -0.25,"Down","No")),"No"))
  
  dat=rbind(dat,submarker)
}
dat<-dat[-1,]
write.csv(dat,"3-Celltype_mutantvsWT.csv",quote=F,row.names=F)

# // cellular composition
count<-table(data$Groups,data$Celltype)
percentage<-apply(count,2,function(x){x/rowSums(count)})
p<-percentage[c("mutant","WT"),]
logfc<-percentage["mutant",]/percentage["WT",]
mege<-rbind(p,abs(log2(logfc)))
rownames(mege)<- c("mutant","WT","logfc(abs)")
write.csv(file="CelltypeLogfc.csv",mege)
pdf("3-CelltypeCompositionBarplot.pdf")
dtype<-as.data.frame(table(data@meta.data[,c('Celltype','Groups')]))
ggplot(dtype, aes(fill=Celltype, y=Freq, x=Groups)) +scale_fill_manual(values=Cellcolors)+
  geom_bar(position="fill", stat="identity")+theme_classic()
ggplot(dtype, aes(fill=Celltype, y=Freq, x=Groups)) +
  geom_bar(position="stack", stat="identity")+ scale_fill_manual(values=Cellcolors)+theme_classic()
dev.off()
# =================== 4. colors score ========================================
gene20<-list(c("AT1G29920","AT1G29910","AT1G29930","AT2G34430","AT2G34420",
 "AT2G05100","AT2G05070","AT3G27690","AT5G54270","AT5G01530","AT3G08940","AT2G40100",
 "AT4G10340","AT1G15820","AT3G54890","AT3G61470","AT1G61520","AT3G47470","AT1G45474","AT1G19150"))  # // light harvest genes
DefaultAssay(object) <- "RNA"

pbmc <- AddModuleScore(
  object = object,
  features = gene,
  ctrl = 100, #默认值是100
  assay="RNA",
  name = 'Lightharvest'
)

scoreDis<-read.table(file="calculatedWTAndmutantRatio.xls",sep="\t",header=T)  # // calculatedWTAndmutantRatio.xls calculated by calculateRangeByScore.pl
scoreDis<-data.frame(Bin=rownames(scoreDis),score=rep(scoreDis$AverageScore.1000.,2),type=c(rep("WT",40),rep("mutant",40)),value=c(scoreDis$wildRatio,scoreDis$mutantRatio))
scoreDis$Bin<-factor(scoreDis$Bin,levels=c(paste("Bin",seq(0,39,1),sep="")))
pdf("4-harvestlightw1000Stacked.pdf",height=8,width=10)
ggplot(scoreDis, aes( x = Bin,y=100 * value,fill = type))+
  geom_col(position = 'stack', width = 0.6)+
  #geom_bar(position = "stack", stat = "identity", width = 0.6)
  theme_bw()+   
  scale_fill_manual(values= Groupcolors)+  
  scale_y_continuous(expand = c(0,0))+
  theme(
    text=element_text(size=12),
    plot.title = element_text(hjust = 0.5,vjust = 0.5),
    axis.text.y=element_text(size=12,color = "black"),
    axis.text.x=element_text(size=12,  color = "black",angle = 45, hjust = 0.5,vjust = 0.5),
    legend.title=element_text(size=12),
    legend.text=element_text(size=12)

  )+
  theme(panel.grid = element_blank(), 
        panel.background = element_rect(color = 'black', fill = 'transparent'),
        )+
  guides(fill=guide_legend(keywidth = 1, keyheight = 1))
dev.off()
# =================== 5. Gene-gene correlations by single-cell RNA sequencing ========================================
# correlations were aaclcualted mainly based on the R packages "scCorr"
# // step1 run cluster
Attest<-subset(x = Atdata, downsample = 2000) # //randomly selected subsets cells (Too many cells can affect running speed)
dat   <- as.matrix(Attest@assays$RNA@data)
dat.t <- t(dat);
tsne  <- Rtsne(dat.t, dims = 2, perplexity=30, verbose=TRUE,num_threads =40,
               max_iter = 2000);      # Consider setting partial_pca=TRUE for large matrices
out   <- tsne$Y;
row.names(out) <- colnames(dat);
outfile        <- paste0("Atdo_tsne30_2000.txt")
write.table(out, outfile, sep="\t", quote=F)

#// step2 GCcluster
options(stringsAsFactors = F);

scale.v <- function(v, v.start=v.start, v.end=v.end){
  v.start + (v - min(v)) * (v.end-v.start)/(max(v) - min(v));
}
print(prefix)
dat  <- read.table("./Atdo_tsne30_2000.txt");
## number of cluster
for(k in c(seq(20,500,10))){
  #k     <- 50;
  ## scale region is from v1 to v2
  cell.number <- nrow(dat);
  v <- 112.65840 + 0.01799 * cell.number
  c1    <- scale.v(dat[,1], -v, v)
  c2    <- scale.v(dat[,2], -v, v)
  dat2  <- data.frame(v1=c1, v2=c2, row.names=row.names(dat))
  ## do GCluster
  clu.o <- GCluster(dat2, k=k);

  clu.i <- clu.o$membership;
  out.c <- data.frame(dat, cluster=clu.i);
  outf  <- paste0("clu_", k, ".txt");
  write.table(out.c, outf, quote=F, sep="\t");
}

# step 3 merge cluster
options(stringsAsFactors = F);
## all cluster sets: from 20 to 1000
clu.n <- c(seq(20,500,10));
clu.n <- sort(unique(clu.n))
out.s <- NULL;
for (n in clu.n){
  inff  <- paste0("./clu_", n, ".txt");
  dat   <- read.table(inff, header=T);
  out.s <- cbind(out.s, dat[,3]);
}

row.names(out.s) <- row.names(dat);
colnames(out.s)  <- paste0("clu", clu.n);
write.table(out.s, "02clust_table_raw.txt", quote=F, sep="\t");

## check the cluster consistency
for (i in 1:ncol(out.s)){
  for (j in i:ncol(out.s)){
    if (i==j){
      next;
    }
    clu.s <- paste0(out.s[,i], "_", out.s[,j]);
    c1.n  <- length(unique(out.s[,i]));
    c2.n  <- length(unique(out.s[,j]));
    c.n   <- max(c1.n, c2.n);
    tmp   <- table(clu.s);
    if (c.n == length(tmp)){
    } else {
      cat(i, " ", j, "\n");
    }
  }
}
# step4 correlations
gene<-read.table(file=geneFile, sep="\t",header=T)
gene20<-list(c("AT1G29920","AT1G29910","AT1G29930","AT2G34430","AT2G34420","AT2G05100","AT2G05070","AT3G27690","AT5G54270","AT5G01530","AT3G08940","AT2G40100","AT4G10340","AT1G15820","AT3G54890","AT3G61470","AT1G61520","AT3G47470","AT1G45474","AT1G19150"))
DefaultAssay(Atdata) <- "RNA"
Atdata <- AddModuleScore(
  object = Atdata,
  features = gene20,
  ctrl = 100, #默认值是100
  assay="RNA",
  name = 'Lightharvest'
)
lightharvestScore<- FetchData(object=Atdata, vars="Lightharvest1")

cellsSet <- rownames(diffClusterMaxtrix)
AtdatasubSet <- subset(Atdata, cells=cellsSet)

dat   <- as.matrix(AtdatasubSet@assays$RNA@data)
ObjcellsSet <- colnames(dat)
LhsuubsetScore <- lightharvestScore[ObjcellsSet,]
dat2 <-rbind(LhsuubsetScore,dat)
j=20
for(i in 1:49){
g_label <- diffClusterMaxtrix[,i]
names(g_label) =colnames(dat2);
Avexpr    <- get_value(dat2, g_label);

outf<-paste(prefix,j,".xls",sep="")
j =j + 10
out.s <- NULL;
for (i in 1:nrow(gene)){
  n1 <- gene[i, 1]; # gene: gene pairs from the pathways
  n1<-gsub("\\s+", "", n1)
  n2 <- gene[i, 2];
  n2<-gsub("\\s+", "", n2)
  print(n1);print(n2)
  i1 <- which(rownames(Avexpr)==n1);
  i2 <- which(rownames(Avexpr)==n2);
  v1 <- Avexpr[i1,]                                        # // dat 表达矩阵
  v2 <- Avexpr[i2,]
  v1.v <- var(v1);
  v2.v <- var(v2);
  if (v1.v==0 | v2.v==0){
    next;
  }
  out <- summary(glm(v1~v2));
  out2  <- cor(v1, v2);
  out3  <- cor(v1, v2, method = "spearman");
  out.s <- rbind(out.s, c(n1, n2, out2, out3, out$coef[2,3], out$coef[2,4]));
}
p.adjust.M <- p.adjust.methods[c(4,7)];
p.adj   <- sapply(p.adjust.M, function(meth) p.adjust(out.s[,6], meth));
out.p   <- cbind(out.s, p.adj);
out.p   <- out.p[order(as.numeric(out.p[,6])),]
colnames(out.p) <- c("gene1", "gene2", "pearson.R", "spearman.R", "t.value", "Pvalue", "bonferroni", "fdr")
write.table(out.p, outf, sep="\t", row.names=F, quote=F)
}

# =================== 6.1 sub type analysis ========================================
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
saveRDS(data,file="./6-scc_integrated_48_0.8_leaf.rds")
p2<-FeaturePlot(data, features = gene,pt.size=0.1,reduction="umap")
data$cluster<-"Mc1"
data$cluster[data$seurat_clusters %in% c(0,2,6,12)]<-"Mc2"
data$cluster[data$seurat_clusters %in% c(4,8)]<-"Mc3"
data$cluster[data$seurat_clusters %in% c(3,11)]<-"Mc4"
data$cluster[data$seurat_clusters %in% c(9)]<-"Mc5"
data$cluster[data$seurat_clusters %in% c(15)]<-"Mc6"

DimPlot(data,label = TRUE,cols= Mccolors,reduction = "umap")
markers <- FindAllMarkers(data, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
markers <-subset(markers, p_val_adj<0.05)
top10 <- markers  %>%  group_by(cluster)  %>%  top_n(n = 10, wt = avg_log2FC)
write.csv(markers,"6-cluster_markers.csv",quote=F,row.names=F)
markers <- FindAllMarkers(data, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
markers <-subset(markers, p_val_adj<0.05)
markers<-read.csv("./6-cluster_markers.csv",header=T)
dat<-AverageExpression(data)
dat<-dat$RNA
anno_col<-data.frame(orig.ident=colnames(dat))
rownames(anno_col)<-factor(colnames(dat),levels=c("Mc1","Mc2", "Mc3", "Mc4", "Mc5", "Mc6"))
ann_colors = list(
  orig.ident = c("Mc1"="#EC5C57","Mc2"="#E59DC5", "Mc3"="#1E4129", "Mc4"="#AC2E83", "Mc5"="#BD966A", "Mc6"="#8D529B")
)
dat<-dat[,c("Mc1","Mc2","Mc3","Mc4","Mc5","Mc6")]

pdf("6-leafcluster-markersgene.pdf",height=20,width=10)
pheatmap(dat[markers$gene,],scale="row",show_rownames=F,gaps_col=c(1,2,3,4,5),annotation_colors = ann_colors,annotation_col = anno_col,cluster_rows = FALSE,cluster_cols = FALSE,color=colorRampPalette(c("#594867","#219186","#e5ea36"))(100),border_color=NA,angle_col=45)
dev.off()

bp <- compareCluster(gene~cluster,data=markers,fun="enrichGO", OrgDb="org.At.tair.db",keyType="TAIR", ont= "BP")
bp <-simplify(bp,cutoff=0.7,by="p.adjust",select_fun=min)
pdf(paste("6-At-AnnoMc_cluster-Go_BP_",top,".pdf",sep=""),width=20,height=20)
dotplot(bp,showCategory=top)
dev.off()
mf <- compareCluster(gene~cluster,data=markers,fun="enrichGO", OrgDb="org.At.tair.db",keyType="TAIR" ,ont= "MF")
mf <-simplify(mf,cutoff=0.7,by="p.adjust",select_fun=min)
pdf(paste("6-At-AnnoMc_cluster-Go_MF_",top,".pdf",sep=""),width=20,height=20)
dotplot(mf,showCategory=top)
dev.off()
cc <- compareCluster(gene~cluster,data=markers,fun="enrichGO", OrgDb="org.At.tair.db",keyType="TAIR" ,ont= "CC")
cc <-simplify(cc,cutoff=0.7,by="p.adjust",select_fun=min)
pdf(paste("6-At-AnnoMc_cluster-Go_CC_",top,".pdf",sep=""),width=20,height=20)
dotplot(cc,showCategory=top)
dev.off()
path<-rbind(kegg_formula_res@compareClusterResult,bp@compareClusterResult,cc@compareClusterResult, mf@compareClusterResult)
write.csv(path,"6-At-AnnoMc_cluster-path.csv",quote=F)

# // Cellular composition
pdf("6-MesophyllCompositionBarplot.pdf")
dtype<-as.data.frame(table(data@meta.data[,c('cluster','Groups')]))
ggplot(dtype, aes(fill=cluster, y=Freq, x=Groups)) +scale_fill_manual(values=Mccolors)+
  geom_bar(position="fill", stat="identity")+theme_classic()
ggplot(dtype, aes(fill=cluster, y=Freq, x=Groups)) +
  geom_bar(position="stack", stat="identity")+ scale_fill_manual(values=Mccolors)+theme_classic()
dev.off()

# Degs analysis in Mesophyll cell subtypes 
Idents(data)<-"cluster"
celltype = as.character(unique(data$cluster))
dat = data.frame(gene="",logFC=0,adj.P=0,Celltype="",State="")
for (i in celltype) {
  subdata=subset(data,idents=i)
  Idents(subdata)<-"orig.ident"
  marker<-FindMarkers(subdata,ident.1=c("mutant1","mutant2"), ident.2=c("WT1", "WT2"),min.pct=0.25,logfc.threshold=0.25)
  submarker = data.frame(gene=row.names(marker),logFC=marker$avg_log2FC,adj.P=marker$p_val_adj,Celltype=i,State=ifelse(marker$p_val_adj<0.05,ifelse(marker$avg_log2FC>0.25,"Up",ifelse(marker$avg_log2FC< -0.25,"Down","No")),"No"))
  dat=rbind(dat,submarker)
}
dat<-dat[-1,]
write.csv(dat,"6-annoClustermutantvsWT.csv",quote=F,row.names=F)

# // Enrichment of DEGs in Mesophyll cell subtypes
datdeg <- read.csv(file="./6-annoClustermutantvsWT.csv", header=T, stringsAsFactors=F)
deg=subset(datdeg,State!="No")
celltype<-unique(data$cluster)
dat = data.frame(celltype="",ID="",Description="",GeneRatio="",BgRatio="",pvalue=0,p.adjust=0,qvalue=0,geneID="",Count=0)
for (i in celltype) {
  ddeg<-subset(deg,Celltype==i)
  x<- ddeg$gene
  ego_BP  <- enrichGO(x, OrgDb = "org.At.tair.db", keyType="TAIR",pvalueCutoff = 0.05 ,ont="BP")
  if (!is.null(ego_BP)) {
    subType<-data.frame(celltype=i,as.data.frame(ego_BP@result))
  }
  dat=rbind(dat,subType)
}
dat<-dat[-1,]
write.csv(dat,"6-subType_DEGsmutantVSWTDegEnrichment.csv",quote=F,row.names=F)

# // Colors distribution in mesophyll cells
DefaultAssay(leafdata) <- "RNA"
seuratdata <- AddModuleScore(
  object = leafdata,
  features = gene20,   # // gene20 lighthargest genes
  ctrl = 100, #默认值是100
  assay="RNA",
  name = 'Lightharvest2'
)
p1<-FeatureDimPlot(
  srt = leafdata, features = "Lightharvest21",
  reduction = "UMAP", theme_use = "theme_blank", split.by = "Groups",
  ncol=2
p2 <- FeatureStatPlot(leafdata, stat.by = "Lightharvest21", group.by = "cluster", add_stat ="mean")+
  scale_fill_manual(values=c("#1c4b2e", "#e6ccd5","#4780a0","#e0d2a1","#93bfd5","#c1b258"))
p3 <- FeatureStatPlot(leafdata, stat.by = "Lightharvest21", group.by = "cluster", split.by = "Groups",add_stat ="mean")

# ===================== 6.2 Pseudotime analysis on subtypes of leaf ================
data<- readRDS("./leafAnnocluster.RDS")
markers<-read.csv("6-cluster_markers.csv",header=T)
subdata<-subset(data, cluster %in% c("Mc1","Mc2","Mc3"))
markers<-subset(markers, cluster %in% c("Mc1","Mc2","Mc3"))
gene<-unique(markers$gene)
subdata<-subset(x=subdata,downsample=10000)
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
save(cds,file="6-PseudotimeClusterMarkersGenes.rData")

# // Pseudotime Analysis
# Distribution of  mesophyll cells Mc1, Mc2 and Mc3.
plotdf2=as.data.frame(t(cds@reducedDimS))
colnames(plotdf2)=c("component1","component2")
plotdf2$Pseudotime=cds$Pseudotime
plotdf2$cluster<-pbmc$cluster[rownames(plotdf2)]
plotdf2$AT2G25080<-AT2G25080Expr[rownames(plotdf2),]
library(viridis)
library(ggridges)
library(RColorBrewer)
pdf("6-PseudotimeDestributionByCluster.pdf")
ggplot(plotdf2, aes(x=Pseudotime,y=cluster,fill = stat(x))) +
  geom_density_ridges_gradient(scale=1) +
  geom_vline(xintercept = c(5,10),linetype=2)+
  scale_fill_gradientn(name="Pseudotime",colors = colorRampPalette(rev(brewer.pal(10, "Spectral")))(99))+
  scale_y_discrete("")+
  theme_minimal()+
  theme(
    panel.grid = element_blank()
  )
dev.off()

pdf("6-Pseudotime-monocle.pdf")
plot_cell_trajectory(cds,show_branch_points=T,show_tree=F)
plot_cell_trajectory(cds,color_by="Pseudotime",show_branch_points=F,show_tree=F)
plot_cell_trajectory(cds,color_by="Groups",show_branch_points=F,show_tree=F)+scale_color_startrek(alpha=0.6)
plot_cell_trajectory(cds,color_by="cluster",show_branch_points=F,show_tree=F)+scale_color_manual( values=c(pal_d3()(6)[4], pal_d3()(6)[5], pal_d3()(6)[6]))
dev.off()

BEAM_res <- BEAM(cds, branch_point = 1, cores = 1)
BEAM_res <- BEAM_res[order(BEAM_res$qval),]
BEAM_res <- BEAM_res[,c("fd", "pval", "qval")]

b<-plot_genes_branched_heatmap(cds[row.names(subset(BEAM_res,
                                          qval < 1e-4)),],
                                          branch_point = 1,
                                          num_clusters = 3,
                                          cores = 1,
                                          use_gene_short_name = T,
                                          show_rownames = F,
                                          return_heatmap=T)
c1<-subset(b$annotation_row,Cluster==1)
c2<-subset(b$annotation_row,Cluster==2)
c3<-subset(b$annotation_row,Cluster==3)
c1<-as.character(row.names(c1))
genelist<-unique(c1)
c1_go<-enrichGO(genelist,OrgDb = "org.At.tair.db", keyType="TAIR",ont = "BP",pAdjustMethod = "BH",pvalueCutoff = 0.05)
write.table(c1_go,"c1_go.xls",quote=F,sep="\t")
c2<-as.character(row.names(c2))
genelist<-unique(c2)
c2_go<-enrichGO(genelist,OrgDb = "org.At.tair.db", keyType="TAIR",ont = "BP",pAdjustMethod = "BH",pvalueCutoff = 0.05)
write.table(c2_go,"c2_go.xls",quote=F,sep="\t")
c3<-as.character(row.names(c3))
genelist<-unique(c3)
c3_go<-enrichGO(genelist,OrgDb = "org.At.tair.db", keyType="TAIR", ont = "BP",pAdjustMethod = "BH",pvalueCutoff = 0.05)
write.table(c3_go,"c3_go.xls",quote=F,sep="\t")

pdf("6-Psedo-time_ByGroups.pdf",width=12,height=8)
pData(cds)$Groups<-factor(pData(cds)$Groups,levels=c("WT","mutant"))
plot_cell_trajectory(cds,color_by="Groups")+facet_wrap (~Groups)+scale_color_manual(values=c("#395659","#958A67"))
dev.off()

# // expression levels of eight genes related to light harvest along the Pseudotime
gene<-c("AT3G44680","AT3G10985","AT4G02380","AT1G62300","AT4G01250","AT4G13250","AT4G22920","AT5G13800")
pdf("6-plot_genes_in_pseudotime_eightGenes.pdf",height=6,width=8)
plot_genes_in_pseudotime(Agingcds[gene,], ncol =4,color_by = "cluster")+scale_colour_manual(values=c("#EC5C57","#E59DC5", "#1E4129"))
dev.off()




