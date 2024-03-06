# // used for producing QC codes

print("setwd(\"./\")\n")
samples=["WT1","WT2","mutant1","mutant2"]
for sample in samples:
    test='''
    Sys.time()
    {0} <- Read10X(data.dir = "./01.QC/cellRenger/{0}/outs/filtered_feature_bc_matrix")
    Sys.time()

    {0} <- CreateSeuratObject(counts ={0},  project = "{0}", min.cells = 3, min.features = 200)
    {0}[["percent.mito"]] <- PercentageFeatureSet({0}, pattern = "^mt-")   # Arabidopsis thaliana Gene ID conversion gene symbol
    VlnPlot( {0}, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3 )
    {0} <- subset({0}, subset = nFeature_RNA > 200 & nFeature_RNA < 4000 & nCount_RNA > 300 & nCount_RNA < 15000 & percent.mito < 15)
    {0} <- NormalizeData({0}, assay = "RNA", normalization.method = "LogNormalize", scale.factor = 10000 )
    {0} <- FindVariableFeatures({0}, selection.method = "vst", nfeatures = 2000, verbose = FALSE)
    # {0} <- SCTransform({0}, return.only.var.genes = FALSE,assay = "RNA", verbose = FALSE)
    {0} <- ScaleData( {0}, features=rownames({0}) )
    {0} <- RunPCA( object={0}, features=VariableFeatures({0}), npcs=20 )
    {0} <- FindNeighbors({0}, dim=1:20)
    {0} <- FindClusters({0}, resolution = 0.3)
    {0} <- RunTSNE( object = {0}, dims=1:20 )
    {0} <- RunUMAP( object = {0}, dims = 1:20 )
    #DimPlot( object = {0}, reduction="umap" )
    saveRDS( {0}, file = "./{0}.rds" )

    ##  remove doubletfinde
    pcaN <- 20
    train <- {0}
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
    pdf( "{0}.DF.cowplot.pdf", width = 6, height = 5 )
    cowplot::plot_grid(ncol = 2, DimPlot(train, group.by = "orig.ident") + NoAxes(),
                       DimPlot(train, group.by = dropLable) + NoAxes())
    dev.off()
    pdf( "{0}.DF.VlnPlot.pdf", width = 6, height = 5 )
    VlnPlot(train, features = "nFeature_RNA", group.by = dropLable, pt.size = 0.1)
    dev.off()
    train = train[, train@meta.data[, dropLable] == "Singlet"]

    pdf( "{0}.tsne.DF.pdf" )
    DimPlot( object = train, reduction="tsne" )
    dev.off()

    pdf( "{0}.umap.DF.pdf" )
    DimPlot( object = train, reduction="umap", label = T )
    dev.off()

    # save Seurat object
    saveRDS( train, file = "./{0}_DF.rds" )

    '''.format(sample)
    print(test)
