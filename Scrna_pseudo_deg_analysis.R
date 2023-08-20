library(MAST)
library(CHETAH)
library(reshape2)
library(Seurat)
library(patchwork)
library(ggplot2)
library(sctransform)
library(tibble)
library(dplyr)

library(Yeskit)
library(topGO)
library(CHETAH)
library(reshape2)
library(Seurat)
library(patchwork)
library(ggplot2)
library(sctransform)
library(tibble)
library(dplyr)
library(scater)
library(loomR)
library(Seurat)
library(patchwork)
library(EnhancedVolcano)


library(Matrix)
install.packages("Matrix")
setwd('D:/base/urinemicrolaura/scrna/newdata_GSE153892/seuratobjects')
list.files()
a=load("SRR12159588_GSM4658185_Sample_01_healthy.RData")
load("SRR12159592_GSM4658186_Sample_01_tumor.RData")
load("SRR12159600_GSM4658188_Sample_02_tumor.RData")
load("SRR12159608_GSM4658190_Sample_03_tumor.RData")


load("SRR12159588_GSM4658185_Sample_01_healthy.RData")
load("SRR12159596_GSM4658187_Sample_02_healthy.RData")
load("SRR12159604_GSM4658189_Sample_03_healthy.RData")

Sample_01_tumor=SRR_object
Sample_02_tumor=SRR_object
Sample_03_tumor=SRR_object

head(Sample_01_tumor@meta.data)

Sample_01_healthy$orig.ident='Sample_01_healthy'

Sample_01_tumor$orig.ident='tumor'
Sample_02_tumor$orig.ident='tumor'
Sample_03_tumor$orig.ident='tumor'

Sample_01_tumor$orig.ident='health'
Sample_02_tumor$orig.ident='health'
Sample_03_tumor$orig.ident='health'



Sample_01_tumor$orig.ident='Sample_01_tumor'
Sample_02_tumor$orig.ident='Sample_02_tumor'
Sample_03_tumor$orig.ident='Sample_03_tumor'

Sample_01_tumor$group='Sample_01_tumor'
Sample_02_tumor$group='Sample_02_tumor'
Sample_03_tumor$group='Sample_03_tumor'

install.packages('Seurat')

library(Seurat)
library(scCustomize)
object_list <- list(Sample_01_tumor,Sample_02_tumor,Sample_03_tumor)
merged_object <- Merge_Seurat_List(list_seurat = object_list)

merged_seurat <- MergeSeurat(object.list = list(Sample_01_tumor,Sample_02_tumor,Sample_03_tumor), 
                             project = c("Sample_01_tumor", "Sample_02_tumor", "Sample_03_tumor"))
pbmc.combined <- merge(Sample_01_tumor,Sample_02_tumor)
merged_seurat <- merge(Sample_01_tumor, y = Sample_02_tumor, add.cell.ids = c("Sample_01_tumor", "Sample_02_tumor"))
combined_seurat <- merge(Sample_01_tumor, y = c(Sample_02_tumor, Sample_03_tumor), add.cell.ids = c("Sample1", "Sample2", "Sample3"))
merged_seurat
head(combined_seurat@meta.data)

head(Sample_01_tumor@meta.data)

seurat1='Sample_01_tumor'
seurat2=seurat1
seurat3=seurat1
pbmc.combined <- 'dsf'
merged_seurat=''
merged_seurat <- merge(seurat1, y = seurat2, add.cell.ids = c("sample1", "sample2"))

# Merge the resulting object with seurat3
combined_seurat <- merge(merged_seurat, y = seurat3, add.cell.ids = c("sample1", "sample2", "sample3"))

cancer=merged_object

cancer=merge(Sample_01_tumor,Sample_02_tumor,Sample_03_tumor)


Sample_01_healthy$group='healthy'
Sample_01_tumor$group='cancer'

m=combined_seurat
m=Sample_01_healthy
meta12=as.data.frame(m@assays[["RNA"]]@data)

metadata=meta12
meta12=''
metadata[is.na(metadata)] <- 0
rownames(metadata) <- gsub("'", "", rownames(metadata))
rownames(metadata) <- gsub("\\[", "", rownames(metadata))
rownames(metadata) <- gsub("\\]", "", rownames(metadata))
rownames(metadata) <- gsub(" ", "_", rownames(metadata))
rownames(metadata) <- gsub("\\(", "", rownames(metadata))
rownames(metadata) <- gsub("\\)", "", rownames(metadata))
metadata_ids <- rownames(metadata)
grep('9372', metadata_ids)
grep('Pseudomonas', metadata_ids)
grep('Aeromonas', metadata_ids)
metadata_ids[33577]
metadata_ids[33541]
dim(metadata)
head(metadata)
metadata[1:4,1:4]
if (nrow(metadata) > 0) {
  object <- Seurat::AddMetaData(object = m, metadata = t(metadata)[,c("Pseudomonas_taxid_3723")],
                                col.name = metadata_ids[33643])
}
if (nrow(metadata) > 0) {
  object <- Seurat::AddMetaData(object = m, metadata = t(metadata)[,c("Aeromonas_taxid_3318")],
                                col.name = metadata_ids[33541])
}


Sample_01_healthy=object
Sample_01_tumor=object
SRR_object='df'
head(cancer1@meta.data)
cancer1=object
cancer1=Sample_01_healthy
cancer1[["percent.mt"]] <- Seurat::PercentageFeatureSet(cancer1, 
                                                        pattern = "^MT-")

head(cancer1)
tail(cancer1)
cancer1<- SetIdent(cancer1, value = "orig.ident")
unique(cancer1$orig.ident)
#cancer1 <- cancer1[!grepl("taxid", rownames(cancer1)), ]
dev.new()
VlnPlot(cancer1, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

cancer1 <- subset(cancer1, subset = nFeature_RNA < 2000 & nCount_RNA < 10000 & percent.mt < 4)
cancer1 <- subset(cancer1, subset = nFeature_RNA < 700 & nCount_RNA < 2500 & percent.mt < 3)
cancer1 <- subset(cancer1, subset = nFeature_RNA < 1500 & nCount_RNA < 10000 & percent.mt < 3)

#cancer1 <- subset(cancer1, subset = percent.mt < 15)
VlnPlot(cancer1, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
cancer1 <- NormalizeData(cancer1)
cancer1 <- FindVariableFeatures(cancer1, selection.method = "vst", nfeatures = 2000)


# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(cancer1), 10)
all.genes <- rownames(cancer1)
cancer1 <- ScaleData(cancer1, features = all.genes)
cancer1 <- RunPCA(cancer1, features = VariableFeatures(object = cancer1))
print(cancer1[["pca"]], dims = 1:5, nfeatures = 5)

cancer1 <- FindNeighbors(cancer1, dims = 1:10)
cancer1 <- FindClusters(cancer1, resolution = 0.5)
head(Idents(cancer1), 5)
cancer1 <- RunUMAP(cancer1, dims = 1:10)
head(cancer1)
top10 <- head(VariableFeatures(cancer1), 10)
cancer1 <- RunTSNE(object = cancer1,check_duplicates = FALSE)
plot1 <- FeatureScatter(cancer1, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(cancer1, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2
plot1 <- VariableFeaturePlot(cancer1)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2


VizDimLoadings(cancer1, dims = 1:2, reduction = "pca")
DimPlot(cancer1, reduction = "pca")
DimHeatmap(cancer1, dims = 1, balanced = TRUE)
DimPlot(cancer1, reduction = "umap")

scrna.sce=cancer1
scrna.sce <- as.SingleCellExperiment(scrna.sce)

counts_data <- assay(scrna.sce)
tsne_data <- reducedDim(scrna.sce, "TSNE")





load('F:/rna/cancer_SRR17327496/ref.rds')
counts_ref <- assay(reference)
celltypes_ref <- reference$celltypes
reference <- SingleCellExperiment(assays = list(counts = counts_ref),
                                  colData = DataFrame(celltypes = celltypes_ref))

input <- SingleCellExperiment(assays = list(counts = counts_data),
                              reducedDims = SimpleList(TSNE = tsne_data))
#input
#input1 <- CHETAHclassifier(input = input, ref_cells = reference)
#input1 <- CHETAHclassifier(input = input, ref_cells = reference,thresh = 1)
input2 <- CHETAHclassifier(input = input, ref_cells = reference,thresh = 0)


PlotCHETAH(input = input2, interm = TRUE)



celltypes <- input2$celltype_CHETAH
#as.data.frame(table(input1$celltype_CHETAH))
#ClassifyReference(ref_cells = reference)
as.data.frame(table(input2$celltype_CHETAH))
typeof(input$celltype_CHETAH)
cancer1$celltype=celltypes
head(cancer1)

pbmc=cancer1
cancer1=''
#DimPlot(pbmc, cols = "Pseudomonas_taxid_3723", label = TRUE, pt.size = 0.5) + NoLegend()
#DimPlot(pbmc, cols = "Pseudomonas_taxid_3723")
#DimPlot(pbmc,feature1 group.by = 'celltype')
unique(pbmc$celltype)

#DimPlot(pbmc, cols = "Pseudomonas_taxid_3723")

#LabelClusters(plot = plot, id = "ident")

head(pbmc)

pbmc.markers <- FindAllMarkers(pbmc, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
head(pbmc.markers)
pbmc.markers %>%
  group_by(cluster) %>%
  slice_max(n = 2, order_by = avg_log2FC)

VlnPlot(pbmc, features = c("Pseudomonas_taxid_3723"), slot = "counts", log = TRUE)
VlnPlot(pbmc, features = c("Aeromonas_taxid_3318"), slot = "counts", log = TRUE)

#VlnPlot(pbmc, features = c("Pseudomonas_sp._CIP-10_taxid_2892442"), slot = "counts", log = TRUE)
head(pbmc@meta.data)
metadata_ids[33559]
metadata_ids[33560]

DimPlot(object = pbmc, split.by = 'ident')


a=FeaturePlot(pbmc, features ="Pseudomonas_taxid_3723")
a=FeaturePlot(pbmc, features ="Aeromonas_taxid_3318")
FeaturePlot(pbmc, features ="nCount_RNA")
FeaturePlot(pbmc, features ="nFeature_RNA")
b=DimPlot(pbmc, reduction = "umap",group.by = 'celltype')
DimPlot(pbmc, reduction = "umap",group.by = 'seurat_clusters')
a+b

head(pbmc@meta.data)
write.csv(pbmc@meta.data,'cancer_combined_pseudo.csv')
write.csv(pbmc@meta.data,'healthy_combined_pseudo.csv')


pbmc.markers %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC) -> top10
DoHeatmap(pbmc, features = top10$gene) + NoLegend()



library(SingleR)
hpca.ref <- celldex::HumanPrimaryCellAtlasData()
monaco.ref <- celldex::MonacoImmuneData()
dice.ref <- celldex::DatabaseImmuneCellExpressionData()
sce <- as.SingleCellExperiment(DietSeurat(pbmc))




monaco.main <- SingleR(test = sce,assay.type.test = 1,ref = monaco.ref,labels = monaco.ref$label.main)
monaco.fine <- SingleR(test = sce,assay.type.test = 1,ref = monaco.ref,labels = monaco.ref$label.fine)
hpca.main <- SingleR(test = sce,assay.type.test = 1,ref = hpca.ref,labels = hpca.ref$label.main)
hpca.fine <- SingleR(test = sce,assay.type.test = 1,ref = hpca.ref,labels = hpca.ref$label.fine)
dice.main <- SingleR(test = sce,assay.type.test = 1,ref = dice.ref,labels = dice.ref$label.main)
dice.fine <- SingleR(test = sce,assay.type.test = 1,ref = dice.ref,labels = dice.ref$label.fine)

table(hpca.main$pruned.labels)
table(dice.main$pruned.labels)
table(monaco.main$pruned.labels)


table(monaco.fine$pruned.labels)
table(hpca.fine$pruned.labels)
table(dice.fine$pruned.labels)
srat=pbmc

srat@meta.data$monaco.main <- monaco.main$pruned.labels
srat@meta.data$monaco.fine <- monaco.fine$pruned.labels

srat@meta.data$hpca.main   <- hpca.main$pruned.labels
srat@meta.data$dice.main   <- dice.main$pruned.labels

srat@meta.data$hpca.fine   <- hpca.fine$pruned.labels
srat@meta.data$dice.fine   <- dice.fine$pruned.labels


srat@meta.data$pseudoinfect <- ifelse(srat@meta.data$Pseudomonas_taxid_3723 >0, "yes", "no")
srat@meta.data$aeroinfect <- ifelse(srat@meta.data$Aeromonas_taxid_3318 >0, "yes", "no")

dev.new()
#FeaturePlot(srat, features ="pseudoinfect")
#FeaturePlot(srat, features ="aeroinfect")
head(srat@meta.data)

table(srat@meta.data$pseudoinfect)

srat <- SetIdent(srat, value = "monaco.fine")

srat <- SetIdent(srat, value = "pseudoinfect")
srat <- SetIdent(srat, value = "aeroinfect")
head(srat@meta.data)



head(FindMarkers(srat, ident.1 = "yes", ident.2 = "no", min.pct = 0.5))
deg=FindMarkers(srat, ident.1 = "yes", ident.2 = "no", min.pct = 0.5)
write.csv(deg,'deg_sample3_cancer_pseudoinfect.csv')


table(srat@meta.data$pseudoinfect)





head(FindMarkers(srat, ident.1 = "Monocytes, CD14+", ident.2 = "Monocytes, CD16+", min.pct = 0.5))
deg=FindMarkers(srat, ident.1 = "Monocytes, CD14+", ident.2 = "Monocytes, CD16+", min.pct = 0.5)
deg
unique(srat@meta.data$monaco.fine)
EnhancedVolcano(deg,
                lab = rownames(deg),
                x = 'avg_log2FC',
                y = 'p_val',
                pointSize = 2.0,
                labSize = 3.0
                )


EnhancedVolcano(deg,
                lab = rownames(deg),
                x = 'avg_log2FC',
                y = 'p_val_adj',
                title = 'PseudoMonas Yes vs No Volcano',
                #pCutoff = 10e-32,
                FCcutoff = 0.25,
                pointSize = 1.0,
                labSize = 3.0
                ,xlim = c(-1, 1))
write.csv(srat@meta.data,'meta_sample1_tumor_dataclusters_final.csv')
dev.new()
srat <- SetIdent(srat, value = "dice.fine")
#dev.new()
DimPlot(srat)
head(srat)
deng <- sc3(srat, ks = 10, biology = TRUE, n_cores = 1)

VlnPlot(pbmc, features = c("Pseudomonas_taxid_3723"), slot = "counts", log = TRUE)
VlnPlot(srat, features = "Pseudomonas_taxid_3723",group.by = "dice.fine")
VlnPlot(pbmc, features = "Pseudomonas_taxid_3723",group.by = "celltype")
VlnPlot(pbmc, features = "Pseudomonas_taxid_3723",group.by = "celltype")

srat <- SetIdent(srat, value = "monaco.fine")
srat <- SetIdent(srat, value = "hpca.fine")
VlnPlot(srat, features = "Pseudomonas_taxid_3723",group.by = "celltype")
VlnPlot(srat, features = "Aeromonas_taxid_3318",group.by = "monaco.main")
dev.new()
VlnPlot(srat, features = "Pseudomonas_taxid_3723",group.by = "monaco.fine")+ NoLegend()
VlnPlot(srat, features = "Pseudomonas_taxid_3723",group.by = "hpca.main")
VlnPlot(srat, features = "Pseudomonas_taxid_3723",group.by = "hpca.fine")+ NoLegend()
VlnPlot(srat, features = "Pseudomonas_taxid_3723",group.by = "dice.main")
VlnPlot(srat, features = "Pseudomonas_taxid_3723",group.by = "dice.fine")



srat <- SetIdent(srat, value = "dice.fine")
dev.new()
head(FindMarkers(srat, ident.1 = "Monocytes, CD14+", ident.2 = "Monocytes, CD16+", min.pct = 0.5))



write.csv(srat@meta.data,'ab.csv')

table(srat@meta.data$pseudoinfect)

srat@meta.data$pseudmonaco.fine =paste(srat@meta.data$monaco.fine,srat@meta.data$pseudoinfect)
srat@meta.data$pseuddice.fine =paste(srat@meta.data$dice.fine,srat@meta.data$pseudoinfect)
table(srat@meta.data$pseuddice.fine)

save(srat, file = "srat_combinedhealthy.RData")

srat <- SetIdent(srat, value = "pseuddice.fine")
deg=FindMarkers(srat, ident.1 = "T cells, CD4+, Th1 yes", ident.2 = "T cells, CD4+, Th1 no", min.pct = 0.5)
head(deg)











srat@meta.data$pseudhpca.fine =paste(srat@meta.data$hpca.fine,srat@meta.data$pseudoinfect)

table(srat@meta.data$pseudmonaco.fine)
table(srat@meta.data$pseudhpca.fine)

prop.table(table(srat@meta.data$aeromonaco.fine)) * 100
