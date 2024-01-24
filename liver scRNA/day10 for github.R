library(Seurat)
library(dplyr)
library(tidyr)
library(ggplot2)
library(patchwork)
library(writexl)

setwd('')

## Read RAW data
ge <- Read10X_h5(file.path("filtered_feature_bc_matrix.h5"))$`Gene Expression`
ac <- Read10X_h5(file.path("filtered_feature_bc_matrix.h5"))$`Antibody Capture`

dim(ge)
dim(ac)


cb <- intersect(colnames(ge), colnames(ac))
ge <- ge[,cb]
ac <- ac[,cb]
dim(ge)

h <- CreateSeuratObject(counts=ge) %>%
  NormalizeData() %>%
  FindVariableFeatures(selection.method="mean.var.plot") %>%
  ScaleData(features=VariableFeatures(.))

h[["HTO"]] <- CreateAssayObject(counts=ac)
h <- NormalizeData(h,assay="HTO",normalization.method="CLR")
h <- HTODemux(h,assay="HTO",positive.quantile=0.99)
Idents(h) <- "HTO_maxID"
saveRDS(h,"seurat-hashed.rds")

as.data.frame(table(h$HTO_maxID)) %>%
  setNames(c("Hashtag","Cells"))


as.data.frame(table(h$HTO_maxID)) %>%
  setNames(c("Hashtag","Cells"))


as.data.frame(table(h$HTO_classification.global)) %>%
  setNames(c("Type","Cells"))


HTOHeatmap(h,assay="HTO") & theme_minimal()


h$sample <- h$HTO_maxID
obj <- subset(h,subset=HTO_classification.global=="Singlet" & sample != "TS-A-0301")

obj$tissue <- dplyr::recode(obj$sample,`TS-A-0303`="Liver",`TS-A-0302`="Aorta")
Idents(obj) <- obj$tissue
obj

head(obj[[]])
table(obj$sample,obj$tissue)
Idents(obj) <- obj$tissue


## QC

obj <- PercentageFeatureSet(obj, "^mt-", col.name = "percent_mito")
obj <- PercentageFeatureSet(obj, "^Rpl.?", col.name = "percent_ribo")

VlnPlot(obj, features = c("nCount_RNA", "nFeature_RNA", "percent_mito", "percent_ribo"), 
        pt.size = 0, ncol = 4) &
  labs(x=NULL) & 
  theme(legend.position="none",
        axis.text.x=element_text(angle=45, hjust=1))


#Violin plots showing distribution of UMIs (RNA counts), number of genes detected, percentage of mitochondrial reads and percentage of ribosomal reads.

### Cell cycle

#Cell cycle phase is estimated using Seurat.


obj <- CellCycleScoring(obj, s.features=stringr::str_to_title(cc.genes$s.genes),
                        g2m.features=stringr::str_to_title(cc.genes$g2m.genes),
                        set.ident=FALSE)
obj$cc_diff <- obj$S.Score-obj$G2M.Score
table(obj$tissue, obj$Phase)


#Number of cells in various phases of cell cycle.

## Filtering


obj <- subset(obj,features = setdiff(rownames(obj),"Gm42418"))


## Normalisation


obj <- NormalizeData(obj, verbose=FALSE, assay="RNA", normalization.method="LogNormalize")
obj <- FindVariableFeatures(obj, selection.method="vst", nfeatures=2000)
obj <- ScaleData(obj, verbose=FALSE, assay="RNA", vars.to.regress="nCount_RNA")


## Dimensionality reduction

set.seed(100)
obj <- RunPCA(obj, assay="RNA", slot="data", verbose = FALSE, seed.use=100)
obj <- RunUMAP(obj, assay="RNA", reduction="pca", dims=1:20, seed.use=100)

p1 <- FeaturePlot(obj, features = "nFeature_RNA")
p2 <- FeaturePlot(obj, features = "nCount_RNA")
p3 <- FeaturePlot(obj, features = "percent_ribo")
p4 <- FeaturePlot(obj, features = "percent_mito")
p5 <- UMAPPlot(obj, group.by = "Phase")

patchwork::wrap_plots(p1,p2,p3,p4,p5,nrow=2,ncol=3) &
  theme(legend.position="top",
        legend.justification=0.5)


## Clustering


obj <- FindNeighbors(obj, reduction = "pca", dims = 1:20, k.param = 10)
obj <- FindClusters(obj, resolution = 0.25)
UMAPPlot(obj)
SaveH5Seurat(obj, filename = "seurat-liver-filtered-aorta-clustered.h5Seurat")
Convert("seurat-liver-filtered-aorta-clustered.h5Seurat", dest = "h5ad")
rm(list = ls())

#Identify cell clusters

Convert("seurat-liver-filtered-aorta-clustered.h5ad", ".h5seurat")
liver <- LoadH5Seurat("seurat-liver-filtered-aorta-clustered.h5Seurat")
DimPlot(liver)

DotPlot(liver,features = c('Apob','Alb','Mat1a','Cps1','Adgre1',
                           'Timd4','Csf1r','S100a8','S100a9','Clec4e','Cx3cr1',
                           'Ms4a7','H2-Eb1','H2-Ab1','Apoe','H2-Aa','Clec4g',
                           'Igfbp7','Kdr','Fosb','Egr1','Cxcl2','Ccl3',
                           'Il1b','Ccl4','Tnf','Cst3','Crip1','Flt3','Ccl5','Trbc2',
                           'Cd3g','Lck','Cd79a','Ebf1','Pax5'))
DotPlot(liver,features = c('Ldlr','Vldlr','Lrp1','Msr1','Marco','Scarb1',
                            'Cd36','Scarb2','Cd68','Cd14','Colec12'))
DotPlot(liver,features = c('C6','Cd5l','Folr2','Grn','Il18bp','Nop58','Inf2','Sdc3','Socs2','Themis2'))

DotPlot(liver,features = c('Cd163','Clec4f','Cd207','Gfra2','Hmox1','Marco','Slc40a1','Tmem26','Vcam1','Vsig4','Id3','Nr1h3','Spic','Mrc1','Esam'))


markers<-FindAllMarkers(liver,only.pos = TRUE, min.pct = 0.1, logfc.threshold =0.5 )
markers
top10_markers <- markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
top10_markers
DoHeatmap(liver, features = top10_markers$gene) + NoLegend()
