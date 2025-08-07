## Import the packages used for single cell RNA seq analysis 
library(Seurat)
library(SeuratObject)
library(Matrix)
library(dplyr)
library(readr)
library(magrittr)
library(ggplot2)
library(patchwork)
library(pheatmap)
library(RColorBrewer)
library(clusterProfiler)
library(enrichplot)
library(DOSE)
library(ReactomePA)
library(org.Hs.eg.db)
library(AnnotationDbi)
library(pathview)
library(EnhancedVolcano)
library(tidygraph)
library(ggraph)

########################################## DATASET IMPORT ##################################################
## datasets with introns
data_H1 <- '/media/ming/Extra_SSD_4TB1/Ishika_sc/Datasets/Ishika_islets_H1/outs/filtered_feature_bc_matrix'
data_F508D <- '/media/ming/Extra_SSD_4TB1/Ishika_sc/Datasets/Ishika_islets_F508del/outs/filtered_feature_bc_matrix'
data_G542X <- '/media/ming/Extra_SSD_4TB1/Ishika_sc/Datasets/Ishika_islets_G542X/outs/filtered_feature_bc_matrix'

datasets <- list(data_H1, data_F508D, data_G542X)

## datasets without introns
data_H1_noIntron <- '/media/ming/Extra_SSD_4TB1/Ishika_sc/Datasets/Ishika_islets_H1_nointron/outs/filtered_feature_bc_matrix'
data_F508d_noIntron <- '/media/ming/Extra_SSD_4TB1/Ishika_sc/Datasets/Ishika_islets_F508del_nointron/outs/filtered_feature_bc_matrix'
data_G542X_noIntron <- '/media/ming/Extra_SSD_4TB1/Ishika_sc/Datasets/Ishika_islets_G542X_nointron/outs/filtered_feature_bc_matrix'
  
datasets_noIntron <- list(data_H1_noIntron, data_F508d_noIntron, data_G542X_noIntron)
datasets_noIntron_names <- c("H1", "F508d", "G542X")
obj_noIntron_names <- c("obj_H1_noIntron", "obj_F508d_noIntron", "obj_G542X_noIntron")

####################################### QUALITY CONTROL ######################################
out <- Read10X(data_G542X)
object <- CreateSeuratObject(counts = out, project = "G542X", min.cells = 1, min.features = 1)
object[["percent.mt"]] <- PercentageFeatureSet(object, pattern = "^MT")
Vlnplot <- VlnPlot(object, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
Vlnplot
ggsave(filename = "Vlnplot1_G542X.pdf", 
       path = '/media/ming/Extra_SSD_4TB1/Ishika_sc/pilot2_figures', dpi = 300,
       device = "pdf", width = 8, height = 8, units = "in")
plot1 <- FeatureScatter(object, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(object, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
FeaturesScatterPlot <- plot1 + plot2
FeaturesScatterPlot
ggsave(filename = "Vlnplot2_G542X.pdf", 
       path = '/media/ming/Extra_SSD_4TB1/Ishika_sc/pilot2_figures', dpi = 300,
       device = "pdf", width = 8, height = 8, units = "in")

obj_G542X <- object

saveRDS(obj_H1, file = '/media/ming/Extra_SSD_4TB1/Ishika_sc/Datasets/pilot2_H1_raw.rds')
saveRDS(obj_F508D, file = '/media/ming/Extra_SSD_4TB1/Ishika_sc/Datasets/pilot2_F508del_raw.rds')
saveRDS(obj_G542X, file = '/media/ming/Extra_SSD_4TB1/Ishika_sc/Datasets/pilot2_G542X_raw.rds')

obj_H1 <- readRDS('/media/ming/Extra_SSD_4TB1/Ishika_sc/Datasets/pilot2_H1_raw.rds')
obj_F508D <- readRDS('/media/ming/Extra_SSD_4TB1/Ishika_sc/Datasets/pilot2_F508del_raw.rds')
obj_G542X <- readRDS('/media/ming/Extra_SSD_4TB1/Ishika_sc/Datasets/pilot2_G542X_raw.rds')

write.table(obj_H1@meta.data, file = '/media/ming/Extra_SSD_4TB1/Ishika_sc/Tables/pilot2/wIntrons/pilot2_H1_metatdata_raw.csv', sep = ',')
write.table(obj_F508D@meta.data, file = '/media/ming/Extra_SSD_4TB1/Ishika_sc/Tables/pilot2/wIntrons/pilot2_F508d_metatdata_raw.csv', sep = ',')
write.table(obj_G542X@meta.data, file = '/media/ming/Extra_SSD_4TB1/Ishika_sc/Tables/pilot2/wIntrons/pilot2_G542X_metatdata_raw.csv', sep = ',')

# H1: 7267
# F508del: 8899
# G542X: 8847 

i = 1
for (x in datasets_noIntron){
  
  obj_name <- datasets_noIntron_names[i]
  obj_to_save <- obj_noIntron_names[i]
  pathname <- '/media/ming/Extra_SSD_4TB1/Ishika_sc/pilot2_figures_nointron'
  
  out <- Read10X(x)
  object <- CreateSeuratObject(counts = out, project = obj_name, min.cells = 1, min.features = 1)
  object[["percent.mt"]] <- PercentageFeatureSet(object, pattern = "^MT")
  Vlnplot <- VlnPlot(object, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
  Vlnplot

  ggsave(file = paste0("VlnPlot_",obj_name,".pdf"), 
         path = pathname,
         dpi = 300,
         device = "pdf",
         width = 8, 
         height = 8, 
         units = "in")
  
  plot1 <- FeatureScatter(object, feature1 = "nCount_RNA", feature2 = "percent.mt")
  plot2 <- FeatureScatter(object, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
  
  FeaturesScatterPlot <- plot1 + plot2
  FeaturesScatterPlot
  
  ggsave(file = paste0("VlnPlot2_",obj_name,".pdf"), 
         path = pathname, 
         dpi = 300,
         device = "pdf", 
         width = 8, 
         height = 8, 
         units = "in")
  
  assign(obj_to_save, get("object"))
  
  i = i + 1
  
}

saveRDS(obj_H1_noIntron, file = '/media/ming/Extra_SSD_4TB1/Ishika_sc/Datasets/pilot2_H1_noIntron_raw.rds')
saveRDS(obj_F508d_noIntron, file = '/media/ming/Extra_SSD_4TB1/Ishika_sc/Datasets/pilot2_F508del_noIntron_raw.rds')
saveRDS(obj_G542X_noIntron, file = '/media/ming/Extra_SSD_4TB1/Ishika_sc/Datasets/pilot2_G542X_noIntron_raw.rds')

# H1: 7404
# F508del: 8899
# G542X: 9104

####################################### SUBSETTING ######################################
hist(object@meta.data$percent.mt, breaks = 100, xlim = c(0,100)) # 7.5
hist(obj_H1@meta.data$nCount_RNA, breaks = 1000, xlim = c(0,20000)) # > 2000
hist(object@meta.data$nFeature_RNA, breaks = 1000, xlim = c(0,8000)) # 200 - 10000

obj_H1@meta.data %>%
  filter(percent.mt < 7.5) %>%
  count() # 6836, 94.7%

obj_F508D@meta.data %>%
  filter(percent.mt < 7.5) %>%
  count() # 8222, 92.4%

obj_G542X@meta.data %>%
  filter(percent.mt < 7.5) %>%
  count() # 7724, 87.3%

obj_H1s <- subset(obj_H1, subset = nFeature_RNA > 200 & percent.mt < 7.5 & nFeature_RNA < 10000 & nCount_RNA > 2000)
obj_F508ds <- subset(obj_F508D, subset = nFeature_RNA > 200 & percent.mt < 7.5 & nFeature_RNA < 10000 & nCount_RNA > 2000)
obj_G542Xs <- subset(obj_G542X, subset = nFeature_RNA > 200 & percent.mt < 7.5 & nFeature_RNA < 10000 & nCount_RNA > 2000)

obj_H1s_nointron <- subset(obj_H1_noIntron, subset = nFeature_RNA > 200 & percent.mt < 7.5 & nFeature_RNA < 10000 & nCount_RNA > 2000)
obj_F508ds_nointron <- subset(obj_F508d_noIntron, subset = nFeature_RNA > 200 & percent.mt < 7.5 & nFeature_RNA < 10000 & nCount_RNA > 2000)
obj_G542Xs_nointron <- subset(obj_G542X_noIntron, subset = nFeature_RNA > 200 & percent.mt < 7.5 & nFeature_RNA < 10000 & nCount_RNA > 2000)

## overlaps 
H1_overlap <- intersect(rownames(obj_H1s@meta.data), rownames(obj_H1s_nointron@meta.data))
F508d_overlap <- intersect(rownames(obj_F508ds@meta.data), rownames(obj_F508ds_nointron@meta.data))
G542X_overlap <- intersect(rownames(obj_G542Xs@meta.data), rownames(obj_G542Xs_nointron@meta.data))

H1_o1 <- subset(obj_H1s, cells = H1_overlap)
H1_o2 <- subset(obj_H1s_nointron, cells = H1_overlap)
H1_diff <- (H1_o1@meta.data$nCount_RNA - H1_o2@meta.data$nCount_RNA)/H1_o1@meta.data$nCount_RNA
mean(H1_diff)
sd(H1_diff)

F508d_o1 <- subset(obj_F508ds, cells = F508d_overlap)
F508d_o2 <- subset(obj_F508ds_nointron, cells = F508d_overlap)
F508d_diff <- (F508d_o1@meta.data$nCount_RNA - F508d_o2@meta.data$nCount_RNA)/F508d_o1@meta.data$nCount_RNA
mean(F508d_diff)
sd(F508d_diff)

G542X_o1 <- subset(obj_G542Xs, cells = G542X_overlap)
G542X_o2 <- subset(obj_G542Xs_nointron, cells = G542X_overlap)
G542X_diff <- (G542X_o1@meta.data$nCount_RNA - G542X_o2@meta.data$nCount_RNA)/G542X_o1@meta.data$nCount_RNA
mean(G542X_diff)
sd(G542X_diff)

####################################### DATASET INTEGRATION ############################################
objs <- list(obj_H1s, obj_F508ds, obj_G542Xs)
# 6276, 7192, 6051

obj_H1s <- NormalizeData(obj_H1s)
obj_H1s <- FindVariableFeatures(obj_H1s, selection.method = "vst", nfeatures = 2000)

obj_F508ds <- NormalizeData(obj_F508ds)
obj_F508ds <- FindVariableFeatures(obj_F508ds, selection.method = "vst", nfeatures = 2000)

obj_G542Xs <- NormalizeData(obj_G542Xs)
obj_G542Xs <- FindVariableFeatures(obj_G542Xs, selection.method = "vst", nfeatures = 2000)

## because we saw some duplicated cell names
intersect(rownames(obj_H1s@meta.data), rownames(obj_F508ds@meta.data))
intersect(rownames(obj_H1s@meta.data), rownames(obj_G542Xs@meta.data))
intersect(rownames(obj_F508ds@meta.data), rownames(obj_G542Xs@meta.data))

barcodes_to_rm1 <- intersect(rownames(obj_H1s@meta.data), rownames(obj_F508ds@meta.data))
barcodes_to_rm2 <- intersect(rownames(obj_H1s@meta.data), rownames(obj_G542Xs@meta.data))
barcodes_to_rm3 <- intersect(rownames(obj_F508ds@meta.data), rownames(obj_G542Xs@meta.data))

obj_H1s <- subset(obj_H1s, cells = setdiff(Cells(obj_H1s), barcodes_to_rm1))
obj_H1s <- subset(obj_H1s, cells = setdiff(Cells(obj_H1s), barcodes_to_rm2))
obj_F508ds <- subset(obj_F508ds, cells = setdiff(Cells(obj_F508ds), barcodes_to_rm3))

# select features that are repeatedly variable across datasets for integration
features <- SelectIntegrationFeatures(object.list = objs)
anchors <- FindIntegrationAnchors(object.list = objs, anchor.features = features)
data_int <- IntegrateData(anchorset = anchors)

data_int <- ScaleData(data_int, verbose = FALSE)
data_int <- RunPCA(data_int, npcs = 30, verbose = FALSE)
data_int <- FindNeighbors(data_int, reduction = "pca", dims = 1:30)
data_int <- FindClusters(data_int, resolution = 0.5)
data_int <- RunUMAP(data_int, reduction = "pca", dims = 1:30, seed.use = 1)
data_int <- RunTSNE(data_int, reduction = "pca", dims = 1:30, seed.use = 1)

write.table(data_int@meta.data, file = '/media/ming/Extra_SSD_4TB1/Ishika_sc/pilot2_data_int_metadata.csv', sep = ',')

# without introns

# 6071 7192 5461
obj_H1s_nointron <- NormalizeData(obj_H1s_nointron)
obj_H1s_nointron <- FindVariableFeatures(obj_H1s_nointron, selection.method = "vst", nfeatures = 2000)

obj_F508ds_nointron <- NormalizeData(obj_F508ds_nointron)
obj_F508ds_nointron <- FindVariableFeatures(obj_F508ds_nointron, selection.method = "vst", nfeatures = 2000)

obj_G542Xs_nointron<- NormalizeData(obj_G542Xs_nointron)
obj_G542Xs_nointron <- FindVariableFeatures(obj_G542Xs_nointron, selection.method = "vst", nfeatures = 2000)

## because we saw some duplicated cell names
intersect(rownames(obj_H1s_nointron@meta.data), rownames(obj_F508ds_nointron@meta.data))
intersect(rownames(obj_H1s_nointron@meta.data), rownames(obj_G542Xs_nointron@meta.data))
intersect(rownames(obj_F508ds_nointron@meta.data), rownames(obj_G542Xs_nointron@meta.data))

barcodes_to_rm1 <- intersect(rownames(obj_H1s_nointron@meta.data), rownames(obj_F508ds_nointron@meta.data))
barcodes_to_rm2 <- intersect(rownames(obj_H1s_nointron@meta.data), rownames(obj_G542Xs_nointron@meta.data))
barcodes_to_rm3 <- intersect(rownames(obj_F508ds_nointron@meta.data), rownames(obj_G542Xs_nointron@meta.data))

obj_H1s_nointron <- subset(obj_H1s_nointron, cells = setdiff(Cells(obj_H1s_nointron), barcodes_to_rm1))
obj_H1s_nointron <- subset(obj_H1s_nointron, cells = setdiff(Cells(obj_H1s_nointron), barcodes_to_rm2))
obj_F508ds_nointron <- subset(obj_F508ds_nointron, cells = setdiff(Cells(obj_F508ds_nointron), barcodes_to_rm3))
objs_nointron <- list(obj_H1s_nointron, obj_F508ds_nointron, obj_G542Xs_nointron)
# select features that are repeatedly variable across datasets for integration
features <- SelectIntegrationFeatures(object.list = objs_nointron)
anchors <- FindIntegrationAnchors(object.list = objs_nointron, anchor.features = features)
data_int_nointron <- IntegrateData(anchorset = anchors)

data_int_nointron <- ScaleData(data_int_nointron, verbose = FALSE)
data_int_nointron <- RunPCA(data_int_nointron, npcs = 30, verbose = FALSE)
data_int_nointron <- FindNeighbors(data_int_nointron, reduction = "pca", dims = 1:30)
data_int_nointron <- FindClusters(data_int_nointron, resolution = 0.5)
data_int_nointron <- RunUMAP(data_int_nointron, reduction = "pca", dims = 1:30, seed.use = 1)
data_int_nointron <- RunTSNE(data_int_nointron, reduction = "pca", dims = 1:30, seed.use = 1)

write.table(data_int_nointron@meta.data, file = '/media/ming/Extra_SSD_4TB1/Ishika_sc/pilot2_data_int_nointron_metadata.csv', sep = ',')
saveRDS(data_int_nointron, file = '/media/ming/Extra_SSD_4TB1/Ishika_sc/Datasets/pilot2_int_data_nointron_unannotated.rds')

############################### INTEGRATED DATASETS IMPORTING ############################
# Import the integrated datasets with and without introns 
data_int_ni <- readRDS('/media/ming/Extra_SSD_4TB1/Ishika_sc/Datasets/pilot2_int_data_nointron_unannotated.rds')
data_int <- readRDS('/media/ming/Extra_SSD_4TB1/Ishika_sc/Datasets/pilot2_int_data_annotated.rds')

############################### CELL ANNOTATION ###########################################
# lists of signiture gene markers for cell subpopulation annotation, not all will be remained 
features_prealpha <- c("INS","GCG")
features_alpha <- c("GCG","DPP4","INS","ARX", "IRX1", "IRX2")
features_beta <- c("PCSK1","INS","NKX6-1","IAPP")
features_mature_beta <- c("INS","G6PC2","MAFA","SIX3")
features_EC <- c("TPH1","DDC","LMX1A","SLC18A1","FEV")
features_delta <- c("SST","HHEX","LEPR")
features_epsilon <- c("GHRL","ASCL1", "PHGR1")
features_gamma <- c("PPY", "AQP3", "ID2")
features_Endocrine <- "CHGA"
features_Exocrine <- c("KRT19","KRT7")
features_neuroendocrine <- c("GAP43","RTN1","CNTNAP2")
features_pancPRO <- c("PDX1","NKX6-1","NEUROG3","NKX2-2","SOX4")
features_stress <- c("DDIT3","CIB1")
features_proliferative <- c("TOP2A","CENPF","PCNA","MKI67")

feature_list <- list(features_prealpha,features_alpha,features_beta,features_mature_beta,features_EC,features_delta,features_epsilon,
                     features_gamma,features_Endocrine,features_Exocrine,features_neuroendocrine,features_pancPRO,features_stress,features_proliferative)

feature_list_nam <- c("features_prealpha","features_alpha","features_beta","features_mature_beta","features_EC","features_delta","features_epsilon",
                      " features_gamma","features_Endocrine","features_Exocrine","features_neuroendocrine","features_pancPRO","features_stress",
                      "features_proliferative")

# Output the feature plot for cell annotation
i = 1
for (x in feature_list){
  
  FeaturePlot(data_int, features = x, pt.size = 1, slot = 'data', reduction = 'tsne')
  
  nam <- paste0('/media/ming/Extra_SSD_4TB1/Ishika_sc/pilot2_figures/',feature_list_nam[i],"_nci_p1_p2.pdf")
  i = i + 1
  
  ggsave(filename = nam,
         device = "pdf", width = 8, height = 8, units = "in") 
}

# cluster-based annotation 
data_int <- RenameIdents(
    object = data_int,
    `0` = 'alpha-2',
    `1` = 'beta',
    `2` = 'alpha-1',
    `3` = 'alpha-3',
    `4` = 'beta-like-LC',
    `5` = 'EC1',
    `6` = 'alpha-like-2-LC',
    `7` = 'EC2',
    `8` = 'proliferative',
    `9` = 'Unknown-1-LC',
    `10` = 'delta',
    `11` = 'Unknown-2-nH1',
    `12` = 'alpha-like3-nG',
    `13` = 'exocrine',
    `14` = 'alpha-like-1-LC'
)

data_int <- RenameIdents(
  object = data_int,
  "alpha-1" = "Alpha",
  "alpha-2" = "Alpha",
  "alpha-3" = "Alpha",
  "alpha-like-1-LC" = "LowCount",
  "alpha-like3-nG" = "Alpha",
  "alpha-like-2-LC" = "LowCount",
  "beta-like-LC" = "LowCount",
  "Unknown-1-LC" = "LowCount",
  "Unknown-2-nH1" = "MSC",
  "EC1" = "EC",
  "EC2" = "EC",
  "delta" = "Delta",
  "beta" = "Beta",
  "exocrine" = "Exocrine",
  "proliferative" = "Proliferative"
  
)

############################################# CFTR READ COUNT #########################################
# Count CFTR transcript 
indMyGene<-"CFTR"
CountMyGeneMyCell <- GetAssayData(obj_F508d_p2s, slot = "counts")[indMyGene,]
write.table(CountMyGeneMyCell, sep=',', file = '/media/ming/Extra_SSD_4TB1/Ishika_sc/CFTR_counts_F508d_p2_ori.csv')
CountMyGeneMyCell <- GetAssayData(obj_G542X_p2s, slot = "counts")[indMyGene,]
write.table(CountMyGeneMyCell, sep=',', file = '/media/ming/Extra_SSD_4TB1/Ishika_sc/CFTR_counts_G542X_p2_ori.csv')
CountMyGeneMyCell <- GetAssayData(obj_H1_p2s, slot = "counts")[indMyGene,]
write.table(CountMyGeneMyCell, sep=',', file = '/media/ming/Extra_SSD_4TB1/Ishika_sc/CFTR_counts_H1_p2_ori.csv')

CountMyGeneMyCell <- GetAssayData(obj_F508d_p2s_nci, slot = "counts")[indMyGene,]
write.table(CountMyGeneMyCell, sep=',', file = '/media/ming/Extra_SSD_4TB1/Ishika_sc/CFTR_counts_F508d_p2_nci.csv')
CountMyGeneMyCell <- GetAssayData(obj_G542X_p2s_nci, slot = "counts")[indMyGene,]
write.table(CountMyGeneMyCell, sep=',', file = '/media/ming/Extra_SSD_4TB1/Ishika_sc/CFTR_counts_G542X_p2_nci.csv')
CountMyGeneMyCell <- GetAssayData(obj_H1_p2s_nci, slot = "counts")[indMyGene,]
write.table(CountMyGeneMyCell, sep=',', file = '/media/ming/Extra_SSD_4TB1/Ishika_sc/CFTR_counts_H1_p2_nci.csv')

## manually re-organized the matrix
CFTR_H1_p2    <- read_csv('/media/ming/Extra_SSD_4TB1/Ishika_sc/Tables/H1_p2_ori_metadata_wCFTRgenotyping.csv')
CFTR_F508d_p2 <- read_csv('/media/ming/Extra_SSD_4TB1/Ishika_sc/Tables/F508d_p2_ori_metadata_wCFTRgenotyping.csv')
CFTR_G542X_p2 <- read_csv('/media/ming/Extra_SSD_4TB1/Ishika_sc/Tables/G542X_p2_ori_metadata_wCFTRgenotyping.csv') 

CFTR_int <- rbind(CFTR_H1_p2, CFTR_F508d_p2)
CFTR_int <- rbind(CFTR_int, CFTR_G542X_p2)
CFTR_int_unique <- CFTR_int[!duplicated(CFTR_int$bc), ]

data_int@meta.data  <- data_int@meta.data %>%
  mutate(CFTR_nointron = 0, CFTR_ori = 0, bc = rownames(data_int@meta.data))

data_int@meta.data    <- data_int@meta.data %>%
  left_join(CFTR_int_unique %>% select(bc, CFTR_ori), by = "bc") %>%
  mutate(CFTR_ori = coalesce(CFTR_ori.y, CFTR_ori.x)) %>%
  select(-CFTR_ori.x, -CFTR_ori.y) %>%
  arrange(bc) 
rownames(data_int@meta.data) <- data_int@meta.data$bc

data_int@meta.data    <- data_int@meta.data %>%
  left_join(CFTR_int_unique %>% select(bc, CFTR_nointron), by = "bc") %>%
  mutate(CFTR_nointron = coalesce(CFTR_nointron.y, CFTR_nointron.x)) %>%
  select(-CFTR_nointron.x, -CFTR_nointron.y) %>%
  arrange(bc) 
rownames(data_int@meta.data) <- data_int@meta.data$bc

data_int@meta.data  <- data_int@meta.data %>%
  mutate(
         CFTR_nointron_normalized = log(CFTR_nointron/nCount_RNA*10000+1), 
         CFTR_ori_normalized = log(CFTR_ori/nCount_RNA*10000+1))

data_int@meta.data <- data_int@meta.data %>%
  arrange(orig.ident,bc)

FeaturePlot(data_int, reduction = 'tsne', 
            features = "CFTR_ori_normalized",
            split.by = "orig.ident",
            label = TRUE,
            pt.size = 1,
            cols = c("#cccccc","red"), 
            order = TRUE,
            min.cutoff = 0,
            keep.scale = "all")

+
  scale_colour_gradient2(low = 'lightgrey', mid = 'blue',high = 'red', midpoint = 0.2)


# Reorder the graphs to match the cell order
data_int@graphs$integrated_nn <- as.Graph(data_int@graphs$integrated_nn[colnames(data_int), colnames(data_int)])
data_int@graphs$integrated_snn <- as.Graph(data_int@graphs$integrated_snn[colnames(data_int), colnames(data_int)])

saveRDS(data_int, file = '/media/ming/Extra_SSD_4TB1/Ishika_sc/Datasets/pilot2_int_annotated_CFTR.rds')


## Visualization

data_int <- subset(data_int, idents = c("Alpha", "Beta","EC","Delta","Exocrine","Proliferative","Unknown"))
Idents(data_int) <- factor(Idents(data_int), levels = c("Alpha", "Beta", "Delta","EC","Exocrine","Proliferative","Unknown"))

cluster_colors <- c(
  'Alpha' = '#F3766E',
  'Beta' = '#1CBDC2',
  'EC' = '#53B448',
  'Delta' = '#9589C1',
  'Exocrine' = '#C39B2D',
  'Proliferative' = '#2BB892',
  'Unknown' = '#D66DAB' #MSCs
)

clone_colors = c(
  "F508D" = '#EE3524',
  "G542X" = "#3D58A7",
  "H1" = "#31B44A"
)

DimPlot(data_ints, reduction = 'tsne', pt.size = 1, cols = cluster_colors)

## Heatmap
selected_genes_heatmap <-  c("GCG","IRX1","INS","DLK1","SST","HHEX","TPH1","SLC18A1","SPP1","KRT19","MKI67","PCNA","VIM","COL3A1")

avg_exp <- AverageExpression(data_ints, assays = "RNA")
avg_matrix <- avg_exp$RNA  # gene-by-cluster matrix
avg_matrix <- avg_matrix[selected_genes_heatmap, ]

pheatmap(avg_matrix, scale = "row", 
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(100))

## ghi CFTR+ (no intron) cells visualization 
data_ints_CFTR <- subset(data_ints, subset = CFTR_nointron > 0)

DimPlot(
  object = data_ints_CFTR,
  reduction = 'tsne',
  cols = cluster_colors,
  #label = TRUE,
  repel = TRUE,
  pt.size = 1,
  split.by = 'orig.ident'
)

data_ints_CFTR@meta.data %>%
  filter(CFTR_nointron >0) %>%
  group_by(orig.ident) %>%
  count()

data_ints@meta.data %>%
  group_by(orig.ident) %>%
  count()


DimPlot(subset(data_ints, idents = "Alpha"), 
        reduction = 'tsne', 
        pt.size = 0.5, 
        group.by = "orig.ident",
        cols = clone_colors,
        shuffle = TRUE)

p<-DimPlot(subset(data_ints, idents = "Beta"), 
        reduction = 'tsne', 
        pt.size = 0.5, 
        group.by = "orig.ident",
        cols = clone_colors)
selected_cells <- CellSelector(plot = p)

DimPlot(subset(data_ints, idents = "Beta"), 
        reduction = 'tsne', 
        pt.size = 0.5, 
        group.by = "orig.ident",
        cells = selected_cells,
        cols = clone_colors,
        shuffle = TRUE)

########################################## DEG ###################################
DefaultAssay(data_ints) <- "RNA"
data_ints <- JoinLayers(data_ints, reorder = TRUE)
data.markers.alpha.G2H <- FindMarkers(subset(data_ints, idents = "Alpha"),
                            ident.1 = "G542X",
                            ident.2 = "H1",
                            group.by = 'orig.ident',
                            logfc.threshold = 0.5, # note this 
                            test.use = "MAST")

data.markers.alpha.F2H <- FindMarkers(subset(data_ints, idents = "Alpha"),
                                    ident.1 = "F508D",
                                    ident.2 = "H1",
                                    group.by = 'orig.ident',
                                    logfc.threshold = 0.5,
                                    test.use = "MAST")

data.markers.alpha.G2F <- FindMarkers(subset(data_ints, idents = "Alpha"),
                                      ident.1 = "G542X",
                                      ident.2 = "F508D",
                                      group.by = 'orig.ident',
                                      logfc.threshold = 0.5,
                                      test.use = "MAST")

data.markers.beta.G2H <- FindMarkers(subset(data_ints, idents = "Beta"),
                                      ident.1 = "G542X",
                                      ident.2 = "H1",
                                      group.by = 'orig.ident',
                                      logfc.threshold = 0.5,
                                      test.use = "MAST")

data.markers.beta.F2H <- FindMarkers(subset(data_ints, idents = "Beta"),
                                      ident.1 = "F508D",
                                      ident.2 = "H1",
                                      group.by = 'orig.ident',
                                      logfc.threshold = 0.5,
                                      test.use = "MAST")

data.markers.beta.G2F <- FindMarkers(subset(data_ints, idents = "Beta"),
                                      ident.1 = "G542X",
                                      ident.2 = "F508D",
                                      group.by = 'orig.ident',
                                      logfc.threshold = 0.5,
                                      test.use = "MAST")

# for GSEA later
data.markers.alpha.G2H0 <- FindMarkers(subset(data_ints, idents = "Alpha"),
                                      ident.1 = "G542X",
                                      ident.2 = "H1",
                                      group.by = 'orig.ident',
                                      logfc.threshold = 0, # note this 
                                      test.use = "MAST")

write.table(data.markers.alpha.G2H, file = '/media/ming/Extra_SSD_4TB1/Ishika_sc/DEG/DEG_G2H_alpha_v2.csv', sep = ',')
write.table(data.markers.alpha.F2H, file = '/media/ming/Extra_SSD_4TB1/Ishika_sc/DEG/DEG_F2H_alpha_v2.csv', sep = ',')
write.table(data.markers.alpha.G2F, file = '/media/ming/Extra_SSD_4TB1/Ishika_sc/DEG/DEG_G2F_alpha_v2.csv', sep = ',')
write.table(data.markers.beta.G2H, file = '/media/ming/Extra_SSD_4TB1/Ishika_sc/DEG/DEG_G2H_beta_v2.csv', sep = ',')
write.table(data.markers.beta.F2H, file = '/media/ming/Extra_SSD_4TB1/Ishika_sc/DEG/DEG_F2H_beta_v2.csv', sep = ',')
write.table(data.markers.beta.G2F, file = '/media/ming/Extra_SSD_4TB1/Ishika_sc/DEG/DEG_G2F_beta_v2.csv', sep = ',')

########################################## DEG VOLCANO PLOT ########################################## 

# gene of interest to label
go_term <- c("GO:0030198","GO:0030073","GO:0009888","GO:0070091","GO:0005179","GO:0000280","GO:0000278","GO:0001653","GO:0007190")

# Get Entrez Gene IDs associated with selected GO term
gene_ids <- AnnotationDbi::select(
  org.Hs.eg.db,
  keys = go_term,
  columns = c("SYMBOL", "ENTREZID"),
  keytype = "GO"
)

volcano.alpha.F2H <- EnhancedVolcano(data.markers.alpha.F2H,
                x = "avg_log2FC",
                y = "p_val_adj",
                selectLab = intersect(
                  rownames(data.markers.alpha.F2H)[
                    data.markers.alpha.F2H$p_val_adj < 0.05 & abs(data.markers.alpha.F2H$avg_log2FC) > 1
                ],   
                  gene_ids$SYMBOL),
                lab = rownames(data.markers.alpha.F2H),
                pCutoff = 5e-02,
                FCcutoff = 1,
                pointSize = 1,
                labSize = 5,
                #labCol = c('#E64B35'),
                col = c('grey30', 'forestgreen', 'royalblue', '#E64B35'),
                colGradientBreaks = 0,
                arrowheads = FALSE,
                gridlines.minor = FALSE,
                gridlines.major = FALSE,
                drawConnectors = TRUE)

EnhancedVolcano(data.markers.alpha.G2H,
                x = "avg_log2FC",
                y = "p_val_adj",
                selectLab = intersect(
                  rownames(data.markers.alpha.G2H)[
                    data.markers.alpha.G2H$p_val_adj < 0.05 & abs(data.markers.alpha.G2H$avg_log2FC) > 1
                  ],   
                  gene_ids$SYMBOL),
                lab = rownames(data.markers.alpha.G2H),
                pCutoff = 5e-02,
                FCcutoff = 1,
                pointSize = 1,
                labSize = 5,
                #labCol = c('#E64B35'),
                col = c('grey30', 'forestgreen', 'royalblue', '#E64B35'),
                colGradientBreaks = 0,
                arrowheads = FALSE,
                gridlines.minor = FALSE,
                gridlines.major = FALSE,
                drawConnectors = TRUE)


EnhancedVolcano(data.markers.alpha.G2F,
                x = "avg_log2FC",
                y = "p_val_adj",
                selectLab = intersect(
                  rownames(data.markers.alpha.G2F)[
                    data.markers.alpha.G2F$p_val_adj < 0.05 & abs(data.markers.alpha.G2F$avg_log2FC) > 1
                  ],   
                  gene_ids$SYMBOL),
                lab = rownames(data.markers.alpha.G2F),
                pCutoff = 5e-02,
                FCcutoff = 1,
                pointSize = 1,
                labSize = 5,
                #labCol = c('#E64B35'),
                col = c('grey30', 'forestgreen', 'royalblue', '#E64B35'),
                colGradientBreaks = 0,
                arrowheads = FALSE,
                gridlines.minor = FALSE,
                gridlines.major = FALSE,
                drawConnectors = TRUE)

EnhancedVolcano(data.markers.beta.F2H,
                x = "avg_log2FC",
                y = "p_val_adj",
                selectLab = intersect(
                  rownames(data.markers.beta.F2H)[
                    data.markers.beta.F2H$p_val_adj < 0.05 & abs(data.markers.beta.F2H$avg_log2FC) > 1
                  ],   
                  gene_ids$SYMBOL),
                lab = rownames(data.markers.beta.F2H),
                pCutoff = 5e-02,
                FCcutoff = 1,
                pointSize = 1,
                labSize = 5,
                #labCol = c('#E64B35'),
                col = c('grey30', 'forestgreen', 'royalblue', '#E64B35'),
                colGradientBreaks = 0,
                arrowheads = FALSE,
                drawConnectors = TRUE)

EnhancedVolcano(data.markers.beta.G2H,
                x = "avg_log2FC",
                y = "p_val_adj",
                selectLab = intersect(
                  rownames(data.markers.beta.G2H)[
                    data.markers.beta.G2H$p_val_adj < 0.05 & abs(data.markers.beta.G2H$avg_log2FC) > 1
                  ],   
                  gene_ids$SYMBOL),
                lab = rownames(data.markers.beta.G2H),
                pCutoff = 5e-02,
                FCcutoff = 1,
                pointSize = 1,
                labSize = 5,
                #labCol = c('#E64B35'),
                col = c('grey30', 'forestgreen', 'royalblue', '#E64B35'),
                colGradientBreaks = 0,
                gridlines.minor = FALSE,
                gridlines.major = FALSE,
                arrowheads = FALSE,
                drawConnectors = TRUE)

EnhancedVolcano(data.markers.beta.G2F,
                x = "avg_log2FC",
                y = "p_val_adj",
                selectLab = intersect(
                  rownames(data.markers.beta.G2F)[
                    data.markers.beta.G2F$p_val_adj < 0.05 & abs(data.markers.beta.G2F$avg_log2FC) > 1
                  ],   
                  gene_ids$SYMBOL),
                lab = rownames(data.markers.beta.G2F),
                pCutoff = 5e-02,
                FCcutoff = 1,
                pointSize = 1,
                labSize = 5,
                #labCol = c('#E64B35'),
                col = c('grey30', 'forestgreen', 'royalblue', '#E64B35'),
                colGradientBreaks = 0,
                gridlines.minor = FALSE,
                gridlines.major = FALSE,
                arrowheads = FALSE,
                drawConnectors = TRUE)

########################################## GSEA ########################################## 

gene_list <- setNames(data.markers.alpha.G2H0$avg_log2FC, rownames(data.markers.alpha.G2H0))
ranked_genes <- sort(gene_list, decreasing = TRUE)

gene_ids <- bitr(names(ranked_genes), fromType = "SYMBOL", 
                 toType = "ENTREZID", OrgDb = org.Hs.eg.db)
entrez_ids <- gene_ids$ENTREZID
names(entrez_ids) <- gene_ids$SYMBOL

ranked_genes_entrez <- ranked_genes[names(ranked_genes) %in% gene_ids$SYMBOL]
names(ranked_genes_entrez) <- entrez_ids[names(ranked_genes_entrez)]
ranked_genes_entrez <- sort(ranked_genes_entrez, decreasing = TRUE)

gsea_res <- gseGO(
  geneList     = ranked_genes_entrez,
  OrgDb        = org.Hs.eg.db,
  ont          = "ALL",
  keyType      = "ENTREZID",
  pvalueCutoff = 0.05,
  verbose      = FALSE
)

gsid <- "GO:0002790" 

genes_in_set <- as.character(gsea_res@geneSets[[gsid]])

top_ranked_entrez <- names(sort(gsea_res@geneList, decreasing = TRUE)[1:400])
top_entrez <- intersect(top_ranked_entrez, genes_in_set)

entrez2symbol <- setNames(gene_ids$SYMBOL, gene_ids$ENTREZID)
top_symbols <- entrez2symbol[top_entrez]

df_label <- data.frame(
  gene = top_symbols,
  x = match(top_entrez, names(gsea_res@geneList)),
  y = gsea_res@geneList[top_entrez]
)
df_label <- df_label[!is.na(df_label$gene), ]

p <- gseaplot2(gsea_res, geneSetID = gsid, 
               title = "GO: Peptide Secretion", base_size = 12)


########################################## GO AND PATHWAY ANALYSIS ########################################## 

## GO
#### Replace the Seurat object each time for entire analysis
data.markers.f <- data.markers.alpha.F2H[data.markers.alpha.F2H$p_val_adj < 0.05,]
#data.markers.f <- data.markers.alpha.G2H[data.markers.alpha.G2H$p_val_adj < 0.05,]
#data.markers.f <- data.markers.alpha.G2F[data.markers.alpha.G2F$p_val_adj < 0.05,]
#data.markers.f <- data.markers.beta.F2H[data.markers.beta.F2H$p_val_adj < 0.05,]
#data.markers.f <- data.markers.beta.G2H[data.markers.beta.G2H$p_val_adj < 0.05,]
#data.markers.f <- data.markers.beta.G2F[data.markers.beta.G2F$p_val_adj < 0.05,]

data.markers_up <- data.markers.f[data.markers.f$avg_log2FC > 1,]
data.markers_down <- data.markers.f[data.markers.f$avg_log2FC < -1,]

GO_up <- clusterProfiler::enrichGO(rownames(data.markers_up), 
                                   "org.Hs.eg.db", 
                                   keyType = "SYMBOL", 
                                   ont = "ALL") 

GO_down <- clusterProfiler::enrichGO(rownames(data.markers_down), 
                                     "org.Hs.eg.db", 
                                     minGSSize = 5,
                                     keyType = "SYMBOL", 
                                     ont = "ALL") 

Dotplot_up <- dotplot(GO_up, showCategory=15)
Emapplot_up <- enrichplot::emapplot(enrichplot::pairwise_termsim(GO_up),
                                    showCategory = 15, cex_label_category = 1)
GO_plot_up <- Dotplot_up + Emapplot_up

go_df <- as.data.frame(GO_up)

## Pathway visulization

# Extract and reshape gene-GO relationships
gene_pathways <- go_df %>%
select(ID, geneID,Description) %>%
separate_rows(geneID, sep = "/") %>%
rename(pathway_id = ID, pathway_desc = Description, gene_id = geneID)

selected_ids <- c(
  "GO:0009914",
  "GO:0043434",
  "GO:0048167",
  "GO:0016049",
  "GO:0099177",
  "GO:0031960",
  "GO:0006979",
  "GO:0062023",
  "GO:0033860",
  "GO:0016055",
  "GO:0046883",
  "GO:0048568",
  "GO:0030198",
  "GO:0099150",
  "GO:0062023",
  "GO:0098982",
  "GO:0062023",
  "GO:0070371",
  "GO:0016049",
  "GO:0099106",
  "GO:0015833",
  "GO:0005178",
  "GO:0031668",
  "GO:0051592",
  "GO:0017080",
  "GO:0046883"
)  

filtered_pathways <- gene_pathways %>%
  filter(pathway_id %in% selected_ids)

gene_pathways <- filtered_pathways

pathway_gene_list <- gene_pathways %>%
  group_by(pathway_desc) %>%
  summarise(genes = list(unique(gene_id)))

sim_matrix <- outer(pathway_gene_list$genes, pathway_gene_list$genes, Vectorize(function(a, b) {
  length(intersect(a, b)) / length(union(a, b))
}))

groups <- igraph::graph_from_adjacency_matrix(sim_matrix >= 0.5 & sim_matrix < 1, mode = "undirected") %>%
  components() %>%
  .$membership

pathway_gene_list$merged_pathway <- paste0("PathwayGroup_", groups)

gene_pathways_merged <- gene_pathways %>%
  left_join(pathway_gene_list %>% select(pathway_desc, merged_pathway), by = "pathway_desc")

edges <- gene_pathways_merged %>%
  mutate(from = merged_pathway, to = gene_id) %>%
  select(from, to)

graph <- tbl_graph(edges = edges, directed = FALSE) %>%
  mutate(
    type = ifelse(name %in% edges$from, "pathway", "gene"),
    node_size = ifelse(type == "pathway", centrality_degree(), 2)
  )

ggraph(graph, layout = "fr") +
  geom_edge_link(alpha = 0.2) +
  geom_node_point(aes(size = node_size, color = type)) +
  geom_node_text(aes(label = name), repel = TRUE, size = 3) +
  scale_size_continuous(range = c(3, 10)) +
  scale_color_manual(values = c("pathway" = "firebrick", "gene" = "grey40")) +
  theme_void()