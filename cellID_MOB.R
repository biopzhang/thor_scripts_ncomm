library(CelliD)
library(tidyverse) # general purpose library for data handling
library(ggpubr) 

# read data and convert them to SeuratObject
ref_mat <- read.csv("GSE121891_OB_6_runs_processed_seurat_filter_cell.dge.csv", header = TRUE, row.names = 1)
metadata <- read.table("GSE121891_Figure_2_metadata.txt", header = TRUE, row.names = 1, sep = "\t")
ref_meta <- data.frame(FinalIds = metadata[,"FinalIds"], row.names = rownames(metadata))

target_mat <- read.csv("MOB_gene_exp.csv", header = TRUE, row.names = 1)
target_meta <- read.csv("MOB_meta.csv", header = TRUE, row.names = 1)
rownames(target_meta) <- gsub("-", ".", rownames(target_meta))

ref_sce <- CreateSeuratObject(counts = ref_mat, project = "ref_MOB", min.cells = 5, meta.data = ref_meta)
target_sce <- CreateSeuratObject(counts = target_mat, project = "thor_MOB", min.cells = 5, meta.data = target_meta)

# because the data was already normalized so we need to put in so the tools can know it already normalized
ref_sce[["RNA"]]$data <- as.matrix(ref_mat)
target_sce[["RNA"]]$data <- as.matrix(target_mat)



# Find variable features (highly variable genes)
ref_sce <- FindVariableFeatures(ref_sce)
# scale
ref_sce <- ScaleData(ref_sce, features = rownames(ref_sce))

# scale
target_sce <- ScaleData(target_sce, features = rownames(target_sce))

# dimensions
ref_sce <- RunMCA(ref_sce, nmcs = 50)
ref_sce <- RunPCA(ref_sce, features = rownames(ref_sce))
ref_sce <- RunUMAP(ref_sce, dims = 1:30)

target_sce <- RunMCA(target_sce, nmcs = 50)
target_sce <- RunPCA(target_sce, features = rownames(target_sce))
target_sce <- RunUMAP(target_sce, dims = 1:30)


# ref
ref_sce_cell_gs <- GetCellGeneSet(ref_sce, dims = 1:50)

HGT_ref_sce_cell_gs <- RunCellHGT(target_sce, pathways = ref_sce_cell_gs, dims = 1:50)



ref_sce_cell_gs_match <- rownames(HGT_ref_sce_cell_gs)[apply(HGT_ref_sce_cell_gs, 2, which.max)]
ref_sce_cell_gs_prediction <- ref_sce$FinalIds[ref_sce_cell_gs_match]
ref_sce_cell_gs_prediction_signif <- ifelse(apply(HGT_ref_sce_cell_gs, 2, max)>2, yes = ref_sce_cell_gs_prediction, "unassigned")
target_sce$ref_sce_cell_gs_prediction <- ref_sce_cell_gs_prediction_signif

write.table(ref_sce_cell_gs_prediction_signif, file='cell_type_predict_result.txt', sep=',',row.names = T,col.names = T)