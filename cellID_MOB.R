# =============================================================================
# Cell Type Prediction using CelliD
# Reference: MOB (Main Olfactory Bulb) data analysis
# =============================================================================

# Load required libraries
library(CelliD)
library(tidyverse)  # General purpose library for data handling
library(ggpubr)

# =============================================================================
# Data Loading and Preprocessing
# =============================================================================

# Read reference data (Gene expression matrix and metadata from GEO)
ref_mat <- read.csv("GSE121891_OB_6_runs_processed_seurat_filter_cell.dge.csv", 
                    header = TRUE, 
                    row.names = 1)
metadata <- read.table("GSE121891_Figure_2_metadata.txt", 
                       header = TRUE, 
                       row.names = 1, 
                       sep = "\t")
ref_meta <- data.frame(FinalIds = metadata[, "FinalIds"], 
                       row.names = rownames(metadata))

# Read target data (Gene expression matrix and metadata from Thor)
target_mat <- read.csv("MOB_gene_exp.csv", 
                       header = TRUE, 
                       row.names = 1)
target_meta <- read.csv("MOB_meta.csv", 
                        header = TRUE, 
                        row.names = 1)
rownames(target_meta) <- gsub("-", ".", rownames(target_meta))

# =============================================================================
# Create Seurat Objects
# =============================================================================

# Create Seurat objects for reference and target data
ref_sce <- CreateSeuratObject(counts = ref_mat, 
                             project = "ref_MOB", 
                             min.cells = 5, 
                             meta.data = ref_meta)
target_sce <- CreateSeuratObject(counts = target_mat, 
                                project = "thor_MOB", 
                                min.cells = 5, 
                                meta.data = target_meta)

# Set normalized data (since data was already normalized)
ref_sce[["RNA"]]$data <- as.matrix(ref_mat)
target_sce[["RNA"]]$data <- as.matrix(target_mat)

# =============================================================================
# Feature Selection and Scaling
# =============================================================================

# Find variable features (highly variable genes) for reference data
ref_sce <- FindVariableFeatures(ref_sce)

# Scale data for both reference and target
ref_sce <- ScaleData(ref_sce, features = rownames(ref_sce))
target_sce <- ScaleData(target_sce, features = rownames(target_sce))

# =============================================================================
# Dimensionality Reduction
# =============================================================================

# Reference data dimensionality reduction
ref_sce <- RunMCA(ref_sce, nmcs = 50)
ref_sce <- RunPCA(ref_sce, features = rownames(ref_sce))
ref_sce <- RunUMAP(ref_sce, dims = 1:30)

# Target data dimensionality reduction
target_sce <- RunMCA(target_sce, nmcs = 50)
target_sce <- RunPCA(target_sce, features = rownames(target_sce))
target_sce <- RunUMAP(target_sce, dims = 1:30)

# =============================================================================
# Cell Type Prediction using CelliD
# =============================================================================

# Get cell gene sets from reference data
ref_sce_cell_gs <- GetCellGeneSet(ref_sce, dims = 1:50)

# Run cell HGT (Hypergeometric Test) for cell type prediction
HGT_ref_sce_cell_gs <- RunCellHGT(target_sce, 
                                  pathways = ref_sce_cell_gs, 
                                  dims = 1:50)

# =============================================================================
# Prediction Results Processing
# =============================================================================

# Get predictions based on maximum HGT scores
ref_sce_cell_gs_match <- rownames(HGT_ref_sce_cell_gs)[apply(HGT_ref_sce_cell_gs, 2, which.max)]
ref_sce_cell_gs_prediction <- ref_sce$FinalIds[ref_sce_cell_gs_match]

# Filter predictions based on significance threshold (>2)
ref_sce_cell_gs_prediction_signif <- ifelse(apply(HGT_ref_sce_cell_gs, 2, max) > 2, 
                                           yes = ref_sce_cell_gs_prediction, 
                                           "unassigned")

# Add predictions to target Seurat object
target_sce$ref_sce_cell_gs_prediction <- ref_sce_cell_gs_prediction_signif

# =============================================================================
# Save Results
# =============================================================================

# Write prediction results to file
write.table(ref_sce_cell_gs_prediction_signif, 
            file = 'cell_type_predict_result.txt', 
            sep = ',', 
            row.names = TRUE, 
            col.names = TRUE)