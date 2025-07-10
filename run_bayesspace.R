# =============================================================================
# BayesSpace Spatial Transcriptomics Analysis
# Spatial clustering and enhancement of spatial transcriptomics data
# =============================================================================

# Install required packages if not already installed
if (!requireNamespace("BayesSpace", quietly = TRUE)) {
    if (!requireNamespace("BiocManager", quietly = TRUE)) {
        install.packages("BiocManager")
    }
    BiocManager::install("BayesSpace")
}

# Load required libraries
library(BayesSpace)
library(ggplot2)

# Load custom functions
source("load.R")
source("sce2h5ad.R")

# =============================================================================
# Main Function: runBayesSpace
# =============================================================================

runBayesSpace <- function(row_data_file, col_data_file, counts_file, output_h5ad_file) {
    
    # Print input arguments for debugging
    cat("Input parameters:\n")
    cat("  row_data_file:", row_data_file, "\n")
    cat("  col_data_file:", col_data_file, "\n")
    cat("  counts_file:", counts_file, "\n")
    cat("  output_h5ad_file:", output_h5ad_file, "\n")
    
    # =============================================================================
    # Data Loading and Preprocessing
    # =============================================================================
    
    # Create SingleCellExperiment object from input files
    sce <- create_sce_from_files(row_data_file, col_data_file, counts_file)
    
    # Set random seed for reproducibility
    set.seed(100)
    
    # Spatial preprocessing
    # Note: counts are not yet log-normalized
    cat("Running spatial preprocessing...\n")
    sce <- spatialPreprocess(sce, 
                            platform = "ST", 
                            n.PCs = 10, 
                            n.HVGs = 1000, 
                            log.normalize = TRUE)
    
    # =============================================================================
    # Parameter Tuning
    # =============================================================================
    
    # Tune the number of clusters (q)
    cat("Tuning number of clusters...\n")
    sce <- qTune(sce, 
                 qs = seq(2, 20), 
                 platform = "ST", 
                 d = 10)
    
    # Generate and save q-plot
    cat("Generating q-plot...\n")
    plot <- qPlot(sce)
    ggplot2::ggsave("qplot.png", plot)
    
    # =============================================================================
    # Spatial Clustering
    # =============================================================================
    
    cat("Running spatial clustering...\n")
    sce <- spatialCluster(sce, 
                         q = 10, 
                         platform = "ST", 
                         d = 10, 
                         init.method = "mclust", 
                         model = "t", 
                         gamma = 3, 
                         nrep = 10000, 
                         burn.in = 100, 
                         save.chain = FALSE)
    
    # =============================================================================
    # Spatial Enhancement
    # =============================================================================
    
    cat("Enhancing spatial features...\n")
    sce.enhanced <- spatialEnhance(sce, 
                                  q = 10, 
                                  platform = "ST", 
                                  d = 10, 
                                  model = "t", 
                                  gamma = 3, 
                                  jitter_prior = 0.3, 
                                  jitter_scale = 3.5, 
                                  nrep = 50000, 
                                  burn.in = 100, 
                                  save.chain = FALSE)
    
    # =============================================================================
    # Feature Enhancement
    # =============================================================================
    
    cat("Getting enhanced features...\n")
    # Uncomment the following lines to enhance specific markers:
    # markers <- c("g_1", "g_500", "g_1000")
    # sce.enhanced <- enhanceFeatures(sce.enhanced, sce, feature_names = markers, nrounds = 0)
    
    # Enhance all features
    sce.enhanced <- enhanceFeatures(sce.enhanced, sce, nrounds = 0)
    
    # =============================================================================
    # Save Results
    # =============================================================================
    
    # Uncomment to save the enhanced SCE object as RDS:
    # saveRDS(sce.enhanced, "sce_enhanced.rds")
    
    # Convert and save as H5AD file
    cat("Converting to H5AD format...\n")
    convert_sce_to_h5ad(sce.enhanced, output_h5ad_file)
    
    cat("Analysis completed successfully!\n")
}

# =============================================================================
# Command Line Interface
# =============================================================================

# Read arguments from command line
args <- commandArgs(trailingOnly = TRUE)

# Check if correct number of arguments provided
if (length(args) != 4) {
    stop("Usage: Rscript run_bayesspace.R <row_data_file> <col_data_file> <counts_file> <output_h5ad_file>")
}

# Run the main function with provided arguments
runBayesSpace(args[1], args[2], args[3], args[4])
