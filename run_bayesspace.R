# Install the Bioconductor package manager, if necessary
if (!requireNamespace("BayesSpace", quietly = TRUE)) {
    if (!requireNamespace("BiocManager", quietly = TRUE)) {
        install.packages("BiocManager")
    }
    BiocManager::install("BayesSpace")
}

library(BayesSpace)
library(ggplot2)

source("load.R")

source("sce2h5ad.R")

# Define the main function
runBayesSpace <- function(row_data_file, col_data_file, counts_file, output_h5ad_file) {

    # print the arguments
    print(paste("row_data_file:", row_data_file))
    print(paste("col_data_file:", col_data_file))
    print(paste("counts_file:", counts_file))
    print(paste("output_h5ad_file:", output_h5ad_file))

    sce <- create_sce_from_files(row_data_file, col_data_file, counts_file)
    set.seed(100)

    # Note the counts are not yet log-normalized
    sce <- spatialPreprocess(sce, platform="ST", n.PCs=10, n.HVGs=1000, log.normalize=TRUE)
    sce <- qTune(sce, qs=seq(2, 20), platform="ST", d=10)
    # qPlot() returns a ggplot object.
    plot <- qPlot(sce)
    # save the plot
    ggplot2::ggsave("qplot.png", plot)



    sce <- spatialCluster(sce, q=10, platform="ST", d=10, init.method="mclust", model="t", gamma=3, nrep=10000, burn.in=100, save.chain=FALSE)
    
    print("Enhancing spatial features...")
    sce.enhanced <- spatialEnhance(sce, q=10, platform="ST", d=10, model="t", gamma=3, jitter_prior=0.3, jitter_scale=3.5, nrep=50000, burn.in=100, save.chain=FALSE)

    print("Get the enhanced features...")
    #markers <- c("g_1", "g_500", "g_1000")
    #sce.enhanced <- enhanceFeatures(sce.enhanced, sce, feature_names=markers, nrounds=0)
    sce.enhanced <- enhanceFeatures(sce.enhanced, sce, nrounds=0)

    # save the enhanced SCE object, 
    #saveRDS(sce.enhanced, "sce_enhanced.rds")
    
    # Convert the enhanced SCE object to an H5AD file
    convert_sce_to_h5ad(sce.enhanced, output_h5ad_file)
}

# read arguments from command line 
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 4) {
    stop("Usage: Rscript run.R <row_data_file> <col_data_file> <counts_file> <output_h5ad_file>")
}
runBayesSpace(args[1], args[2], args[3], args[4])
