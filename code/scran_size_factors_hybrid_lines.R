# This creates size factors using scran normalization. This is step 1 in the normalization process.

# CLEAR ALL OBJECTS FROM WORKSPACE
rm(list=ls())


# SET WORKING DIRECTORY WHERE FILES ARE STORED
setwd("/project2/gilad/awchen55/differentialDispersion/analysis")

# LOAD PACKAGES
library(reticulate)                # for loading python
reticulate::py_config()            # load python 3.7
library(SingleCellExperiment)      # dimension reduciton and size factors for cells
library(scran)                     # normalization 
library(tictoc)                    # timing
pd <- import('pandas')             # import pandas from python
np <- import('numpy')              # import numpy from python

# LOAD DATA

# # random sample of 10,000 individuals data
# obs_data <- pd$read_pickle("/project2/gilad/awchen55/ebqtl/simulated_null_analysis/adata_obs_sim_data.pkl")

# # genes with labels
# var_data_1 <- pd$read_pickle("/project2/gilad/awchen55/ebqtl/simulated_null_analysis/adata_var_sim_data.pkl")
# var_data <- head(var_data_1,10000)


# subset_test
sim_data <- pd$read_csv('/project2/gilad/awchen55/differentialDispersion/data/human_ASE_subset.csv', index_col=0).T



# CREATE SCE OBJECT
sce <- SingleCellExperiment(
    assays      = list(counts = t(sim_data)) )#,
#     colData     = obs_data,
#     rowData     = var_data
#)
print("SCE object made")
# SCRAN NORMALIZATION
#     tic('running scran')
#     tic('clustering with quickCluster')
clusters <- quickCluster(sce)
print("quickcluster")

sce$normalization_clusters <- clusters
print("norm_cluster")
#     toc()
#     tic('computing size factors')
sce <- computeSumFactors(sce, clusters=clusters)
print("compute sum factor")
#     toc()
#     toc()
#     summary(sizeFactors(sce))
save_output <- paste("/project2/gilad/awchen55/differentialDispersion/data/human_ASE_subset_scran_size_factors".Rds", sep = "")
saveRDS(sce, save_output)






