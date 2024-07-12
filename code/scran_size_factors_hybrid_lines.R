# This creates size factors using scran normalization. This is step 1 in the normalization process.

# CLEAR ALL OBJECTS FROM WORKSPACE
rm(list=ls())


# SET WORKING DIRECTORY WHERE FILES ARE STORED
setwd("/project2/gilad/awchen55/differentialDispersion/analysis")

# LOAD PACKAGES
library(reticulate)                # for loading python
Sys.setenv(RETICULATE_PYTHON = "/project2/gilad/kenneth/miniconda3/envs/scvi_pytorch/bin/python") # set python environment
reticulate::py_config()            # load python 3.7
library(SingleCellExperiment)      # dimension reduciton and size factors for cells
library(scran)                     # normalization
library(tictoc)                    # timing
pd <- import('pandas')             # import pandas from python
np <- import('numpy')              # import numpy from python
anndata <- import('anndata')       # import anndata from python

# LOAD DATA
replicate <- "Rep1"

# load human hybrid data
human_path <- paste("r'/project2/gilad/awchen55/differentialDispersion/data/hybrid_lines_raw_data/human.ASE.", replicate, ".h5ad'", sep = "")
adata_human_ase <- anndata$read_h5ad(human_path)

# load chimp hybrid data
chimp_path <- paste("r'/project2/gilad/awchen55/differentialDispersion/data/hybrid_lines_raw_data/chimp.ASE.", replicate, ".h5ad'", sep = "")
adata_chimp_ase <- anndata$read_h5ad(chimp_path)


# CREATE SCE OBJECT
sce_human <- SingleCellExperiment(
    assays      = list(counts = t(adata_human_ase$X$toarray())) ,
     colData     = adata_human_ase$obs,
     rowData     = adata_human_ase$var
)

sce_chimp <- SingleCellExperiment(
    assays      = list(counts = t(adata_chimp_ase$X$toarray())) ,
     colData     = adata_chimp_ase$obs,
     rowData     = adata_chimp_ase$var
)

print("SCE object made")
# SCRAN NORMALIZATION
#     tic('running scran')
#     tic('clustering with quickCluster')
clusters_human <- quickCluster(sce_human)
clusters_chimp <- quickCluster(sce_chimp)

print("quickcluster")

sce_human$normalization_clusters <- clusters_human
sce_chimp$normalization_clusters <- clusters_chimp

print("norm_cluster")
#     toc()
#     tic('computing size factors')
sce_human <- computeSumFactors(sce_human, clusters=clusters_human)
sce_chimp <- computeSumFactors(sce_chimp, clusters=clusters_chimp)
print("compute sum factor")
#     toc()
#     toc()
#     summary(sizeFactors(sce))
save_output_human <- paste("/project2/gilad/awchen55/differentialDispersion/data/hybrid_lines_scran_normalized_data/human_ASE_", replicate, "_scran_size_factors.Rds", sep = "")
save_output_chimp <- paste("/project2/gilad/awchen55/differentialDispersion/data/hybrid_lines_scran_normalized_data/chimp_ASE_", replicate, "_scran_size_factors.Rds", sep = "")
saveRDS(sce_chimp, save_output)






