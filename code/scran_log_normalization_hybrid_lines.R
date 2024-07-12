# This performs log normalization. This is step 2 in the scran normalization process.


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
library(MASS)                      # save matrix file
pd <- import('pandas')             # import pandas from python
np <- import('numpy')              # import numpy from python



# LOAD SCE CLUSTERS
replicate <- "Rep1"

load_output_human <- paste("/project2/gilad/awchen55/differentialDispersion/data/hybrid_lines_scran_normalized_data/human_ASE_", replicate, "_scran_size_factors.Rds", sep = "")
sce_human = readRDS(load_output_human)

load_output_chimp <- paste("/project2/gilad/awchen55/differentialDispersion/data/hybrid_lines_scran_normalized_data/chimp_ASE_", replicate, "_scran_size_factors.Rds", sep = "")
sce_chimp = readRDS(load_output_chimp)

# LOG NORMALIZE DATA
log_norm_human <- logNormCounts(sce_human)
log_norm_data_human <- logcounts(log_norm_human)
log_norm_data_df_human <- as.data.frame(log_norm_data_human)
save_output_human <-
paste("/project2/gilad/awchen55/differentialDispersion/data/hybrid_lines_scran_normalized_data/human_ASE_", replicate, "scran_log_normalized_data.csv",sep ="")
write.csv(log_norm_data_df_human, file=save_output_human, row.names = TRUE)

log_norm_chimp <- logNormCounts(sce_chimp)
log_norm_data_chimp <- logcounts(log_norm_chimp)
log_norm_data_df_chimp <- as.data.frame(log_norm_data_chimp)
save_output_chimp <-
paste("/project2/gilad/awchen55/differentialDispersion/data/hybrid_lines_scran_normalized_data/chimp_ASE_", replicate, "scran_log_normalized_data.csv",sep ="")
write.csv(log_norm_data_df_chimp, file=save_output_chimp, row.names = TRUE)


# NORMALIZED DATA
# norm_1 <- logNormCounts(sce, log=FALSE)
# norm_data <- normcounts(norm_1)
# norm_data_df <- as.data.frame(norm_data)
# write.csv(norm_data_df, file="/project2/gilad/awchen55/ebqtl/simulated_poisson_null_analysis/sim_scDesign_normalized_data_Cardiomyocytes.csv", row.names = TRUE)

