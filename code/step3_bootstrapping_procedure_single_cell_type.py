import numpy as np
import pandas as pd
from scipy import stats
import scipy.stats.mstats as mstats
from statsmodels.nonparametric.smoothers_lowess import lowess
from statsmodels.stats.multitest import fdrcorrection
import scipy
import time
import pickle
import random

import statsmodels.stats.api as sms

######### FUNCTIONS
test_statistic = lambda x: np.std(x)/np.mean(x)

######## LOAD DATA

# paths and names
path = r'/project/gilad/brendan/dispersion/pilot/cHDC_data/cellranger_cluster-mode_trial/analysis/datasets/'
input_name_res = path + "simulation4000_metrics_residuals_1.25_cpm.csv"
input_name_metrics = path + "simulation4000_metrics_1.25_cpm.csv"
output_name = path + "simulation4000_bootstrap_results_1.25_cpm.csv"


# expression data
expression_data = pd.read_csv(path + "lane_a_card_raw_count_matrix.csv")
expression_data = np.transpose(expression_data.set_index("Unnamed: 0"))
expression_data = np.array(expression_data)

#cell_data = np.array(simulated_data)


# load residual calculations
nb_res = pd.read_csv(input_name_res).set_index("metrics")
nb_res = nb_res.drop(nb_res.columns[0], axis=1)

# load metric data
nb_metrics = pd.read_csv(input_name_metrics).set_index("metrics")
nb_metrics = nb_metrics.drop(nb_metrics.columns[0], axis=1)


gene_list = [f"Gene_{i:d}" for i in range(10000)]
#no_nan_gene_list = list( set(nb_res.columns) & set(cell_data))


# get cell data and gene list from non na values
#cell_df_t = cell_data


######## RUN CODE
cell_st = time.time()

# Resampling
samples = []
for _ in range(5000):
    x = np.random.choice(expression_data.index, size=len(expression_data), replace=True)
    cv_samp = np.array(test_statistic(expression_data.loc[x,:]))
    samples.append(cv_samp)

samp_df = pd.DataFrame(samples)
samp_df.columns = gene_list 
samp_df = samp_df.fillna(0)

# Estimates, stderr and 95% CI
est = np.mean(samp_df)
stderr = np.std(samp_df)
lwr = np.percentile(samp_df,[2.5],axis=0)
upr = np.percentile(samp_df,[97.5],axis=0)

boot_val =np.vstack((lwr,upr,est,stderr)).T
ci_df = pd.DataFrame(boot_val,columns=["lwr","upr","est","stderr"],index = gene_list)



# Get CV, UMI, Cell count, LOWESS estimates
merged = nb_metrics.T.merge(nb_res.T,left_index=True, right_index=True)[["cv","umi","mean","nonzero_cell_cnt","lowess"]]

# Log of UMI
merged["umi_log"] = np.log(merged["umi"])

# Mean UMI
#merged["mean_umi"] = merged["umi"]/merged["cell_cnt"]

# Add jackknife estimates
merged = merged.merge(ci_df, left_index=True, right_index=True)

# Z score
merged["z_score"] = (merged["lowess"] - merged["est"])/merged["stderr"]

# p-value
merged["p_value"] = scipy.stats.norm.sf(abs(merged["z_score"]))*2

# p_value_adj
merged["p_value_adj"] = fdrcorrection(merged["p_value"])[1]
cell_fn = time.time()

######## OUTPUT  
print("Cell Time:", cell_fn-cell_st)
merged.to_csv(output_name)   

