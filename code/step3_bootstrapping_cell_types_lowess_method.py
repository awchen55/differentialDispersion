
######### PACKAGES
import numpy as np
import pandas as pd
from scipy import stats
import scipy.stats.mstats as mstats
from statsmodels.stats.multitest import fdrcorrection
import scipy
import time
import random
import anndata


######### FUNCTIONS
test_statistic = lambda x: np.std(x)/np.mean(x)

######## LOAD DATA


# inputs
replicate = 'Rep1'
cell_type = "Cardiomyocytes"
cell_type_name = "Cardiomyocytes"
species = "chimp"
path = r'/project2/gilad/awchen55/differentialDispersion/data/hybrid_lines_cpm_normalized/'
input_name = path + species +'_ase_cpm_' + replicate
output_name = path + species +'_ase_cpm_' + replicate + '_' + cell_type_name
cell_label_name = 'labels'

#  anndata object
adata_human_ase = anndata.read_h5ad(r'/project2/gilad/awchen55/differentialDispersion/data/hybrid_lines_raw_data/' + species + '.ASE.' + replicate +'.h5ad')
sc.pp.normalize_total(adata_human_ase, target_sum=1e6)
adata_clusters = adata_human_ase.obs.reset_index()
adata_clusters.rename(columns={'index': 'cell'}, inplace=True)

# get expression data for cell type
idx = adata_clusters.loc[adata_clusters["labels"] == str(cell_type)].index
expression_data = adata_human_ase.X.toarray()[idx,:]


# get lowess, cv, umi data
lowess = pd.read_csv(input_name + "_lowess_values.csv", index_col=0)
cv_data = pd.read_csv(input_name + '_cv.csv', index_col=0)
umi_data = pd.read_csv(input_name + '_umi_data.csv', index_col=0)
cell_type_counts = pd.read_csv(input_name + '_cell_type_counts.csv', index_col=0)
cell_counts =cell_type_counts.set_index(cell_label_name)

# residual value and get gene list
resid_cv = pd.read_csv(input_name + "_residuals.csv", index_col=0)
gene_list = list(resid_cv.index)

# get expression data for cell type and gene list from non na values
cell_df = pd.DataFrame(expression_data.T)
cell_df.index = gene_list
no_nan_gene_list = list(resid_cv[cell_type].dropna().index)
cell_df_t = cell_df.loc[no_nan_gene_list,:].T


######## RUN CODE
cell_st = time.time()

# Resampling
samples = []
for _ in range(5000):
    x = np.random.choice(cell_df_t.index, size=len(cell_df_t), replace=True)
    cv_samp = np.array(test_statistic(cell_df_t.iloc[x,:]))
    samples.append(cv_samp)

samp_df = pd.DataFrame(samples)
samp_df.columns = no_nan_gene_list 
samp_df = samp_df.fillna(0)

# Estimates, stderr and 95% CI
est = np.mean(samp_df)
stderr = np.std(samp_df)
lwr = np.percentile(samp_df,[2.5],axis=0)
upr = np.percentile(samp_df,[97.5],axis=0)

boot_val =np.vstack((lwr,upr,est,stderr)).T
ci_df = pd.DataFrame(boot_val,columns=["lwr","upr","est","stderr"],index = no_nan_gene_list)


# Create merged data set with p value calculations
cv = cv_data.set_index(cell_label_name).T

# Add cv calculation
merged = ci_df.merge(cv[cell_type],left_index=True,right_index=True, how="inner")
merged = merged.rename(columns={cell_type: 'cv'})

# Add UMI sum
merged=merged.merge(umi_data.set_index(cell_label_name).T[cell_type],left_index=True,right_index=True, how="inner")
merged = merged.rename(columns={cell_type: 'umi'})

# Add cell counts
merged=merged.merge(cell_counts.T[str(cell_type + '_nonzero')],left_index=True,right_index=True, how="inner")
merged = merged.rename(columns={str(cell_type + '_nonzero'): 'cell_count'})

# Add LOWESS estimate
merged=merged.merge(lowess[cell_type],left_index=True,right_index=True, how="inner")
merged = merged.rename(columns={cell_type: 'lowess'})

# Log of UMI
merged["umi_log"] = np.log(merged["umi"])

# Mean UMI
merged["mean_umi"] = merged["umi"]/merged["cell_count"]

# Z score
merged["z_score"] = (merged["lowess"] - merged["est"])/merged["stderr"]

# p-value
merged["p_value"] = scipy.stats.norm.sf(abs(merged["z_score"]))*2

# p_value_adj
merged["p_value_adj"] = fdrcorrection(merged["p_value"])[1]
cell_fn = time.time()

######## OUTPUT  
print(cell_type)
print("Cell Time:", cell_fn-cell_st)
merged.to_csv(output_name +'_dispersion.csv')   

