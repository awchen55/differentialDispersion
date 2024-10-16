
######### PACKAGES
import numpy as np
import pandas as pd
from statsmodels.nonparametric.smoothers_lowess import lowess


# paths and names
path = r"/project2/gilad/awchen55/differentialDispersion/data/simulations/"
input_name = path + "simulation4000_metrics_1.25_cpm.csv"
output_name = path + "simulation4000_metrics_residuals_1.25_cpm.csv"



######## LOAD DATA

metrics = pd.read_csv(input_name).set_index("metrics")
metrics = metrics.drop(metrics.columns[0], axis=1)



######## RUN CODE

# pivot data and merge CV and UMI
cv_umi_mod = metrics.loc[["cv","umi"]].T
cv_umi_mod_sub = cv_umi_mod[cv_umi_mod["umi"]>=10]
gene_list = list(cv_umi_mod_sub.index)

# initialize dataframe outputs
resid_df = []

# get array of CV and UMI separtely
cv_sub = np.array(cv_umi_mod_sub["cv"].T)
umi_sub = np.log(np.array(cv_umi_mod_sub["umi"].T))

# LOWESS SORTED
lowess_sort = lowess(cv_sub,umi_sub, return_sorted=True,frac=0.01)

# LOWESS values unsorted
lowess_unsort = lowess(cv_sub,umi_sub, return_sorted=False, frac=0.01)

# Residuals calculation
resid = cv_sub - lowess_unsort



# Covert values into dataframe
resid_df.append(resid)
resid_df.append(lowess_unsort)
resid_df.append(lowess_sort[:,0])
resid_df.append(lowess_sort[:,1])

resid_df = pd.DataFrame(resid_df)
resid_df.columns = gene_list
resid_df.insert(loc=0, column="metrics", value=list(["residuals","lowess","lowess_sort_x","lowess_sort_y"]))


######## OUTPUT
resid_df.to_csv(output_name)

print("lowess calculated")