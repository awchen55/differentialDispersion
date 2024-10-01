######### PACKAGES
import numpy as np
import pandas as pd
from statsmodels.nonparametric.smoothers_lowess import lowess

######## LOAD DATA

#replicate
replicate = 'Rep1'
cell_label_name = 'labels'
species = "human"
path = r'/project2/gilad/awchen55/differentialDispersion/data/hybrid_lines_cpm_normalized/'
input_name = path + species + '_ase_cpm_' + replicate
output_name = path + species + '_ase_cpm_' + replicate 


var_data = pd.read_csv(input_name + '_variance.csv', index_col=0)
std_data = pd.read_csv(input_name + '_std.csv', index_col=0)
mean_data = pd.read_csv(input_name + '_mean.csv', index_col=0)
cv_data = pd.read_csv(input_name + '_cv.csv', index_col=0)
cell_type_counts = pd.read_csv(input_name + '_cell_type_counts.csv', index_col=0)
umi_data = pd.read_csv(input_name + '_umi_data.csv', index_col=0)

# get CV and UMI values
cv_mod = cv_data.set_index("labels").T
umi_mod = umi_data.set_index("labels").T
cv_umi_mod = cv_mod.merge(umi_mod,left_index=True, right_index=True)

# initialize residual and lowess dataframes
resid_df = umi_mod
resid_df["Cardiomyocytes"]=0
resid_df = pd.DataFrame(resid_df["Cardiomyocytes"])

lowess_df = umi_mod
lowess_df["Cardiomyocytes"]=0
lowess_df = pd.DataFrame(lowess_df["Cardiomyocytes"])

# get the list of cell types
unique_cells = list(set(mean_data["labels"]))
lowess_dict = {}
cv_umi_sub_dict = {}
lowess_resid = []
lowess_val = []

# create the lowess and residual values
for cell in unique_cells:
    
    # get CV and UMI values in columns and filter for UMI >=10
    cv_col = str(cell) + "_x"
    umi_col = str(cell) + "_y"
    cv_umi_mod_sub = cv_umi_mod[[cv_col,umi_col]]
    cv_umi_mod_sub = cv_umi_mod_sub[cv_umi_mod_sub[umi_col]>=10]
    
    # get vectors of CV and log(UMI) for regression
    cv_sub = np.array(cv_umi_mod_sub[cv_col])
    umi_sub = np.log(np.array(cv_umi_mod_sub[umi_col]))
    
    # get index values
    idx = cv_umi_mod_sub.index
    
    # apply lowess regression
    lowess_dict[str(cell)] = lowess(cv_sub,umi_sub, return_sorted=True, frac=0.1)
    cv_umi_sub_dict[str(cell)] = cv_umi_mod_sub
    
    # unsorted lowess regression
    lowess_unsort = lowess(cv_sub,umi_sub, return_sorted=False, frac=0.1)
    lowess_val.append(lowess_unsort)
    lowess_values = lowess_unsort
    resid = cv_sub - lowess_unsort
    resid_df_sub = pd.DataFrame(resid).set_index(idx).rename(columns={0: str(cell)})
    lowess_df_sub = pd.DataFrame(lowess_values).set_index(idx).rename(columns={0: str(cell)})
    
    # integrate into final dataframe
    lowess_df = lowess_df.merge(lowess_df_sub,how="left",left_index = True, right_index = True)
    resid_df = resid_df.merge(resid_df_sub,how="left",left_index = True, right_index = True)

# rename columns
del resid_df["Cardiomyocytes_x"]
del lowess_df["Cardiomyocytes_x"]
resid_df = resid_df.rename(columns={"Cardiomyocytes_y": "Cardiomyocytes"})
lowess_df = lowess_df.rename(columns={"Cardiomyocytes_y": "Cardiomyocytes"})


######## OUTPUT
resid_df.to_csv(output_name + '_residuals.csv')
lowess_df.to_csv(output_name + '_lowess_values.csv')


