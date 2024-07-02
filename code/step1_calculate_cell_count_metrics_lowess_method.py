
######### PACKAGES
import numpy as np
import pandas as pd
import scanpy as sc
import anndata
import scipy.sparse as ss
import scipy.special as sp
import scipy.stats as st
from scipy import stats
import scipy.stats.mstats as mstats
from statsmodels.nonparametric.smoothers_lowess import lowess
import scipy
import pickle
import random

######### FUNCTIONS
def mean_func(array):
    mean = sum(array)/len(array[:,0])
    return mean

def var_func(array):
    '''
    calculate variance. takes sparse matrix as input.
    '''
    var = mean_func(array**2) - mean_func(array)**2
    return var


def cv_func(array):
    '''
    calculate coeff of variation. takes sparse matrix as input.
    '''
    var = var_func(array)
    mean = mean_func(array)
    cv = np.sqrt(var)/mean
    return cv

def var_mean_of_expression(adata, annotation_title):
    '''
    Calculates variance, mean, coefficient of variation, and n counts of cell type. 
    Takes in AnnData object, cluster mappings, and column name of cell annotation.
    '''
    
    # get obs data
    adata_clusters = adata.obs.reset_index()
    adata_clusters.rename(columns={'index': 'cell'}, inplace=True)
    
    # get expression data as an array
    expr = adata.X.toarray()
    
    # create subsets based on cell type and store as dictionary
    unqique_cell_types = set(adata_clusters[annotation_title])
    dict_subsets = {}
    for cell in unqique_cell_types:
        if pd.isna(cell):
            df_curr = adata_clusters.loc[adata_clusters[annotation_title].isna()]
        else:
            df_curr = adata_clusters.loc[adata_clusters[annotation_title] == str(cell)]
        dict_subsets[str(cell)] = df_curr
    
    # get varaince, mean, coeff of variation, and n counts per cell types per genes
    cv_num = []
    n_num = []
    n_names = []
    var_num = []
    mean_num = []
    umi_num = []
    for cell in unqique_cell_types:
        # get subset indices that correspond to dictionary subset indices per cell type
        subset = pd.DataFrame(expr).iloc[dict_subsets[str(cell)].index,:]
        # convert back to csr matrix array
        subset_csr = ss.csr_matrix(subset).toarray()
        # calculate cv
        cv_num.append(cv_func(subset_csr))
        # n count
        #n_num.append(len(subset))
        # variance
        var_num.append(var_func(subset_csr))
        # mean
        mean_num.append(mean_func(subset_csr))
        # zero and non-zero counts
        zero_count = np.count_nonzero(subset==0,axis=0)
        nonzero_count = np.count_nonzero(subset,axis=0)
        total_count = zero_count + nonzero_count
        n_num.append(zero_count)
        n_num.append(nonzero_count)
        n_num.append(total_count)
        n_names.append(str(cell)+"_zero")
        n_names.append(str(cell)+"_nonzero")
        n_names.append(str(cell)+"_total")
        
        # umi counts
        umi_num.append(np.sum(subset,axis=0))
        
        
    
    # variance data and std data
    var_data = pd.DataFrame(var_num)
    std_data = pd.DataFrame(np.sqrt(var_num))
    var_data.columns = adata.var_names
    var_data.insert(loc=0, column=annotation_title, value=list(unqique_cell_types))
    std_data.columns = adata.var_names
    std_data.insert(loc=0, column=annotation_title, value=list(unqique_cell_types))
  
    # mean data
    mean_data = pd.DataFrame(mean_num)
    mean_data.columns = adata.var_names
    mean_data.insert(loc=0, column=annotation_title, value=list(unqique_cell_types))
    
    # coeff of variation data
    cv_data = pd.DataFrame(cv_num)
    cv_data.columns = adata.var_names
    cv_data.insert(loc=0, column=annotation_title, value=list(unqique_cell_types))
    
    # n counts per cell type per gene
    cell_type_counts = pd.DataFrame(n_num)
    cell_type_counts.columns = adata.var_names
    cell_type_counts.insert(loc=0, column=annotation_title, value=list(n_names))
    
    # UMI counts
    umi_data = pd.DataFrame(umi_num)
    umi_data.columns = adata.var_names
    umi_data.insert(loc=0, column=annotation_title, value=list(unqique_cell_types))
    
    return {'var_data': var_data, 'std_data':std_data, 'mean_data':mean_data, 'cv_data':cv_data, 'cell_type_counts': cell_type_counts, "umi_data": umi_data}
    
    

######## LOAD DATA

# cell data
adata_human_ase = anndata.read_h5ad(r'/project2/gilad/awchen55/differentialDispersion/analysis/human.ASE.h5ad')


######## RUN CODE
results = var_mean_of_expression(adata_human_ase,"labels")

######## OUTPUT
results['var_data'].to_pickle(r'/project2/gilad/awchen55/differentialDispersion/analysis/cell_count_metrics_by_cell_type/human_ase_variance.pkl')
results['std_data'].to_pickle(r'/project2/gilad/awchen55/differentialDispersion/analysis/cell_count_metrics_by_cell_type/human_ase_std.pkl')
results['mean_data'].to_pickle(r'/project2/gilad/awchen55/differentialDispersion/analysis/cell_count_metrics_by_cell_type/human_ase_mean.pkl')
results['cv_data'].to_pickle(r'/project2/gilad/awchen55/differentialDispersion/analysis/cell_count_metrics_by_cell_type/human_ase_cv.pkl')
results['cell_type_counts'].to_pickle(r'/project2/gilad/awchen55/differentialDispersion/analysis/cell_count_metrics_by_cell_type/human_ase_cell_type_counts.pkl')
results['umi_data'].to_pickle(r'/project2/gilad/awchen55/differentialDispersion/analysis/cell_count_metrics_by_cell_type/human_ase_umi_data.pkl')

