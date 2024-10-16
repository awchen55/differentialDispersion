
######### PACKAGES
import numpy as np
import pandas as pd
import anndata
import scanpy as sc
import scipy.sparse as ss

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

def fano_func(array):
    '''
    calculate fano factor. takes sparse matrix as input.
    '''
    var = var_func(array)
    mean = mean_func(array)
    fano = var/mean
    return fano

def var_mean_of_expression(adata,gene_list):
    '''
    Calculates variance, mean, coefficient of variation, and n counts of cell type. 
    Takes in AnnData object, cluster mappings, and column name of cell annotation.
    '''
 
    # get varaince, mean, coeff of variation, and n counts per cell types per genes
    metrics = []

    # calculate cv
    metrics.append(cv_func(adata))

    # variance
    metrics.append(var_func(adata))

    # std
    metrics.append(np.sqrt(var_func(adata)))

    # mean
    metrics.append(mean_func(adata))
    
    # non-zero counts
    nonzero_count = np.count_nonzero(adata,axis=0)
    metrics.append(nonzero_count)

    # umi counts
    metrics.append(np.sum(adata,axis=0))

    # calculate fano
    metrics.append(fano_func(adata))

    # results
    res = pd.DataFrame(metrics)
    res.columns = gene_list
    res.insert(loc=0, column="metrics", value=list(["cv","var","std","mean","nonzero_cell_cnt","umi","fano"]))
    
    return res
    
    

######## LOAD DATA

# paths and names
path = r'/project/gilad/brendan/dispersion/pilot/cHDC_data/cellranger_cluster-mode_trial/analysis/datasets/'
output_name = path + "OUTPUT_NAME.csv"

expression_data = pd.read_csv(path + "lane_a_card_raw_count_matrix.csv")
expression_data = np.transpose(expression_data.set_index("Unnamed: 0"))
expression_data = np.array(expression_data)

gene_list = list(file.columns)

######## RUN CODE
results = var_mean_of_expression(expression_data,gene_list)

######## OUTPUT
results.to_csv(output_name)

