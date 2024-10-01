
######### PACKAGES
import numpy as np
import pandas as pd
import anndata
import scanpy as sc
import scipy.sparse as ss

######### FUNCTIONS
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
    res.insert(loc=0, column="metrics", value=list(["cv","var","std","mean","cell_cnt","umi","fano"]))
    
    return res
    
    

######## LOAD DATA

# paths and names
path = r"/project2/gilad/awchen55/differentialDispersion/data/simulations/"
output_name = path + "simulation_metrics_1.25.csv"

simulated_data = pd.read_csv(path + "Y_1.25_cpm.csv",sep=" ", header=None)
sim_data2 = np.array(simulated_data)

gene_list = [f"Gene_{i:d}" for i in range(10000)]

######## RUN CODE
results = var_mean_of_expression(sim_data2,gene_list)

######## OUTPUT
results.to_csv(output_name)

