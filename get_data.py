import pandas as pd
import numpy as np
from sklearn import preprocessing
import os

# change wd
os.chdir("/Users/tomdoyle/Documents/University/Southampton/Course/BIOL6068-Research_Project/Python")

# genepy matrix file in relation to wd
genepy_matrix_file = "../Data/GenePy/GENEPY_JUNE22_CADDCUTOFF15.matrix"

# loeuf excel file in relation to wd
loeuf_excel_file = "../Data/LOEUF/gnomad2.xlsx"

# patient phenotype table
patient_phenotype_file = '../Data/ibd_phe.txt'

# import genepy matix as df. 'Sampleid' is column name for the sample ids set to be row index.
genepy_df = pd.read_table(genepy_matrix_file, 
                          index_col = 'Samid')

# read patient phenotype into df
patient_phenotype = pd.read_table(patient_phenotype_file, index_col=0)


#########
# LOEUF #
#########

# read LOEUF scores into df
LOEUF_df = pd.read_excel(loeuf_excel_file, 
                         sheet_name="Sheet1")

# create series and dict of loeuf scores upper
LOEUF_series = pd.Series(LOEUF_df["oe_lof_upper"].values, 
                         index = LOEUF_df["gene"])

LOEUF_dict = LOEUF_series.to_dict()


# normalise genepy across genes using Min max scale, every value to be between 0 and 1.
# DEPRECIATED
def norm_genepy():
    '''
    Function to normalise scores for each gene.
    '''
    return genepy_df.apply(preprocessing.minmax_scale, 0)


# binerise series (required for binarise_genepy)
def binarise_series(gene_series, percent=95):
    '''
    Function to binarise series of patients for each gene based 
    on the percentile of scores for that gene across cohort.
    
    See 'binarise_genepy' to binarise Genepy matrix.
    
    Parameters
    ----------
    
    gene_series : array_like, e.g. pd.Series
        Input array or object that can be converted to an array. Needs to be numeric.
    
    percent : array_like of float, default = 95
        Percentile or sequence of percentiles to compute, which must be between 0 and 100 inclusive. See np.percentile().
        
    Returns
    -------
    
    List of 'binarised' scores for input array that are either greater than (1) or less than/equal to (0) the percentile. 
    
    Applied to 'binarise_genepy()' function to be applied to an entire df.
    '''
    new_series = []
    # percentile calculated using numpy.percentile
    percentile = np.percentile(gene_series,
                               percent)
    for row in gene_series:
        if row > percentile:
            new_series.append(1)
        else:
            new_series.append(0)
    return new_series



def binarise_genepy(percent=95):
    '''
    Function to binarise Genepy scores.
    
    See above for function 'binarise_series'.
    
    'binarise_series()' applied to entire Genepy df with apply function.
    
    Parameters
    ----------
    percent : array_like of float, default = 95
        Percentile or sequence of percentiles to compute, which must be between 0 and 100 inclusive. See np.percentile().
    
    Returns
    -------
    
    Genepy DataFrame with matching index and column names, but with binarised scores for given percentile as described in 'binarise_series'
    
    '''
    return(genepy_df.apply(binarise_series,
                           0,
                           percent=percent))
    
# removes 'nan' values
#LOEUF_dict = {k:v if not np.isnan(v) else 1 for k,v in
#              LOEUF_dict.items() }
# DEPRECIATED: LOEUF score does not filter out 'noise'
def loeuf_score(df):
    list_LOEUF_scores = []
    # check that gene name in df is in the LOEUF dict
    for i in genepy_df.columns:
        if i in LOEUF_dict:
            # if it is, append the dictionary value
            list_LOEUF_scores.append(LOEUF_dict[i])
        else:
            # if not, append NaN
            list_LOEUF_scores.append(np.nan)
            
    # change list back into series with df.column names     
    LOEUF_score_series = pd.Series(list_LOEUF_scores, 
                                   index = df.columns)
    # Fill any NaN values so that it can be devided 
    filled_series = LOEUF_score_series.fillna(1)
    return filled_series

'''
# gets a list of loeuf scores for genepy matrix
def get_loeuf_score(df_for_loeuf):
'''
    # function to get loeuf scores for df containing gene names as column headings
'''
    loeuf_score_list = []
    for i in df_for_loeuf.columns:
        if i in LOEUF_dict.keys():
            if LOEUF_dict[i] == str('NaN'):
                loeuf_score_list.append(1)
            else:
                loeuf_score_list.append(LOEUF_dict[i])
        else:
            loeuf_score_list.append(1)
    return(loeuf_score_list)
'''
    
# weight normalised scores by LOEUF
# DEPRECIETED
def divide_df_by_loeuf(df):
    '''
    Function to divide df rows by loeuf score for each gene
    '''
    return df/loeuf_score(df)



def get_diagnosis(diagnosis):
    '''
    Function to check diagnosis of patients. 
    
    Returns a series of boolean values for each patient to describe if their diagnosis matches the input diagnosis
    
    Parameters
    ----------
    
    diagnosis : {'CD', 'UC', 'UIBD', 'UC_CD', 'Not_IBD'}
        IBD diagnosis values as assigned in patient_phenotype df. 
        
    Returns
    -------
    
    Series of patient ID's as index with boolean value, True if it matches input 'diagnosis' or False if not. 
    
    
    '''
    return patient_phenotype.loc[:,'Diagnosis'] == diagnosis
    

# get patient diagnosis subsets
def get_diagnosis_df(df, diagnosis):
    '''
    Function to return subset of df based on clincal diagnosis. Diagnosis = string object e.g. 'CD'. 
    
    Parameters
    ----------
    
    df : pd.DataFrame 
        Genepy_df to be subset into diagnosis
        
    diagnosis : {'CD', 'UC', 'UIBD', 'UC_CD', 'Not_IBD'}
        IBD diagnosis values as assigned in patient_phenotype df.
    
    Returns
    -------
    df['diagnosis' subset] : pd.Dataframe
        Returns a subset of the original df based on diagnosis and matching True index.
    
    '''
    index_value = patient_phenotype.iloc[np.where(get_diagnosis(diagnosis))[0]].index
    return df.loc[index_value]


def remove_ensembl_IDs(df):
    '''
    Function to remove ensembl IDs from df, e.g. _ESB0000338
    
    Parameters
    ----------
    df : pd.DataFrame
        Genepy_df to remove ensemble IDs in format geneID_ESB0000338 to return geneID only. 
    
    Returns
    -------
    pd.DataFrame
        Returns matching df but with ensembl IDs removed.
    
    
    '''
    return df.rename(columns = lambda colname: colname.split("_")[0])


'''
# create subset df
CD_subset = get_diagnosis_df(genepy_norm_loeuf, 'CD')
UC_subset = get_diagnosis_df(genepy_norm_loeuf, 'UC')
IBDU_subset = get_diagnosis_df(genepy_norm_loeuf, 'IBDU')
NOT_IBD_subset = get_diagnosis_df(genepy_norm_loeuf, 'NOT_IBD')
'''

# binerise genepy
genepy_bin_95 = remove_ensembl_IDs(binarise_genepy(percent=95))
genepy_bin_97_5 = remove_ensembl_IDs(binarise_genepy(percent=97.5))
genepy_bin_99 = remove_ensembl_IDs(binarise_genepy(percent=99))

# get subsets
# 95th percentile
CD_bin_95 = get_diagnosis_df(genepy_bin_95, 'CD')
UC_bin_95 = get_diagnosis_df(genepy_bin_95, 'UC')
IBDU_bin_95 = get_diagnosis_df(genepy_bin_95, 'IBDU')
UC_CD_bin_95 = get_diagnosis_df(genepy_bin_95, 'UC_CD')
NOT_IBD_95 = get_diagnosis_df(genepy_bin_95, 'NOT_IBD')

# 97.5th percentile
CD_bin_97_5 = get_diagnosis_df(genepy_bin_97_5, 'CD')
UC_bin_97_5 = get_diagnosis_df(genepy_bin_97_5, 'UC')
IBDU_bin_97_5 = get_diagnosis_df(genepy_bin_97_5, 'IBDU')
UC_CD_bin_97_5 = get_diagnosis_df(genepy_bin_97_5, 'UC_CD')
NOT_IBD_97_5 = get_diagnosis_df(genepy_bin_97_5, 'NOT_IBD')

# 99th percentile
CD_bin_99 = get_diagnosis_df(genepy_bin_99, 'CD')
UC_bin_99 = get_diagnosis_df(genepy_bin_99, 'UC')
IBDU_bin_99 = get_diagnosis_df(genepy_bin_99, 'IBDU')
UC_CD_bin_99 = get_diagnosis_df(genepy_bin_99, 'UC_CD')
NOT_IBD_99 = get_diagnosis_df(genepy_bin_99, 'NOT_IBD')


def get_top_bin_sum(matrix, head=30):
    '''
    Returns a series of the sum of patients with disease in that gene.
    
    Parameters
    ----------
    
    Returns
    -------
    
    '''
    return matrix.sum().sort_values(ascending=False).head(head)


def get_as_percentage(matrix, cohort):
    '''
    returns as a percentage of cohort, either 'CD' or 'UC'
    
    Parameters
    ----------
    
    Returns
    -------
    
    '''
    top_series = get_top_bin_sum(matrix)
    if cohort == 'CD':
        return((top_series/681)*100)
    elif cohort == 'UC':
        return((top_series/368)*100)
    else:
        return(print("incorrect cohort value, please use 'CD' or 'UC'"))
    return


def get_matching_UC(CD_matrix, UC_matrix, as_percentage=False):
    '''
    return matching top 30 UC genes
    
    Parameters
    ----------
    
    Returns
    -------
    
    '''
    UC_matching = UC_matrix[get_top_bin_sum(CD_matrix).index].sum()
    if as_percentage is False:
        return(UC_matching)
    else:
        return((UC_matching/368)*100)

    
def get_matching_as_dataframe(CD_matrix, UC_matrix):
    '''
    return df with top 30 genes
    
    Parameters
    ----------
    
    
    Returns
    -------
    
    '''
    CD_series = get_top_bin_sum(CD_matrix).rename("CD_Cohort")
    UC_series = get_matching_UC(CD_matrix,UC_matrix).rename("UC_Cohort")
    return(pd.merge(CD_series, UC_series,left_index=True,right_index=True))


def get_matching_as_df_percentage(CD_matrix, UC_matrix):
    '''
    return df with top 30 genes as percentage 
    
    Parameters
    ----------
    
    
    Returns
    -------
    
    '''
    CD_series = get_as_percentage(CD_matrix, 
                                  "CD").rename("CD_Cohort")
    UC_series = get_matching_UC(CD_matrix,
                                UC_matrix,
                                as_percentage=True).rename("UC_Cohort")
    return(pd.merge(CD_series, 
                    UC_series,
                    left_index=True,
                    right_index=True))