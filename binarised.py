import get_data
import pandas as pd
'''
DEPRECIATED:
MOVED TO GET_DATA.PY
'''
# binerise genepy
genepy_bin_95 = get_data.remove_ensembl_IDs(get_data.binarise_genepy(percent=95))
genepy_bin_97_5 = get_data.remove_ensembl_IDs(get_data.binarise_genepy(percent=97.5))
genepy_bin_99 = get_data.remove_ensembl_IDs(get_data.binarise_genepy(percent=99))

# get subsets
# 95th percentile
CD_bin_95 = get_data.get_diagnosis_df(genepy_bin_95, 'CD')
UC_bin_95 = get_data.get_diagnosis_df(genepy_bin_95, 'UC')
IBDU_bin_95 = get_data.get_diagnosis_df(genepy_bin_95, 'IBDU')
UC_CD_bin_95 = get_data.get_diagnosis_df(genepy_bin_95, 'UC_CD')
NOT_IBD_95 = get_data.get_diagnosis_df(genepy_bin_95, 'NOT_IBD')

# 97.5th percentile
CD_bin_97_5 = get_data.get_diagnosis_df(genepy_bin_97_5, 'CD')
UC_bin_97_5 = get_data.get_diagnosis_df(genepy_bin_97_5, 'UC')
IBDU_bin_97_5 = get_data.get_diagnosis_df(genepy_bin_97_5, 'IBDU')
UC_CD_bin_97_5 = get_data.get_diagnosis_df(genepy_bin_97_5, 'UC_CD')
NOT_IBD_97_5 = get_data.get_diagnosis_df(genepy_bin_97_5, 'NOT_IBD')

# 99th percentile
CD_bin_99 = get_data.get_diagnosis_df(genepy_bin_99, 'CD')
UC_bin_99 = get_data.get_diagnosis_df(genepy_bin_99, 'UC')
IBDU_bin_99 = get_data.get_diagnosis_df(genepy_bin_99, 'IBDU')
UC_CD_bin_99 = get_data.get_diagnosis_df(genepy_bin_99, 'UC_CD')
NOT_IBD_99 = get_data.get_diagnosis_df(genepy_bin_99, 'NOT_IBD')


def get_top_bin_sum(matrix, head=30):
    '''
    Returns a series 
    
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