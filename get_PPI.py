import pandas as pd
# import get string ids
import get_patient


def get_protein_interactions(df, patientID):
    '''
    function that returns the StringDB API response for PPI data.
    
    Parameters
    ----------
    
    
    Returns
    -------
    '''
    import requests
    string_api_url = "https://version-11-5.string-db.org/api"
    output_format = "tsv-no-header"
    method = "network"

    ## Construct URL
    request_url = "/".join([string_api_url, output_format, method])

    # set confience_score between 0-1
    species_id = 9606

    ## Set parameters
    my_genes = [gene.split("_")[0] for gene in get_patient.get_genes_above_zero(df, patientID)]
    
    
    params = {
        "identifiers" : "%0d".join(my_genes), # your protein
        "species" : species_id, # species NCBI identifier 
        "add_nodes" : 0,
        "caller_identity" : "Research_Project", # your app name
        "show_query_node_labels" : 1
    }

    """ other params
"required_score" :	# threshold of significance to include a interaction, a number between 0 and 1000 (default depends on the network)
    "network_type" :	# network type: functional (default), physical
    "add_nodes"	: # adds a number of proteins with to the network based on their confidence score
    "show_query_node_labels" :	# when available use submitted names in the preferredName column when (0 or 1) (default:0)
    """

    ## Call STRING
    response = requests.post(request_url, data=params)
    return response


def get_PPI_df(df, patientID):
    '''
    Function pulls down PPI data from StringDB. Output is a df.
    
    Parameters
    ----------
    
    
    Returns
    -------
    '''
    response = get_protein_interactions(df, patientID)
    patient_df = pd.DataFrame(columns = ['stringId_A',
                                         'stringId_B',
                                         'preferredName_A', 
                                         'preferredName_B',
                                         'ncbiTaxonId',
                                         'score',
                                         'nscore',
                                         'fscore',
                                         'pscore',
                                         'ascore',
                                         'escore',
                                         'dscore', 
                                         'tscore'])
  
    
    for line in response.text.strip().split("\n"):
        # seperates each line into list
        l = line.strip().split("\t")
        patient_df.loc[len(patient_df)] = l
        
    patient_df = patient_df.drop_duplicates(subset=['preferredName_A',
                                                    'preferredName_B'],
                                            keep='first')
    return patient_df


def get_PPI_with_confidence(df, patientID, confidence_score = 0.6):
    '''
    Function to return PPI data to specific confidence score.
    Parameters
    ----------
    
    
    Returns
    -------
    
    '''
    response = get_protein_interactions(df, patientID)
    patient_df = pd.DataFrame(columns = ['stringId_A', 
                                         'stringId_B',
                                         'preferredName_A',
                                         'preferredName_B',
                                         'ncbiTaxonId',
                                         'score',
                                         'nscore',
                                         'fscore',
                                         'pscore',
                                         'ascore',
                                         'escore',
                                         'dscore',
                                         'tscore'])
    
    for line in response.text.strip().split("\n"):
        l = line.strip().split("\t")
        
        ## filter the interaction according to combined score
        combined_score = float(l[5])
        
        if combined_score > confidence_score:
            patient_df.loc[len(patient_df)] = l
    patient_df = patient_df.drop_duplicates(subset=['preferredName_A',
                                                    'preferredName_B'],
                                            keep='first')
    return(patient_df)


def get_physical_protein_interactions(df, patientID):
    '''
    function returns the StringDB API response for PPI data.
    Parameters
    ----------
    
    
    Returns
    -------
    '''
    import requests
    string_api_url = "https://version-11-5.string-db.org/api"
    output_format = "tsv-no-header"
    method = "network"

    ## Construct URL
    request_url = "/".join([string_api_url, output_format, method])

    # set confience_score between 0-1
    species_id = 9606

    ## Set parameters
    my_genes = [gene.split("_")[0] for gene in get_patient.get_genes_above_zero(df, patientID)]
    
    params = {
        "identifiers" : "%0d".join(my_genes), # your protein
        "species" : species_id, # species NCBI identifier 
        "network_type" : "physical",
        "show_query_node_labels" : 1,
        "caller_identity" : "Research_Project" # your app name
    }


    ## Call STRING
    response = requests.post(request_url, data=params)
    return response


def get_physical_PPI_df(df, patientID):
    '''
    Function pulls down PPI data from StringDB. Output is a df.
    Parameters
    ----------
    
    
    Returns
    -------
    '''
    response = get_physical_protein_interactions(df, patientID)
    patient_df = pd.DataFrame(columns = ['stringId_A',
                                         'stringId_B',
                                         'preferredName_A', 
                                         'preferredName_B',
                                         'ncbiTaxonId',
                                         'score',
                                         'nscore',
                                         'fscore',
                                         'pscore',
                                         'ascore',
                                         'escore',
                                         'dscore', 
                                         'tscore'])
  
    
    for line in response.text.strip().split("\n"):
        # seperates each line into list
        l = line.strip().split("\t")

        # adds each line to the data frame 
        patient_df.loc[len(patient_df)] = l
        
    patient_df = patient_df.drop_duplicates(subset=['preferredName_A',
                                                    'preferredName_B'],
                                            keep='first')
    return patient_df