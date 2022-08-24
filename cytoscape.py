import os
import sys
import pandas as pd
import py4cytoscape as py4
import get_data
import get_patient
import get_PPI



# open cytoscape, check its connected
dir(py4)
py4.cytoscape_ping()
py4.cytoscape_version_info()


def create_cytoscape_network(df, patientID, method = 'network'):
    '''
    Function to pull down PPI data from string and import into cyctoscape.
    DEPRECIATED: Use network file/NetworkX instead.
    '''
    if method == 'network':
        network = get_PPI.get_PPI_df(df, patientID)
    elif method == 'interactions':
        network = get_PPI.get_all_PPI_df(df, patientID)
    else:
        print("Incorrect method selected")
    
    edge_data = {'source': network["preferredName_A"],
                 'target': network["preferredName_B"],
                 'String_Score': network["score"]
                }
    edges = pd.DataFrame(data=edge_data, columns=['source', 
                                                  'target', 
                                                  'String_Score'])
    return py4.create_network_from_data_frames(edges=edges,
                                               title='Patient_Network '+str(patientID),
                                               collection="Patient_Collection")


def create_network_from_tab_file(patient_name, directory):
    '''
    Function to return cytoscape network from filename and directory.
    '''
    # iterate over files in that directory
    file_list = os.listdir(directory)
    filename = patient_name+".txt"
    py4.sandbox_send_to(directory+filename)
    py4.import_network_from_tabular_file(file=filename,
                                         first_row_as_column_names=False,
                                         start_load_row=1, 
                                         column_type_list='s,t,i', 
                                         delimiters='\t')
    

def create_network_in_cytoscape(i):
    '''
    Function to create network 
    '''
    patient_net = get_PPI.get_PPI_df(CD_bin_97_5, i)
    edge_data = {'source': patient_net["preferredName_A"],
                 'target': patient_net["preferredName_B"],
                 'String_Score': patient_net["score"]
                }
    edges = pd.DataFrame(data=edge_data, columns=['source', 
                                                  'target', 
                                                  'String_Score'])
    py4.create_network_from_data_frames(edges=edges, 
                                        title='Patient_Network'+str(i),
                                        collection="Patient_Collection"+str(i))
    return(print(f"Patient {i} completed"))


def view_networkx_in_cytoscape(graph):
    '''
    Function to visualise NetworkX graph in Cytoscape.
    Requires Cytoscape application to be open.
    
    Parameters
    ----------
    graph : NetworkX graph
        a NetworkX graph object
    
    Returns
    -------
    Visualisation in Cytocape of NetworkX graph
    '''
    return py4.create_network_from_networkx(graph)

# TODO: add node 
# TODO: list comprehension 
def create_multigraph_network_in_cytoscape(multigraph):
    list4 = []
    for u, v, data in multigraph.edges.data():
        list4.append([u, v, data['name'], data['weight']])
    
    df = pd.DataFrame(list4)
    df['4'] = 'CD'
    df['5'] = 'pp'
    df.columns = ['source', 
                  'target',
                  'patientID', 
                  'weight',
                  'diagnosis',
                  'interaction']
    return py4.create_network_from_data_frames(nodes= pd.unique(df[['source', 'target']].values.ravel('K')), edges=df)
    