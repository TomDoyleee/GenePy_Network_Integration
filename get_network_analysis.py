import os
import pandas as pd
import numpy as np
import networkx as nx
import get_patient
import get_PPI
from get_data import UC_bin_97_5, CD_bin_97_5, CD_bin_99, UC_bin_99, CD_bin_95, UC_bin_95
import seaborn as sns
import matplotlib.pyplot as plt


def create_multiple_graph(directory, score=0.7):
    '''
    Function to generate NetworkX multigraph from directory of edge lists for individual graphs. 
    '''
    G = nx.MultiGraph(subgroup = directory.split("/")[1])
    filelist = [f for f in os.listdir(directory) if not f.startswith('.')]
    if "CD" in directory:
        if "99" in directory:
            p_df = CD_bin_99
        if "97_5" in directory:
            p_df = CD_bin_97_5
        if "95" in directory:
            p_df = CD_bin_95 
    else: 
        if "99" in directory:
            p_df = UC_bin_99
        if "97_5" in directory:
            p_df = UC_bin_97_5
        if "95" in directory:
            p_df = UC_bin_95
    for file in filelist:
        patient_name = file.split(".")[0]
        G.add_nodes_from([gene.split("_")[0] for gene in get_patient.get_genes_above_zero(p_df,patient_name)], 
                         name = patient_name)
        patient_edges = nx.read_weighted_edgelist(directory+file,
                                                  delimiter = "\t",
                                                  encoding='utf-8').edges
        G.add_weighted_edges_from([(u,v,e) for u,v,e in patient_edges.data('weight') if e >= score],
                                  name = patient_name)
    return G

'''
def create_multiple_graph(directory):

    G = nx.MultiGraph()
    patient = {}
    list2 = [f for f in os.listdir(directory) if not f.startswith('.')]
    for file in list2:
        patient_edge = nx.to_edgelist(nx.read_edgelist(directory+file, 
                                                       data =[('Weight', 
                                                               float)]))
        G.add_edges_from(patient_edge, **patient)
    return G
'''

def create_graph_list_from_directory(directory):
    '''
    Function to generate a list of NetworkX graphs from directory of edge lists for each graphs.
    '''
    
    Glist = []
    filelist = [f for f in os.listdir(directory) if not f.startswith('.')]
    if "CD" in directory:
        if "99" in directory:
            p_df = CD_bin_99
        if "97_5" in directory:
            p_df = CD_bin_97_5
        if "95" in directory:
            p_df = CD_bin_95 
    else: 
        if "99" in directory:
            p_df = UC_bin_99
        if "97_5" in directory:
            p_df = UC_bin_97_5
        if "95" in directory:
            p_df = UC_bin_95
    for file in filelist:
        patient_name = file.split(".")[0]
        G = nx.Graph(nx.read_weighted_edgelist(directory+file,
                                               delimiter = "\t",
                                               encoding='utf-8'),
                     name = patient_name)
        G.add_nodes_from([gene.split("_")[0] for gene in get_patient.get_genes_above_zero(p_df,patient_name)])
        Glist.append(G)
    return Glist


## slow
def create_graph_list_from_df(df):
    Glist = []
    for i in df.index:
        try:
            ppi_df = get_PPI.get_PPI_df(df, i)[["preferredName_A",
                                                "preferredName_B",
                                                "score"]]
            G = nx.Graph(nx.from_pandas_edgelist(ppi_df,
                                                 source="preferredName_A",
                                                 target="preferredName_B",
                                                 edge_attr="score"), 
                         name=i)
            G.add_nodes_from([gene.split("_")[0] for gene in get_patient.get_genes_above_zero(df,i)])
            Glist.append(G)
        except ValueError:
            print(f"patient {str(i)} has >2000 disease genes")
            continue 
            
    return Glist


def get_patient_graph(patient_name, graph_list):
    '''
    Function to return an individual patient graph from graph list.
    
    '''
    for graph in graph_list:
        if graph.name == patient_name:
            return graph
        else:
            continue
    return


def get_total_degree(directory):
    '''
    DEPRECIATED
    '''
    multi = create_multiple_graph(directory=directory)
    dic = dict(multi.degree(list(multi.nodes)))
    return pd.Series(data=dic).sort_values(ascending=False)


def get_node_degree_df(graph_list):
    '''
    Function to generate a df of node degrees for each patient.
    '''
    df_degree = pd.DataFrame(dtype='int64')
    for graph in graph_list:
        GSeries = pd.Series(data = [val for (node, val) in graph.degree()],
                            index = [node for (node, val) in graph.degree()],
                            name = graph.name,
                            dtype = 'int64')
        df_degree = pd.concat([df_degree, GSeries], axis=1)
    df_degree = df_degree.fillna(0).T
    return df_degree


# select edges based on string score i.e. "weight" 
def get_edges_with_weight(graph, edge_weight=0.7):
    '''
    Function to return graph with edges greater than the specified edge weight.
    '''
    edge_list = []
    for u, v, weight in graph.edges.data("weight"):
        if weight < edge_weight:
            edge_list.append((u,v,weight))
        else:
            continue
    graph = graph.remove_edges_from(edge_list)
    return graph


def get_edges_with_weight_in_graph_list(graph_list, edge_weight=0.7):
    '''
    Function to apply get_edges_with_weight function to graph list using list comprehension.
    '''
    [get_edges_with_weight(graph, edge_weight) for graph in graph_list]
    
    
def get_top_genes(data, head=50):
    '''
    Get top genes from a dataframe or series, e.g. node degree df.
    '''
    return data.sum().sort_values(ascending=False).head(head).index


def create_boxplot(data, subset, name, save=False):
    '''
    Function to generate a boxplot from a dataframe. Requires a subset index to limit number of entries.
    '''
    sns.set(rc = {'figure.figsize':(40,30)})
    degree_boxplot = data[subset].boxplot(whis=(0, 100))
    degree_boxplot.tick_params(axis='x',
                               labelrotation=90)
    degree_boxplot.tick_params(labelsize=20)
    degree_boxplot.set_title('Boxplot for distrubution of node degree in ' + name,
               size=40)
    if save == False:
        pass
    else: 
        degree_boxplot.figure.savefig(name+".png")
    return plt.show(degree_boxplot)


def get_missing_patients(directory):
    '''
    Function to return missing patient networks in directory.
    '''
    if "CD" in directory:
        if "99" in directory:
            p_df = CD_bin_99
        if "97_5" in directory:
            p_df = CD_bin_97_5
        if "95" in directory:
            p_df = CD_bin_95 
    else: 
        if "99" in directory:
            p_df = UC_bin_99
        if "97_5" in directory:
            p_df = UC_bin_97_5
        if "95" in directory:
            p_df = UC_bin_95
    patient_network_list = [filename.split('.')[0] for filename in os.listdir(directory)]
    return[patient for patient in p_df.index if patient not in patient_network_list]


def get_no_of_disease_genes(list_of_patients, reference_df):
    '''
    Function to return dictionary of patients with missing networks and the number of disease genes.
    '''
    dictionary = {}
    for i in list_of_patients:
        dictionary[i] = reference_df.loc[i,:].sum()
    return dictionary


def get_no_of_disease_genes_for_missing_patient_networks(directory):
    '''
    Function to combine above 2 functions. 
    '''
    if "CD" in directory:
        if "99" in directory:
            p_df = CD_bin_99
        elif "97_5" in directory:
            p_df = CD_bin_97_5
        elif "95" in directory:
            p_df = CD_bin_95 
    elif "UC" in directory: 
        if "99" in directory:
            p_df = UC_bin_99
        elif "97_5" in directory:
            p_df = UC_bin_97_5
        elif "95" in directory:
            p_df = UC_bin_95
    return get_number_of_disease_genes(get_missing_patients(directory),
                                       p_df)


def get_total_degree_from_multigraph(multigraph):
    dic = dict(multigraph.degree(list(multigraph.nodes)))
    return pd.Series(data=dic).sort_values(ascending=False)


def get_top_node_degree_and_patients(multigraph, reference_df, head=30):
    node_degree_series = get_total_degree_from_multigraph(multigraph).head(head)
    bin_df_sum_top_nd = reference_df.sum()[node_degree_series.index]
    patient_degree_df = pd.concat([node_degree_series, bin_df_sum_top_nd], axis=1)
    patient_degree_df.columns = ["Total_Node_degree", "No_of_Patients"]
    return patient_degree_df


def count_patient_edges(multigraph, gene=None):
    edge_df = [[edge[0], edge[1], multigraph.number_of_edges(edge[0], edge[1])] for edge in multigraph.edges]
    patient_edge_df = pd.DataFrame(data=edge_df, 
                                   columns=['source', 'target', 'Patients']).drop_duplicates().sort_values(by = 'Patients', ascending=False).reset_index(drop=True)
    if gene == None:
        return patient_edge_df
    else:
        gene_edge1 = patient_edge_df[patient_edge_df['source']==gene]
        gene_edge2 = patient_edge_df[patient_edge_df['target']==gene]
        gene_edge_df = pd.concat([gene_edge1, gene_edge2])
        return gene_edge_df

    
def get_only_patient_names_from_multigraph(multigraph, gene):
    return [[dicItems[1]['name'] for dicItems in multigraph.get_edge_data(str(line[0]),str(line[1])).items()] for line in count_patient_edges(multigraph, gene).itertuples(index=False,name=None)]
    

def get_patient_names_from_multigraph(multigraph, gene = None):
    patient_list = [[dicItems[1]['name'] for dicItems in multigraph.get_edge_data(str(line[0]),str(line[1])).items()] for line in count_patient_edges(multigraph, gene).itertuples(index=False,name=None)]
    df = count_patient_edges(multigraph, gene)
    df['patients_names'] = patient_list
    return df