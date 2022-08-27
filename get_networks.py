import get_network_analysis
import networkx as nx
import matplotlib.pyplot as plt
import numpy as np 
import seaborn as sns
import get_data
from get_data import genepy_df, genepy_bin_95, genepy_bin_99, IBDU_bin_99, IBDU_bin_95, UC_CD_bin_99, UC_CD_bin_95, NOT_IBD_95, NOT_IBD_99, CD_bin_95, CD_bin_97_5, CD_bin_99, UC_bin_95, UC_bin_97_5, UC_bin_99

from get_network_analysis import get_edges_with_weight_in_graph_list
## DEPRECIATED
# Graph lists
CD_95_GL = get_network_analysis.create_graph_list_from_directory("patient_networks/CD_95_network_3/")
UC_95_GL = get_network_analysis.create_graph_list_from_directory("patient_networks/UC_95_network_3/")

CD_99_GL = get_network_analysis.create_graph_list_from_directory("patient_networks/CD_99_network_3/")
UC_99_GL = get_network_analysis.create_graph_list_from_directory("patient_networks/UC_99_network_3/")

# 
get_edges_with_weight_in_graph_list(CD_95_GL)
get_edges_with_weight_in_graph_list(UC_95_GL)

get_edges_with_weight_in_graph_list(CD_99_GL)
get_edges_with_weight_in_graph_list(UC_99_GL)


# Multigraphs
CD_95_MG = get_network_analysis.create_multiple_graph("patient_networks/CD_95_network_3/")
UC_95_MG = get_network_analysis.create_multiple_graph("patient_networks/UC_95_network_3/")

CD_99_MG = get_network_analysis.create_multiple_graph("patient_networks/CD_99_network_3/")
UC_99_MG = get_network_analysis.create_multiple_graph("patient_networks/UC_99_network_3/")


# Node degree df
from get_network_analysis import get_node_degree_df

CD_95_ND = get_node_degree_df(CD_95_GL)
UC_95_ND = get_node_degree_df(UC_95_GL)

CD_99_ND = get_node_degree_df(CD_99_GL)
UC_99_ND = get_node_degree_df(UC_99_GL)