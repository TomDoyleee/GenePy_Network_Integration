import os
import pandas as pd
import numpy as np
import get_data 
import get_patient
import get_PPI
from get_data import genepy_df
import networkx as nx
import get_network_analysis
from get_data import UC_bin_97_5, CD_bin_97_5, CD_bin_99, UC_bin_99, CD_bin_95, UC_bin_95
import concurrent.futures
from multiprocessing import Pool

binarised_dfs = [CD_bin_95,
                 UC_bin_95,
                 CD_bin_97_5,
                 UC_bin_97_5,
                 CD_bin_99,
                 UC_bin_99
                ]

folder_names = ['CD_bin_95_network',
                'UC_bin_95_network',
                'CD_bin_97_5_network',
                'UC_bin_97_5_network',
                'CD_bin_99_network',
                'UC_bin_99_network'
                ]

def get_patient_networks(df, i, folder_name):
    patient_net = get_PPI.get_PPI_df(df, i)
    edge_data ={'source': patient_net["preferredName_A"],
                'target': patient_net["preferredName_B"],
                'String_Score': patient_net["score"]
               }
    edges = pd.DataFrame(data=edge_data,
                         columns =['source',
                                   'target',
                                   'String_Score'])
    edges.to_csv(f"patient_networks_full_test/{folder_name}/{df.index[i]}.txt",
                 sep="\t",
                 header=False,
                 index=False,
                 encoding="utf-8")
    return(print(f"Patient {df.index[i]} completed"))


def get_patient_networks_patient_name(df, i, folder_name):
    patient_net = get_PPI.get_PPI_df(df, i)
    edge_data ={'source': patient_net["preferredName_A"],
                'target': patient_net["preferredName_B"],
                'String_Score': patient_net["score"]
               }
    edges = pd.DataFrame(data=edge_data,
                         columns =['source',
                                   'target',
                                   'String_Score'])
    edges.to_csv(f"patient_networks_full_test/{folder_name}/{i}.txt",
                 sep="\t",
                 header=False,
                 index=False,
                 encoding="utf-8")
    return(print(f"Patient {i} completed"))


def main(binarised_df, folder_name):
    with concurrent.futures.ProcessPoolExecutor() as executor:
        patients = range(len(binarised_df))
        executor.map(get_patient_networks, binarised_df, patients, folder_name)
        print(f'Waiting for {folder_name} to build networks..')
    print(f'{folder_name} complete')



def main(binarised_df, folder_name):
    with concurrent.futures.ProcessPoolExecutor() as executor:
        patients = range(len(binarised_df))
        executor.map(get_patient_networks, binarised_df, patients, folder_name)
        print(f'Waiting for {folder_name} to build networks..')
    print(f'{folder_name} complete')


def make_directories():
    os.mkdir("patient_networks_full_test/")
    for folder in folder_names:
        try:
            os.mkdir("patient_networks_full_test/"+folder)
        except OSError as error:
            print(error)
            continue
            
def make_network_files():
    folder_index = 0
    for df in binarised_dfs:
        main(df, "patient_networks_full_test/"+folder_names[folder_index]+"/")
        
        # wait here for the result to be available before continuing
        folder_index += 1

def check_missing_patients():
    for folder in folder_names:
        df_index = folder_names.index(folder)
        missing = get_network_analysis.get_missing_patients("patient_networks_full_test/"+folder)
        print("Missing patient networks:")
        print(missing)
        for patient in missing:
            get_patient_networks_patient_name(binarised_dfs[df_index], patient, folder)
    
if __name__ == '__main__':
    make_directories()
    make_network_files()
    check_missing_patients()
   