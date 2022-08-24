import os
import pandas as pd
import get_PPI
import concurrent.futures
from multiprocessing import Pool
from binerised import UC_bin_95



def create_patient_networks(df):
    for i in range(len(df)):
        patient_net = get_PPI.get_PPI_df(df, i)
        edge_data = {'source': patient_net["preferredName_A"],
                     'target': patient_net["preferredName_B"],
                     'String_Score': patient_net["score"]
                    }
        edges = pd.DataFrame(data=edge_data, columns=['source', 'target', 'String_Score'])
        edges.to_csv("patient_networks/patient_"+str(i)+".txt",
                     sep="\t",
                     header=False,
                     index=False)
        return


def fast(i, df):
    patient_net = get_PPI.get_PPI_df(df, i)
    edge_data = {'source': patient_net["preferredName_A"],
                 'target': patient_net["preferredName_B"],
                 'String_Score': patient_net["score"]
                }
    edges = pd.DataFrame(data=edge_data, columns=['source', 
                                                  'target',
                                                  'String_Score'])
    edges.to_csv("patient_networks/UC_95_network_3/"+str(df.index[i])+".txt",
                 sep="\t",
                 header=False,
                 index=False,
                 encoding="utf-8")
    return(print(f"Patient {df.index[i]} completed"))

def fast2(i, df):
    patient_net = get_PPI.get_physical_PPI_df(df, i)
    edge_data = {'source': patient_net["preferredName_A"],
                 'target': patient_net["preferredName_B"],
                 'String_Score': patient_net["score"]
                }
    edges = pd.DataFrame(data=edge_data, columns=['source', 
                                                  'target',
                                                  'String_Score'])
    edges.to_csv("patient_networks/UC_95_network_physical/"+str(df.index[i])+".txt",
                 sep="\t",
                 header=False,
                 index=False,
                 encoding="utf-8")
    return(print(f"Patient {df.index[i]} completed"))
    
def main(data):
    with concurrent.futures.ProcessPoolExecutor() as executor:
        patients = range(len(data))
        results = [executor.submit(fast2, patient, data) for patient in patients]

if __name__ == '__main__':
    main(UC_bin_95)
    

    '''
    pool = Pool(os.cpu_count())
    pool.map(fast, range(len(CD_bin_97_5))
    '''