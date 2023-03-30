# import 
import numpy as np
import pandas as pd
import scanorama

import sys, os
sys.path.append('./scr')

from config import config
from utils import getNodes, getEdges, get_graph
import utils
import Model
import Results_Analysis

import scipy
from scipy import sparse
import pickle
import scipy.linalg

import openpyxl

from matplotlib import pyplot as plt

import torch
from torch_geometric.data import Data, DataLoader
device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')

# Configuration: para. that frequently used to adjust model
DATA_NAME = 'graph_data_BvsD_SubCell_dgi.csv' # choose dataset # graph_data_YvsD_baseline
NEWDATA_FOLDER = 'BvsD_SubCell_dgi_10cluster/' # choose directory correpsonding to the dataset
# TODO: generate folder automatically

config['data_path'] = config['ROOT_DATA_PATH'] + DATA_NAME
config['generated_data_path'] = config['ROOT_OUTPUT_DATA_PATH'] + NEWDATA_FOLDER
config['args_result_path'] = config['ROOT_RESULTS_PATH'] + NEWDATA_FOLDER
config['args_n_clusters'] = 10

config['args_num_epoch'] = 20000 # epoch
config['threshold'] = 0.38 # threshold for edges 


# Get nodes features
Node_features = getNodes(config)


# get edges
getEdges(config) # Edge data will be generated. to do: use 'if' to choose whether generate data


# Get graph data
X_data, adj = utils.get_GraphData(config)


import torch
device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')


X_data, adj = utils.get_GraphData(config)


# Trainning 
if config['args_DGI']:
    print("-----------Deep Graph Infomax-------------")
    data_list = utils.get_graph(adj, X_data)
    data_loader = DataLoader(data_list, batch_size= config['batch_size'])
    np.random.seed(0)
    DGI_model = Model.train_DGI(config, data_loader=data_loader)
    for data in data_loader:
        data.to(device)
        X_embedding, _, _ = DGI_model[0](data)
        X_embedding = X_embedding.cpu().detach().numpy()
        X_embedding_filename =  config['generated_data_path'] + 'lambdaI' + str(config['args_lambda_I']) + '_epoch' + str(config['args_num_epoch']) + '_Embed_X.npy'
        np.save(X_embedding_filename, X_embedding)


def plot(DGI_model):
    all_loss_list = DGI_model[2]
    plt.plot(range(len(all_loss_list)), all_loss_list, 'b')  
    plt.savefig(config['args_result_path'] + 'lossplot.pdf')
    df_all_loss_list = pd.DataFrame(all_loss_list)
    df_all_loss_list.to_csv(config['args_result_path'] +'/loss.csv')
    np.savetxt(config['args_result_path'] +'/loss.txt', np.array(all_loss_list), fmt='%.4f', delimiter='\t') 

plot(DGI_model)


for data in data_loader:
    data.to(device)
    X_embedding, _, _ = DGI_model[0](data)
    X_embedding = X_embedding.cpu().detach().numpy()
    X_embedding_filename =  config['generated_data_path'] + 'lambdaI' + str(config['args_lambda_I']) + '_epoch' + str(config['args_num_epoch']) + '_Embed_X.npy'
    np.save(X_embedding_filename, X_embedding)


if config['args_cluster']:
    print("-----------Clustering-------------")
    X_embedding_filename =  config['generated_data_path'] +'lambdaI' + str(config['args_lambda_I']) + '_epoch' + str(config['args_num_epoch']) + '_Embed_X.npy'
    X_embedding = np.load(X_embedding_filename)
    if config['args_PCA']:
        X_embedding = Results_Analysis.PCA_process(X_embedding, nps=30)


# continue the function in above cell. 
if config['args_cluster']:   
    print('Shape of data to cluster:', X_embedding.shape)
    cluster_labels, score = Results_Analysis.Kmeans_cluster(X_embedding, config['args_n_clusters']) 
    config['args_n_clusters'] = cluster_labels.max()+1
    all_data = cluster_labels
    np.savetxt(config['args_result_path'] + '/types.txt', np.array(all_data), fmt='%3d', delimiter='\t')


# save clustering labels and print SCI score
np.savetxt(config['args_result_path'] + '/types.txt', np.array(all_data), fmt='%3d', delimiter='\t')
print('SCI score (Clustering quality) is:', score)