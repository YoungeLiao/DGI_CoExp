import numpy as np
import pandas as pd
import os

from sklearn.metrics.pairwise import euclidean_distances

import config

import matplotlib.pyplot as plt
from numpy import zeros,arange

def getNodes(config):
    all_data = pd.read_csv(config['data_path'])
    all_features = all_data.values[:, 3: ] 
    all_features[all_features == ''] = 0.0
    all_features = all_features.astype(np.float64) 
    
    all_features = np.swapaxes(all_features, 0, 1) 
    config['num_feature'] = all_features.shape[0]
    print('feature shape: ', all_features.shape)
    generated_data_path = config['generated_data_path']
    np.save( generated_data_path + 'features.npy', all_features)
    return all_features

def getEdges(config):
    all_data = pd.read_csv(config['data_path'])
    all_edges = all_data.values[:, :3] 
    np.save( config['generated_data_path'] + 'edges.npy', all_edges)
    edge_data = all_edges
    X = np.array(edge_data[:, 1]).astype(np.float64)
    Y = np.array(edge_data[:, 2]).astype(np.float64)
    num_otu = len(X)
    distance_list = []
    distance_matrix = np.zeros((num_otu, num_otu))
    for i in range(num_otu):
        for j in range(num_otu):
            if i!=j:
                dist = np.linalg.norm(np.array([X[i], Y[i]])-np.array([X[j], Y[j]]))
                distance_list.append(dist)
                distance_matrix[i,j] = dist
    distance_array = np.array(distance_matrix)

    thre = config['threshold']
    for threshold in [thre]:
        num_big = np.where(distance_array<threshold)[0].shape[0]
        print ('Threshold:',threshold, '\n', 'Links number:', num_big, '\n', 'Average Links:', str(num_big/(num_otu)))
        distance_matrix_threshold_I_list = []
        distance_matrix_threshold_W_list = []
        
        distance_matrix_threshold_I = np.zeros(distance_matrix.shape)
        distance_matrix_threshold_W = np.zeros(distance_matrix.shape)
        for i in range(distance_matrix_threshold_I.shape[0]):
            for j in range(distance_matrix_threshold_I.shape[1]):
                if distance_matrix[i,j] <= threshold and distance_matrix[i,j] > 0:
                    distance_matrix_threshold_I[i,j] = 1
                    distance_matrix_threshold_W[i,j] = distance_matrix[i,j]

        distance_matrix_threshold_I_N = np.float32(distance_matrix_threshold_I) 
        from scipy import sparse
        import pickle
        distance_matrix_threshold_I_N_crs = sparse.csr_matrix(distance_matrix_threshold_I_N)
        with open(config['generated_data_path'] + 'Adjacent', 'wb') as fp:
            pickle.dump(distance_matrix_threshold_I_N_crs, fp)

def get_GraphData(config):
    import pickle
    from scipy import sparse
    with open(config['generated_data_path'] + 'Adjacent', 'rb') as fp: 
        adj_0 = pickle.load(fp) 

    X_data = np.load(config['generated_data_path'] + 'features.npy').transpose()

    num_points = X_data.shape[0]
    adj_I = np.eye(num_points) 
    adj_I = sparse.csr_matrix(adj_I) 
    adj = (1-config['args_lambda_I'])*adj_0 + config['args_lambda_I']*adj_I 

    num_cell = X_data.shape[0]
    num_feature = X_data.shape[1]
    print('Adj:', adj.shape, 'Edges:', len(adj.data))
    print('X:', X_data.shape)
    n_clusters = config['args_n_clusters']
    return X_data, adj

def get_graph(adj, X):
    # create sparse matrix
    import torch
    from torch_geometric.data import Data, DataLoader
    row_col = []
    edge_weight = []
    rows, cols = adj.nonzero() 
    edge_nums = adj.getnnz() 
    for i in range(edge_nums):
        row_col.append([rows[i], cols[i]])
        edge_weight.append(adj.data[i]) 
    edge_index = torch.tensor(np.array(row_col), dtype=torch.long).T
    edge_attr = torch.tensor(np.array(edge_weight), dtype=torch.float)

    graph_bags = []
    graph = Data(x=torch.tensor(X, dtype=torch.float), edge_index=edge_index, edge_attr=edge_attr)  
    graph_bags.append(graph)
    return graph_bags