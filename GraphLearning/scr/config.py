# Here is configuration setting.

config = {
    'batch_size': 1,
    'ROOT_DATA_PATH' : './data/', # root data pathway
    'ROOT_OUTPUT_DATA_PATH' : './output/generated_data/',
    'ROOT_RESULTS_PATH': './output/results/', 
    'data_path':'./data/', # './data/graph_data_AllDEGs.csv',
    'generated_data_path': './output/generated_data/', # to do
    # 'edge_data_path': 'generated_data_path/' + 'edges.npy', # to do
    'args_DGI': 1,
    'args_lambda_I': 0.8,
    'args_n_clusters': 24, # number of clusters desired
    'args_cluster': 10, # run cluster or not 
    'args_num_epoch': 10000,
    'args_hidden': 256,
    'args_load': 0,
    'args_PCA': 1,
    'args_draw_map': 1,
    'threshold' : 1,
    'num_feature': 9,
    'args_embedding_data_path': './output/generated_data/', # './output/generated_data/',
    # 'args_data_path': './output/generated_data/', # './output/generated_data/',
    # 'args_model_path': './output/generated_data/', # './output/generated_data/',
    'args_result_path': './output/results/',
    'args_diff_gene': '1'
}

# config: {'INPUT_SHAPE': (450, 10),
#           'ENCODER_SHAPE': [512, 256],
#           'DECODER_SHAPE': [256, 512],
#           'ACTIVATION': 'relu',
#           'LAST_ACTIVATION': 'softmax',
#           'DROPOUT': 0,
#           'LATENT_SCALE': 5,
#           'OPTIMIZER': 'Adam',
#           'BATCH_SIZE': 64,
#           'EPOCHS': 300,
#           'STEPS_PER_EPOCH': None,
#           'VALIDATION_SPLIT': 0.2,
#           'VALIDATION_STEPS': 10,
#           'LATENT_OFFSET': 10,
#           'DECODER_BIAS': 'last',
#           'DECODER_REGULARIZER': 'var_l1',
#           'DECODER_REGULARIZER_INITIAL': 0.0001,
#           'DECODER_RELU_THRESH': 0,
#           'BASE_LOSS': 'mse',
#           'DECODER_BN': False,
#           'CB_MONITOR': 'val_loss',
#           'CB_LR_USE': True,
#           'CB_LR_FACTOR': 0.2,
#           'CB_LR_PATIENCE': 15,
#           'CB_LR_MIN_DELTA': 1e-8,
#           'CB_ES_USE': True,
#           'CB_ES_PATIENCE': 30,
#           'CB_ES_MIN_DELTA': 0,
#           'MULTI_GPU': False
#           }