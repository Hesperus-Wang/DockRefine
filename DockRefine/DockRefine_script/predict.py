#!/usr/bin/env python
import pickle
from gnn import gnn
import numpy as np
import utils
import torch.nn as nn
import torch
import time
import os
from sklearn.metrics import roc_auc_score
import argparse
import time
from torch.utils.data import DataLoader
from dataset import MolDataset, collate_fn, DTISampler

scriptdir=os.path.dirname(os.path.abspath(__file__))
rootdir=os.path.dirname(scriptdir)

parser = argparse.ArgumentParser()
parser.add_argument("--lr", help="learning rate", type=float, default = 0.0001)
parser.add_argument("--epoch", help="epoch", type=int, default = 30)
parser.add_argument("--ngpu", help="number of gpu", type=int, default = 1)
parser.add_argument("--batch_size", help="batch_size", type=int, default = 32)
parser.add_argument("--num_workers", help="number of workers", type=int, default = 7)
parser.add_argument("--n_graph_layer", help="number of GNN layer", type=int, default = 4)
parser.add_argument("--d_graph_layer", help="dimension of GNN layer", type=int, default = 140)
parser.add_argument("--n_FC_layer", help="number of FC layer", type=int, default = 4)
parser.add_argument("--d_FC_layer", help="dimension of FC layer", type=int, default = 128)
parser.add_argument("--total_key", help="file path of dude data", type=str, default='data/')
parser.add_argument("--save_dir", help="save directory of model parameter", type=str, default = './save/')
parser.add_argument("--initial_mu", help="initial value of mu", type=float, default = 4.0)
parser.add_argument("--initial_dev", help="initial value of dev", type=float, default = 1.0)
parser.add_argument("--dropout_rate", help="dropout_rate", type=float, default = 0.0)
parser.add_argument("--train_keys", help="train keys", type=str, default='keys/train_keys.pkl')
parser.add_argument("--test_keys", help="test keys", type=str, default='keys/test_keys.pkl')
parser.add_argument("--predict_keys", help="predict keys", type=str, default='keys/predict_keys.pkl')
args = parser.parse_args()

model = gnn(args)
model.load_state_dict(torch.load(scriptdir + '/save/save_428.pt', map_location=torch.device('cpu'), ),False)
model.eval()

with open (args.predict_keys, 'rb') as fp:
    predict_keys = pickle.load(fp)
predict_dataset = MolDataset(predict_keys, args.total_key)
predict_dataloader = DataLoader(predict_dataset, 1, \
     shuffle=False, num_workers = args.num_workers, collate_fn=collate_fn)
#################
# with open (args.test_keys, 'rb') as fp:
#     test_keys = pickle.load(fp)
# test_dataset = MolDataset(test_keys, args.dude_data_fpath)
# test_dataloader = DataLoader(test_dataset, 1, \
#      shuffle=False, num_workers = args.num_workers, collate_fn=collate_fn)
##################
device = torch.device( "cpu")
# for i_batch, sample in enumerate(predict_dataloader):
for i_batch, sample in enumerate(predict_dataloader):
    # for i_batch, sample in enumerate(trn_loader):
    model.zero_grad()
    H, A1, A2, Y, V, keys = sample
    H, A1, A2, Y, V = H.to(device), A1.to(device), A2.to(device), \
                      Y.to(device), V.to(device)

    # train neural network
    pred = model.train_model((H, A1, A2, V))
    pred1 = pred.detach().numpy()
    if pred1 >= 0.5:
        pred1 =1
    else:
        pred1 =0
    # print(pred,keys)
    A =[]
    A.append(keys)
    A.append(pred1)
    print(A)

    with open(rootdir+'/GAT_result.txt','a+')as f:
        f.write(str(A)+'\n')



