from torch.utils.data import Dataset
from torch.utils.data.sampler import Sampler
import utils
import numpy as np
import torch
import random
from rdkit import Chem
from rdkit.Chem import AllChem
from scipy.spatial import distance_matrix
from rdkit.Chem.rdmolops import GetAdjacencyMatrix

random.seed(0)

def get_atom_feature(m, is_ligand=True):
    n = m.GetNumAtoms()
    H = []
    for i in range(n):
        H.append(utils.atom_feature(m, i, None, None))
    H = np.array(H)        
    if is_ligand:
        H = np.concatenate([H, np.zeros((n,28))], 1)
    else:
        H = np.concatenate([np.zeros((n,28)), H], 1)
    return H        

class MolDataset(Dataset):

    def __init__(self, keys, data_dir):
        self.keys = keys
        self.data_dir = data_dir

    def __len__(self):
        return len(self.keys)

    def __getitem__(self, idx):
        #idx = 0
        key = self.keys[idx]
        with open(self.data_dir + '/' + key, 'r') as f:
        # with open( key, 'r') as f:
            l,p = f.readline().split(',')

            # m1 = Chem.MolFromMolFile(l,sanitize=False)  # 将pdb表示分子转换为mol对象
            m1 = Chem.SDMolSupplier(str(l),sanitize=False)
            m2 = Chem.AllChem.MolFromPDBFile(p,sanitize=False)
            # print(type(m1), type(m1))
            # # # m1, m2 = pickle.load(f)
            n1 = m1[0].GetNumAtoms()
            n2 = m2.GetNumAtoms()
            c1 = m1[0].GetConformers()[0]
            d1 = np.array(c1.GetPositions())
            adj1 = GetAdjacencyMatrix(m1[0])+np.eye(n1)
            H1 = get_atom_feature(m1[0], True)

                #prepare protein
                # n2 = m2.GetNumAtoms()
            c2 = m2.GetConformers()[0]
            d2 = np.array(c2.GetPositions())
            adj2 = GetAdjacencyMatrix(m2)+np.eye(n2)
            H2 = get_atom_feature(m2, False)
                # H2 = get_atom_feature(m2, True)

                #aggregation
            H = np.concatenate([H1, H2], 0)
            agg_adj1 = np.zeros((n1+n2, n1+n2))
            agg_adj1[:n1, :n1] = adj1
            agg_adj1[n1:, n1:] = adj2
            agg_adj2 = np.copy(agg_adj1)
            dm = distance_matrix(d1,d2)
            agg_adj2[:n1,n1:] = np.copy(dm)
            agg_adj2[n1:,:n1] = np.copy(np.transpose(dm))

                #node indice for aggregation
            valid = np.zeros((n1+n2,))
            valid[:n1] = 1

                #pIC50 to class

                # Y = 1 if 'CHEMBL' in key else 0
            Y = 0 if 'inactive' in key else 1
                #if n1+n2 > 300 : return None
            sample = {
                          'H':H, \
                          'A1': agg_adj1, \
                          'A2': agg_adj2, \
                          'Y': Y, \
                          'V': valid, \
                          'key': key, \
                          }

            return sample

class DTISampler(Sampler):

    def __init__(self, weights, num_samples, replacement=True):
        weights = np.array(weights)/np.sum(weights)
        self.weights = weights
        self.num_samples = num_samples
        self.replacement = replacement
    
    def __iter__(self):
        #return iter(torch.multinomial(self.weights, self.num_samples, self.replacement).tolist())
        retval = np.random.choice(len(self.weights), self.num_samples, replace=self.replacement, p=self.weights) 
        return iter(retval.tolist())

    def __len__(self):
        return self.num_samples

def collate_fn(batch):
    max_natoms = max([len(item['H']) for item in batch if item is not None])
    
    H = np.zeros((len(batch), max_natoms, 56))
    A1 = np.zeros((len(batch), max_natoms, max_natoms))
    A2 = np.zeros((len(batch), max_natoms, max_natoms))
    Y = np.zeros((len(batch),))
    V = np.zeros((len(batch), max_natoms))
    keys = []
    
    for i in range(len(batch)):
        natom = len(batch[i]['H'])
        H[i, :natom] = batch[i]['H']
        A1[i,:natom,:natom] = batch[i]['A1']
        A2[i,:natom,:natom] = batch[i]['A2']
        Y[i] = batch[i]['Y']
        V[i,:natom] = batch[i]['V']
        keys.append(batch[i]['key'])
        # try:
        #     natom = len(batch[i]['H'])       ####################raise problem
        # except:
        #     print(batch)
        # else:
        #     H[i,:natom] = batch[i]['H']
        #     A1[i,:natom,:natom] = batch[i]['A1']
        #     A2[i,:natom,:natom] = batch[i]['A2']
        #     Y[i] = batch[i]['Y']
        #     V[i,:natom] = batch[i]['V']
        #     keys.append(batch[i]['key'])


    H = torch.from_numpy(H).float()
    A1 = torch.from_numpy(A1).float()
    A2 = torch.from_numpy(A2).float()
    Y = torch.from_numpy(Y).float()
    V = torch.from_numpy(V).float()
    
    return H, A1, A2, Y, V, keys

