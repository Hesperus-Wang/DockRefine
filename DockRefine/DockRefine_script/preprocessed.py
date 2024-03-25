#!/usr/bin/env python
import logging
import os
import subprocess
from rdkit.Chem.SaltRemover import SaltRemover
from rdkit.Chem import AllChem as Chem
scriptdir=os.path.dirname(os.path.abspath(__file__))
rootdir=os.path.dirname(scriptdir)
chosen_pdb_dir = rootdir + '/Before_docking_processed/chosen_pdb/'
chosen_mol2_dir = rootdir + '/Before_docking_processed/chosen_mol2/'
preprocessed_pdb_dir = rootdir + '/Before_docking_processed/preprocessed_pdb/'
preprocessed_mol2_dir = rootdir + '/Before_docking_processed/preprocessed_mol2/'
preprocessed_sdf_dir = rootdir + '//Before_docking_processed/preprocessed_sdf/'

#extract protein from pdbbind
def extract_pdb_file(outFolder,out_file_name):
    # outFolder = '/home/jianping/data/wanglin/test/'
    allfile = os.listdir(outFolder)
    # print(allfile)
    for file in allfile:
        dir_path = os.path.join(outFolder, file)
        file_list = os.listdir(dir_path)
        # print(file_list)
        for target in file_list:
            if target.endswith("protein.pdb"):
                ori_file = outFolder + file + '/' + target
                print(ori_file)
                # out_file_name = '/home/jianping/data/wanglin/wl_project/extract_pdb/'
                os.system('scp -r %s %s' % (ori_file, out_file_name))

#extract ligand from pdbbind
def extract_mol2_file(outFolder,out_file_name):
    # outFolder = '/home/jianping/data/wanglin/test/'
    allfile = os.listdir(outFolder)
    # print(allfile)
    for file in allfile:
        dir_path = os.path.join(outFolder, file)
        file_list = os.listdir(dir_path)
        # print(file_list)
        for target in file_list:
            if target.endswith("ligand.mol2"):
                ori_file = outFolder + file + '/' + target
                print(ori_file)
                # out_file_name = '/home/jianping/data/wanglin/wl_project/extract_pdb/'
                os.system('scp -r %s %s' % (ori_file, out_file_name))

def get_file_num(file_dir):
    M = []
    for root, dirs, file_list in os.walk(file_dir):
        for items in file_list:
            a = items.split('.')[0]
            M.append(a)
        return M

def remove_water(m):
    remover = SaltRemover(defnData="[O]")
    mol = Chem.MolFromPDBFile(m)
    return remover.StripMol(mol)

def preprocess_pdb_files(pdb_in_file, pdb_out_file):
    cmd_list = [
        'obabel',
        '-ipdb', pdb_in_file,
        '-opdb',
        '-O', pdb_out_file,
        '-p'
        '-partialcharge', 'gasteiger'
    ]
    # yapf: enable
    cmd_return = subprocess.run(cmd_list, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    output = cmd_return.stdout.decode('utf-8')
    logging.debug(output)

def preprocess_mol2_files(mol2_in_file, mol2_out_file):
    # yapf: disable
    cmd_list = [
        'obabel',
        '-imol2', mol2_in_file,
        '-omol2',
        '-O', mol2_out_file,
        '-p'
        '-partialcharge', 'gasteiger'
    ]
    # yapf: enable
    cmd_return = subprocess.run(cmd_list, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    output = cmd_return.stdout.decode('utf-8')
    # logging.debug(output)
    print(output)

def change_mol2_to_sdf(mol2_in_file, sdf_out_file):
    cmd_list = [
        'obabel',
        '-imol2', mol2_in_file,
        '-osdf',
        '-O', sdf_out_file,
        '-m'
    ]
    # yapf: enable
    cmd_return = subprocess.run(cmd_list, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    output = cmd_return.stdout.decode('utf-8')

L = []
for root, dirs, files in os.walk(chosen_pdb_dir):
    for file in files:
        # print(file)
        if file.split('_')[-1] == '.pdb':
            L.append(os.path.join(root , file))

#remove water for protein
for m in L:
    # print(m)
    mol = remove_water(m)
    Chem.MolToPDBFile(mol,m)

#add hydrogen and gasteiger charge for proteins,then save them as pdbqt files
protein_name_num = get_file_num(chosen_pdb_dir)
for j in range(protein_name_num.__len__()):
    # print(j)
    preprocess_pdb_files(chosen_pdb_dir + protein_name_num[j]+'.pdb', preprocessed_pdb_dir + protein_name_num[j]+'.pdb')

#add hydrogen and gasteiger charge for ligands,then save them as mol2 files
ligand_name_num = get_file_num(chosen_mol2_dir)
for l in range(ligand_name_num.__len__()):
    preprocess_mol2_files(chosen_mol2_dir + ligand_name_num[l] + '.mol2',preprocessed_mol2_dir + ligand_name_num[l] + '.mol2')

mol2_num = get_file_num(preprocessed_mol2_dir)
for l in range(mol2_num.__len__()):
    change_mol2_to_sdf(preprocessed_mol2_dir + mol2_num[l] + '.mol2',preprocessed_sdf_dir + mol2_num[l] + '.sdf')

