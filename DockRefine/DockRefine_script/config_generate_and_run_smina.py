#!/usr/bin/env python
import logging
import os
import subprocess
scriptdir=os.path.dirname(os.path.abspath(__file__))
rootdir=os.path.dirname(scriptdir)
pdb_file_dir = rootdir + '/Before_docking_processed/preprocessed_pdb/'
mol2_file_dir = rootdir + '/Before_docking_processed/preprocessed_mol2/'
config_dir = rootdir + '/Before_docking_processed/config_output/'
docking_result_dir = rootdir + '/2.docking_result/'

list_dir = os.listdir(pdb_file_dir)
pdb_list = []
for item in list_dir:
    num = item.split('.')[0]
    pdb_list.append(num)
print(pdb_list)#find ligands according to this list -> pdb_list
#
for target in pdb_list:
    receptor_path = os.path.join(pdb_file_dir,target+'.pdb')
    # print(receptor_path)
    info= 'receptor='+receptor_path
    f = open(config_dir + target + '_config.txt', 'w')
    f.write(info+'\n')
    f.close()

for lig in os.listdir(mol2_file_dir):
    ligand_path = os.path.join(mol2_file_dir, lig)
    for i in range(len(pdb_list)):
        if pdb_list[i] == lig.split('_')[0] + lig.split('_')[2]:
            # info = 'ligand='+ligand_path+'\n'
            info = 'ligand=' + ligand_path
            target = info.split('/')[-1].split('_')[0]+info.split('/')[-1].split('_')[2]
            f = open(config_dir + target + '_config.txt', 'a')
            f.write(info+'\n')
            f.close()
            break


for i in os.listdir(config_dir):
    out_info = 'out='+docking_result_dir+i.split('_')[0]+'_out.sdf'
    f = open(config_dir + i, 'a')
    f.write(out_info+'\n')
    f.close()


#################################################
def run_smina(config,ligand):
    cmd_list = [
        'smina',
        '--config', config,
        '--autobox_ligand',ligand,
        '--autobox_add 4 '
        '--exhaustiveness=8 '
        '--num_modes=20 '
    ]
    # yapf: enable
    # cmd_return = subprocess.run(cmd_list, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    # output = cmd_return.stdout.decode('utf-8')
    # logging.debug(output)
    line = ' '.join(cmd_list)
    print(line)
    os.system(' '.join(cmd_list))

def split_SDF(file_name):
    file_str_list = []
    index = 0
    with open(file_name,'r+') as f:
        B = file_name.split('/')[-1]
        print(B)
        for line in f:
            if line != '$$$$\n':
                file_str_list.append(line)
                print(file_str_list)
            else:
                index += 1
                with open(rootdir+'/2.1_docking_result_split/'+B.split('.')[0]+'{0}.sdf'.format(index),'w+') as wt:
                    for ds in file_str_list:
                        wt.write(ds)
                file_str_list = []

for j in os.listdir(mol2_file_dir):
    ligand_path = os.path.join(mol2_file_dir, j)
    for e in range(len(pdb_list)):
        if pdb_list[e] == j.split('_')[0] + j.split('_')[2]:
            run_smina(config_dir + j.split('_')[0] + j.split('_')[2]+'_config.txt', mol2_file_dir + j)
# run_smina('/home/jianping/Downloads/config_inactive.txt','/home/jianping/Desktop/ADRB2/6mxt_ligand.mol2')
# run_smina('/home/jianping/Downloads/config_active.txt','/home/jianping/Desktop/ADRB2/6mxt_ligand.mol2')
A = []
for file in os.listdir(docking_result_dir):
    dir_path = os.path.join(docking_result_dir, file)
    A.append(dir_path)
    print(A)
for item in A:
    split_SDF(item)