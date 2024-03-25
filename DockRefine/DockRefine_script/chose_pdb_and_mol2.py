#!/usr/bin/env python
import os
import pandas as pd
import shutil
import time
#############################################
scriptdir=os.path.dirname(os.path.abspath(__file__))
rootdir=os.path.dirname(scriptdir)
file = scriptdir + '/input_pocket/pocket_align.POC'


N = []
H = []
L = []
count = 0
for count,line in enumerate(open(file,'r',encoding = 'utf-8').readlines()):
    count += 1
print(count)   #count how many templete pocket was chosen
############################################
# get protein name which corresponding to templete pocket,then written to a txt files
if os.path.isfile(file):
    with open (file,'r')as f:
        for line in f.readlines():
            data = line.split(',')[0]
            # print(data)
            lig = data.split('.')[0] + '.mol2'
            print(lig)
            L.append(lig)
            N.append(data)
        # print(N)
for item in N:
    code = item.split('_')[0]+item.split('_')[2]
    H.append(code +'.pdb')
    H.append(code)
print(H)
#############################################
p = open(scriptdir+'/input_pocket/pocket_search_result_list.txt','w')
p.write('name'+'\n')
p.write('\n'.join(H))
p.close()

l = open(scriptdir+'/input_pocket/ligand_search_result_list.txt','w')
l.write('name'+'\n')
l.write('\n'.join(L))
l.close()  #the next step need a column name,I add it manually,because every time write-in will cover the column name,so I annotation 35-40
###############################################
#list the lost protein
df = pd.read_csv(scriptdir+'/input_pocket/pocket_search_result_list.txt')
rf = pd.read_csv(scriptdir+'/input_pocket/ligand_search_result_list.txt')
print(len(df))
receptor_file_path = '/media/jianping/839282ac-3ae6-4224-98cb-8c5cf2c2cf3e/wanglin/receptor/'
ligand_file_path = rootdir + '/BioLip/MOL2/'
big = os.listdir(receptor_file_path)
small = df['name'].tolist()
gap = set(small) - set(big)
print(len(gap),gap)
###################################################
receptor_save_path = rootdir + '/Before_docking_processed/chosen_pdb/'
ligand_save_path = rootdir +'/Before_docking_processed/chosen_mol2/'
F = []
start = time.time()
#pocket
for j,filename in enumerate(os.listdir(receptor_file_path)):
    old_dir = os.path.join(receptor_file_path,filename)
    if j%100 ==0:
        print(time.time()-start,j)
    for i in range(len(df)):
        if str(df['name'][i]) == filename:
            new_dir = os.path.join(receptor_save_path,filename)
            shutil.copy(old_dir,new_dir)
#ligand
for s,ligname in enumerate(os.listdir(ligand_file_path)):
    old_dir = os.path.join(ligand_file_path,ligname)
    if s%100 ==0:
        print(time.time()-start,s)
    for a in range(len(rf)):
        if str(rf['name'][a]) == ligname:
            new_dir = os.path.join(ligand_save_path,ligname)
            shutil.copy(old_dir,new_dir)
###############
#         else:
#             if filename not in F:
#                 F.append(filename)
# lost = open('/home/jianping/data/wanglin/wl_project/pocket_search/pocket_comepare_result/pocket_search_lost.txt','w')
# lost.write('\n'.join(F))
# lost.close()
#####################################################

