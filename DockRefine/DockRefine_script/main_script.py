import sys,os
import string
import subprocess
import time

scriptdir=os.path.dirname(os.path.abspath(__file__))
rootdir=os.path.dirname(scriptdir)
BioLip = rootdir + '/BioLip/'  ### update

########################1.data preprocess part

# step 1: pocket search
print ("step 1/7 run pocket_search.py")
command= scriptdir + '/pocket_search.py'
stdout,stderr = subprocess.Popen(';'.join([command]),shell=True).communicate()

# step 2: chosen pocket and ligand
print ("step 2/7 run chose_pdb_and_mol2.py")
command= scriptdir + '/chose_pdb_and_mol2.py'
stdout,stderr = subprocess.Popen(';'.join([command]),shell=True).communicate()

#step 3: preprocessed before docking
print ("step 3/7 run preprocessed.py")
command= scriptdir + '/preprocessed.py'
stdout,stderr = subprocess.Popen(';'.join([command]),shell=True).communicate()
#
# step 4: config generate and run smina
print ("step 4/7 run config_generate_and_run_smina.py")
command= scriptdir + '/config_generate_and_run_smina.py'
stdout,stderr = subprocess.Popen(';'.join([command]),shell=True).communicate()

# step 5: calculate RMSD
print ("step 5/7 run not_align_rmsd.py")
command= scriptdir + '/not_align_rmsd.py'
stdout,stderr = subprocess.Popen(';'.join([command]),shell=True).communicate()

# step 6: generate train and test set
print ("step 6/7 run generate_train_test_set.py")
command= scriptdir + '/generate_train_test_set.py'
stdout,stderr = subprocess.Popen(';'.join([command]),shell=True).communicate()

# step 7: find template pocket_topN
print ("step7/7 run ranking and select the maximum 50 templates")
command= scriptdir + '/pocket_topN.py'
stdout,stderr = subprocess.Popen(';'.join([command]),shell=True).communicate()
command = 'cp ' + scriptdir +'/input_pocket/pocket_topN* ' + scriptdir +'/input_ligand_mol2/'
stdout,stderr = subprocess.Popen(';'.join([command]),shell=True).communicate()

########################2.ligand similarity search part
# step 1: LS-align
print ("step1/2 LS-align")
command1 = 'cd ' + scriptdir +'/input_ligand_mol2/'
total_mol2 = list()
lsalign_index = 0
for mol2_item in (os.listdir(scriptdir +'/input_ligand_mol2/')):
    if ('pocket_topN.mol2' in mol2_item):
        continue
    if '.mol2' in mol2_item :
        total_mol2.append(mol2_item)

for ii in range(len(total_mol2)):
    lsalign_index = lsalign_index + 1
    print ('LSalign ' + total_mol2[ii] + ' ----' + str(lsalign_index) + '/' + str(len(total_mol2)))
    command2 = scriptdir + '/LSalignO pocket_topN.mol2 ' + total_mol2[ii] + ' >' + total_mol2[ii].split('.mol2')[0]+'.lsalign'
    stdout,stderr = subprocess.Popen(';'.join([command1, command2]),shell=True).communicate()

# step 2: LS-score
print ("step2/2 run LS-score")
command= scriptdir + '/LS_score.py'
stdout,stderr = subprocess.Popen(';'.join([command]),shell=True).communicate()

########################3.DTI prediction part
# step 1: generate complex info
print ("step 1/4 run generate_train_test_set.py")
command= scriptdir + '/generate_dir_txt.py'
stdout,stderr = subprocess.Popen(';'.join([command]),shell=True).communicate()

# step 2: generate keys.pkl
print ("step 2/4 run divide.py")
command= scriptdir + '/divide.py'
stdout,stderr = subprocess.Popen(';'.join([command]),shell=True).communicate()

# step 3: train model
print ("step 3/4 run train.py")
cmd_list = ['python train.py','--epoch=500','--batch_size=32','--dropout_rate=0.3','--n_graph_layer=3 ','--total_key= /home/jianping/data/wanglin/DockRefine/data_for_classification/key/total/','--ngpu=0']
line = ' '.join(cmd_list)
print(line)
os.system(' '.join(cmd_list))

# step 4: predict
print ("step 4/4 run predict.py")
cmd_list = ['python predict.py','--batch_size=1']
line = ' '.join(cmd_list)
print(line)
os.system(' '.join(cmd_list))

########################4.ranking part
print ("run duplicates_and_merge_score.py")
command= scriptdir + '/duplicates_and_merge_score.py'
stdout,stderr = subprocess.Popen(';'.join([command]),shell=True).communicate()


