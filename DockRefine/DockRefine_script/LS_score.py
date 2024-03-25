#!/usr/bin/env python

import sys,os
import string
import subprocess
import numpy
from numpy import *
import math

scriptdir=os.path.dirname(os.path.abspath(__file__))
rootdir=os.path.dirname(scriptdir)
###############################################################
def bubble_sort(scorelist):
    count = len(scorelist)
    scoreindex = list()
    for i in range(count):
        scoreindex.append(i)
    
    for i in range(count):
        for j in range(i + 1, count):
            if scorelist[i] < scorelist[j]:
                scorelist[i], scorelist[j] = scorelist[j], scorelist[i]
                scoreindex[i],scoreindex[j] = scoreindex[j],scoreindex[i]
    return scorelist,scoreindex
###############################################################

#########################################################
pocket_list = list()
ligand_list = list()
#########################################################
#pocket score
pocketname = 'pocket_topN.info'

fp_pocket = open(scriptdir + '/input_pocket/' + pocketname)
pocket_reader = fp_pocket.read()
fp_pocket.close()

pocket_score = list()

for line in pocket_reader.splitlines():
    array = line.split(',')
    pocket_list.append(array[0])
    pocket_score.append(float(array[1]))

#########################################################
#ls-score and docking score
ls_score = list()
vina_score = list()
L_score = list()
PL_score = list()
TS_ave = list()
TS_max = list()

for item in (os.listdir(scriptdir +'/input_ligand_mol2/')):
   if('pocket_topN' in item or 'ranking' in item):
       continue
   if ('.lsalign' in item):
        ligand_list.append(item.split('.lsalign')[0])
for ii in range(len(ligand_list)):
    lsalign_name = '.lsalign' 
    #lsalign_name = '_replica-1.lsalign'
    fp_lsalign = open(scriptdir +'/input_ligand_mol2/'+ligand_list[ii]+lsalign_name)
    lsalign_reader = fp_lsalign.read()
    fp_lsalign.close()
    sub_ls_score = list()
    if('PLEASE PROVIDE TEMPLs.mol2 FILE!!!' in lsalign_reader):
        print (ligand_list[ii]+lsalign_name)
        continue
    for lsalign_line in lsalign_reader.splitlines():
        if (len(lsalign_line.split())<10):
            continue
        if('QEURY_NAME' in lsalign_line):
            continue
        if('problem' in lsalign_line):
            print (lsalign_line)
            continue
        array = lsalign_line.split()
        sub_ls_score.append(float(array[2]))
        
    #check if miss some pockets
    if (len(sub_ls_score)<len(pocket_list)):
        missnum = len(pocket_list)-len(sub_ls_score)
        for kk in range(missnum):
            sub_ls_score.append(-1.0)
    ls_score.append(sub_ls_score)
#print (ls_score)

#### testing
for tt in range(len(ligand_list)):
    #print (ligand_list[tt])
    tmp_TS = list()
    tmp_PL = list()
    for mm in range(len(pocket_score)):
        #print pocket_score[mm]
        #print pocket_list[mm], ligand_list[tt], vina_score[tt], pocket_score[mm],ls_score[tt][mm],pocket_score[mm]*ls_score[tt][mm]
        if (ls_score[tt][mm]<0):
            continue

        tmp_TS.append(ls_score[tt][mm])
        tmp_PL.append(pocket_score[mm]*ls_score[tt][mm])

    tmp_ave = numpy.mean(tmp_TS)
    tmp_max = max(tmp_TS)

    tmp_PL_ave = numpy.mean(tmp_PL)
    tmp_PL_max = max(tmp_PL)
    
    TS_ave.append(tmp_ave)
    TS_max.append(tmp_max)

    L_score.append(math.sqrt(0.5*tmp_max+0.5*tmp_ave))
    PL_score.append(math.sqrt(0.5 * tmp_PL_max + 0.5 * tmp_PL_ave))
    #L_score.append(0.5*tmp_max+0.5*tmp_ave + (-0.005)*vina_score[tt])
    
    #print ligand_list[tt],"%.3f " % float(0.5*tmp_max+0.5*tmp_ave+ (-0.1)*vina_score[tt])

sorted_score,sorted_indx = bubble_sort(L_score)
sorted_PL_score,sorted_PL_indx = bubble_sort(PL_score)

fpout = open(rootdir + 'LS_score.result' ,'w')
for ii in range(0,len(ligand_list)):
    labelflag = 0
    fpout.write(ligand_list[sorted_indx[ii]]+'\t'+str("%.4f " % sorted_score[ii])+'\n')
fpout.close()
##########################################################################################
