#!/usr/bin/env python
import sys, os
import string
import subprocess

###############################################################
scriptdir = os.path.dirname(os.path.abspath(__file__))
rootdir = os.path.dirname(scriptdir)
# outputdir = rootdir + '/user_input/results_docking/pocket/'
outputdir = rootdir + '/wl_project/1.pocket_search/input_pocket/'

BioLip = rootdir + '/BioLip/'
# BioLip = '/media/jianping/839282ac-3ae6-4224-98cb-8c5cf2c2cf3e/wanglin/EViS-main/BioLip/'
poc_size = 'POC'


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
                scoreindex[i], scoreindex[j] = scoreindex[j], scoreindex[i]
    return scorelist, scoreindex
###############################################################
# get target pocket
# print (os.listdir(rootdir+'/user_input/results_docking/pocket/')
# for ii in os.listdir(rootdir + '/user_input/results_docking/pocket/'):
# for ii in os.listdir(rootdir + '/user_input/results_docking/pocket/'):
for ii in os.listdir(rootdir+'/1.pocket_search/input_pocket/'):
    # print(ii)
    if ('.poc' in ii):
        pocketname = str(ii)
        print(pocketname)
        break
targetname = pocketname.split('.poc')[0]

fp = open('pocket_num', 'r')
fpreader = fp.read()
fp.close()

pocket_information = list()
flagindex = 0
template_num = len(fpreader.splitlines())
for line in fpreader.splitlines():
    flagindex = flagindex + 1
    array = line.split()
    templateligand = array[0]
    templatepoc = array[0].split('.mol2')[0] + '.pdb'
    command1 = 'cd ' + outputdir
    command2 = scriptdir + '/PPSalign ' + pocketname + ' ' + BioLip + str(
        poc_size) + '/' + templatepoc + ' >align_info' + '_' + poc_size
    stdout, stderr = subprocess.Popen(';'.join([command1, command2]), shell=True).communicate()

    falign = open(outputdir + 'align_info_' + poc_size, 'r')
    align_reader = falign.read()
    falign.close()

    tmp_score = 0.0
    tmp_query_aa = ''
    tmp_templ_aa = ''
    for align_line in align_reader.splitlines():
        if targetname in align_line:
            # print align_line
            tmp_score = float(align_line.split()[2])
            # print flagindex,tmp_score
            continue
        if 'Query AA Index:' in align_line:
            tmp_query_aa = align_line.split(':')[1]
            continue
        if 'Templ AA Index:' in align_line:
            tmp_templ_aa = align_line.split(':')[1]
            break
    print(str(flagindex) + '/' + str(template_num) + '\t' + templatepoc + '\t' + str(tmp_score))
    if (tmp_score >= 0.5):
        pocket_information.append(templatepoc + ',' + str(tmp_score) + ',' + tmp_query_aa + ',' + tmp_templ_aa)
        print(templatepoc + ',' + str(tmp_score) + ',' + tmp_query_aa + ',' + tmp_templ_aa, '#', flagindex)

fpout = open(outputdir + '/pocket_align.' + poc_size, 'w')
for letter in pocket_information:
    fpout.write(letter + '\n')
fpout.close()
