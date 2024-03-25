#!/usr/bin/env python
import os
import random
import shutil

scriptdir=os.path.dirname(os.path.abspath(__file__))
rootdir=os.path.dirname(scriptdir)

def capturerandomfile(filedir,traindir):#generate train set
    files = os.listdir(filedir)
    file_number = len(files)
    rate1 = 0.7
    # rate1 = 0.799
    picknumber_1 = int(file_number * rate1)
    random.seed(42)
    sample_train = random.sample(files,picknumber_1)
    for name in sample_train:
        shutil.copyfile(filedir +name , traindir +name)

def capturerepeatfile(filedir,traindir,testdir):#generate test set
    files = os.listdir(filedir)
    train_list = os.listdir(traindir)
    test_list = os.listdir(testdir)
    for i in files:
        if i not in train_list:
            shutil.copyfile(filedir +i , testdir +i)
        else:
            pass

def generate_new_inact(filedir,new_inact_dir,rate):#generate train set
    files = os.listdir(filedir)
    file_number = len(files)
    picknumber_1 = int(file_number * rate)
    random.seed(42)
    sample_train = random.sample(files,picknumber_1)
    for name in sample_train:
        shutil.copyfile(filedir +name , new_inact_dir +name)

act_filedir = rootdir + '/data_for_classification/act/'
inact_filedir = rootdir +'/data_for_classification/inact/'
new_inact_filedir = rootdir +'/data_for_classification/new_inact/'
train_dir = rootdir +'/data_for_classification/train/'
test_dir = rootdir +'/data_for_classification/test/'
total_dir = rootdir +'/data_for_classification/total/'

#generate_new_inact_set
if len(os.listdir(act_filedir)) < len(os.listdir(inact_filedir)):
   extract_rate = round((len(os.listdir(act_filedir))*2.5)/len(os.listdir(inact_filedir)),2)
   print(extract_rate)
   generate_new_inact(inact_filedir, new_inact_filedir, extract_rate)

capturerandomfile(act_filedir,train_dir)
capturerepeatfile(act_filedir,train_dir,test_dir)
capturerandomfile(new_inact_filedir,train_dir)
capturerepeatfile(new_inact_filedir,train_dir,test_dir)
for i in os.listdir(train_dir):
    shutil.copyfile(train_dir+i,total_dir+i)
for j in os.listdir(test_dir):
    shutil.copyfile(test_dir+j,total_dir+j)
print('done')

