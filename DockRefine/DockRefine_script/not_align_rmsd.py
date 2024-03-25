#!/usr/bin/env python
import os
import math
import shutil
from oddt.toolkits import rdk as toolkit
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem.AllChem import AlignMol

scriptdir=os.path.dirname(os.path.abspath(__file__))
rootdir=os.path.dirname(scriptdir)

def GetBestRMSD(probe, ref, refConfId=-1, probeConfId=-1, maps=None, align=True):
    # When mapping the coordinate of probe will changed!!!
    ref.pos = orginXYZ(ref)
    probe.pos = orginXYZ(probe)
    try:
        name = probe.GetProp("_Name")
    except KeyError as e:
        name = "NaN"
    if not maps:
        matches = ref.GetSubstructMatches(probe, uniquify=False)
        if not matches:
            raise ValueError(
                "mol %s does not match mol %s"
                % (ref.GetProp("_Name"), probe.GetProp("_Name"))
            )
        if len(matches) > 1e6:
            print(
                "{} matches detected for molecule {}, this may lead to a performance slowdown.".format(
                    len(matches), name
                )
            )
        maps = [list(enumerate(match)) for match in matches]
    bestRMSD = 10000.0
    for amap in maps:
        if align:
            rmsd = AlignMol(probe, ref, probeConfId, refConfId, atomMap=amap)
        else:
            rmsd = RMSD_NotAlign(probe, ref, amap)
        bestRMSD = min(bestRMSD, rmsd)
    return bestRMSD

def RMSD_NotAlign(probe, ref, amap):
    rmsd = 0.0
    atomNum = ref.GetNumAtoms() + 0.0
    for (pi, ri) in amap:
        posp = probe.pos[pi]
        posf = ref.pos[ri]
        rmsd += dist_2(posp, posf)
    rmsd = math.sqrt(rmsd / atomNum)
    return rmsd

def orginXYZ(mol):
    mol_pos = {}
    for i in range(0, mol.GetNumAtoms()):
        pos = mol.GetConformer().GetAtomPosition(i)
        mol_pos[i] = pos
    return mol_pos

def dist_2(atoma_xyz, atomb_xyz):
    dis2 = 0.0
    for i, j in zip(atoma_xyz, atomb_xyz):
        dis2 += (i - j) ** 2
    return dis2

raw_file_dir = os.path.abspath(rootdir+'/Before_docking_processed/preprocessed_sdf/')
docking_out_dir = os.path.abspath(rootdir+'/2.1_docking_result_split/')
act = rootdir+'/data_for_classification/act/'
inact = rootdir+'/data_for_classification/inact/'
sdf_dir = rootdir+'/Before_docking_processed/preprocessed_sdf/'

count = 0
process =0
for file in os.listdir(raw_file_dir):
    file_path = os.path.join(raw_file_dir,file)
    ref_str = file.split('_')[0] + file.split('_')[2]
    # ref_str = file.split('_')[1]
    target_list = [os.path.join(docking_out_dir, i) for i in os.listdir(docking_out_dir) if ref_str in i]
    for target in target_list:
###########################################
        rfile = file_path
        pfile = target
###########################################

        ref_fmt = rfile.split(".")[-1]
        ref_oddt = next(toolkit.readfile(ref_fmt, rfile))
        try:
            ref_rdk = Chem.RemoveHs(ref_oddt.Mol)
            # count+=1
        except:
            process+=1
            pass
        else:
            probe_fmt = pfile.split(".")[-1]
            probe_oddt_supp = toolkit.readfile(probe_fmt, pfile)

            data = []
            for i, probe_oddt in enumerate(probe_oddt_supp):
                process += 1
                if probe_oddt is None:
                    name = "NaN"
                    rmsd_notalign = 10000.0
                    rmsd_align = 10000.0
                else:
                    probe_rdk = Chem.RemoveHs(probe_oddt.Mol)
                    try:
                        name = probe_rdk.GetProp("_Name")
                        name = "_".join(name.split())
                    except KeyError as e:
                        name = "NaN"
                    # ref = ref_rdk

                    try:
                        ref = AllChem.AssignBondOrdersFromTemplate(probe_rdk, ref_rdk)
                    except:
                        count+=1
                    else:
                        rmsd_notalign = GetBestRMSD(probe_rdk, ref, align=False)
                        # data.append((name,rmsd_notalign))
                        # print(data)
                        if rmsd_notalign <= 2:
                            # act
                            shutil.copy(target, target.replace(docking_out_dir, act))
                        elif rmsd_notalign >= 4:
                            shutil.copy(target, target.replace(docking_out_dir, inact))
                        else:
                            print(rmsd_notalign)
print(count,'total',process)

for i in os.listdir(act):
    new_name = 'active_' +i
    os.rename(os.path.join(act,i),os.path.join(act,new_name))
for j in os.listdir(inact):
    new_name = 'inactive_' +j
    os.rename(os.path.join(inact, j), os.path.join(inact, new_name))
for r in os.listdir(sdf_dir):
    new_name = 'active_'+r.split('_')[0]+r.split('_')[2]+'_raw.sdf'
    os.rename(os.path.join(sdf_dir, r), os.path.join(sdf_dir, new_name))
for k,filename in enumerate(os.listdir(raw_file_dir)):
    old_dir = os.path.join(raw_file_dir, filename)
    if filename.split('.')[-1] == 'sdf':
        new_dir = os.path.join(act, filename)
        shutil.copy(old_dir, new_dir)



