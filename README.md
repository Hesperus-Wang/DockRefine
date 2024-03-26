# DockRefine
This is the implementation of DockRefine: a post-docking rescoring method based on ligand similarity comparison and Drug-Target interaction
## About DockRefine
The DockRefine rescoring process is based on the ligand similarity comparison method and uses LS-align to search the similarity of compound libraries and give a similarity score. GAT was also used to learn the DTI patterns of target proteins and their corresponding active compounds, and DTI classification predictions were given for the compound library. By lowering the ranking of compounds with the incorrect DTI and increasing the ranking of compounds with the correct DTI to modify the ligand similarity score, the problem of poor specificity of the ligand similarity search method and more false positive results of the head compound can be solved.
## Installation
DockRefine is built on Python3, we recommend using Anaconda environment as enviroment management for the installation of DockRefine and its dependencies. The virtual environment can be created as follows:
```
conda create -n your_environment python==3.7
conda activate your_environment
```
Install ZeroBind dependencies as following:
```
conda install rdkit rdkit
conda install openbabel -c conda-forge
conda install -c conda-forge smina
conda install pytorch==1.9.0 torchvision==0.10.0 torchaudio==0.9.0 cuda-toolkit=10.2 –c pytorch
```
## Executing workflow
DockRefine needs two consecutive steps (please ensure the last step was completed before running the next step):
1. prepare the screening compounds (MOL2 format and SDF format) and pocket ('.POC' format) into the input folder.
put all compounds into (DockRefine/input_ligand_mol2/ and DockRefine/input_ligand_sdf/) subfolder.
put the pocket file (named as "nativepocket.poc") into the (DockRefine/input_pocket/) subfolder.The pocket is like as:
```
POC     nativepocket
ATOM    253  N   VAL A  49       6.937   3.497 -28.457  1.00 34.68           N
ATOM    254  CA  VAL A  49       6.465   2.179 -28.887  1.00 34.68           C
ATOM    255  C   VAL A  49       5.665   2.297 -30.179  1.00 34.68           C
ATOM    256  O   VAL A  49       5.965   1.651 -31.192  1.00 34.68           O
...
TRE
```
