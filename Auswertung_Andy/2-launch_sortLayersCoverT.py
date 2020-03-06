import os
import sys
import subprocess

# Folder which contains pos folders with tif. images
rootdir = r'H:\Experiment_data\B01_AnWi_Invasion_2020-03-04\Platte2_U87\CNTRL'

# code that searches for all pos folders
pos_folder_list = []
print('Searching for pos folders')
for subdir, dirs, files in os.walk(rootdir):
    if "pos" in os.path.split(subdir)[1]:
        print(subdir)
        pos_folder_list.append(subdir)
print('\n -> ' + str(len(pos_folder_list)) + ' position folders found')




# launch other script, pass src_path as argument
for pos in pos_folder_list:
    process = subprocess.Popen(['python','sortLayersCoverT.py',pos], shell=True)
    process.wait()

print('Finished!')