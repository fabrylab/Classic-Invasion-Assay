

"""
Takes a folder with .tif images in them and moves them based on their filenames into newly generated position folders.
This is particularly helpful for the Calculate-drift python programme of Lucas Heublein.
"""
import numpy as np
import os
import glob
import re
import shutil


# Method to create analysis folder
def createFolder(directory):
    try:
        if not os.path.exists(directory):
            os.makedirs(directory)
    except OSError:
        print ('Error: Creating directory. ' +  directory)






# Folder with .tif images in there. All of them will be sorted into new pos-folders
rootdir =r"H:\Experiment_data\B01_AnWi_Invasion_2020-03-04\Platte2_U87\HA"


# Define empty lists
positions = []





# Iterate over folder and search for available positions
for subdir, dirs, files in os.walk(rootdir):
    # get all available positions
    if 'pos' not in subdir:
        # NOTE: this is a hack, you will only get one file per position (just rep 0000), positions are still actually image file names
        positions = [path.split('\\')[-1] for path in glob.glob(subdir + '\*rep0_pos*_x0_y0_modeBF*_z*0.tif')]
        print('Searching folder %s - %d positions found \n ' % (subdir, len(positions)))

        # Extract pos numbers from image file names (=named position)
        # iterate over positions (which are the filenames of the first image of each position)
        for position in positions:
            search1 = re.compile('.*(?P<date>\d{8})-(?P<time>\d{6})_(?P<mic>.*)_rep(?P<rep>\d{1,6})_pos(?P<pos>\d{1,6})_.*')
            image_dict1 = search1.match(position).groupdict()

            # Ceate folder for each position
            print('Creating folder:  '+ rootdir + '/pos' + image_dict1['pos'])
            pos_folder = rootdir + '/pos' + image_dict1['pos']
            createFolder(pos_folder)
    else:
        continue


    #'Verschiebe alle tif. Files innerhalb dieses Ordners in den entsprechenden Pos-Folder'
    # Loop over .tif files within directory
    for filename in os.listdir(rootdir):
        if filename.endswith(".tif"):
            # Check .tif file for position
            search2 = re.compile('.*(?P<date>\d{8})-(?P<time>\d{6})_(?P<mic>.*)_rep(?P<rep>\d{1,6})_pos(?P<pos>\d{1,6})_.*')
            image_dict2 = search2.match(filename).groupdict()

            # Define path of file and pos-corresponding destination path
            destination_path = os.path.join(rootdir, 'pos' + image_dict2['pos'])
            file_path = os.path.join(rootdir,filename)

            # Actual command to move .tif files
            if not os.path.exists(destination_path):
                print("couldn't find position folder for %s"%filename)
                continue

            try:
                shutil.move(file_path, destination_path)
                print('\t Moved file: ' + filename)
            except:
                print('!!! Cannot move file: ' + filename)
                pass


        else:
            continue



print('Finished!')









