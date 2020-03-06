# Author: Andreas Bauer 07.02.2020, Institut für Biophyisk
# Script to find the z_positions of cells.
# Input: A z-stack of fluorescent images, sorted into a clickpoints-database
# This script adds markers to the clickpoints database in the frame of the z-position of the cell. Additionally it
# generates a csv-file with all x,y and z- positions of cells.
# This will also change default display options for tracks.




import clickpoints
import matplotlib.pyplot as plt
import os
import numpy as np
import copy
import sys
from tqdm import tqdm
from skimage.filters import gaussian, threshold_otsu
from skimage.measure import label as measure_label
from skimage.morphology import remove_small_objects
from scipy.ndimage.morphology import binary_fill_holes, binary_dilation, binary_erosion
from skimage.measure import regionprops



def detection_with_all_images():
    #### unused ####
    for i in tqdm(range(db.getImageCount())):
        image = db.getImage(frame=i, layer="modeFluo1")
        im = image.data.astype(float)
        mask, detections = detect_dog(im)

        mask_cdb = mask.astype(np.uint8)
        db.setMask(data=mask_cdb, image=image)
        # setting markers at the weighted centroids of cells in the cdb files
        try:
            db.setMarkers(image=image, x=detections[:, 0], y=detections[:, 1], type='detection_prelim')
        except Exception as e:
            print("Error", e)
        plt.figure();
        plt.imshow(mask)


def detect_dog(img,gauss_1=1,gauss_2=2,threshold="otsu", exclude_close_to_edge=False,threshold_factor=1):
    '''
    Segmentation (=identifying the area of cells). The image is bandpass-filtered (removing large/unsharp objects
    and small objects). Then the cell area is identified by thresholding. You can use otsus method for
    thresholding ("otsu"), a threshodl based on the histogram of pixels ("mean_std") or use a fixed threshold ("absolute").
    You can also increase or decrease all thresholds with a factor (threshold_factor). If you choose "absolute", the
    threshold is set to 1 and you can only change it with the threshold_factor.
    :param img: Np.ndarray; Image, e.g. the maximums projection.
    :param gauss_1: lower size for the bandpass filter
    :param gauss_2: upper size for the bandpass filter
    :param threshold: Method of thresholding. Possible values are "otsu","mean_std" and "absolute".
    :param exclude_close_to_edge: boolean; Choose if cells close to the image edge are ignored. (Probably not necessary)
    :param threshold_factor: Additional factor for the threshold.
    :return:
    '''
    th=None
    img2 = gaussian(img, gauss_1) - gaussian(img, gauss_2)
    if threshold == "otsu":
        th = threshold_otsu(img2)
    if threshold == "mean_std":
        mu, std = np.mean(np.ravel(img2)), np.std(np.ravel(img2), ddof=1)
        th = mu + 5 * std
    if threshold == "absolute":
        th = 1

    mask = img2 >th*threshold_factor
    labeled = measure_label(mask)
    regions = regionprops(labeled, intensity_image=img2)

    detections = []
    for r in regions:
        y, x = r.weighted_centroid # optional filtering all detection close to the image edge
        close_to_edge = not ((75 < x < img.shape[1] - 75) and (75 < y < img.shape[0] - 75))
        if not close_to_edge or not exclude_close_to_edge:
            detections.append((x, y))
        else:
            mask[mask==r.label]=0 # removing label from mask

    detections = np.array(detections)
    return mask,detections


def clean_up_mask(mask,closing_iterations=4,area_factor=1.5):
    '''
    Cleaning up the segmentation of cells by:
    1) Removing small holes. This somwwhat controlled by "closing_iterations" parameter. More
    iterations will fill larger holes, but will ultimately cause wierd object shapes.
    2) Excluding small objects. Objects with a size of mu - area_factor*std
    (mu: average object area, std: standard deviation of the object area) are excluded. You can choose the area_factor;
    a high factor will result in less objects beeing removed.

    :param mask: mask of cells
    :param closing_iterations: number of iterations during a binary_closing operation
    :param area_factor: Factor defining the threshold to exclude small objects. A large area_factor
    allows smaller objects (see above)
    :return:
    '''
    # binary closing
    mask_clean=copy.deepcopy(mask)
    mask_clean = binary_dilation(mask_clean,iterations=closing_iterations)
    mask_clean = binary_erosion(mask_clean, iterations=closing_iterations)
    # filling holes
    mask_clean=binary_fill_holes(mask_clean)


    # excluding small areas
    labeled = measure_label(mask_clean)
    regions = regionprops(labeled)
    areas=[r.area for r in regions]
    mu = np.mean(areas)
    std = np.std(areas, ddof=1)
    mask_clean=remove_small_objects(mask_clean, mu - area_factor*std)

    return mask_clean



def create_z_stack(db,layer):
    '''
    Creating minimum- and maximum-projections from images in database. The Images need to be sorted correctly.
    You have to specify the layer that is used for the projection. This is optimized for minimal RAM-usage
    :param db: clickpoints Database object
    :param layer: string; the layer name.
    :return:
    '''
    im_shape=db.getImage(0).data.shape
    # initialize arrays
    min_proj = np.zeros(im_shape, dtype=np.uint16) +np.inf
    max_proj = np.zeros(im_shape, dtype=np.uint16)
    min_indices = np.zeros(im_shape, dtype=np.uint16)
    max_indices = np.zeros(im_shape, dtype=np.uint16)

    for height,i in tqdm(enumerate(range(db.getImageCount()))):
        shot=db.getImage(frame=i,layer=layer).data.astype(float)
        mask = shot < min_proj
        min_proj[mask] = shot[mask]
        min_indices[mask] = height

        mask = shot > max_proj
        max_proj[mask] = shot[mask]
        max_indices[mask] = height

    return max_indices, min_indices, max_proj, min_proj

def get_max_indices_and_position(mask,max_indices):
    '''
    Estimating the z-position of cells from a segmentation mask. individual objects are identified by labeling, then
    the z-position is calculated by taking the mean of the maximum-indices in the area of each objects. This also
    returns the x-y-positions of cells by calculating the centroid of each object. Additionally it calculates the
    standard deviation of the maximum indices. A large standard is a signe for problems
    :param mask: Boolean-segmentation mask
    :param max_indices: map of maximum indices
    :return:
    '''
    labeled = measure_label(mask)
    regions=regionprops(labeled)
    max_indices_list=[]
    index_variation=[]
    pos_list=[]
    for r in regions:
        max_indices_list.append(np.mean(max_indices[r.coords[:,0],r.coords[:,1]]))
        index_variation.append(np.std(max_indices[r.coords[:,0],r.coords[:,1]]))
        pos_list.append(r.centroid)
    return max_indices_list,index_variation,pos_list

def write_to_db(max_indices_list,pos_list,index_variation,var_flags,layer="modeFluo5",marker_type_name="cells"):
    '''
    Adding the cell positions as markers (of type "track) to the database.
    :param max_indices_list: list of z-positions of cells
    :param pos_list: list of tuples; list of xy-positions of cells
    :param index_variation: list of standard deviations of maximum indices in the cell area
    :param var_flags: list,str; List of strings, that are used as Annotations to the markers in the Database
    This is used as a warning in case of high standard deviation of maximum indices.
    :return:
    '''
    db.setMarkerType(marker_type_name,color="#1fff00",mode=4)
    for ind,pos,var,var_flag in tqdm(zip(max_indices_list,pos_list,index_variation,var_flags)):
        frame=int(np.round(ind))
        new_track = db.setTrack(marker_type_name)
        db.setMarker(frame=frame, layer=layer, x=pos[1], y=pos[0],text=var_flag,track=new_track)

def write_textfile(folder,max_indices_list,pos_list):
    '''
    writing a text file with x,y and z positions of the cell
    :param folder:
    :param max_indices_list:
    :param pos_list:
    :return:
    '''
    xyz_array=np.round(np.array([[x,y,z] for (y,x),z in zip(pos_list,max_indices_list)]),2)
    np.savetxt(os.path.join(folder,"xyz_positions.txt"),xyz_array,fmt='%1.2f',header ="x,y,z")

def flag_variation(index_variation,threshold=2):
    '''
    Identifying problematic cells (cells where the standard deviation of the maximum indices in the cell area is higher
    the the threshold is problematic). This returns a list of annotations, that are added to the markers in the Database
    :param index_variation: list of standard deviations of the maximum indices in the cell area
    :param threshold: threshold of standard deviations of the maximum indices in the cell area that defines which cells are
    problematic
    :return:
    '''
    var_flags=["" if var<threshold else "\nhigh variation (%s)"%str(np.round(var,1)) for var in index_variation ]
    return var_flags

def qq_comparison(set1, set2):
    '''
    Calculating the quantiles used for a qq-plot comparing the data in set1 and set2
    :param set1: 1-D np.ndarray
    :param set2: 1-D np.ndarray
    :return:
    '''
    l = np.max([len(set1), len(set2)])
    percentile_range = np.linspace(0, 100, l)
    p_set1 = []
    p_set2 = []
    for p in percentile_range:
        p_set1.append(np.percentile(set1, p))
        p_set2.append(np.percentile(set2, p))
    p_set1 = np.array(p_set1)
    p_set2 = np.array(p_set2)
    return p_set1, p_set2

def result_comparison():
    '''
    Making a qq-plot to compare the heights of cells, predicted in two data sets. The files for the data sets are
    defined bellow. Currently this reads a text file, with comma-separated integers in the first line.
    :return:
    '''

    with open("/media/user/NG_TRANSFER/Experiment_data/B01_AnWi_Invasion_2020-01-31/U87/Ctrl/pos002/U87_Ctrl_2020-01-31_p02.txt","r") as f:
        heights_clicked=f.readline().strip().split(",")
    with open("/media/user/GINA1-BK/B01_AnWi_Invasion_2020-01-31/U87/Ctrl/pos000/U87_Ctrl_2020-01-31_p00.txt",
              "r") as f:
        heights2 = f.readline().strip().split(",")

    heights_clicked=[int(x) for x in heights_clicked]
    heights2 = [int(x) for x in heights2]
    print("numbe of cells clicked:",len(heights_clicked))
    print("numbe of cells automatically detected:", len(heights2))
    heights_clicked_qq,heights2_qq=qq_comparison(heights_clicked, heights2)
    plt.figure()
    plt.plot(heights_clicked_qq,heights2_qq)
    plt.xlabel("quatiles clicked")
    plt.ylabel("quatiles automatic")
    plt.plot([0,100],[0,100],color="C1")


import re
if __name__ == "__main__":
    # opening a database file, enter the path to the cdb-file here
    for dir,subdir,files in os.walk(r"H:\Experiment_data\B01_AnWi_Invasion_2020-03-04"):
        print(dir)
        if re.search("pos\d{2}",dir):
           #db_path=r"H:\Experiment_data\B01_AnWi_Invasion_2020-03-04\Platte2_U87\ALG\pos03\sorted.cdb"
            db_path=os.path.join(dir,"sorted.cdb")
            db=clickpoints.DataFile(db_path,"r")
            # generating mninium, maximum projections and corresponding index-maps
            max_indices, min_indices, max_proj, min_proj=create_z_stack(db,layer="modeFluo5")
            # finding the area of cell (nuclei?) by using the maximums projection
            mask, detections = detect_dog(max_proj,threshold="otsu", exclude_close_to_edge=False,threshold_factor=1)
            # filling small holes in objects and excluding small objects
            mask_clean = clean_up_mask(mask,closing_iterations=4,area_factor=1)
            # calculating z-position of cell by taking the mean of maximum-indices in the area of the cell.
            max_indices_list,index_variation,pos_list=get_max_indices_and_position(mask_clean,max_indices)
            # identifying cells where the maximum-indices have a high standard deviation, these could be problematic and are
            # annotated
            var_flags=flag_variation(index_variation,threshold=2)
            # adding markers (as tracks) to the database
            write_to_db(max_indices_list,pos_list,index_variation,var_flags,layer="modeFluo5",marker_type_name="cell_in_focus")
            # writing a text file with x,y,z positions of the cells
            write_textfile(os.path.split(db_path)[0],max_indices_list,pos_list)

            # trying to set display options for tracks (does this work?)
            try:
                db.setOption("tracking_show_trailing", 300)
                db.setOption("tracking_show_leading", 300)
            except:
                pass
            # closing the database object
            db.db.close()






#plt.figure();plt.imshow(max_proj)
#plt.figure();plt.imshow(min_proj)
#plt.figure();plt.imshow(max_indices)