# -*- coding: utf-8 -*-
"""
Created on Fri Jan  5 14:46:11 2024

@author: danli
"""
import re
import os
import shutil
import glob

os.chdir('D:/Maestria/MER/SOTON/clases/deepSeaEco/assigment2/video2')
dir_list=os.listdir('D:/Maestria/MER/SOTON/clases/deepSeaEco/assigment2/video2');

output_folder = "extracted_fig"  # Name of the new folder
os.makedirs(output_folder, exist_ok=True)  # Create the folder if it doesn't exist

img_files = glob.glob('f*.jpg')

for ii in range(0,len(img_files)+1,12):
    filename = img_files[ii]
    source_file = os.path.join(os.getcwd(), filename)  # Use os.getcwd() to get the current working directory
    new_filename = f"GENIA_{filename[6:11]}_DEEP.jpg"  # New filename format
    destination_file = os.path.join(output_folder, new_filename)
    shutil.copyfile(source_file, destination_file)  # Copy and rename the image
   
    
