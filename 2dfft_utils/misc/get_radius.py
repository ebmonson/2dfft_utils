'''
    get_radius.py
    Purpose: Automates the process of finding the radii of simulation images. Uses numpy and astropy.
    Inputs:  all_centers.txt (the output from get_center.py) and all FITS images in the working directory
    Outputs: r_max.txt, a list of the galaxies' radii in the same order as all_centers.txt
    Notes: Line 28 in the Radius function removes the scales superimposed on the image from the calculation.
           You may need to modify the slice indices, depending on whether or not your image has scales and 
           where they are located.
           It is best to run this script from a directory containing a complete set of simulation images (0.000 - 3.000 Gyr).
           Otherwise, the entries in r_max.txt will not line up with those in all_centers.txt
           It is not strictly necessary to loop over the entire image; assuming symmetry and looping over a smaller part 
           would increase program speed.
    Author: Erik Monson
    
    Last modified May 29, 2015
'''
#!/usr/bin/env/ python
import numpy as np
from astropy.io import fits
import glob
import math
import sys

# Begin subroutine Radius
def Radius(filename, CENTER_X, CENTER_Y):
    img_data = fits.getdata(filename, ignore_missing_end = True) #Grab the pixel values from the image
    X_SIZE,Y_SIZE = img_data.shape
    img_data[0:80,0:600] = img_data[520:600,0:600] = 0 #Get rid of the scales at the top and bottom of the image
    radius = 1
    for i in range(Y_SIZE):
        for j in range(X_SIZE):
            if img_data[i][j] != 0: #If the pixel at (j,i) isn't black ...
                distance = int(math.ceil(math.sqrt((j-CENTER_X)**2+(i-CENTER_Y)**2))) #find the distance from the center to the pixel at (j,i) ...
                if distance > radius: #and if it's greater than the current record ...
                    radius = distance #set the distance as the new record.
    return radius
# End subroutine

# Begin main program
imgs = glob.glob("*.fit") #Grab all the FITS files in the working directory
if len(imgs) == 0:
    print "There are no FITS files (.fit) in this directory\n"
    sys.exit()
names = np.loadtxt("all_centers.txt", dtype ='S', usecols = (0,)) #Grab the base names of the images to compare with the files in the directory
centers = np.loadtxt("all_centers.txt",usecols = (1,2)) #Grab the x and y coords for the center from all_centers.txt
outfile = open("r_max.txt","w")

for i in range(len(names)):
    if (names[i]+".fit") in imgs: #Make sure that the image is in the working directory
        if len(names) > 1:
            center_x = centers[i][0]
            center_y = centers[i][1]
        else:
            center_x = centers[0]
            center_y = centers[1]
        radius = Radius(names[i]+".fit",center_x,center_y)
        outfile.write(str(radius))
        outfile.write('\n')
        print names[i] + ' ' + str(radius) #This line can be commented
outfile.close()
# End main program
