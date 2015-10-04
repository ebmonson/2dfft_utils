'''
    fits_to_text.py
    Purpose: Converts FITS files to .txt.
    Input: FITS image files. The string in line 13 can be modified to reflect your needs and file
    naming conventions.
    Output: .txt files with the pixel values for the corresponding image and the same basename.
'''

import glob
from pyraf import iraf

#Get filenames of all the suitable fits files in working directory
fits_list_all = glob.glob("*.fit") + glob.glob("*.fits")

if (len(fits_list_all) == 0):
    print "ERROR: No fits files matching search string were found.\nCheck directory or search string."
else:
    #Convert the .fit files to .txt
    for i in range(0,len(fits_list_all)):
        if ".fits" in fits_list_all[i]:
            basename = fits_list_all[i][0:(len(fits_list_all[i])-5)]
        else:
            basename = fits_list_all[i][0:(len(fits_list_all[i])-4)]
        iraf.wtextimage(fits_list_all[i],basename+".txt",header='no',pixels='yes')
#End program
