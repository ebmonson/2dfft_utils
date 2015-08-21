'''
    sim_prep.py
    A simple script using SimImage.py to prepare an entire sim batch at once.
    
    Erik Monson
    Last edited August 21, 2015
'''

#/usr/bin/env python
import glob
from SimImage import SimImage

#------Begin program------
all_fits = glob.glob("*.fit")
sim_name = str(raw_input("Enter the name of this simulation (e.g. a1_lores): "))

scripter_file = open("scripter_input-"+sim_name+".txt",'w')
scripter_file.write('\n')

for filename in all_fits:
    image = SimImage(filename)
    if not image.GetCenter():
        image.GetRadius()
        print image.ImgStats()
        if not image.Crop():
            crop_image = SimImage(image.name[:-4]+"_crop.fit")
            crop_image.ConvertText()
            del crop_image
            scripter_file.write(image.ScriptString()+'\n')
    del image

scripter_file.write('\n')
scripter_file.close()
#-----End program-----
