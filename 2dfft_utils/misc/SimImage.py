'''
    SimImage.py
    A python class interface that collects many of the tasks needed to prepare simulated galaxy images
    for 2DFFT and pitch angle measurement in one place.
    
    To use this class, include this import statement (without quotes) in the header of your python script 
    (or type it into your python shell of choice):
        "from SimImage import SimImage"
    
    Methods:    SimImage(filename, radius = None, center = (None,None)) - Constructor; can be passed a
                    radius and a center tuple.
                GetCenter() - Finds the image center with an IRAF (Pyraf) task.
                GetRadius() - Finds the galaxy's radial extent using Astropy. Requires a defined center.
                Crop() - Crops the image with an IRAF (Pyraf) task. Requires a defined radius and center.
                ConvertText() - Converts the image to a .txt file with an IRAF (Pyraf) task.
                ScriptString() - Returns a string formatted for the 2DFFT scripter's input.
                ImgStats() - Returns a string containing the name, center, and radius of the image.
    
    Methods include exception handling where relevant and return a Linux standard exit status (0 = success, 1 = error) 
    where relevant. However, IRAF converts some errors into warnings which can't be caught by try-except. This is 
    common when the name of an input file is entered incorrectly, so it is best to double-check filenames.
    
    This class was designed and tested with a very specific set of disk-on simulation images in mind, so 
    modification may be necessary to extend it to other simulations or real galaxy images.
    
    Several methods are modified or copied from programs written by J.E. Berlanga Medina.
    
    The methods GetCenter, Crop, and ConvertText use Pyraf; GetRadius uses Astropy.
    
    Class written by E. Monson
    Last edited June 22, 2015.
    
TODO - make GetCenter work with filenames outside our convention.
'''

from pyraf import iraf
from astropy.io import fits
import numpy as np

class SimImage():
    
    #Constructor expects only the name of the image file, but can be passed
    #the radius and center as well.
    def __init__(self, filename, radius = None, center = (None,None)):
        self.name = filename
        self.radius = radius
        self.center = center
    
    #GetCenter uses the IRAF imcntr task to find the center of the image.
    #Based on the size of the image, the first guess for center may need to be changed.
    #cboxsize should be an odd number.
    #See get_center.py (J. Berlanga) for the original implementation and documentation.
    def GetCenter(self):
        try:
            center_string = iraf.imcntr(self.name,300,300,cboxsize=31, Stdout=1)
            x_coord = int(round(float(center_string[0][19:26])))
            y_coord = int(round(float(center_string[0][33:40])))
            self.center = (x_coord,y_coord)
            return 0
        except IndexError:
            print "IndexError: Image not found or incompatible."
            return 1
        except iraf.IrafError as e:
            print "IRAF Error: " + e
            return 1

    #GetRadius uses the fits module from astropy to pull the pixel values of the image.
    #Based on the size of the image and whether it has overlaid scales, the indices
    #in the fourth line may need to be changed.
    #Note that GetRadius requires that self.center be defined, e.g. by calling GetCenter().
    #See get_radius.py (E. Monson) for the orignal implementation and documentation.
    def GetRadius(self):
        try:
            img_data = fits.getdata(self.name, ignore_missing_end = True)
            X_SIZE,Y_SIZE = img_data.shape
            img_data[0:80,0:600] = img_data[520:600,0:600] = 0
            radius = 1
            for i in range(Y_SIZE):
                for j in range(X_SIZE):
                    if img_data[i][j] != 0:
                        distance = int(np.ceil(np.sqrt((j-self.center[0])**2+(i-self.center[1])**2)))
                        if distance > radius:
                            radius = distance
            self.radius = radius
            status = 0
        except TypeError:
            print "TypeError: Image center not defined."
            status = 1
        except IOError:
            print "IOError: Image not found or incompatible."
            status = 1
        finally:
            return status
    
    #Crop uses the IRAF imcopy task to crop a galaxy image to a square filled by the galaxy.
    #Creates a new file name_crop.fit in the current directory.
    #See auto_crop_fits.py (J. Berlanga) for the original implementation and documentation.
    def Crop(self):
        try:
            iraf.imcopy(self.name+"["+str(int(self.center[0]-self.radius))+":"+str(int(self.center[0]+self.radius))+","+str(int(self.center[1]-self.radius))+":"+str(int(self.center[1]+self.radius))+"]",self.name[:-4]+"_crop.fit")
            status = 0
        except TypeError:
            print "TypeError: Image center or radius not defined."
            status = 1
        except iraf.IrafError as e:
            print "IRAF Error: " + e
            status = 1
        finally:
            return status
                
    #ConvertText uses the IRAF wtextimage task to convert images (usually cropped images) to text files.
    #Creates a new file, name.txt, in the current directory.
    #ConvertText is particularly susceptible to returning a misleading exit status due to the way that IRAF
    #handles errors - if given an invalid filename, its exit status will remain 0 even though the process does not
    #complete.
    def ConvertText(self):
        try:
            iraf.wtextimage(self.name,self.name[:-4]+".txt",header='no',pixels='yes')
            return 0
        except iraf.IrafError as e:
            print "IRAF Error: " + e
            return 1

    #ScriptString returns a string formatted for the 2DFFT scripter input.
    def ScriptString(self):
        return self.name[0:8]+"_crop.txt"+","+self.name[0:8]+","+str(int(self.radius-1))

    #ImgStats returns a string containing the image keyword (i.e., the snapshot date),
    #the center coordinates, and the radius, separated by tabs.
    def ImgStats(self):
        return self.name[0:8]+"\t"+str(self.center)+"\t"+str(self.radius)
