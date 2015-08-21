'''
    spiral_overlay.py
    
    DESCRIPTION:        Overlay spirals onto a FITS image.
    
    USAGE:              python spiral_overlay.py <filename>
    
    INPUTS:             Filename.
                        FITS filename can be passed as a command line argument,
                        but the program will prompt for it if it isn't.
    
    OUTPUT:             Interactive plot in python window.
    
    NOTES:              Uses pylab (including matplotlib and numpy) and astropy.
                        All necessary libraries included in Ureka package.
    
                        Eqns used:
    
                            r = ae^(b*theta)
    
                            x = r*cos(theta)
                            y = r*sin(theta)
    
                            tan(phi) = b
                            phi = pitch angle
    
                        Edit theta_max, title, origin of the plot (imshow)
                        & the colormap (imshow(...,cmap=..) as needed.
                        Colormaps schemes available can be seen at:
                        http://matplotlib.org/examples/color/colormaps_reference.html
    
                        Recommended cmap options are included as comments.
                    
                        Rotation angle and pitch can be changed while the program is running
                        using the sliders below the plot. The reset button initializes the
                        slider to default values (0.0 for rotation, 25.0 for pitch).
                        
                        Scale can be swapped between logarithmic and linear modes with the radio
                        buttons on the left side of the plot. Defaults to linear scale.
                        
                        Chirality (corresponding to the sign of the pitch angle: CCW = +, CW = -)
                        can be changed with the second set of radio buttons on the left side of
                        the plot. Defaults to CCW.
                        
                        Arm number (corresponding to mode) can be changed with the radio buttons
                        on the right side of the plot. Defaults to 2 arms.
                    
    REVISION HISTORY:   Written by J.E. Berlanga Medina, modified by E. Monson.
                        Last edited August 21, 2015.
'''

#-----------------------------Program begins----------------------------
#!/usr/bin/python
from pylab import *
import numpy as np
from astropy.io import fits
from matplotlib.colors import LogNorm
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider, Button, RadioButtons
import sys

#-----------------------------Command line input-----------------------------
if len(sys.argv) != 1:
    if sys.argv[1][-4:] == ".fit":
        gal_file = sys.argv[1]
    else:
        print "Usage: python spiral_overlay.py <filename>"
        print "<filename> must be a .fit file."
        sys.exit()
else:
    gal_file = str(raw_input("\nEnter filename for galaxy image.\n> "))
#endif

try:
    #pull the pixel values and save them as a uint8 numpy array
    galaxy_image = fits.getdata(gal_file)
except IOError:
    print "\nIOError: File not found. Check wd or filename and try again"
    sys.exit()
#-----------------------------Setup for plot-----------------------------
fig, ax = plt.subplots()
plt.subplots_adjust(left=0.15,bottom=0.25) #set placement of plot
plt.axis([0,int(len(galaxy_image)),0,int(len(galaxy_image))]) #initial axis bounds set by this array
#title("Image: "+gal_file+"  |  Pitch: "+str(pitch_angle)+" deg.")
fig.text(0.52,0.92,"Image: "+gal_file, fontsize='14',ha='center')

#-----------------------------Spiral arms-----------------------------
#pitch_rad = pitch_angle*np.pi/180
#b = np.tan(abs(pitch_rad))

#Find the diameter of the galaxy system
#gal_pixel_diameter = float(len(galaxy_image))

#Find the radius of the galaxy image
#global gal_pixel_radius = gal_pixel_diameter/2

#Define the maximum angle that the arms will wind through
theta_max = 3*np.pi

#Find the scale factor a by assuming a max radius slightly larger than the actual image radius
#global a = (1.15*gal_pixel_radius)/np.exp(b*theta_max)

#Convert the rotation angle to radians
#rotation_angle_rad = float(rotation_angle)*np.pi/180

#Create an array for theta values. For small pitch angles (<15 degrees), it may be a good
#   idea to increase the radians to more than 6.
#theta = np.arange(0,theta_max,0.10)

def SpiralPlot(galaxy_image,arm_number,pitch_angle,rotation_angle = 0.0,colorscale_option = 'n', chirality = 'CCW'):

    if str(colorscale_option)=='n':
        ax.imshow(galaxy_image,cmap='gray',origin='lower')
    elif str(colorscale_option) == 'y':
        # seismic & Greys recommended for images faint in log scale.
        #ax.imshow(galaxy_image,cmap='Greys',norm=LogNorm(),origin='lower')
        #ax.imshow(galaxy_image,cmap='seismic',norm=LogNorm(),origin='lower')
        #ax.imshow(galaxy_image,cmap='YlOrBr',norm=LogNorm(),origin='lower')
        ax.imshow(galaxy_image,cmap='gist_yarg',norm=LogNorm(),origin='lower')
    #endif

    pitch_rad = pitch_angle*np.pi/180
    b = np.tan(abs(pitch_rad))

    gal_pixel_diameter = float(len(galaxy_image))
    gal_pixel_radius = gal_pixel_diameter/2

    a = (1.15*gal_pixel_radius)/np.exp(b*theta_max)
    rotation_angle_rad = float(rotation_angle)*np.pi/180
    theta = np.arange(0,theta_max,0.10)

    # Plot each spiral arm.
    for i in range(0,arm_number):
        #Add the rotation angle to the phase
        phase=i*2*np.pi/arm_number + rotation_angle_rad
        #Create arrays (implicit numpy arrays) for the x- & y-values using the theta array.
        if chirality == 'CCW':
            x = (-a*np.cos(theta+phase)*np.exp(b*theta) + gal_pixel_radius)
        elif chirality == 'CW':
            x = (a*np.cos(theta+phase)*np.exp(b*theta) + gal_pixel_radius)
        #End of if statements.
        y = (a*np.sin(theta+phase)*np.exp(b*theta) + gal_pixel_radius)
        #Plot the parametrically defined spiral
        ax.plot(x,y)
        #End of loop
#End subroutine

#These are default values only
rotation_angle = 0.0
colorscale_option = 'n'
pitch_angle = 25.0
arm_number = 2
if pitch_angle > 0:
    chirality = 'CCW'
    chir_flag = 0
elif pitch_angle < 0:
    chirality = 'CW'
    chir_flag = 1
#endif

SpiralPlot(galaxy_image,arm_number,pitch_angle,rotation_angle,colorscale_option)

#-----------------------------Basic GUI elements-----------------------------
#Define a consistent UI color
axcolor = '#FFFFAD'

#Define a new "axis" for the slider
rot_ax = plt.axes([0.20, 0.12, 0.65, 0.03], axisbg=axcolor)

#Slider ranges from -180..180, initialized at the user's value (default 0.0).
rot_slider = Slider(rot_ax, 'Rotation', -180.0, 180.0, valinit=rotation_angle)

pitch_ax = plt.axes([0.20, 0.08, 0.65,0.03], axisbg=axcolor)
pitch_slider = Slider(pitch_ax, 'Pitch', 0.001, 89.125, valinit=pitch_angle)

#Subroutine update
#   Define a set of actions to take when the slider is modified
#   Note that since the image and all the arms have to be redrawn, this function is fairly
#   slow. As such, it's best to click on the slider rather than click and drag.
#   Changing the rotation can also make the image jump, but the built-in cursor can fix that.
def Update(val):
    #Clear the axes
    ax.cla()
    global rotation_angle #This is probably bad style
    rotation_angle = rot_slider.val
    global pitch_angle
    pitch_angle = pitch_slider.val
    SpiralPlot(galaxy_image,arm_number,pitch_angle,rotation_angle,colorscale_option,chirality)
    #Force the buffer to flip
    fig.canvas.draw_idle()
#End subroutine

rot_slider.on_changed(Update)
pitch_slider.on_changed(Update)

#Define a new "axis" for the reset button
reset_ax = plt.axes([0.47, 0.03, 0.1, 0.04])
reset_button = Button(reset_ax, 'Reset', color=axcolor, hovercolor='0.975')

#Subroutine reset
#   Define actions to take when button is pressed
def Reset(event):
    #Reset sliders to initial values
    rot_slider.reset()
    pitch_slider.reset()
#End subroutine

reset_button.on_clicked(Reset)

#Initialize a set of radio buttons for colorscale
log_ax = plt.axes([0.025, 0.60, 0.15, 0.15], axisbg=axcolor)
log_radio = RadioButtons(log_ax, ('linear','log'), active=0)
fig.text(0.10,0.76,"Scale", fontsize='12',ha='center')

#Subroutine ScaleChange
#   Allow the user to swap between linear/log scales on the fly
def ScaleChange(label):
    global colorscale_option
    if label == 'linear':
        colorscale_option = 'n'
    elif label == 'log':
        colorscale_option = 'y'
    ax.cla()
    SpiralPlot(galaxy_image,arm_number,pitch_angle,rotation_angle,colorscale_option,chirality)
    fig.canvas.draw_idle()
#End subroutine

log_radio.on_clicked(ScaleChange)

#Initialize a set of radio buttons for chirality
chir_ax = plt.axes([0.025, 0.40, 0.15, 0.15], axisbg=axcolor)
chir_radio = RadioButtons(chir_ax, ('CCW','CW'), active=chir_flag)
fig.text(0.10,0.56,"Chirality", fontsize='12',ha='center')

#Subroutine ChirChange
#   Allow the user to change the overlay's chirality on the fly
def ChirChange(label):
    global chirality
    chirality = label
    ax.cla()
    SpiralPlot(galaxy_image,arm_number,pitch_angle,rotation_angle,colorscale_option,chirality)
    fig.canvas.draw_idle()
#End subroutine

chir_radio.on_clicked(ChirChange)

#Initialize a set of radio buttons for the number of arms
arm_ax = plt.axes([0.8255, 0.50, 0.12, 0.25], axisbg=axcolor)
arm_radio = RadioButtons(arm_ax, ('1','2','3','4','5','6'), active=1)
fig.text(0.8855,0.76,"Arm Number", fontsize='12',ha='center')

#Subroutine ArmChange
#   Allow the user to change the number of overlaid arms on the fly
def ArmChange(label):
    global arm_number
    arm_number = int(label)
    ax.cla()
    SpiralPlot(galaxy_image,arm_number,pitch_angle,rotation_angle,colorscale_option,chirality)
    fig.canvas.draw_idle()
#End subroutine

arm_radio.on_clicked(ArmChange)


plt.show()

#-----------------------------End Program-----------------------------
