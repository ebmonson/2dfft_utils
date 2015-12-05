'''
    spiral_overlay.py
    
    DESCRIPTION:        Overlay spirals onto a FITS image.
    
    USAGE:              python spiral_overlay.py <filename>
    
    INPUTS:             Filename.
                        FITS filename can be passed as a command line argument,
                        but the program will prompt for it if it isn't.
    
    OUTPUT:             Interactive plot in MatPlotLib window.
    
    NOTES:              Uses pylab (including matplotlib and numpy) and astropy.
                        All necessary libraries included in Ureka package.
    
                        Eqns used:
    
                            r = ae^(b*theta)
    
                            x = r*cos(theta)
                            y = r*sin(theta)
    
                            tan(phi) = b
                            phi = pitch angle
    
                        
                        It is recommended that you use cropped images with this program. Image
                        and spiral plotting slows down dramatically as image size increases.
                    
                        Rotation angle and pitch can be changed while the program is running
                        using the sliders below the plot. The reset button initializes the
                        slider to default values (0.0 for rotation, 20.0 for pitch).
                        
                        Scale can be swapped between logarithmic and linear modes with the radio
                        buttons on the left side of the plot. There are two colormaps for the lograthmic
                        scale: Seismic and Greyscale. Defaults to linear scale.
                        
                        Chirality (corresponding to the sign of the pitch angle: CCW = -, CW = +)
                        can be changed with the second set of radio buttons on the right side of
                        the plot. Defaults to positive chirality.
                        
                        Arm number (corresponding to mode) can be changed with the radio buttons
                        on the right side of the plot. Defaults to 2 arms.
                        
                        Default values can be changed in the block of constants below.
'''

#-----------------------------Program begins----------------------------
#!/usr/bin/python
from pylab import *
import numpy as np
import glob
from astropy.io import fits
from matplotlib.colors import LogNorm
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider, Button, RadioButtons
import sys

#-----------------------------Default/initial values-----------------------------
DEFAULT_ARMS = 2    #Between 1 and 6, inclusive
DEFAULT_PITCH = 20.0  #Between about 0.0 and 45.0 degrees.
DEFAULT_ROTATION = 0.0 #Between -180.0 and 180.0, inclusive
DEFAULT_COLORSCALE = 'n' #n is linear, yg is logarithmic grey scale, ys is logarithmic seismic scale
#-----------------------------Command line input-----------------------------
if len(sys.argv) != 1:
    if sys.argv[1][-4:] == ".fit" or sys.argv[1][-5:] == ".fits":
        gal_file = sys.argv[1]
    else:
        print "Usage: python spiral_overlay.py <filename>"
        print "<filename> must be a FITS file."
        sys.exit()
    #endif
else:
    print "Available images: "
    for filename in (glob.glob("*.fit") + glob.glob("*.fits")):
        print filename
    #End of loop
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
plt.subplots_adjust(left= 0.15,bottom= 0.25) #set placement of plot
plt.tick_params(axis= 'both',which= 'both',bottom= 'off',top= 'off',left= 'off',right= 'off',labelbottom= 'off',labelleft= 'off')
plt.axis([0, int(len(galaxy_image)), 0, int(len(galaxy_image))]) #initial axis bounds set by this array
fig.text(0.52, 0.92, "Image: "+gal_file, fontsize= '14', ha= 'center')

#-----------------------------Spiral arms-----------------------------

#Define the maximum angle that the arms will wind through
theta_max = 3*np.pi

def SpiralPlot(galaxy_image,arm_number,pitch_angle,rotation_angle = 0.0,colorscale_option = 'n', chirality = 'CCW', line_color = ''):

    if str(colorscale_option)=='n':
        ax.imshow(galaxy_image,cmap='gray',origin='lower')
    elif str(colorscale_option) == 'yg':
        ax.imshow(galaxy_image,cmap='Greys',norm=LogNorm(),origin='lower')
    elif str(colorscale_option) == 'ys':
        ax.imshow(galaxy_image,cmap='seismic',norm=LogNorm(),origin='lower')
    #endif

    #Convert the pitch angle to radians
    pitch_rad = pitch_angle*np.pi/180
    b = np.tan(abs(pitch_rad))

    #Find the diameter and radius of the galaxy image
    gal_pixel_diameter = float(len(galaxy_image))
    gal_pixel_radius = gal_pixel_diameter/2

    #Find the scale factor a by assuming a max radius slightly larger than the actual image radius
    a = (1.15*gal_pixel_radius)/np.exp(b*theta_max)
    
    #Convert the rotation angle to radians
    rotation_angle_rad = float(rotation_angle)*np.pi/180
    
    #Create an array for theta values. For small pitch angles (<15 degrees), it may be a good
    #   idea to increase the radians to more than 6.
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
        l, = ax.plot(x,y,linewidth=2.0)
        if line_color == 'red':
            plt.setp(l,color='red')
        elif line_color == 'black':
            plt.setp(l,color='black')
        elif line_color == 'white':
            plt.setp(l,color='white')
    #End of loop
#End subroutine

#Initialize the overlay variables with the default values
rotation_angle = DEFAULT_ROTATION
colorscale_option = DEFAULT_COLORSCALE
pitch_angle = DEFAULT_PITCH
arm_number = DEFAULT_ARMS
line_color = ''
#The flags here control which button starts out pressed
if pitch_angle >= 0:
    chirality = 'CCW'
    chir_flag = 0
elif pitch_angle < 0:
    chirality = 'CW'    #If the pitch is negative, flip the chirality and drop the sign
    chir_flag = 1
    pitch_angle = np.abs(pitch_angle)
#endif
if colorscale_option == 'n':
    colorscale_flag = 0
elif colorscale_option == 'yg':
    colorscale_flag = 1
elif colorscale_option == 'ys':
    colorscale_flag = 2
#endif

SpiralPlot(galaxy_image,arm_number,pitch_angle,rotation_angle,colorscale_option,chirality,line_color)

#-----------------------------Basic GUI elements-----------------------------
#Define a consistent UI color
axcolor = 'white'

#Define a new "axis" for the slider
rot_ax = plt.axes([0.20, 0.12, 0.65, 0.03], axisbg=axcolor)

#Slider ranges from -180..180, initialized at the user's value (default 0.0).
rot_slider = Slider(rot_ax, 'Rotation angle', -180.0, 180.0, valinit=rotation_angle)

pitch_ax = plt.axes([0.20, 0.08, 0.65,0.03], axisbg=axcolor)
pitch_slider = Slider(pitch_ax, 'Pitch angle', 0.001, 45.00, valinit=pitch_angle)

#Subroutine update
#   Define a set of actions to take when the slider is modified
#   Note that since the image and all the arms have to be redrawn, this function is fairly
#   slow. As such, it's best to click on the slider rather than click and drag.
#   Changing the rotation can also make the image jump, but the built-in cursor can fix that.
def Update(val):
    #Clear the axes
    ax.cla()
    global rotation_angle
    rotation_angle = rot_slider.val
    global pitch_angle
    pitch_angle = pitch_slider.val
    SpiralPlot(galaxy_image,arm_number,pitch_angle,rotation_angle,colorscale_option,chirality,line_color)
    #Force the buffer to flip
    fig.canvas.draw_idle()
#End subroutine

rot_slider.on_changed(Update)
pitch_slider.on_changed(Update)

#Define a new "axis" for the reset button
reset_ax = plt.axes([0.47, 0.03, 0.15, 0.04])
reset_button = Button(reset_ax, 'Reset sliders', color=axcolor, hovercolor='lightsteelblue')

#Subroutine reset
#   Define actions to take when  resetbutton is pressed
def Reset(event):
    #Reset sliders to initial values
    rot_slider.reset()
    pitch_slider.reset()
#End subroutine

reset_button.on_clicked(Reset)

#Define a new axis for the save button
save_ax = plt.axes([0.80, 0.03, 0.1, 0.04])
save_button = Button(save_ax, 'Save', color=axcolor, hovercolor='lightsteelblue')

#Subroutine save
#   Define actions to take when save button is pressed
def Save(event):
    #Find the extent of the subplot that shows the image
    extent = ax.get_window_extent().transformed(fig.dpi_scale_trans.inverted())
    fig.savefig(gal_file+'_'+str(pitch_angle)+'degrees'+'.png', bbox_inches=extent)
#End subroutine

save_button.on_clicked(Save)

#Initialize a set of radio buttons for colorscale
log_ax = plt.axes([0.025, 0.60, 0.17, 0.15], axisbg=axcolor)
log_radio = RadioButtons(log_ax, ('linear','log grey', 'log seismic'), active=colorscale_flag)
fig.text(0.10,0.76,"Colorscale", fontsize='12',ha='center')

#Subroutine ScaleChange
#   Allow the user to swap between linear/log scales on the fly
def ScaleChange(label):
    global colorscale_option
    if label == 'linear':
        colorscale_option = 'n'
    elif label == 'log grey':
        colorscale_option = 'yg'
    elif label == 'log seismic':
        colorscale_option = 'ys'
    ax.cla()
    SpiralPlot(galaxy_image,arm_number,pitch_angle,rotation_angle,colorscale_option,chirality,line_color)
    fig.canvas.draw_idle()
#End subroutine

log_radio.on_clicked(ScaleChange)

#Initialize a set of radio buttons for chirality
chir_ax = plt.axes([0.8155, 0.30, 0.15, 0.15], axisbg=axcolor)
chir_radio = RadioButtons(chir_ax, ('CCW (-)','CW (+)'), active=chir_flag)
fig.text(0.8855,0.46,"Chirality", fontsize='12',ha='center')

#Subroutine ChirChange
#   Allow the user to change the overlay's chirality on the fly
def ChirChange(label):
    global chirality
    if label == 'CCW (-)':
        chirality = 'CCW'
    elif label == 'CW (+)':
        chirality = 'CW'
    ax.cla()
    SpiralPlot(galaxy_image,arm_number,pitch_angle,rotation_angle,colorscale_option,chirality,line_color)
    fig.canvas.draw_idle()
#End subroutine

chir_radio.on_clicked(ChirChange)

#Initialize a set of radio buttons for the number of arms
arm_ax = plt.axes([0.8255, 0.50, 0.12, 0.25], axisbg=axcolor)
arm_radio = RadioButtons(arm_ax, ('1','2','3','4','5','6'), active=(DEFAULT_ARMS-1))
fig.text(0.8855,0.76,"Number of Arms", fontsize='12',ha='center')

#Subroutine ArmChange
#   Allow the user to change the number of overlaid arms on the fly
def ArmChange(label):
    global arm_number
    arm_number = int(label)
    ax.cla()
    SpiralPlot(galaxy_image,arm_number,pitch_angle,rotation_angle,colorscale_option,chirality,line_color)
    fig.canvas.draw_idle()
#End subroutine

arm_radio.on_clicked(ArmChange)

#Initialize a set of radio buttons for the line color
color_ax = plt.axes([0.032, 0.30, 0.12, 0.25], axisbg=axcolor)
color_radio = RadioButtons(color_ax, ('default','red','black','white'), active=0)
fig.text(0.10,0.56,"Line Color", fontsize='12',ha='center')

#Subroutine ColorChange
#   Allow the user to change the color of overlaid arms on the fly
def ColorChange(label):
    global line_color
    line_color = label
    ax.cla()
    SpiralPlot(galaxy_image,arm_number,pitch_angle,rotation_angle,colorscale_option,chirality,line_color)
    fig.canvas.draw_idle()
#End subroutine

color_radio.on_clicked(ColorChange)

plt.show()

#-----------------------------End Program-----------------------------
