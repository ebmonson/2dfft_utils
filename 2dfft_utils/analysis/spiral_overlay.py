'''
    spiral_overlay.py
    
    DESCRIPTION:        Overlay spirals onto a FITS image.
    
    USAGE:              python spiral_overlay.py <filename>
    
    INPUTS:             Pitch angle, number of arms, rotation angle.
                        FITS filename can be passed as a command line argument,
                        but the program will prompt for it if it isn't.
    
    OUTPUT:             Plot in python window.
    
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
                    
                        Rotation angle can be changed while the program is running using
                        the slider below the plot. The reset button initializes the slider
                        and rotation to the user-defined rotation (defaulting to 0).
                    
    REVISION HISTORY:   Written by J.E. Berlanga Medina, modified by E. Monson.
                        Last edited July 10, 2015.
'''

#-----------------------------Program begins----------------------------
#!/usr/bin/python
from pylab import *
import numpy as np
from astropy.io import fits
from matplotlib.colors import LogNorm
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider, Button
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

flag = True
while flag:
    try:
        pitch_angle = float(raw_input("\nEnter a pitch angle in degrees; include the sign.\n> "))
        flag = False
    except ValueError:
        print "\nValueError: Cannot convert input to float. Try again."
#end loop

flag = True
while flag:
    try:
        arm_number = int(raw_input("\nEnter the number of spiral arms.\n> "))
        flag = False
    except ValueError:
        print "\nValueError: Cannot convert input to int. Try again."
#end loop

print "\n"
print "Enter a rotation angle in degrees if applicable.  Hit Return or Enter if not."
print "NOTE: Positive rotation corresponds to chirality."
print "	e.g., positive rotation is CCW if spiral arms wind CCW."
rotation_angle = raw_input("> ")

if len(rotation_angle) == 0:
    rotation_angle = 0.0
else:
    flag = True
#endif

    while flag:
        try:
            rotation_angle = float(rotation_angle)
            flag = False
        except ValueError:
            print "\nValueError: Could not convert input to float.\n"
            rotation_angle = raw_input("Please enter a valid rotation angle.\n> ")
    #end loop

#Pull the rotation angle into the bounds of the slider
if rotation_angle > 180.0:
    rotation_angle = rotation_angle % 180
elif rotation_angle < -180.0:
    rotation_angle = rotation_angle % -180
#endif

print("\n")
print("Indicate use of a logarithmic color scale by typing 'y' ")
print(" or use the default linear scale by typing 'n' (no quotes).")
colorscale_option = raw_input("> ")

while str(colorscale_option)!='y' and str(colorscale_option)!='n':
    print("Please enter exactly one letter, either 'y' or 'n'.")
    colorscale_option = raw_input("> ")
#end loop

#-----------------------------Setup for plot-----------------------------
fig, ax = plt.subplots()
plt.subplots_adjust(left=0.15,bottom=0.25) #set placement of plot
plt.axis([0,int(len(galaxy_image)),0,int(len(galaxy_image))]) #initial axis bounds set by this array
#title("Image: "+gal_file+"  |  Pitch: "+str(pitch_angle)+" deg.")
fig.text(0.27,0.92,"Image: "+gal_file+"  |  Pitch: "+str(pitch_angle)+" deg.",fontsize=14)

#-----------------------------Show the FITS image-----------------------------
if str(colorscale_option)=='n':
    imshow(galaxy_image,cmap='gray',origin='lower')
elif str(colorscale_option) == 'y':
    # seismic & Greys recommended for images faint in log scale.
    #imshow(galaxy_image,cmap='Greys',norm=LogNorm(),origin='lower')
    #imshow(galaxy_image,cmap='seismic',norm=LogNorm(),origin='lower')
    #imshow(galaxy_image,cmap='YlOrBr',norm=LogNorm(),origin='lower')
    imshow(galaxy_image,cmap='gist_yarg',norm=LogNorm(),origin='lower')
#endif

#-----------------------------Spiral arms-----------------------------
pitch_rad = pitch_angle*np.pi/180
b = np.tan(abs(pitch_rad))

#Find the diameter of the galaxy system
gal_pixel_diameter = float(len(galaxy_image))

#Find the radius of the galaxy image
gal_pixel_radius = gal_pixel_diameter/2

#Define the maximum angle that the arms will wind through
theta_max = 3*np.pi

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
    if pitch_angle > 0:
        x = (-a*np.cos(theta+phase)*np.exp(b*theta) + gal_pixel_radius)
    if pitch_angle < 0:
        x = (a*np.cos(theta+phase)*np.exp(b*theta) + gal_pixel_radius)
    #End of if statements.
    y = (a*np.sin(theta+phase)*np.exp(b*theta) + gal_pixel_radius)
    #Plot the parametrically defined spiral
    plt.plot(x,y)
#End of loop

#-----------------------------Basic GUI elements-----------------------------
#Define a consistent UI color
axcolor = 'lightgoldenrodyellow'

#Define a new "axis" for the slider
axrot = plt.axes([0.20, 0.12, 0.65, 0.03], axisbg=axcolor)

#Slider ranges from -180..180, initialized at the user's value (default 0.0).
srot = Slider(axrot, 'Rotation', -180.0, 180.0, valinit=rotation_angle)

#Subroutine update
#   Define a set of actions to take when the slider is modified
#   Note that since the image and all the arms have to be redrawn, this function is fairly
#   slow. As such, it's best to click on the slider rather than click and drag.
#   Changing the rotation can also make the image jump, but the built-in cursor can fix that.
def update(val):
    
    #Clear the axes
    ax.cla()

    #Restore the image with appropriate cmap options
    if str(colorscale_option)=='n':
        ax.imshow(galaxy_image,cmap='gray',origin='lower')
    elif str(colorscale_option)=='y':
        # seismic & Greys recommended for images faint in log scale.
        #imshow(galaxy_image,cmap='Greys',norm=LogNorm(),origin='lower')
        #imshow(galaxy_image,cmap='seismic',norm=LogNorm(),origin='lower')
        #imshow(galaxy_image,cmap='YlOrBr',norm=LogNorm(),origin='lower')
        ax.imshow(galaxy_image,cmap='gist_yarg',norm=LogNorm(),origin='lower')
    #endif
    
    #Pull the new rotation angle from the slider
    rotation_angle = srot.val
    rotation_angle_rad = float(rotation_angle)*np.pi/180
    
    #Re-plot each spiral arm.
    for i in range(0,arm_number):
        phase=i*2*np.pi/arm_number + rotation_angle_rad
        #Create arrays for the x- & y-values using the theta array.
        if pitch_angle > 0:
            x = (-a*np.cos(theta+phase)*np.exp(b*theta) + gal_pixel_radius)
        if pitch_angle < 0:
            x = (a*np.cos(theta+phase)*np.exp(b*theta) + gal_pixel_radius)
        #End of if statements.
        y = (a*np.sin(theta+phase)*np.exp(b*theta) + gal_pixel_radius)
        #Plot the parametrically defined spiral
        ax.plot(x,y)
    #End of loop
    
    #Force the buffer to flip
    fig.canvas.draw_idle()
#End subroutine

srot.on_changed(update)

#Define a new "axis" for the reset button
resetax = plt.axes([0.47, 0.04, 0.1, 0.04])
button = Button(resetax, 'Reset', color=axcolor, hovercolor='0.975')

#Subroutine reset
#   Define actions to take when button is pressed
def reset(event):
    #Reset slider to initial value
    srot.reset()
#End subroutine

button.on_clicked(reset)

plt.show()

#-----------------------------End Program-----------------------------
