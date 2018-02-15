
# coding: utf-8

# In[ ]:

# In Class Lab 
# Feb 15 2018


# In[ ]:

# import modules
import numpy as np
import astropy.units as u
from astropy.constants import G

# import plotting modules
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
get_ipython().magic(u'matplotlib inline')

# my modules
from ReadFile import Read
from CenterOfMass import CenterOfMass
from MassProfile import MassProfile

# yt 
import yt
from yt.units import kiloparsec, Msun, Gyr


# In[ ]:

# Create a COM of object for MW Disk Using Code from Assignment 4
# using highest res files - HR

COMD = CenterOfMass("MW_HR_000.txt",2)


# In[ ]:

# Compute COM of MW using disk particles
COMP = COMD.COM_P(0.1, 4.0)
# In assignment 6 you are asked to modify the CenterOfMass.py file so that 
# in addition to a 'delta' (tolerance) there is also a volume decrement value.  
# dividing the volume by 4.0 actually works better than dividing by half (RMAX/4 instead)

# store COM Velocity
COMV = COMD.COM_V(COMP[0],COMP[1],COMP[2])


# ## (PART 1) Edit Below:  
# ### Define the remaining data properties: yD, zD, vxD, vyD, vzD 

# In[ ]:

# Determine positions and velocities of disk particles relative to COM motion
xD = COMD.x


# ## (PART 2) Edit Below: 
# ### What position coordinates should you plot against each other to create a face-on representation of a galaxy's disk? Use the hist2d function in matplotlib to test this. Start with bins=100.

# In[ ]:

# MW Disk Density, face on

fig = plt.figure(figsize=(13,10))
ax = plt.subplot(111)

#### PLACE FUNCTION HERE ####
plt.colorbar()

# Add axis labels
plt.xlabel('? (kpc)', fontsize=22)
plt.ylabel('? (kpc)', fontsize=22)

#set axis limits
plt.ylim(-30,30)
plt.xlim(-30,30)

#adjust tick label font size
label_size = 22
matplotlib.rcParams['xtick.labelsize'] = label_size 
matplotlib.rcParams['ytick.labelsize'] = label_size

# Save to a file
ax.set_rasterized(True)


# ## (PART 3) Edit above:
# 
# ### 1. Add the norm=LogNorm() argument to the hist2d function (Note: You imported this module at the top of the script). Your colorbar should now show logarithmic units. What does the colorbar represent?
# 
# ### 2. Increase the number of bins used in hist2d to get a higher resolution plot.
# ## ----------------------------------------------------------------------------------------------------------

#  ## (PART 4) Edit Below:
#  
#  ### Now let's make a density plot that shows the edge-on MW disk. Using hist2d again, which position coordinates should you plot against each other to get an edge-on view of the disk? 
#  
#  ### Use the code below to create this plot. You'll want to use a similar number of bins and the norm=LogNorm() argument again.
# 

# In[ ]:

# MW Disk Density, edge on

fig = plt.figure(figsize=(18,5))
ax = plt.subplot(111)

#### PLACE FUNCTION HERE ####
plt.colorbar()

# Add axis labels
plt.xlabel('? (kpc)', fontsize=22)
plt.ylabel('? (kpc)', fontsize=22)

#set axis limits
plt.ylim(-3,3)
plt.xlim(-30,30)

#adjust tick label font size
label_size = 22
matplotlib.rcParams['xtick.labelsize'] = label_size 
matplotlib.rcParams['ytick.labelsize'] = label_size


# Save to a file
ax.set_rasterized(True)


# You'll see the disk is fairly thin with a scale height < 1 kpc. Notice the X-shape at the center of the disk as well. This is due to the interaction between disk stars and the bar. Bulge particles are not shown here. 

# ## (PART 5) Edit Above: 
# 
# ### Does the densest part of the MW disk in your face-on and edge-on plots reside at position=(0,0)? Re-center your particles using COMP. You can do this by changing the quantity COMP back to floats by doing the following to all position particles:
# 
# xD = COMD.x - float(COMP[0]/u.kpc)}
# 
# ### Do the same for velocity using COMV. Be careful with the units.
# 
# ### Recreate the two plots above in the center of mass reference frame.
# ## ----------------------------------------------------------------------------------------------------------

# # *** OPTIONAL *** 
# 
# ### You can change the colorbar used in the hist2d function using the cmap argument. Set cmap to one of the colorbars (i.e. cmap='magma') from this page: https://matplotlib.org/users/colormaps.html
# 
# ### The perceptually uniform sequential color maps are the best options because they are also legible in grayscale.
# ## ----------------------------------------------------------------------------------------------------------

# ## (PART 6) Edit Below:
# 
# ### Now we will create a phase diagram and compare it to the circular velocity assuming spherical symmetry. If you look at the MW's disk along the y-axis, what velocity component will give you the line of sight velocity? Use your intuition from the edge-on disk plot above. 
# 
# ### Use hist2d again to plot the y-position of the disk particles against this velocity component. 
# 
# ### Compare your plot to the plots in InClass5_Slides.pdf. Does your plot match the observational data?

# In[ ]:

# MW Disk Velocity Field edge on.

fig = plt.figure(figsize=(10,10))
ax = plt.subplot(111)

# looking at MW edge on along y axis, what is line of sight velocity component?

#### PLACE FUNCTION HERE ###
plt.colorbar()


# Add axis labels
plt.xlabel('? (kpc)', fontsize=22)
plt.ylabel('Velocity ? (km/s)', fontsize=22)

#set axis limits
#plt.ylim(0,200)
plt.xlim(-30,30)

#adjust tick label font size
label_size = 22
matplotlib.rcParams['xtick.labelsize'] = label_size 
matplotlib.rcParams['ytick.labelsize'] = label_size


# Save to a file
ax.set_rasterized(True)


# In[ ]:

# Create a Mass Profile Object:  Code from Assignment 5
MWV = MassProfile("MW",0)
# note that since CenterOfMass has been modified to haveVolDec input, 
# MassProfile.py must be changed to reflect the additional input argument

# Input array of radii
R = np.arange(0.01,30,0.2)

# Store the total circular velocity for later use 
VCirc = MWV.CircularVelocityTotal(R)


# ## (PART 7) Edit Above:
# 
# ### Now overplot the MW's circular velocity computed with the CircularVelocityTotal() against the radius array. These are already defined as 'VCirc' and 'R' above. 
# 
# ### Notice the axis limits of the y-axis above. Hint: If you have an array called X, then -X will multiply the whole array by -1. 
# 
# ### How does the phase diagram compare to the spherically averaged circular velocity? Why is there a spread in velocities?
# 
# ## ----------------------------------------------------------------------------------------------------------

# ## (PART 8) Edit Below: 
# 
# ### Let's isolate an interesting set of particles and explore their origin. This may be helpful to many of your for your final projects. For example, where do the particles with the highest speeds reside and why do they have such high speeds?
# 
# ### 1. Create a velocity mask that selects all particles from the plot above that have a velocity of 250 km/s or higher. Note that these velocities can be positive or negative.
# 
# ### 2. Using the plotting code above, create an edge-on plot of just these disk particles. How did these disk stars arrive at this position?
# ## ----------------------------------------------------------------------------------------------------------

# # *** OPTIONAL: Using the yt package to make particle plots ***
# 
# Using Anaconda you can install YT on your computer. From the command line type:
# conda install yt
# 
# You can then open an interactive notebook in your working directory (where your data files are located)
# yt notebook
# 
# Which will give a message stating that it needs a password (give it any number) and it will also state : "The notebook is now live at: http: blah" Copy the URL into a web browser and it will open a Jupyter notebook that will access your working directory. Then proceed as normal with python script.
# 
# YT DOCS: http://yt-project.org/docs/dev/visualizing/plots.html#particle-plots http://yt-project.org/doc/examining/generic_particle_data.html?highlight=loading%20generic%20data
# 
# ## *See the solutions for this lab for examples of how to use yt to make the same plots as above.*

# In[ ]:



