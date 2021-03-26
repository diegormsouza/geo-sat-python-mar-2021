# Training: Python and GOES-R Imagery: Script 8 - Custom Colormaps - Enhancing IR Channels
#----------------------------------------------------------------------------------------------------------- 
# Required modules
from netCDF4 import Dataset          # Read / Write NetCDF4 files
import matplotlib.pyplot as plt      # Plotting library
from datetime import datetime        # Basic Dates and time types
from matplotlib import cm            # Colormap handling utilities
import numpy as np                   # Scientific computing with Python
from utilities import loadCPT        # Import the CPT convert function
#----------------------------------------------------------------------------------------------------------- 

# COLORMAP EXAMPLE 1
# Custom colormap joining matplotlib colormaps
# Available matplotlib colormaps: https://matplotlib.org/stable/tutorials/colors/colormaps.html

vmin1 = -80                                                       # Min. value
vmax1 = 40                                                        # Max. value

gray_cmap = cm.get_cmap('gray_r', 120)                            # Read the reversed 'gray' cmap
gray_cmap = gray_cmap(np.linspace(0, 1, 120))                     # Create the array
jet_cmap  = cm.get_cmap('jet_r', 40)                              # Read the reversed 'jet' cmap 
jet_cmap  = jet_cmap(np.linspace(0, 1, 40))                       # Create the array
gray_cmap[:40, :] = jet_cmap                                      # Join both cmaps arrays
my_cmap1 = cm.colors.ListedColormap(gray_cmap)                    # Create the custom colormap

#-----------------------------------------------------------------------------------------------------------

# COLORMAP EXAMPLE 2 
# Creating the INPE DISSM IR colormap
# Online color picker: https://imagecolorpicker.com/

vmin2 = -80                                                       # Min. value
vmax2 = 40                                                        # Max. value

gray_cmap = cm.get_cmap('gray_r', 120)                            # Read the reversed 'gray' cmap
gray_cmap = gray_cmap(np.linspace(0, 1, 120))                     # Create the array
colors = ["#ffa0ff", "#0806ff", "#3bcfff", "#feff65", "#ff7516"]  # Custom colors
my_colors = cm.colors.ListedColormap(colors)                      # Create a custom colormap
my_colors = my_colors(np.linspace(0, 1, 50))                      # Create the array
gray_cmap[:50, :] = my_colors                                     # Join both cmaps arrays
my_cmap2 = cm.colors.ListedColormap(gray_cmap)                    # Create the custom colormap 

#-----------------------------------------------------------------------------------------------------------

# COLORMAP EXAMPLE 3
# Creating a linear colormap
# Online color picker: https://imagecolorpicker.com/

colors = ["#bc8462", "#ae656f", "#a44a79", "#962e97", "#6158c5", "#2b8ffb", "#5fcdff", "#94fff0", "#a5ff94", "#fff88c", "#ffbf52", "#ec7b27", "#b84827", "#a1333d", "#bd5478", "#cc6a99", "#d982b8"]
my_cmap3 = cm.colors.LinearSegmentedColormap.from_list("", colors)# Create a custom linear colormap
vmin3 = -80                                                       # Min. value
vmax3 = 40                                                        # Max. value

#-----------------------------------------------------------------------------------------------------------

# COLORMAP EXAMPLE 4
# Converts a CPT file to be used in Python
# CPT archive: http://soliton.vm.bytemark.co.uk/pub/cpt-city/

cpt = loadCPT('IR4AVHRR6.cpt')                                    # Load the CPT file   
my_cmap4 = cm.colors.LinearSegmentedColormap('cpt', cpt)          # Create a custom linear colormap
vmin4 = -103.0                                                    # Min. value
vmax4 = 84.0                                                      # Max. value

#-----------------------------------------------------------------------------------------------------------
# Open the GOES-R image
# Download files at this link: http://home.chpc.utah.edu/~u0553130/Brian_Blaylock/cgi-bin/goes16_download.cgi
file = Dataset("OR_ABI-L2-CMIPF-M6C13_G16_s20191981200396_e20191981210116_c20191981210189.nc")
 
# Get the pixel values
data = file.variables['CMI'][:] - 273.15

#----------------------------------------------------------------------------------------------------------- 

# Choose the plot size (width x height, in inches)
fig, axs = plt.subplots(2,2, figsize=(10,10)) # 2 rows x 2 columns
 
#-----------------------------------------------------------------------------------------------------------
# Plot 1 (first row, first column)

# Plot the image
img1 = axs[0,0].imshow(data, vmin=vmin1, vmax=vmax1, origin='upper', cmap=my_cmap1)
 
# Add a colorbar
plt.colorbar(img1, extend='both', orientation='vertical', pad=0.05, fraction=0.05, ax=axs[0,0])
 
# Add a title
axs[0,0].set_title('Example 1 - Joining 2 cmaps', fontweight='bold', fontsize=10, loc='center')
#-----------------------------------------------------------------------------------------------------------
# Plot 2 (first row, second column)

# Plot the image
img2 = axs[0,1].imshow(data, vmin=vmin2, vmax=vmax2, origin='upper', cmap=my_cmap2)
 
# Add a colorbar
plt.colorbar(img2, label='Brightness Temperatures (°C)', extend='both', orientation='vertical', pad=0.05, fraction=0.05, ax=axs[0,1])
 
# Add a title
axs[0,1].set_title('Example 2 - Custom + cmap', fontweight='bold', fontsize=10, loc='center')
#-----------------------------------------------------------------------------------------------------------
# Plot 3 (second row, first column)

# Plot the image
img3 = axs[1,0].imshow(data, vmin=vmin3, vmax=vmax3, origin='upper', cmap=my_cmap3)
 
# Add a colorbar
plt.colorbar(img3, extend='both', orientation='vertical', pad=0.05, fraction=0.05, ax=axs[1,0])
 
# Add a title
axs[1,0].set_title('Example 3 - Only custom colors', fontweight='bold', fontsize=10, loc='center')
#-----------------------------------------------------------------------------------------------------------
# Plot 4 (second row, second column)

# Plot the image
img4 = axs[1,1].imshow(data, vmin=vmin4, vmax=vmax4, origin='upper', cmap=my_cmap4)
 
# Add a colorbar
plt.colorbar(img4, label='Brightness Temperatures (°C)', extend='both', orientation='vertical', pad=0.05, fraction=0.05, ax=axs[1,1])
 
# Add a title
axs[1,1].set_title('Example 4 - CPT file', fontweight='bold', fontsize=10, loc='center')
#-----------------------------------------------------------------------------------------------------------
# Save the image
plt.savefig('Image_08.png')
 
# Show the image
plt.show()
