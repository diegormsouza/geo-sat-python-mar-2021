# Training: Python and GOES-R Imagery: Script 7 - RGB Creation
#----------------------------------------------------------------------------------------------------------- 
# Required modules
from netCDF4 import Dataset                # Read / Write NetCDF4 files
import matplotlib.pyplot as plt            # Plotting library
from datetime import datetime              # Basic Dates and time types
import cartopy, cartopy.crs as ccrs        # Plot maps
import numpy as np                         # Import the Numpy packag
#-----------------------------------------------------------------------------------------------------------
# Open the GOES-R image
# Download files at this link: http://home.chpc.utah.edu/~u0553130/Brian_Blaylock/cgi-bin/goes16_download.cgi
file1 = Dataset("OR_ABI-L2-CMIPF-M6C08_G16_s20191981200396_e20191981210104_c20191981210182.nc")
# Get the pixel values, and convert to Celsius
data1 = file1.variables['CMI'][:,:][::4 ,::4] - 273.15

# Open the GOES-R image
file2 = Dataset("OR_ABI-L2-CMIPF-M6C10_G16_s20191981200396_e20191981210116_c20191981210188.nc")
# Get the pixel values, and convert to Celsius
data2 = file2.variables['CMI'][:,:][::4 ,::4] - 273.15

# Open the GOES-R image
file3 = Dataset("OR_ABI-L2-CMIPF-M6C12_G16_s20191981200396_e20191981210111_c20191981210185.nc")
# Get the pixel values, and convert to Celsius
data3 = file3.variables['CMI'][:,:][::4 ,::4] - 273.15

# Open the GOES-R image
file4 = Dataset("OR_ABI-L2-CMIPF-M6C13_G16_s20191981200396_e20191981210116_c20191981210189.nc")
# Get the pixel values, and convert to Celsius
data4 = file4.variables['CMI'][:,:][::4 ,::4] - 273.15
#-----------------------------------------------------------------------------------------------------------
# RGB Quick Guide: http://rammb.cira.colostate.edu/training/visit/quick_guides/QuickGuide_GOESR_AirMassRGB_final.pdf 

# RGB Components
R = data1 - data2
G = data3 - data4
B = data1
 
# Minimuns, Maximuns and Gamma
Rmin = -26.2
Rmax = 0.6

Gmin = -43.2
Gmax = 6.7

Bmin = -29.25
Bmax = -64.65

R[R > Rmax] = Rmax
R[R < Rmin] = Rmin

G[G > Gmax] = Gmax
G[G < Gmin] = Gmin

B[B < Bmax] = Bmax # Inverted!
B[B > Bmin] = Bmin # Inverted!

gamma_R = 1
gamma_G = 1
gamma_B = 1
 
# Normalize the data
R = ((R - Rmin) / (Rmax - Rmin)) ** (1/gamma_R)
G = ((G - Gmin) / (Gmax - Gmin)) ** (1/gamma_G)
B = ((B - Bmin) / (Bmax - Bmin)) ** (1/gamma_B) 

# Create the RGB
RGB = np.stack([R, G, B], axis=2)

# Eliminate values outside the globe
mask = (RGB == [R[0,0],G[0,0],B[0,0]]).all(axis=2)
RGB[mask] = np.nan
#-----------------------------------------------------------------------------------------------------------
# Choose the plot size (width x height, in inches)
plt.figure(figsize=(7,7)) 
 
# Use the Geostationary projection in cartopy
ax = plt.axes(projection=ccrs.Geostationary(central_longitude=-75.0, satellite_height=35786023.0))
img_extent = (-5434894.67527,5434894.67527,-5434894.67527,5434894.67527)

# Add coastlines, borders and gridlines
ax.coastlines(resolution='10m', color='white', linewidth=0.8)
ax.add_feature(cartopy.feature.BORDERS, edgecolor='white', linewidth=0.8)
ax.gridlines(color='white', alpha=0.5, linestyle='--', linewidth=0.5)
 
# Plot the image
img = ax.imshow(RGB, origin='upper', extent=img_extent)

# Extract the date
date = (datetime.strptime(file1.time_coverage_start, '%Y-%m-%dT%H:%M:%S.%fZ'))
	
# Add a title
plt.title('GOES-16 Airmass RGB ' + date.strftime('%Y-%m-%d %H:%M') + ' UTC', fontweight='bold', fontsize=10, loc='left')
plt.title('Full Disk', fontsize=10, loc='right')
#----------------------------------------------------------------------------------------------------------- 
# Save the image
plt.savefig('Image_07.png')
 
# Show the image
plt.show()
