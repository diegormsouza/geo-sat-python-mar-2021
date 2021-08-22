# Training: Python and GOES-R Imagery: Script 13 - Cropping the Full Disk and Creating an RGB
#-----------------------------------------------------------------------------------------------------------
# Required modules
from netCDF4 import Dataset              # Read / Write NetCDF4 files
import matplotlib.pyplot as plt          # Plotting library
from datetime import datetime            # Basic Dates and time types
import cartopy, cartopy.crs as ccrs      # Plot maps
import numpy as np                       # Scientific computing with Python
import os                                # Miscellaneous operating system interfaces
from utilities import download_CMI       # Our own utilities
from utilities import geo2grid, convertExtent2GOESProjection      # Our own utilities
#-----------------------------------------------------------------------------------------------------------
# Input and output directories
input = "Samples"; os.makedirs(input, exist_ok=True)
output = "Output"; os.makedirs(output, exist_ok=True)

# AMAZON repository information 
# https://noaa-goes16.s3.amazonaws.com/index.html
bucket_name = 'noaa-goes16'
product_name = 'ABI-L2-CMIPF'
yyyymmddhhmn = '202102181800'

# Desired extent
extent = [-64.0, -36.0, -40.0, -15.0] # Min lon, Max lon, Min lat, Max lat

#-----------------------------------------------------------------------------------------------------------
# Download the necessary bands from AWS
file_ch13 = download_CMI(yyyymmddhhmn, 13, input)
file_ch02 = download_CMI(yyyymmddhhmn, 2, input)
file_ch05 = download_CMI(yyyymmddhhmn, 5, input)

#-----------------------------------------------------------------------------------------------------------
# Open the GOES-R images
file_ch13 = Dataset(f'{input}/{file_ch13}.nc')
file_ch02 = Dataset(f'{input}/{file_ch02}.nc')
file_ch05 = Dataset(f'{input}/{file_ch05}.nc')
#-----------------------------------------------------------------------------------------------------------                   
# Convert lat/lon to grid-coordinates
lly, llx = geo2grid(extent[1], extent[0], file_ch13)
ury, urx = geo2grid(extent[3], extent[2], file_ch13)

# Get the pixel values
data_ch13 = file_ch13.variables['CMI'][ury:lly, llx:urx] - 273.15  
#-----------------------------------------------------------------------------------------------------------
# Convert lat/lon to grid-coordinates
lly, llx = geo2grid(extent[1], extent[0], file_ch02)
ury, urx = geo2grid(extent[3], extent[2], file_ch02)

# Get the pixel values
data_ch02 = file_ch02.variables['CMI'][ury:lly, llx:urx][::4 ,::4] 
#-----------------------------------------------------------------------------------------------------------
# Convert lat/lon to grid-coordinates
lly, llx = geo2grid(extent[1], extent[0], file_ch05)
ury, urx = geo2grid(extent[3], extent[2], file_ch05)

# Get the pixel values
data_ch05 = file_ch05.variables['CMI'][ury:lly, llx:urx][::2 ,::2] 
#-----------------------------------------------------------------------------------------------------------
# Make the arrays equal size
cordX = np.shape(data_ch02)[0], np.shape(data_ch05)[0], np.shape(data_ch13)[0]
cordY = np.shape(data_ch02)[1], np.shape(data_ch05)[1], np.shape(data_ch13)[1]

minvalX = np.array(cordX).min()
minvalY = np.array(cordY).min()

data_ch02 = data_ch02[0:minvalX, 0:minvalY]
data_ch05 = data_ch05[0:minvalX, 0:minvalY]
data_ch13 = data_ch13[0:minvalX, 0:minvalY]
#-----------------------------------------------------------------------------------------------------------
# Compute data-extent in GOES projection-coordinates
img_extent = convertExtent2GOESProjection(extent)              
#-----------------------------------------------------------------------------------------------------------
# RGB Quick Guide: http://rammb.cira.colostate.edu/training/visit/quick_guides/QuickGuide_DayCloudConvectionRGB_final.pdf 

# RGB Components
R = data_ch13
G = data_ch02
B = data_ch05

# Minimuns and Maximuns
Rmin = -53.5
Rmax = 7.5

Gmin = 0.0
Gmax = 0.78

Bmin = 0.01
Bmax = 0.59

R[R<Rmin] = Rmin
R[R>Rmax] = Rmax

G[G<Gmin] = Gmin
G[G>Gmax] = Gmax

B[B<Bmin] = Bmin
B[B>Bmax] = Bmax

# Choose the gamma
gamma = 1

# Normalize the data
R = ((R - Rmax) / (Rmin - Rmax)) ** (1/gamma)
G = ((G - Gmin) / (Gmax - Gmin)) ** (1/gamma)
B = ((B - Bmin) / (Bmax - Bmin)) ** (1/gamma) 

# Create the RGB
RGB = np.stack([R, G, B], axis=2)
#-----------------------------------------------------------------------------------------------------------
# Choose the plot size (width x height, in inches)
plt.figure(figsize=(10,10))

# Use the Geostationary projection in cartopy
ax = plt.axes(projection=ccrs.Geostationary(central_longitude=-75.0, satellite_height=35786023.0))

# Plot the image
img = ax.imshow(RGB, origin='upper', extent=img_extent)

# Add coastlines, borders and gridlines
ax.coastlines(resolution='10m', color='white', linewidth=0.8)
ax.add_feature(cartopy.feature.BORDERS, edgecolor='white', linewidth=0.5)
ax.gridlines(color='white', alpha=0.5, linestyle='--', linewidth=0.5)

# Extract the date
date = (datetime.strptime(file_ch13.time_coverage_start, '%Y-%m-%dT%H:%M:%S.%fZ'))

# Add a title
plt.title('GOES-16 DCP RGB ' + date.strftime('%Y-%m-%d %H:%M') + ' UTC', fontweight='bold', fontsize=10, loc='left')
plt.title('Reg.: ' + str(extent) , fontsize=10, loc='right')
#-----------------------------------------------------------------------------------------------------------
# Save the image
plt.savefig(f'{output}/' + 'DCP_RGB.png', bbox_inches='tight', pad_inches=0, dpi=300)
            
# Show the image
plt.show()
