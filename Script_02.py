# Training: Python and GOES-R Imagery: Script 2 - Basic Operation / Colorbar / Title / Date
#----------------------------------------------------------------------------------------------------------- 
# Required modules
from netCDF4 import Dataset      # Read / Write NetCDF4 files
import matplotlib.pyplot as plt  # Plotting library
from datetime import datetime    # Basic Dates and time types
#----------------------------------------------------------------------------------------------------------- 
# Open the GOES-R image
# Download files at this link: http://home.chpc.utah.edu/~u0553130/Brian_Blaylock/cgi-bin/goes16_download.cgi
file = Dataset("OR_ABI-L2-CMIPF-M6C13_G16_s20191981200396_e20191981210116_c20191981210189.nc")
 
# Get the pixel values
data = file.variables['CMI'][:] - 273.15
#-----------------------------------------------------------------------------------------------------------
# Choose the plot size (width x height, in inches)
plt.figure(figsize=(7,7))
 
# Plot the image
plt.imshow(data, vmin=-80, vmax=40, cmap='jet')
 
# Add a colorbar
plt.colorbar(label='Brightness Temperature (°C)', extend='both', orientation='horizontal', pad=0.05, fraction=0.05)

# Extract the date
date = (datetime.strptime(file.time_coverage_start, '%Y-%m-%dT%H:%M:%S.%fZ'))
	
# Add a title
plt.title('GOES-16 Band 13 ' + date.strftime('%Y-%m-%d %H:%M') + ' UTC', fontweight='bold', fontsize=10, loc='left')
plt.title('Full Disk', fontsize=10, loc='right')
#----------------------------------------------------------------------------------------------------------- 
# Save the image
plt.savefig('Image_02.png')
 
# Show the image
plt.show()

















# Algumas ideias para modificar:
# O colormap
# O título
# Posição da colorbar (vertical)
# Colocar nan
#data[data < 23] = np.nan 
# Topos de nuvens mais frios
#data[data > -30] = np.nan
# Apenas dados entre dois valores
#import numpy as np
#data[np.logical_or(data < 5, data > 15)] = np.nan
 

