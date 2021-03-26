# Training: Python and GOES-R Imagery: Script 6 - ABI and GLM 
#----------------------------------------------------------------------------------------------------------- 
# Required modules
from netCDF4 import Dataset          # Read / Write NetCDF4 files
import matplotlib.pyplot as plt      # Plotting library
from datetime import datetime        # Basic Dates and time types
import cartopy, cartopy.crs as ccrs  # Plot maps
#----------------------------------------------------------------------------------------------------------- 
# Open the GOES-R image
# Download files at this link: http://home.chpc.utah.edu/~u0553130/Brian_Blaylock/cgi-bin/goes16_download.cgi
file = Dataset("OR_ABI-L2-CMIPF-M6C13_G16_s20191981200396_e20191981210116_c20191981210189.nc")
 
# Get the pixel values
data = file.variables['CMI'][:] - 273.15
 
# Choose the plot size (width x height, in inches)
plt.figure(figsize=(7,7))
 
# Use the Geostationary projection in cartopy
ax = plt.axes(projection=ccrs.Geostationary(central_longitude=-75.0, satellite_height=35786023.0))
img_extent = (-5434894.67527,5434894.67527,-5434894.67527,5434894.67527)
 
# Add coastlines, borders and gridlines
ax.coastlines(resolution='10m', color='white', linewidth=0.8)
ax.add_feature(cartopy.feature.BORDERS, edgecolor='white', linewidth=0.5)
ax.gridlines(color='white', alpha=0.5, linestyle='--', linewidth=0.5)
 
# Plot the image
img = ax.imshow(data, vmin=-80, vmax=40, origin='upper', extent=img_extent, cmap='Greys')
#-----------------------------------------------------------------------------------------------------------
# Open the GLM file
# Download files at this link: http://home.chpc.utah.edu/~u0553130/Brian_Blaylock/cgi-bin/goes16_download.cgi
fileGLM = Dataset("OR_GLM-L2-LCFA_G16_s20191981200000_e20191981200200_c20191981200224.nc") 
 
e_lats = fileGLM.variables['event_lat'][:]
e_lons = fileGLM.variables['event_lon'][:]

g_lats = fileGLM.variables['group_lat'][:]
g_lons = fileGLM.variables['group_lon'][:]

f_lats = fileGLM.variables['flash_lat'][:]
f_lons = fileGLM.variables['flash_lon'][:] 

img2 = ax.plot(e_lons,e_lats,'.r', markersize=10, transform=ccrs.PlateCarree(), alpha=0.01)
img3 = ax.plot(g_lons,g_lats,'.y', markersize=5, transform=ccrs.PlateCarree(), alpha=0.5)
img4 = ax.plot(f_lons,f_lats,'.g', markersize=2.5, transform=ccrs.PlateCarree(), alpha=1)
#-----------------------------------------------------------------------------------------------------------
# Add a colorbar
plt.colorbar(img, label='Brightness Temperatures (°C)', extend='both', orientation='horizontal', pad=0.05, fraction=0.05)
 
# Extract the date
date = (datetime.strptime(file.time_coverage_start, '%Y-%m-%dT%H:%M:%S.%fZ'))
	
# Add a title
plt.title('GOES-16 Band 13 and GLM ' + date.strftime('%Y-%m-%d %H:%M') + ' UTC', fontweight='bold', fontsize=10, loc='left')
plt.title('Full Disk', fontsize=10, loc='right')
#----------------------------------------------------------------------------------------------------------- 
# Save the image
plt.savefig('Image_06.png')
 
# Show the image
plt.show()
