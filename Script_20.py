# Training: Python and GOES-R Imagery: Script 20 - GLM Density (hotspots)
#-----------------------------------------------------------------------------------------------------------
# Required modules
from netCDF4 import Dataset                           # Read / Write NetCDF4 files
import matplotlib.pyplot as plt                       # Plotting library
import cartopy, cartopy.crs as ccrs                   # Plot maps
from datetime import timedelta, date, datetime        # Basic Dates and time types
import cartopy, cartopy.crs as ccrs                   # Plot maps
import os                                             # Miscellaneous operating system interfaces
from osgeo import gdal                                # Python bindings for GDAL
import numpy as np                                    # Scientific computing with Python
from matplotlib import cm                             # Colormap handling utilities
from utilities import download_CMI, download_GLM      # Our function for download
from utilities import reproject                       # Our function for reproject
from scipy.ndimage.filters import gaussian_filter     # To make a heatmap
gdal.PushErrorHandler('CPLQuietErrorHandler')         # Ignore GDAL warnings
#-----------------------------------------------------------------------------------------------------------
# Input and output directories
input = "Samples"; os.makedirs(input, exist_ok=True)
output = "Output"; os.makedirs(output, exist_ok=True)

# Desired extent
loni = -74
lonf = -34.1
lati = -34.8
latf = 5.5
extent = [loni, lonf, lati, latf] # Min lon, Max lon, Min lat, Max lat

# AMAZON repository information 
# https://noaa-goes16.s3.amazonaws.com/index.html
bucket_name = 'noaa-goes16'
product_name = 'ABI-L2-RRQPEF'
yyyymmddhhmn = '202102081800'

#-----------------------------------------------------------------------------------------------------------
# Get the Band 13 Data

# Download the file
file_ir = download_CMI(yyyymmddhhmn, 13, input)
#-----------------------------------------------------------------------------------------------------------
# Variable
var = 'CMI'

# Open the file
img = gdal.Open(f'NETCDF:{input}/{file_ir}.nc:' + var)

# Data Quality Flag (DQF)
dqf = gdal.Open(f'NETCDF:{input}/{file_ir}.nc:DQF')

# Read the header metadata
metadata = img.GetMetadata()
scale = float(metadata.get(var + '#scale_factor'))
offset = float(metadata.get(var + '#add_offset'))
undef = float(metadata.get(var + '#_FillValue'))
dtime = metadata.get('NC_GLOBAL#time_coverage_start')

# Load the data
ds_cmi = img.ReadAsArray(0, 0, img.RasterXSize, img.RasterYSize).astype(float)
ds_dqf = dqf.ReadAsArray(0, 0, dqf.RasterXSize, dqf.RasterYSize).astype(float)

# Apply the scale, offset and convert to celsius
ds_cmi = (ds_cmi * scale + offset) - 273.15

# Apply NaN's where the quality flag is greater than 1
ds_cmi[ds_dqf > 1] = np.nan

# Reproject the file
filename_ret = f'{output}/IR_{yyyymmddhhmn}.nc'
reproject(filename_ret, img, ds_cmi, extent, undef)

# Open the reprojected GOES-R image
file = Dataset(filename_ret)

# Get the pixel values
data = file.variables['Band1'][:]
#-----------------------------------------------------------------------------------------------------------
# Get the GLM Data

# Initialize arrays for latitude, longitude, and event energy
lats = np.array([])
lons = np.array([])
#-----------------------------------------------------------------------------------------------------------
# Initial time and date
yyyy = datetime.strptime(yyyymmddhhmn, '%Y%m%d%H%M').strftime('%Y')
mm = datetime.strptime(yyyymmddhhmn, '%Y%m%d%H%M').strftime('%m')
dd = datetime.strptime(yyyymmddhhmn, '%Y%m%d%H%M').strftime('%d')
hh = datetime.strptime(yyyymmddhhmn, '%Y%m%d%H%M').strftime('%H')
mn = datetime.strptime(yyyymmddhhmn, '%Y%m%d%H%M').strftime('%M')

date_ini = str(datetime(int(yyyy),int(mm),int(dd),int(hh),int(mn)) - timedelta(minutes=10))
date_end = str(datetime(int(yyyy),int(mm),int(dd),int(hh),int(mn)))

# GLM accumulation loop
while (date_ini <= date_end):
 
    # Date structure
    yyyymmddhhmnss = datetime.strptime(date_ini, '%Y-%m-%d %H:%M:%S').strftime('%Y%m%d%H%M%S')

    # Download the file
    file_glm = download_GLM(yyyymmddhhmnss, input)

    # Read the file
    glm = Dataset(f'{input}/{file_glm}.nc')

    # Append lats / longs / event energies
    lats = np.append(lats, glm.variables['event_lat'][:])
    lons = np.append(lons, glm.variables['event_lon'][:])
    
    # Increment the date_ini
    date_ini = str(datetime.strptime(date_ini, '%Y-%m-%d %H:%M:%S') + timedelta(seconds=20))
#-----------------------------------------------------------------------------------------------------------
# Stack and transpose the lat lons
values = np.vstack((lons,lats)).T

values = values[(values[:,0] >= loni) & (values[:,0] <= lonf)]
values = values[(values[:,1] >= lati) & (values[:,1] <= latf)]

from scipy.ndimage.filters import gaussian_filter
heatmap, x, y = np.histogram2d(values[:,0], values[:,1], bins=(120,120), range=[[loni,lonf],[lati,latf]])
heatmap = gaussian_filter(heatmap,1)

print(heatmap.min(), heatmap.max(), heatmap.shape)

#-----------------------------------------------------------------------------------------------------------
# Choose the plot size (width x height, in inches)
plt.figure(figsize=(10,10))

# Use the Geostationary projection in cartopy
ax = plt.axes(projection=ccrs.PlateCarree())

# Define the image extent
ax.set_extent([extent[0], extent[2], extent[1], extent[3]], ccrs.PlateCarree())

# Define the data extent
img_extent = [extent[0], extent[2], extent[1], extent[3]]

# Plot the image
img = ax.imshow(data, vmin=-50, vmax=80, cmap='gray_r', origin='upper', extent=img_extent)

# Plot the GLM Data
glm = ax.imshow(heatmap.T, vmin=0, vmax=600, origin="lower", cmap="hot", interpolation="nearest", extent=img_extent, alpha=0.6)

# Add coastlines, borders and gridlines
ax.coastlines(resolution='10m', color='white', linewidth=0.8, zorder=3)
ax.add_feature(cartopy.feature.BORDERS, edgecolor='white', linewidth=0.5, zorder=4)
gl = ax.gridlines(crs=ccrs.PlateCarree(), color='gray', alpha=1.0, linestyle='--', linewidth=0.25, xlocs=np.arange(-180, 180, 5), ylocs=np.arange(-90, 90, 5), draw_labels=True)
gl.top_labels = False
gl.right_labels = False

# Add the img colorbar
#plt.colorbar(img, label='Brightness Temperatures (°C)', extend='both', orientation='horizontal', pad=0.05, fraction=0.05)

# Add the glm colorbar
plt.colorbar(glm, label='GLM Density', extend='max', orientation='vertical', pad=0.05, fraction=0.05)

# Add a title
plt.title('GOES-16 IR+GLM', fontweight='bold', fontsize=10, loc='left')
date_ini = str(datetime(int(yyyy),int(mm),int(dd),int(hh),int(mn)) - timedelta(minutes=10))
date_end = str(datetime(int(yyyy),int(mm),int(dd),int(hh),int(mn)))
plt.title(str(date_ini) + " - " + str(date_end), fontsize="10", loc="right")

#-----------------------------------------------------------------------------------------------------------
# Save the image
plt.savefig(f'{output}/GLM_DENS_HOT_{yyyymmddhhmn}.png', bbox_inches='tight', pad_inches=0, dpi=300)

# Show the image
plt.show()
